"""lib/gff.py
Shared GFF / gffutils helpers used by metric plugins.

Centralises:
  * GFF3 ID normalisation (strip namespace prefix, strip version suffix)
    so manifest IDs and GFF IDs reconcile across Ensembl/GENCODE/RefSeq.
  * gffutils SQLite DB open-or-build with mtime check.
  * Transcript index built from 'features that are Parents of exon
    features' — handles mRNA / transcript / primary_transcript / lnc_RNA /
    ncRNA / etc. without enumeration.
  * Manifest reader.
"""
from __future__ import annotations

import csv
import logging
import re
from pathlib import Path

log = logging.getLogger('lib.gff')


# ---------------------------------------------------------------------------
# ID normalisation
# ---------------------------------------------------------------------------

# Known namespace prefixes seen in GFF3 ID attributes:
#   Ensembl/GENCODE: 'transcript:ENST...', 'gene:ENSG...'
#   NCBI/RefSeq:     'rna-NM_...', 'gene-FOO' (also 'rna-XM_', 'rna-NR_', etc.)
# Conservative list — only patterns that uniquely identify a namespace
# prefix and won't false-strip anything that could appear inside a real ID.
KNOWN_PREFIXES: tuple[str, ...] = ('transcript:', 'gene:', 'rna-', 'gene-')

# Trailing version suffix: '.<digits>' at end of string (e.g. ENST00000000001.7)
_VERSION_RE = re.compile(r'\.\d+$')


def strip_prefix(s: str) -> str:
    """Strip a known GFF3 ID namespace prefix; return s unchanged if none match."""
    s_lower = s.lower()
    for pfx in KNOWN_PREFIXES:
        if s_lower.startswith(pfx):
            return s[len(pfx):]
    return s


def strip_version(s: str) -> str:
    """Strip a trailing '.<digits>' version suffix.
    Only strips dot-then-pure-digits-to-end so e.g. 'ENST00000000001.7' -> 'ENST00000000001'
    but 'gene.alpha' is preserved.
    """
    return _VERSION_RE.sub('', s)


def normalise_id(raw: str) -> str:
    """Maximal normalisation: strip namespace prefix AND trailing version."""
    return strip_version(strip_prefix(raw))


# ---------------------------------------------------------------------------
# gffutils DB handling
# ---------------------------------------------------------------------------

def open_or_build_db(gff_path: Path, db_path: Path):
    """Open an existing gffutils DB or build it fresh.

    Rebuild if the GFF is newer than the DB (defensive: handles a manual
    re-run of 01_extract.py that produced a new canonical.gff).

    gffutils is imported lazily so plugins that don't need a DB aren't
    forced to declare the dependency.
    """
    try:
        import gffutils
    except ImportError:
        log.error("gffutils not installed (pip install gffutils).")
        raise

    needs_build = (
        not db_path.exists()
        or db_path.stat().st_mtime < gff_path.stat().st_mtime
    )
    if needs_build:
        if db_path.exists():
            log.info(f"GFF newer than DB; rebuilding {db_path.name}")
            db_path.unlink()
        else:
            log.info(f"Building gffutils DB: {db_path}")
        # Use gffutils defaults for ID assignment. Custom id_spec remapping
        # of mRNA IDs would break Parent= references on child features
        # (exons / CDS / UTR), since those refer to the *original* IDs in
        # the GFF. ID mismatches between manifest and GFF (version, prefix)
        # are handled in lookup_transcript instead.
        db = gffutils.create_db(
            str(gff_path),
            dbfn=str(db_path),
            force=False,
            keep_order=False,
            merge_strategy='create_unique',
            sort_attribute_values=False,
        )
    else:
        log.debug(f"Loading existing DB: {db_path}")
        db = gffutils.FeatureDB(str(db_path))
    return db


# ---------------------------------------------------------------------------
# Transcript index
# ---------------------------------------------------------------------------

def _index_keys_for_feature(feat) -> list[str]:
    """Generate every plausible lookup key for a transcript-like feature.

    Covers ID forms (raw, prefix-stripped, fully normalised), the
    transcript_id attribute (same three forms), and a reconstructed
    versioned form using a separate transcript_version / version attribute
    if present (Ensembl convention where ID is unversioned but a separate
    attribute carries the version).
    """
    keys = []
    fid = feat.id
    keys.extend([fid, strip_prefix(fid), normalise_id(fid)])

    for tid in feat.attributes.get('transcript_id', []):
        keys.extend([tid, strip_prefix(tid), normalise_id(tid)])
        # Reconstruct 'BARE.VERSION' if version stored separately
        bare = normalise_id(tid)
        for vattr in ('transcript_version', 'version'):
            for v in feat.attributes.get(vattr, []):
                keys.append(f"{bare}.{v}")

    # Deduplicate while preserving order; drop empties
    seen = set()
    out = []
    for k in keys:
        if k and k not in seen:
            seen.add(k)
            out.append(k)
    return out


def build_transcript_index(db) -> dict:
    """Map every plausible ID form to the gffutils feature.id of its transcript.

    Discovers transcripts as 'features that are Parents of exon features'.
    This catches mRNA, transcript, primary_transcript, lnc_RNA, ncRNA, etc.
    without needing to enumerate feature types up front.

    First-write-wins on collisions: if two features share a normalised ID,
    the first one indexed keeps the slot. Real annotations don't collide
    here once the canonical.gff filter has trimmed to one transcript per gene.
    """
    parent_ids = set()
    for exon in db.features_of_type('exon'):
        for pid in exon.attributes.get('Parent', []):
            parent_ids.add(pid)

    index = {}
    feat_types = set()
    for pid in parent_ids:
        try:
            feat = db[pid]
        except Exception:
            continue
        feat_types.add(feat.featuretype)
        for key in _index_keys_for_feature(feat):
            index.setdefault(key, feat.id)

    log.debug(f"Indexed {len(parent_ids)} transcript-like features "
              f"(types: {sorted(feat_types)}) under {len(index)} keys")
    return index


def lookup_transcript(db, manifest_tx_id: str, index: dict):
    """Look up a transcript by trying several forms of the manifest ID."""
    for key in (manifest_tx_id,
                strip_prefix(manifest_tx_id),
                normalise_id(manifest_tx_id)):
        feat_id = index.get(key)
        if feat_id is not None:
            try:
                return db[feat_id]
            except Exception:
                pass
    return None


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def read_manifest_transcripts(manifest_tsv: Path) -> list[str]:
    """Read manifest.tsv; return unique transcript_ids preserving first-seen order."""
    seen = set()
    ordered = []
    with open(manifest_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'transcript_id' not in (reader.fieldnames or []):
            raise ValueError(
                f"manifest.tsv has no 'transcript_id' column; got {reader.fieldnames}")
        for row in reader:
            tid = row['transcript_id']
            if tid and tid not in seen:
                seen.add(tid)
                ordered.append(tid)
    return ordered

# ---------------------------------------------------------------------------
# Fasta header parsing
# ---------------------------------------------------------------------------

def split_composite_fasta_id(
    record_id: str,
    region: str | None = None,
) -> tuple[str, str, str] | None:
    """Split a composite FASTA header into (gene_id, transcript_id, region).

    01_extract.py writes FASTA headers in the form
    '<gene_id>_<transcript_id>_<region>' (e.g.
    'ENSG00000137075_ENST00000259605.11_CDS'). This helper extracts the
    three components.

    If `region` is provided, the function strips that exact suffix. This
    is the recommended path for plugin code: each plugin knows which
    region's FASTA it's reading, so passing it as a hint is robust to
    multi-token region names (notably 'UTR_pair').

    Without a region hint, falls back to assuming the region is the final
    '_'-delimited token, with everything between the first '_' and the
    final '_' as the transcript_id. This handles unknown regions but
    misparses 'UTR_pair'-style names.

    Returns None if the header doesn't conform to the composite pattern
    (e.g. no underscores, or fewer than the expected parts).
    """
    if region is not None:
        suffix = '_' + region
        if not record_id.endswith(suffix):
            return None
        rest = record_id[:-len(suffix)]
        if '_' not in rest:
            return None
        gene, _, transcript = rest.partition('_')
        if not gene or not transcript:
            return None
        return gene, transcript, region

    parts = record_id.split('_')
    if len(parts) < 3:
        return None
    return parts[0], '_'.join(parts[1:-1]), parts[-1]
