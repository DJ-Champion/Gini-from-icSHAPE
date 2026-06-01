"""lib/gini_core.py

Core, bedtools-independent logic for the Gini-from-icSHAPE pipeline:

  * population Gini (matches the source paper; no n-1 correction)
  * GFF region extraction (exon / CDS / five_prime_UTR / three_prime_UTR
    read verbatim, projected GFF 1-based-inclusive -> BED 0-based half-open)
  * parsing/expanding `bedtools intersect -wo` output into per-base score
    vectors, with overlap-bp clipping and [0,1] score validation

Everything here is pure Python and unit-testable without bedtools installed.
The driver (run_gini.py) handles sorting, the bedtools call, and orchestration.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

log = logging.getLogger('gini.core')

# Region names used throughout. 'full' is built from exon features; the
# other three from the correspondingly-named GFF feature types.
REGIONS = ('full', '5UTR', 'CDS', '3UTR')

# Map region name -> GFF feature type(s) whose union defines it.
_REGION_FEATURETYPES = {
    'full': ('exon',),
    'CDS': ('CDS',),
    '5UTR': ('five_prime_UTR',),
    '3UTR': ('three_prime_UTR',),
}

# Minimum valid per-base observations per (transcript, region) to emit a
# Gini. Matches the paper's icSHAPE.computeGINI valid_cutoff > 15.
VALID_CUTOFF = 15


# ---------------------------------------------------------------------------
# Gini
# ---------------------------------------------------------------------------

def population_gini(scores: list[float]) -> float:
    """Population Gini coefficient of a list of non-negative scores.

    Uses the sorted-cumsum form:
        sort ascending; n = len; S = cumulative sum
        G = (n + 1 - 2 * sum(S) / S[-1]) / n

    No n-1 sample correction (matches the source paper). Caller is
    responsible for the n > VALID_CUTOFF and zero-sum filters; this
    function assumes a non-empty vector with a positive sum and will
    raise ZeroDivisionError otherwise (a deliberate fail-loud, since the
    driver must have filtered those cases already).

    Verified anchors (see tests):
        [0, 1]     -> 0.5
        [0, 0, 1]  -> 2/3
        [1, 1, 1]  -> 0.0
    """
    s = sorted(scores)
    n = len(s)
    cum = 0.0
    cumsum_of_cumsum = 0.0
    for x in s:
        cum += x
        cumsum_of_cumsum += cum
    total = cum  # S[-1]
    return (n + 1 - 2 * cumsum_of_cumsum / total) / n


def population_gini_weighted(pairs: list[tuple[float, int]]) -> float:
    """Population Gini on a run-length-encoded multiset of scores.

    `pairs` is a list of (value, count): the same scores as
    population_gini([...]) would take, but with repeats collapsed to a
    count. This is mathematically identical to expanding each (value,
    count) into `count` copies and calling population_gini — but it never
    instantiates the repeats, so memory scales with the number of distinct
    RASP intervals per (transcript, region), not the number of bases.

    Derivation: within a run of `c` copies of value v, sitting on a prefix
    of summed-value P, the cumulative sums S_i are P+v, P+2v, ..., P+cv, so
    the run contributes  c*P + v*c*(c+1)/2  to sum(S). Total count n and
    grand total S[-1] are the weighted sums.

    Same contract as population_gini: assumes total > 0 (caller filters
    zero-sum) and n > 0.
    """
    ps = sorted(pairs)
    n = 0
    total = 0.0
    for v, c in ps:
        n += c
        total += v * c
    cumsum_of_cumsum = 0.0
    prefix = 0.0
    for v, c in ps:
        cumsum_of_cumsum += c * prefix + v * c * (c + 1) / 2.0
        prefix += v * c
    return (n + 1 - 2 * cumsum_of_cumsum / total) / n


# ---------------------------------------------------------------------------
# GFF region extraction
# ---------------------------------------------------------------------------

@dataclass
class RegionBedRecord:
    """One BED interval for a (transcript, region). 0-based half-open."""
    chrom: str
    start: int          # BED start (0-based)
    end: int            # BED end (exclusive)
    transcript_id: str
    region: str
    strand: str


def _gff_to_bed_interval(feat) -> tuple[int, int]:
    """GFF 1-based inclusive [start, end] -> BED 0-based half-open [start, end)."""
    # gffutils exposes .start (1-based inclusive) and .end (1-based inclusive).
    return feat.start - 1, feat.end


def extract_regions_for_transcript(
    db,
    transcript_feat,
) -> tuple[list[RegionBedRecord], list[str]]:
    """Return (records, missing_regions) for one transcript.

    `records` are BED intervals across all annotated regions for this
    transcript. `missing_regions` lists region names that produced no
    feature (e.g. CDS/UTRs on a non-coding transcript, or a coding
    transcript lacking an annotated UTR) so the driver can log them.

    'full' is computed independently from exon features (not as the union
    of the sub-regions), so it is robust to annotation quirks where
    UTR + CDS != exon union.
    """
    tid = transcript_feat.id
    strand = transcript_feat.strand
    records: list[RegionBedRecord] = []
    missing: list[str] = []

    for region in REGIONS:
        ftypes = _REGION_FEATURETYPES[region]
        feats = []
        for ft in ftypes:
            feats.extend(
                db.children(transcript_feat, featuretype=ft, order_by='start')
            )
        if not feats:
            missing.append(region)
            continue
        for f in feats:
            bstart, bend = _gff_to_bed_interval(f)
            if bend <= bstart:
                # Defensive: skip degenerate intervals rather than emit them.
                continue
            records.append(RegionBedRecord(
                chrom=f.seqid,
                start=bstart,
                end=bend,
                transcript_id=tid,
                region=region,
                strand=strand,
            ))
    return records, missing


def write_region_bed(records: Iterable[RegionBedRecord], path: Path) -> int:
    """Write region records to a 6-column BED. Returns row count.

    Column layout (so bedtools -s finds strand in col 6, and the
    name column carries transcript|region for downstream keying):
        chrom  start  end  transcript_id|region  .  strand
    """
    n = 0
    with open(path, 'w') as fh:
        for r in records:
            name = f"{r.transcript_id}|{r.region}"
            fh.write(f"{r.chrom}\t{r.start}\t{r.end}\t{name}\t.\t{r.strand}\n")
            n += 1
    return n


# ---------------------------------------------------------------------------
# Intersect output parsing / per-base expansion
# ---------------------------------------------------------------------------

@dataclass
class ScoreAccumulator:
    """Run-length per-base score observations keyed by (transcript_id, region).

    Instead of materialising `count` copies of each score (which blows up
    memory for large RASP files — a single intersect row can cover dozens
    of bases, and there are millions of rows), we hold a {score: total_count}
    map per key. This is exactly the multiset the per-base vector would be,
    just run-length encoded. population_gini_weighted consumes it directly.
    """
    counts: dict[tuple[str, str], dict[float, int]] = field(default_factory=dict)
    out_of_range_dropped: int = 0

    def add(self, key: tuple[str, str], score: float, count: int) -> None:
        m = self.counts.setdefault(key, {})
        m[score] = m.get(score, 0) + count

    def merge_from(self, other: 'ScoreAccumulator') -> None:
        """Fold another accumulator into this one (used if combining regions)."""
        for key, m in other.counts.items():
            dst = self.counts.setdefault(key, {})
            for score, c in m.items():
                dst[score] = dst.get(score, 0) + c
        self.out_of_range_dropped += other.out_of_range_dropped


def parse_intersect_line(
    line: str,
    acc: ScoreAccumulator,
    region_bed_cols: int = 6,
    rasp_cols: int = 6,
) -> None:
    """Parse one `bedtools intersect -a region -b rasp -s -wo` line into acc.

    With `-a` having `region_bed_cols` columns and `-b` having `rasp_cols`
    columns, each `-wo` line is:

        [region 6 cols] [rasp 6 cols] [overlap_bp]

    so the RASP score is field (region_bed_cols + 4) and the overlap base
    count is the final field. The region name (col 4 of -a) is
    'transcript_id|region'. The overlap_bp is the *clipped* overlap width
    (already the in-region base count), so it is used directly as the
    run-length count — clip-and-expand in one step, with no expansion.

    Scores outside [0, 1] (or unparseable) are dropped; their base count is
    added to acc.out_of_range_dropped for the caller to log.
    """
    line = line.rstrip('\n')
    if not line:
        return
    f = line.split('\t')
    name_idx = 3
    score_idx = region_bed_cols + 4
    overlap_idx = region_bed_cols + rasp_cols
    raw_score = f[score_idx]
    overlap_bp = int(f[overlap_idx])
    try:
        score = float(raw_score)
    except ValueError:
        acc.out_of_range_dropped += overlap_bp
        return
    if not (0.0 <= score <= 1.0):
        acc.out_of_range_dropped += overlap_bp
        return
    name = f[name_idx]
    if '|' in name:
        tid, region = name.split('|', 1)
    else:
        tid, region = name, 'full'
    acc.add((tid, region), score, overlap_bp)


def parse_intersect_wo(
    lines: Iterable[str],
    region_bed_cols: int = 6,
    rasp_cols: int = 6,
) -> ScoreAccumulator:
    """Iterable convenience wrapper around parse_intersect_line.

    Builds and returns a fresh ScoreAccumulator. The driver streams via
    parse_intersect_line directly to avoid holding all lines at once; this
    form is kept for tests and small in-memory uses.
    """
    acc = ScoreAccumulator()
    for line in lines:
        parse_intersect_line(line, acc, region_bed_cols, rasp_cols)
    return acc


# ---------------------------------------------------------------------------
# Aggregation -> Gini rows
# ---------------------------------------------------------------------------

@dataclass
class GiniRow:
    transcript_id: str
    gene_id: str
    species: str
    compartment: str
    region: str
    gini: float
    n_valid: int


@dataclass
class DropoutRow:
    transcript_id: str
    species: str
    compartment: str
    region: str
    reason: str
    n_valid: int


def gini_rows_from_accumulator(
    acc: ScoreAccumulator,
    species: str,
    compartment: str,
    tx_to_gene: dict[str, str],
) -> tuple[list[GiniRow], list[DropoutRow]]:
    """Apply filters and compute Gini per (transcript, region).

    Filters, in order:
      * n_valid <= VALID_CUTOFF        -> dropout reason 'n_valid<=15'
      * sum(scores) == 0 (zero-sum)    -> dropout reason 'zero_sum'
    Survivors get a population Gini and a GiniRow.

    `tx_to_gene` maps the transcript_id (as written into the region BED,
    i.e. the gffutils feature id) to an unversioned gene_id for the output.
    """
    gini_rows: list[GiniRow] = []
    dropouts: list[DropoutRow] = []

    for (tid, region), counts in acc.counts.items():
        gene_id = tx_to_gene.get(tid, '')
        n = sum(counts.values())
        if n <= VALID_CUTOFF:
            dropouts.append(DropoutRow(
                tid, species, compartment, region, 'n_valid<=15', n))
            continue
        total = sum(v * c for v, c in counts.items())
        if total == 0:
            dropouts.append(DropoutRow(
                tid, species, compartment, region, 'zero_sum', n))
            continue
        g = population_gini_weighted(list(counts.items()))
        gini_rows.append(GiniRow(
            transcript_id=tid,
            gene_id=gene_id,
            species=species,
            compartment=compartment,
            region=region,
            gini=g,
            n_valid=n,
        ))
    return gini_rows, dropouts
