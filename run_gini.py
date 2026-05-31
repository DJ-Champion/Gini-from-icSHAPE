#!/usr/bin/env python3
"""run_gini.py

Driver for the Gini-from-icSHAPE-reactivity pipeline.

Per species: build (cached) a gffutils DB from the reduced one-transcript-
per-gene GFF, enumerate transcripts via exon parents, write one region BED
per region (full / 5UTR / CDS / 3UTR). Per compartment: sort-and-cache the
RASP score BED, then for each region run

    bedtools intersect -a <region>.bed -b <rasp.sorted>.bed -s -wo

parse the overlap output into per-base score vectors (overlap-bp clipped),
filter (n_valid > 15, drop zero-sum), and compute population Gini.

Outputs:
  * <out>/gini_long.tsv   columns: transcript_id gene_id species
                          compartment region gini n_valid
  * <out>/dropouts.tsv    columns: transcript_id species compartment
                          region reason n_valid

No SLURM, no half-life join (done downstream in R). bedtools is shelled out
to via subprocess; the only Python dependency is gffutils.

Usage:
    python run_gini.py --config config.yaml
or with explicit flags (see argparse below). A YAML config is the tidy path;
flags are there for quick one-offs.
"""
from __future__ import annotations

import argparse
import csv
import logging
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

# Make the sibling lib/ importable whether run from repo root or bin/.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib import gff as gfflib            # noqa: E402  (DJ's helper)
from lib import gini_core as core        # noqa: E402

log = logging.getLogger('gini.driver')


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

@dataclass
class SpeciesConfig:
    name: str                       # 'human' / 'mouse'
    gff: Path                       # reduced one-tx-per-gene GFF
    compartments: dict[str, Path]   # compartment name -> RASP score BED


def load_config(path: Path) -> list[SpeciesConfig]:
    """Load a small YAML config. Structure:

    species:
      human:
        gff: /path/MANE.reduced.gff
        compartments:
          nucleoplasm: /path/icSHAPE_hg38_nucl_vivo-score.bed
          cytoplasm:   /path/icSHAPE_hg38_cy_invivo-score.bed
      mouse:
        gff: /path/vM25.reduced.gff3
        compartments:
          nucleoplasm: /path/icSHAPE_mm10_nucl_vivo-score.bed
          cytoplasm:   /path/icSHAPE_mm10_cy_vivo-score.bed
    """
    import yaml
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    out = []
    for sp_name, sp in raw['species'].items():
        out.append(SpeciesConfig(
            name=sp_name,
            gff=Path(sp['gff']),
            compartments={k: Path(v) for k, v in sp['compartments'].items()},
        ))
    return out


# ---------------------------------------------------------------------------
# Region BED construction (once per species)
# ---------------------------------------------------------------------------

def build_region_beds(
    sp: SpeciesConfig,
    work_dir: Path,
) -> tuple[dict[str, Path], dict[str, str], list[core.DropoutRow]]:
    """Build one BED per region for this species.

    Returns:
      region_bed_paths : region name -> BED path
      tx_to_gene       : transcript feature id -> unversioned gene id
      missing_dropouts : DropoutRow(reason='no_feature') for absent regions
                         (logged once per species; compartment left blank).
    """
    db_path = work_dir / f"{sp.name}.gffutils.db"
    db = gfflib.open_or_build_db(sp.gff, db_path)

    # Enumerate transcripts the same way the helper indexes them: features
    # that are parents of exons. This is the transcript universe.
    parent_ids = set()
    for exon in db.features_of_type('exon'):
        for pid in exon.attributes.get('Parent', []):
            parent_ids.add(pid)

    records_by_region: dict[str, list[core.RegionBedRecord]] = {
        r: [] for r in core.REGIONS
    }
    tx_to_gene: dict[str, str] = {}
    missing_dropouts: list[core.DropoutRow] = []

    n_tx = 0
    for pid in sorted(parent_ids):
        try:
            tx = db[pid]
        except Exception:
            continue
        n_tx += 1
        # gene_id: prefer the attribute, fall back to Parent; strip version.
        gene_raw = (tx.attributes.get('gene_id', [None])[0]
                    or tx.attributes.get('Parent', [None])[0]
                    or '')
        tx_to_gene[tx.id] = gfflib.strip_version(gfflib.strip_prefix(gene_raw))

        recs, missing = core.extract_regions_for_transcript(db, tx)
        for r in recs:
            records_by_region[r.region].append(r)
        for region in missing:
            missing_dropouts.append(core.DropoutRow(
                transcript_id=tx.id, species=sp.name, compartment='',
                region=region, reason='no_feature', n_valid=0))

    region_bed_paths: dict[str, Path] = {}
    for region, recs in records_by_region.items():
        # Sort the region BED so the later bedtools intersect is valid.
        recs.sort(key=lambda r: (r.chrom, r.start, r.end))
        p = work_dir / f"{sp.name}.{region}.bed"
        n = core.write_region_bed(recs, p)
        region_bed_paths[region] = p
        log.info(f"[{sp.name}] region {region}: {n} intervals -> {p.name}")

    log.info(f"[{sp.name}] {n_tx} transcripts; "
             f"{len(missing_dropouts)} missing-region dropouts")
    return region_bed_paths, tx_to_gene, missing_dropouts


# ---------------------------------------------------------------------------
# RASP sorting (once per file) + intersect (per region x compartment)
# ---------------------------------------------------------------------------

def sort_rasp(rasp_path: Path, work_dir: Path) -> Path:
    """Coordinate-sort a RASP BED to a cached copy; skip if already current."""
    out = work_dir / f"{rasp_path.stem}.sorted.bed"
    if out.exists() and out.stat().st_mtime >= rasp_path.stat().st_mtime:
        log.debug(f"sorted RASP cache hit: {out.name}")
        return out
    log.info(f"sorting RASP {rasp_path.name} -> {out.name}")
    with open(out, 'w') as fh:
        # LC_ALL=C for byte-wise sort matching bedtools' expectation.
        subprocess.run(
            ['sort', '-k1,1', '-k2,2n', str(rasp_path)],
            stdout=fh, check=True, env={'LC_ALL': 'C', 'PATH': _path()},
        )
    return out


def run_intersect(region_bed: Path, rasp_sorted: Path) -> list[str]:
    """bedtools intersect -a region -b rasp -s -wo; return stdout lines.

    -s : require same strand (essential; antisense bleed otherwise)
    -wo: write A, B, and the number of overlapping bases (the clip width)
    """
    proc = subprocess.run(
        ['bedtools', 'intersect',
         '-a', str(region_bed), '-b', str(rasp_sorted),
         '-s', '-wo'],
        capture_output=True, text=True, check=True,
    )
    if proc.stdout:
        return proc.stdout.splitlines()
    return []


def _path() -> str:
    import os
    return os.environ.get('PATH', '/usr/bin:/bin')


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_outputs(
    gini_rows: list[core.GiniRow],
    dropouts: list[core.DropoutRow],
    out_dir: Path,
) -> None:
    gini_path = out_dir / 'gini_long.tsv'
    with open(gini_path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['transcript_id', 'gene_id', 'species',
                    'compartment', 'region', 'gini', 'n_valid'])
        for r in gini_rows:
            w.writerow([r.transcript_id, r.gene_id, r.species,
                        r.compartment, r.region, f"{r.gini:.6f}", r.n_valid])
    log.info(f"wrote {len(gini_rows)} Gini rows -> {gini_path}")

    drop_path = out_dir / 'dropouts.tsv'
    with open(drop_path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['transcript_id', 'species', 'compartment',
                    'region', 'reason', 'n_valid'])
        for d in dropouts:
            w.writerow([d.transcript_id, d.species, d.compartment,
                        d.region, d.reason, d.n_valid])
    log.info(f"wrote {len(dropouts)} dropout rows -> {drop_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', type=Path, help='YAML config (see module docstring)')
    ap.add_argument('--work-dir', type=Path, default=Path('runs/work'),
                    help='cache dir for DBs, region BEDs, sorted RASP')
    ap.add_argument('--out-dir', type=Path, default=Path('runs/out'),
                    help='output dir for gini_long.tsv and dropouts.tsv')
    ap.add_argument('-v', '--verbose', action='store_true')
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
    )

    if not args.config:
        ap.error('--config is required')

    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    species_configs = load_config(args.config)

    all_gini: list[core.GiniRow] = []
    all_drop: list[core.DropoutRow] = []

    for sp in species_configs:
        region_beds, tx_to_gene, missing_drop = build_region_beds(sp, args.work_dir)
        all_drop.extend(missing_drop)

        for comp_name, rasp_path in sp.compartments.items():
            rasp_sorted = sort_rasp(rasp_path, args.work_dir)
            acc_total = core.ScoreAccumulator()

            # One intersect per region; merge into a per-compartment view
            # then filter+gini. (Kept per-region so a region's parse can't
            # leak into another; ScoreAccumulator keys on (tid, region).)
            for region, region_bed in region_beds.items():
                lines = run_intersect(region_bed, rasp_sorted)
                acc = core.parse_intersect_wo(lines)
                for key, vec in acc.vectors.items():
                    acc_total.vectors.setdefault(key, []).extend(vec)
                acc_total.out_of_range_dropped += acc.out_of_range_dropped

            if acc_total.out_of_range_dropped:
                log.warning(f"[{sp.name}/{comp_name}] dropped "
                            f"{acc_total.out_of_range_dropped} out-of-range "
                            f"score-bases")

            gini_rows, drops = core.gini_rows_from_accumulator(
                acc_total, sp.name, comp_name, tx_to_gene)
            all_gini.extend(gini_rows)
            all_drop.extend(drops)
            log.info(f"[{sp.name}/{comp_name}] {len(gini_rows)} gini rows, "
                     f"{len(drops)} dropouts")

    write_outputs(all_gini, all_drop, args.out_dir)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
