"""tests/test_gini_core.py — pytest suite for the bedtools-independent core.

Run: pytest -q   (from the gini_pipeline/ root)

Covers:
  * population Gini anchors from validation_history.md
  * GFF 1-based-inclusive -> BED 0-based half-open projection
  * intersect -wo parsing with overlap-bp expansion (the clip-and-expand)
  * [0,1] score validation (log-and-drop)
  * n_valid<=15 and zero-sum filters
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib import gini_core as core


# --- Gini anchors (from validation_history.md) -----------------------------

def test_gini_zero_one():
    assert math.isclose(core.population_gini([0, 1]), 0.5, rel_tol=1e-12)

def test_gini_zero_zero_one():
    assert math.isclose(core.population_gini([0, 0, 1]), 2/3, rel_tol=1e-12)

def test_gini_perfect_equality():
    assert math.isclose(core.population_gini([1, 1, 1]), 0.0, abs_tol=1e-12)

def test_gini_scale_invariant():
    # Gini is scale-invariant: multiplying all scores by k leaves it unchanged.
    a = core.population_gini([0.1, 0.2, 0.7])
    b = core.population_gini([1.0, 2.0, 7.0])
    assert math.isclose(a, b, rel_tol=1e-12)


# --- coordinate projection -------------------------------------------------

class _FakeFeat:
    """Minimal stand-in for a gffutils feature for projection tests."""
    def __init__(self, seqid, start, end, strand, fid):
        self.seqid = seqid; self.start = start; self.end = end
        self.strand = strand; self.id = fid

def test_gff_to_bed_projection():
    # GFF 1-based inclusive [100, 200] -> BED 0-based half-open [99, 200).
    f = _FakeFeat('chr1', 100, 200, '+', 'x')
    start, end = core._gff_to_bed_interval(f)
    assert (start, end) == (99, 200)
    assert end - start == 101   # 1-based inclusive width is end-start+1 = 101


def test_gini_weighted_matches_expanded():
    # The run-length weighted form must equal the expanded list form.
    import random
    rng = random.Random(0)
    for _ in range(50):
        vals = [round(rng.random(), 3) for _ in range(rng.randint(16, 40))]
        # build counts map
        counts = {}
        for v in vals:
            counts[v] = counts.get(v, 0) + 1
        a = core.population_gini(vals)
        b = core.population_gini_weighted(list(counts.items()))
        assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)

def test_weighted_anchors():
    assert math.isclose(core.population_gini_weighted([(0, 1), (1, 1)]), 0.5, rel_tol=1e-12)
    assert math.isclose(core.population_gini_weighted([(0, 2), (1, 1)]), 2/3, rel_tol=1e-12)
    assert math.isclose(core.population_gini_weighted([(1, 3)]), 0.0, abs_tol=1e-12)


# --- intersect parsing / expansion / clipping ------------------------------

def _wo_line(tid, region, chrom, rs, re_, rasp_s, rasp_e, score, overlap,
             strand='+'):
    """Build one `bedtools intersect -a region -b rasp -s -wo` line."""
    a = [chrom, str(rs), str(re_), f"{tid}|{region}", '.', strand]
    b = [chrom, str(rasp_s), str(rasp_e), 'NA', str(score), strand]
    return '\t'.join(a + b + [str(overlap)])

def test_expand_uses_overlap_not_full_interval():
    # RASP interval spans 4 bases but only overlaps the region by 2.
    # The -wo overlap column (2) must drive the count, not the RASP width.
    line = _wo_line('ENST1', 'CDS', 'chr1', 100, 102, 98, 102, 0.5, 2)
    acc = core.parse_intersect_wo([line])
    assert acc.counts[('ENST1', 'CDS')] == {0.5: 2}

def test_expand_multiple_rows_accumulate():
    lines = [
        _wo_line('ENST1', 'full', 'chr1', 0, 10, 0, 3, 0.2, 3),
        _wo_line('ENST1', 'full', 'chr1', 0, 10, 3, 5, 0.8, 2),
    ]
    acc = core.parse_intersect_wo(lines)
    assert acc.counts[('ENST1', 'full')] == {0.2: 3, 0.8: 2}

def test_out_of_range_scores_dropped():
    lines = [
        _wo_line('ENST1', 'CDS', 'chr1', 0, 10, 0, 1, 1.5, 1),   # >1, drop
        _wo_line('ENST1', 'CDS', 'chr1', 0, 10, 1, 2, -0.1, 1),  # <0, drop
        _wo_line('ENST1', 'CDS', 'chr1', 0, 10, 2, 3, 0.4, 1),   # ok
    ]
    acc = core.parse_intersect_wo(lines)
    assert acc.counts[('ENST1', 'CDS')] == {0.4: 1}
    assert acc.out_of_range_dropped == 2


# --- filters ---------------------------------------------------------------

def test_filter_n_valid_cutoff():
    acc = core.ScoreAccumulator()
    acc.add(('ENST_small', 'CDS'), 0.5, 15)   # exactly 15 -> dropped (need >15)
    # ENST_ok: 16 obs with variation so it's not zero-sum
    acc.add(('ENST_ok', 'CDS'), 0.1, 8)
    acc.add(('ENST_ok', 'CDS'), 0.9, 8)
    gini, drops = core.gini_rows_from_accumulator(
        acc, 'human', 'nucleoplasm', {'ENST_small': 'G1', 'ENST_ok': 'G2'})
    kept = {r.transcript_id for r in gini}
    dropped = {(d.transcript_id, d.reason) for d in drops}
    assert 'ENST_ok' in kept
    assert ('ENST_small', 'n_valid<=15') in dropped

def test_filter_zero_sum():
    acc = core.ScoreAccumulator()
    acc.add(('ENST_zero', 'full'), 0.0, 20)   # 20 obs, all zero -> zero_sum
    gini, drops = core.gini_rows_from_accumulator(
        acc, 'mouse', 'cytoplasm', {'ENST_zero': 'G'})
    assert gini == []
    assert drops[0].reason == 'zero_sum'
    assert drops[0].n_valid == 20
