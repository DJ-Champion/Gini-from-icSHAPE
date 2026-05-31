# Gini-from-icSHAPE pipeline

Computes per-transcript population Gini coefficients of icSHAPE reactivity
for human (HEK293) and mouse (mESC), in two compartments (nucleoplasm,
cytoplasm), at four region granularities (full, 5UTR, CDS, 3UTR). Output is
a long-format TSV ready for the half-life join in R downstream.

## Layout

```
gini_pipeline/
  run_gini.py            driver (orchestration only)
  config.yaml            paths to the 2 GFFs + 4 RASP files
  lib/
    gff.py               DJ's GFF/gffutils helper (unmodified)
    gini_core.py         Gini, region extraction, intersect parsing
  tests/
    test_gini_core.py    pytest suite (no bedtools needed)
```

## Run

```bash
conda activate yaml_env          # gffutils + bedtools on PATH
python run_gini.py --config config.yaml \
    --work-dir runs/work --out-dir runs/out
```

Outputs:
- `runs/out/gini_long.tsv` — `transcript_id gene_id species compartment region gini n_valid`
- `runs/out/dropouts.tsv` — `transcript_id species compartment region reason n_valid`

Dropout reasons: `no_feature` (region not annotated for that transcript,
e.g. non-coding has no CDS/UTR), `n_valid<=15` (too few valid positions),
`zero_sum` (all-zero vector, Gini undefined).

## Tests

```bash
python -m pytest tests/ -q
```

Includes the Gini anchors from `validation_history.md`
(`[0,1]`→0.5, `[0,0,1]`→2/3, `[1,1,1]`→0) and the clip-and-expand check.

## Design decisions (and why)

- **Regions read verbatim from GFF features** (`exon`, `CDS`,
  `five_prime_UTR`, `three_prime_UTR`), not derived by slicing exons against
  CDS bounds. The reduced GFFs already carry these features, so the
  stop-codon-in-CDS convention comes for free. Verified on Cdh2 and TSPAN6:
  `full` (from exon union) equals `5UTR+CDS+3UTR` to the base.
- **`full` built independently from exon features**, not as the union of the
  sub-regions, so it is correct for non-coding transcripts (full only) and
  robust to annotation quirks.
- **`bedtools intersect -s -wo`**, not `-wa -wb`. The `-wo` overlap-bp column
  is the clipped overlap width, so expanding to that many per-base
  observations clips and expands in one step — a RASP row straddling a
  region boundary contributes only its in-region bases. `-s` is required
  (antisense reactivity bleed otherwise).
- **subprocess to bedtools**, not pybedtools — keeps the dependency surface
  to gffutils only.
- **Population Gini, no n-1 correction** — matches the source paper.
- **RASP files sorted-and-cached once** (`sort -k1,1 -k2,2n`, `LC_ALL=C`).
- **Scores validated to [0,1]**; out-of-range bases dropped and counted in
  a per-compartment warning.
- **Transcript universe = exon parents in the reduced GFF.** No separate
  manifest; the GFF is already one-transcript-per-gene from the gene list.

## Not in scope

- Half-life join and correlation (done in R downstream).
- liftOver of mouse coordinates — deliberately avoided; the entire mouse
  computation stays in mm10 with vM25, keyed by ENSMUSG for the downstream
  join.
