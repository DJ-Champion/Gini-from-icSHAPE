# Data sources

Data files are not tracked in this repository (see `.gitignore`). This file records where each input came from, and how the Gini scores were checked.

## icSHAPE RNA structure probing data

**Source data**
- In vivo icSHAPE reactivity scores for nucleoplasm and cytoplasm fractions:
  - HEK293 (human, hg38)
  - mouse ES cells (mm10)
- From Sun *et al.* 2019, "RNA structure maps across mammalian cellular compartments", *Nat Struct Mol Biol* 26:322–330. doi:[10.1038/s41594-019-0200-7](https://doi.org/10.1038/s41594-019-0200-7)
- Obtained as genome-coordinate `.bed` files from the RASP database, 2021 release:
  - Li *et al.* 2021, "RASP: an atlas of transcriptome-wide RNA secondary structure probing data", *Nucleic Acids Res* 49:D183–D191. doi:[10.1093/nar/gkaa880](https://doi.org/10.1093/nar/gkaa880)
  - Site: rasp.zhanglab.net (now offline)
- Downloaded in 2024 (files dated 3 Sep 2024). Originals are archived in `probingFiles/RASP/`.

| Species | Compartment | File (in `config.yaml`) |
|---|---|---|
| human | nucleoplasm | `human/icSHAPE_hg38_nucl_vivo-score.bed` |
| human | cytoplasm | `human/icSHAPE_hg38_cy_invivo-score.bed` |
| mouse | nucleoplasm | `mouse/icSHAPE_mm10_nucl_vivo-score.bed` |
| mouse | cytoplasm | `mouse/icSHAPE_mm10_cy_vivo-score.bed` |

**Annotations**: one canonical transcript per gene.
- Human: MANE (`human/canonical_human.gff`).
- Mouse: GENCODE vM25, mm10 (`mouse/canonical_mouse_m25.gff`).
- Mouse must use an mm10 annotation to match the mm10 scores. The older `mouse/canonical_mouse.gff` is mm39 and must not be used.

## Gini scores

Computed by `run_gini.py` per transcript, for the full transcript, 5′UTR, CDS and 3′UTR:
- **Formula:** population Gini over all scored nucleotides. Missing (NULL) positions are excluded.
- **Minimum:** more than 15 scored positions per region.
- **Output:** `gini_long.tsv` (current run: `gini_long_memfix.tsv`).

## Validation (October 2026)

- **Exact reproduction:** an independent re-implementation, written separately in numpy, reproduced all 58,838 Gini values exactly.
- **Same formula as the authors:** the formula is algebraically identical to the one in Sun *et al.*'s lab code (IPyRSSA / icSHAPE-pipe `calc_gini`); maximum difference on test vectors was 9 × 10⁻¹⁶.
- **Correct coordinates:** a per-nucleotide alignment against Sun *et al.*'s own transcript-coordinate scores (~48,000 nucleotides) confirmed the genome-to-transcript mapping. There are no off-by-one or strand errors.
- **Reproduces the paper:** applied to Sun *et al.*'s GEO-deposited scores (GSE117840) with the same HEK293 half-lives they used (Schueler *et al.* 2014), the method reproduces their published structure–half-life correlation: r = −0.233 vs −0.23 (nucleoplasm, Fig 2c).

## Note on data versions

Two processed versions of Sun *et al.*'s icSHAPE scores are publicly available:
1. **RASP 2021 `.bed` files (used here).** A separate processing of the same experiment from the scores Sun *et al.* deposited in GEO.
   - Highly correlated with GEO (r ≈ 0.7–0.9 per transcript), but not identical.
   - Scores more positions.
2. **GEO-deposited scores (GSE117840).** Also available in genome coordinates as RASP v2 "raw" files ([rasp2.zhanglab.net](https://rasp2.zhanglab.net/); *Nucleic Acids Res* 2025, 53:D211).

Results were checked on both versions (RASP v2 raw config: `RASP/config_rasp_v2.yaml`), and the conclusions agree.

RASP v2 "imputed" files, where missing values are filled in by a deep-learning model, were **not** used. They should not be used for Gini, because filling in missing positions changes the score distribution.

## Half-life and expression data (used downstream)

- **Primary half-life measure:** PC1 consensus half-life from Agarwal & Kelley 2022, "The genetic and biochemical determinants of mRNA degradation rates in mammals", *Genome Biol* 23:245. doi:[10.1186/s13059-022-02811-x](https://doi.org/10.1186/s13059-022-02811-x) (Additional file 3).
- **Validation only:** HEK293 half-lives from Schueler *et al.* 2014, *Genome Biol* 15:R15, taken from the "Dieterich_4sU_HEK293" columns of Agarwal & Kelley Additional file 2.
- **HEK293 expression:** Human Protein Atlas v25.1, cell line RNA data (`rna_celline.tsv`), filtered to HEK293 (`hpa_hek293/`).
