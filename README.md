# bulk-de

Which genes are expressed differently between groups of bulk RNA-seq samples (tumour vs naïve,
osteoarthritic vs healthy cartilage, male vs female), and which pathways do the changes point
to? This step of the bulk RNA-seq workflows fits DESeq2's model to gene counts for the design
you give, tests each comparison, and ranks genes for gene set enrichment.

- **For a biologist:** give it a count table and a sheet saying which group each sample is in.
  It writes, per comparison, a table of every gene with its fold change and adjusted p-value,
  a volcano plot and an MA plot, the gene sets (MSigDB Hallmark, GO biological process) whose
  genes move together, and a PCA showing how the samples group. A block term (donor, sex,
  batch, study) can be added to the design so its differences are not mistaken for the effect.
- **For a computational analyst:** PyDESeq2 (DESeq2's negative binomial GLM in Python) on an
  additive design of categorical sample-sheet columns, Wald tests for `column:level:reference`
  contrasts, Benjamini-Hochberg with independent filtering, apeglm-type LFC shrinkage, VST +
  PCA (DESeq2's `plotPCA` recipe), and preranked GSEA on the Wald statistic with GSEApy, gene
  sets from MSigDB 2026.2 baked into the image. It reads rnaseq-quant's output or count tables
  as deposited in data repositories (the common case for RE-JOIN knee/spine tissue datasets).

```
A. rnaseq-quant (+ index from rnaseq-index) -> bulk-de -> bulk-report
B.                         deposited counts -> bulk-de -> bulk-report
```

Light step: runs as a container or a Lambda function (same image).

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `counts` | input (required) | count table(s) + a sample sheet | See *Input layouts*. Without links everything is found anywhere under `INPUT_DIR` (symlinks followed). |
| `de` | output | `<LABEL>/...` (below) | Everything under a folder named after `LABEL`, so several comparisons can be linked into one report. |

### Input layouts

The first layout found is used (in this order); the others are ignored and listed in `de.json`.

| Layout | Recognised by | Notes |
|---|---|---|
| rnaseq-quant | `gene_counts.csv` (`gene_id, gene_name`, one column per sample) | length-scaled counts, as conventions.md describes |
| featureCounts | header `Geneid Chr Start End Strand Length ...` | sample = BAM file name without `Aligned.sortedByCoord.out.bam` etc. |
| genes x samples table | any other CSV/TSV/XLSX with >= 100 rows and >= 2 numeric columns | first column = gene id (Ensembl, Entrez or symbol), a text column = gene names, numeric columns = samples; a title row above the header is skipped (as in many deposited spreadsheets); `Length`/`Chr`/`Start`/... columns are not samples; a sample column named twice keeps its first copy (noted); duplicated gene ids are summed (noted) |
| STAR | `*ReadsPerGene.out.tab`, one per sample | column 3 or 4 when it holds >= 3x the other (stranded), else column 2 (unstranded); noted |
| RSEM | `*.genes.results`, one per sample | `expected_count`, rounded |
| two-column | one file per sample, gene and count (htseq-count; Discover 427's `*_Genes_ReadCount.txt`) | `__no_feature` etc. dropped |

Per-sample files are named after the nearest `sam-`/`pool-` folder above them (SPARC layout),
else after the file name without its tool suffix. Non-integer counts are rounded (as tximport
does for estimated counts), and the rounding is noted. Only genes present in every sample's
file are kept.

**Sample metadata.** Any CSV/TSV/XLSX that is not the count table can describe the samples.
Column names become snake_case (`OA grade` -> `oa_grade`). The sheet and column whose values name
the most count samples is used (compared without case and without a `sam-`/`sub-`/`pool-`
prefix); `SAMPLE_SHEET` picks a file by name. When that sheet has `subject_id` and another sheet
is keyed by `subject_id` (SPARC `subjects.xlsx`), the subject columns are joined. When several
rows match one library (SPARC `samples.xlsx` lists the animals of a pooled library), a column is
kept when they agree and left empty otherwise. `condition` falls back to SPARC's
`subject_experimental_group`, then `sample_experimental_group`. **Only the DESIGN columns are
written to the outputs**, so clinical fields in a sheet never spread.

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `LABEL` | string | *(empty)* | Output folder (`[a-z0-9-]+`). Empty: the last DESIGN term with `_` as `-`. |
| `DESIGN` | string | `~ condition` | Additive formula of sheet columns, variable of interest last: `~ condition`, `~ sex + condition`, `~ donor + tissue`, `~ study + group`. Every term is a factor (numbers are levels). No interactions; for a sex-specific effect, run once per sex with a sheet listing only those samples. |
| `CONTRASTS` | string | *(empty)* | `column:level:reference`, comma-separated, e.g. `condition:T:N,sex:M:F`. log2 fold change > 0 = higher in `level`. Empty: every level of the last term vs its reference. |
| `REFERENCE_LEVELS` | string | *(empty)* | `column:level,...`: each factor's base level. Unset: the reference of its first contrast, else the first level alphabetically. |
| `SAMPLE_SHEET` | string | *(empty)* | File name of the sheet to use; empty = automatic (above). |
| `MIN_COUNT` | integer | `10` | A gene is tested with >= `MIN_COUNT` counts in >= `MIN_SAMPLES` samples. |
| `MIN_SAMPLES` | integer | `0` | `0` = the size of the smallest group of the last term (edgeR `filterByExpr`'s idea). |
| `FDR` | number | `0.05` | BH-adjusted p below which genes and gene sets count as significant (also DESeq2's independent-filtering alpha). |
| `LFC_SHRINK` | boolean | `true` | Add `log2_fold_change_shrunk` (apeglm-type prior) for contrasts against the factor's base level. Tests are unchanged. |
| `GENE_SETS` | string | `hallmark,go-bp` | `hallmark`, `go-bp`, both, or `none`. |
| `ORGANISM` | string | `auto` | `human`, `mouse`, `rat`, or `auto`: the sheet's `species`/`organism` column, else the Ensembl prefix (`ENSG`, `ENSMUSG`, `ENSRNOG`), else symbol case (upper = human, else mouse; **symbols cannot tell rat from mouse: set `rat`**). |
| `SEED` | integer | `42` | GSEA permutations. |
| `THREADS` | integer | `0` | PyDESeq2, GSEA, BLAS threads; `0` = the cgroup CPU quota (2 on Lambda). |

Samples without a value for every DESIGN term are left out (logged, listed in `design.json`).
Exit codes: `0` success; `1` no count table, fewer than 3 usable samples, < 10 genes pass the
filter, or the DESIGN columns are missing because no sheet matched; `2` bad parameters (formula,
contrast, reference level, a DESIGN column not in the sheet, a term with one level, a design not
of full rank, e.g. `~ sex + group` when the group already holds the sex).

## Outputs (`OUTPUT_DIR/<LABEL>/`)

`<contrast>` is `<column>-<level>-vs-<reference>` (characters outside `[A-Za-z0-9._-]` as `_`).

| File | Content |
|---|---|
| `<contrast>.csv` | One row per tested gene, by p-value: `gene_id`, `gene_name`, `base_mean` (mean normalised count), `log2_fold_change` (MLE), `lfc_se`, `stat` (Wald), `pvalue`, `padj` (BH after independent filtering; empty when filtered or a Cook's outlier), `log2_fold_change_shrunk` (when available) |
| `<contrast>_gsea.csv` | Preranked GSEA by `stat`: `collection`, `term`, `size` (genes of the set in the ranking), `es`, `nes`, `pvalue`, `fdr`, `n_leading_edge`, `leading_edge` (`;`-separated MSigDB symbols) |
| `<contrast>_volcano.png`, `<contrast>_ma.png` | log2 FC vs -log10 p (top 12 named); MA with the shrunk LFC when present |
| `pca.csv`, `pca.png` | Design columns and PC1..PCk of the VST (blind) of the 500 most variable genes; variance explained in `summary.json` |
| `samples.csv` | `sample_id`, design columns, `library_size`, `n_genes_detected`, `size_factor`, `used` |
| `design.json` | Formula, terms, contrasts, reference and all levels, samples used (id + design values), samples left out with the reason, gene filter (`input`, `tested`, `min_count`, `min_samples`) |
| `de.json` | What was found: layout, count files, genes, samples, sheet / matching column / subjects sheet, organism and how it was decided, gene-set release and genes mapped, skipped files, notes |
| `summary.json` | Per contrast: samples per side, genes tested, n up/down at FDR, top 10 up/down, per collection sets tested / significant and the top 10; PCA variance; notes |
| `provenance.json` | App, version, tools, parameters, resolved design (tokens), gene-set release with checksums, inputs, timings, peak memory, per-sample status, `citations` (`{doi, pmid, text}`) |

## Methods

Following Love et al. 2015 (RNA-seq workflow) and Conesa et al. 2016 (best practices):

1. **Counts** are combined into one integer genes x samples matrix (above). No normalisation is
   applied to them: DESeq2 models raw counts.
2. **Filter:** genes with >= `MIN_COUNT` counts in >= `MIN_SAMPLES` samples (default the smallest
   group). This mainly saves time and stabilises the dispersion trend; DESeq2's independent
   filtering then optimises the number of rejections at `FDR`.
3. **Model:** PyDESeq2 (Muzellec et al. 2023), a reimplementation of DESeq2 (Love et al. 2014):
   median-of-ratios size factors, gene-wise dispersions shrunk to a parametric trend, negative
   binomial GLM on the design, Cook's-distance outlier refitting when a group has >= 7 samples,
   Wald test per contrast, BH (Benjamini & Hochberg 1995). Columns and levels are encoded as
   tokens before reaching the formula parser, so any level name works.
4. **Shrinkage:** `lfc_shrink` with an apeglm-type heavy-tailed prior (Zhu et al. 2019), for
   ranking and plotting by effect size; possible for coefficients of the model, i.e. contrasts
   against the factor's base level (others get none, logged).
5. **PCA:** variance-stabilising transformation (Anders & Huber 2010), blind to the design,
   of the 500 most variable genes (DESeq2's `plotPCA`).
6. **GSEA:** genes ranked by the Wald statistic (sign = direction, magnitude = evidence; the
   standard preranked input), mapped to MSigDB symbols with MSigDB's own chip files (Ensembl
   id, else symbol/alias; rat -> mouse orthologs); duplicates keep the largest |stat|. GSEA
   (Subramanian et al. 2005) with GSEApy's Rust core (Fang et al. 2023), 1,000 gene-set
   permutations, sets of 15-500 genes, weight 1, FDR per collection. Collections: MSigDB Hallmark
   (Liberzon et al. 2015; mouse MH, Castanza et al. 2023) and GO biological process (Ashburner
   et al. 2000; GO Consortium 2023) as MSigDB C5/M5 GO:BP.

**Deviations from the reference workflow and why.** No interaction terms or continuous
covariates (keeps one clear path; covers the designs of the first users). The PCA uses
PyDESeq2's VST rather than `rlog` (equivalent for > 30 samples, faster). GSEApy is called
through its Rust core (`gseapy.gse.prerank_rs`) rather than `gseapy.prerank`, which keeps every
gene set's running-sum vector as Python lists (3.4 GB for GO:BP on 13 k genes; 1.2 GB this
way): same algorithm, same results up to permutation noise (|ΔNES| <= 0.07 on Discover 427).
Gene-set permutation (preranked GSEA) is the only option with few samples; it ignores
gene-gene correlation, so its p-values are optimistic for sets of co-regulated genes.

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python      # pytest on synthetic fixtures (no Docker)
make run PYTHON=.venv/bin/python       # image (amd64) over the RE-JOIN-shaped fixture -> data/output/
make run INPUT=/path/to/counts         # or over a folder of your own (set DESIGN etc. in dev.env)
```

Locally without Docker, build the gene sets once (`python genesets/build.py data/genesets`) and
set `GENESETS_DIR=data/genesets`. The tests cover: the RE-JOIN shape (title row above the header,
versioned Ensembl ids, separate `sample_info.xlsx` with `Sample ID`, `Tissue`, `OA grade`,
`Donor`, `Age`) with `~ donor + tissue` (all 32 planted genes found, the planted gene set top
of GSEA, PC1 = tissue, `Age` never in the outputs); a second factor and a reversed reference;
the SPARC pool layout (htseq files in `pool-*` folders, `samples.xlsx` with two animals per
pool joined to `subjects.xlsx`); the `condition` fallback; a user sheet with a missing value;
featureCounts, STAR (strandedness), RSEM, rnaseq-quant precedence; negative counts; a symlink
leaving the input; bad parameters; and the Lambda handler.

## Runtime

`computeTypes: [standard, lambda]`, 2 vCPU / 3008 MB, 15 min timeout. Writes only to
`OUTPUT_DIR` and `/tmp` (checked with `docker run --read-only --tmpfs /tmp --network none`).
Measured on an M-series Mac shared with other builds (load average ~10), so wall times are
upper bounds:

`docker run --platform linux/amd64 --memory 3g --cpus 2 --network none --read-only` (amd64
emulation), peak RSS from `getrusage`, cgroup `memory.peak` (includes page cache):

| Data | Design | Wall time (PyDESeq2 / GSEA) | Peak RSS | cgroup peak |
|---|---|---|---|---|
| Discover 427, authors' counts: 12 libraries x 22,609 genes (13,338 tested), 2 contrasts, Hallmark + GO:BP | `~ sex + condition` | 202 s (54 s / 143 s) | 2,369 MB | 2,803 MB |
| Discover 457, deposited XLSX: 48 samples x 23,451 genes (16,631 tested), 1 contrast, Hallmark + GO:BP | `~ study + group` | 158 s (60 s / 84 s) | 2,248 MB | 2,441 MB |
| Discover 427, same without GSEA (native arm64, SPARC sheets only) | `~ subject_experimental_group` | 18 s | 346 MB | |
| Synthetic RE-JOIN-shaped fixture: 12 x 400 | `~ donor + tissue` | 3 s | 270 MB | |

Memory is dominated by the GO:BP enrichment (all ~3,700 sets of 15-500 genes held by GSEApy's
core) and, for large XLSX tables, by openpyxl. Time grows with samples x genes for PyDESeq2 and
with sets x permutations for GSEA. 3008 MB fits about 50 samples x 25 k genes with both
collections; for hundreds of samples, run on `standard` with more memory.

Image: 463 MB unpacked, 154 MB compressed (no R, no torch; gene sets 9 MB).

## Results on the demo data

**Discover 427** (RE-JOIN; retro-labelled tongue-innervating trigeminal neurons, FACS-sorted,
pooled 2 mice per library; 12 single-end libraries, sex F/M x naïve N / tongue tumour T, 3 pools
each). The authors' per-library `*_Genes_ReadCount.txt` (htseq-style, gene symbols) were read
directly. The SPARC sheets list the mice (`sam-M1-N1` ...) with a `pool id` and their subjects'
`sex` and `subject experimental group` (`Female-Control` ... `Male-Tumor`); **the male pools are
named `N1-N3`, `T1`, `T2`, `MT`, not `MN`/`MT1-3`**, so libraries cannot be classified from the
folder letters alone. The SPARC sheets agree with the pool names (every pool's two mice share sex
and group). The demo sheet `data/work/tg-tumour/de-input/design.csv` (library, sex F/M,
condition N/T) was derived from them; with the SPARC sheets alone the app matches the pools
itself (`pool_id`) and `~ subject_experimental_group` works too. Outputs:
`data/work/tg-tumour/de/de/condition/`.

- `~ sex + condition`, `condition:T:N`: **0 genes at FDR 0.05** (19 at p < 0.01), but GSEA finds
  immune and inflammatory sets up in tumour-bearing mice: Hallmark TNF-α signalling via NF-κB
  (NES 1.94) and allograft rejection (1.77); 23 GO:BP sets, led by MHC class II antigen
  presentation (NES 2.45).
- `sex:M:F`: **30 up and 652 down in males**. Besides the expected Y genes (Ddx3y, Eif2s3y, Uty,
  Kdm5d) and Xist/Tsix, the female-higher genes are immune-cell markers: C1qa (log2 FC -2.8,
  padj 5e-4), C1qc, C1qb, Lyz2, Cd68, Fcer1g, Tyrobp, S100a8 (Csf1r -2.0, padj 0.15).
  **This matches the paper** (Ibrahim et al. 2023): "Naïve female neurons innervating the tongue
  expressed immune cell markers such as Csf1R, C1qa and others, that weren't expressed in
  males". The paper also reports that "very few differentially expressed genes (DEGs)
  overlapped between the sexes" upon tumour growth; an additive model, which estimates one
  tumour effect shared by both sexes, is then expected to find little, as here. For the
  sex-specific tumour response, run once per sex (a sheet listing only those libraries).
- PCA: PC1 holds 68 % of the variance of the 500 most variable genes and is not explained by
  sex or condition (R² 0.38 for sex): pool-to-pool variation dominates, as expected for sorted,
  pooled neurons.

**From rnaseq-quant** (salmon on reads subsampled to 1 M per library, GENCODE M39;
`data/work/tg-tumour/quant/counts/`, output `data/work/tg-tumour/de-quant/de/condition/`): 9,231
genes pass the filter (vs 13,338 at full depth); `sex:M:F` 28 up / 68 down, 45 of them among the
682 found from the authors' counts; Wald statistics correlate with the authors'-count run at
Spearman 0.78 (sex) and 0.49 (condition) over 7,752 shared genes. `condition:T:N`: 1 gene. The
loss is depth (1 M vs ~20 M reads), not method.

**Discover 457** (Newton MD et al., osteoarthritis synovium meta-analysis, doi:10.26275/robc-8zol,
CC BY 4.0), the shape RE-JOIN's knee data take: `primary/RawCounts-Bulk-RNAseq.xlsx` is a
deposited genes x samples workbook with a study banner row above the header, Entrez ids + symbols,
`Description/Chr/Start/End/Strand/Length` columns, and a sample column named twice
(`sub-GSE89408-healthy8`; the second copy is ignored, noted). Sample columns are SPARC subject
ids, matched to `samples.xlsx` (`subject_id`) and joined to `subjects.xlsx`. `~
geo_study_accession_id + sample_experimental_group` (2 GEO studies, 48 samples): **2,977 genes
up and 1,980 down in OA** vs non-OA synovium (top up: CSN1S1, S100A8, SCRG1, S100A9; down:
SLC7A10, AQP7B, G0S2), 28 of 50 Hallmark sets (oxidative phosphorylation, protein secretion, MYC
targets, allograft rejection up) and 411 GO:BP sets at FDR 0.05. Output:
`data/work/oa-synovium-457/de/de/oa/`.

## Limitations

- **Few replicates, little power.** With 3 pools per group, only large, consistent changes reach
  significance; a gene that is not significant is not shown to be unchanged.
- **Additive designs only.** `~ sex + condition` estimates one tumour effect shared by both
  sexes; when the sexes respond differently (as in Discover 427), run each sex separately.
- **Confounding.** A block helps only when it is not confounded with the effect (the app
  refuses designs that are not of full rank, but cannot detect a batch nobody recorded).
- **Pooled libraries** average out differences between animals; replicates are pools.
- **Gene identifiers.** GSEA maps genes by Ensembl id or symbol; Entrez-only tables without a
  symbol column, or genes absent from MSigDB's chip files, are not ranked (the share mapped is
  in `de.json`). Rat genes use mouse orthologs and mouse gene sets.
- **GSEA** with gene-set permutation is optimistic for sets of co-regulated genes; GO sets
  overlap heavily, so many significant terms can reflect the same genes.
- Counts must be raw (not TPM/FPKM/normalised); non-integer values are rounded, which is right
  for estimated counts (salmon, RSEM) but wrong for TPM. Check `de.json` notes.

## Licences

Our code: Apache-2.0 (`LICENSE`). Everything in the image is permissive (PyDESeq2 MIT, GSEApy
BSD-3-Clause, NumPy/SciPy/pandas/scikit-learn BSD-3-Clause, openpyxl MIT, ...); nothing is GPL.
The gene sets are MSigDB 2026.2 Hallmark and GO:BP, **CC BY 4.0** (MSigDB licence terms, checked
2026-10-31; KEGG and BioCarta sets, which have other terms, are not included). DESeq2
(LGPL-3) and apeglm (GPL-2) are methods reimplemented by PyDESeq2, not bundled. See
`THIRD_PARTY.md`.

## Credits

- **PyDESeq2** (Boris Muzellec, Maria Teleńczuk, Vincent Cabeli, Mathieu Andreux; Owkin): the
  DESeq2 model, Wald tests, shrinkage and VST used here.
- **DESeq2** (Michael I. Love, Wolfgang Huber, Simon Anders): the method PyDESeq2 implements;
  **apeglm** (Anqi Zhu, Joseph G. Ibrahim, Michael I. Love): the shrinkage prior; the
  **variance-stabilising transformation** (Simon Anders, Wolfgang Huber).
- **GSEA** (Aravind Subramanian, Pablo Tamayo and colleagues at the Broad Institute): the
  enrichment method; **GSEApy** (Zhuoqing Fang, Xinyuan Liu, Gary Peltz): its implementation.
- **MSigDB** (Broad Institute and UC San Diego; Arthur Liberzon et al. for the Hallmark
  collection, Anna Castanza et al. for the mouse collections): the gene sets and identifier
  mappings; **Gene Ontology Consortium**: the GO annotations behind the GO:BP sets.
- **Benjamini & Hochberg**: the false discovery rate.
- **Demo data:** Discover dataset 427 by Alexei V. Tumanov, Sergey A. Shein, Zhao Lai, Yidong
  Chen and Shivani B. Ruparel (RE-JOIN; CC BY 4.0), described by Ibrahim et al. 2023; the
  authors' per-library gene counts were used. Discover dataset 457 by Michael D. Newton, Hannah
  Swahn, Dana Orange, Joseph B. Lesnak, Theodore J. Price, Anne-Marie Malfait, Rachel E. Miller,
  Martin K. Lotz and Tristan Maerz (CC BY 4.0) was used as a real-world deposited count workbook.

## Citations

Checked on 2026-10-31 against Crossref (title, authors and their number and order, journal,
year, volume, issue, pages) and PubMed (esearch on the DOI, then esummary: title, first author,
author count, year, pages); datasets against DataCite and the Discover API. Full author lists
are in `CITATION.cff`; `scripts/verify_citations.py` rewrites `processor/references.json`.

- Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq
  differential expression analysis. *Bioinformatics* 39(9):btad547 (2023).
  doi:10.1093/bioinformatics/btad547. PMID 37669147.
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq
  data with DESeq2. *Genome Biol* 15(12):550 (2014). doi:10.1186/s13059-014-0550-8. PMID 25516281.
- Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing
  the noise and preserving large differences. *Bioinformatics* 35(12):2084-2092 (2019).
  doi:10.1093/bioinformatics/bty895. PMID 30395178.
- Anders S, Huber W. Differential expression analysis for sequence count data. *Genome Biol*
  11(10):R106 (2010). doi:10.1186/gb-2010-11-10-r106. PMID 20979621.
- Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful
  Approach to Multiple Testing. *J R Stat Soc Series B* 57(1):289-300 (1995).
  doi:10.1111/j.2517-6161.1995.tb02031.x. (Crossref; not in PubMed)
- Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based
  approach for interpreting genome-wide expression profiles. *Proc Natl Acad Sci U S A*
  102(43):15545-15550 (2005). doi:10.1073/pnas.0506580102. PMID 16199517.
- Fang Z, Liu X, Peltz G. GSEApy: a comprehensive package for performing gene set enrichment
  analysis in Python. *Bioinformatics* 39(1):btac757 (2023). doi:10.1093/bioinformatics/btac757.
  PMID 36426870.
- Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular
  Signatures Database Hallmark Gene Set Collection. *Cell Syst* 1(6):417-425 (2015).
  doi:10.1016/j.cels.2015.12.004. PMID 26771021.
- Castanza AS, Recla JM, Eby D, Thorvaldsdóttir H, Bult CJ, Mesirov JP. Extending support for
  mouse data in the Molecular Signatures Database (MSigDB). *Nat Methods* 20(11):1619-1620
  (2023). doi:10.1038/s41592-023-02014-7. PMID 37704782.
- Ashburner M, Ball CA, Blake JA, et al. Gene Ontology: tool for the unification of biology.
  *Nat Genet* 25(1):25-29 (2000). doi:10.1038/75556. PMID 10802651.
- The Gene Ontology Consortium, Aleksander SA, Balhoff J, et al. The Gene Ontology knowledgebase
  in 2023. *Genetics* 224(1):iyad031 (2023). doi:10.1093/genetics/iyad031. PMID 36866529.
- Love MI, Anders S, Kim V, Huber W. RNA-Seq workflow: gene-level exploratory analysis and
  differential expression. *F1000Res* 4:1070 (2015). doi:10.12688/f1000research.7035.1.
  PMID 26674615. (practice followed; not run)
- Conesa A, Madrigal P, Tarazona S, et al. A survey of best practices for RNA-seq data analysis.
  *Genome Biol* 17(1):13 (2016). doi:10.1186/s13059-016-0881-8. PMID 26813401. (practice followed)

Demo data: Tumanov AV, Shein SA, Lai Z, Chen Y, Ruparel SB. Sex-dependent differences in the
genomic profile of lingual sensory neurons in naive and tongue-tumor bearing mice. SPARC Portal
(2025). doi:10.26275/gaap-zlpj. CC BY 4.0. Described by Ibrahim T, Wu P, Wang LJ, et al.
*Sci Rep* 13(1):13117 (2023). doi:10.1038/s41598-023-40380-6. PMID 37573456.

Test data: Newton MD, Swahn H, Orange D, Lesnak JB, Price TJ, Malfait AM, Miller RE, Lotz MK,
Maerz T. Cross-platform transcriptomic data integration of human osteoarthritis synovium. SPARC
Portal (2025). doi:10.26275/robc-8zol. CC BY 4.0. (DataCite and the Discover API)
