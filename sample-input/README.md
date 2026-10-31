# sample-input

bulk-de reads one input, `counts`: a count table (or one file per sample) and a sheet naming the
samples (see the README, *Input layouts*; the contract is in
`workflow-bulk-rnaseq/docs/conventions.md`).

## Synthetic fixtures (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{rejoin,sparc,genesets}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/rejoin with dev.env (~ donor + tissue)
```

- `rejoin/`: the shape of RE-JOIN's deposited knee/spine tables: `counts/knee_spine_counts.csv`
  (a title row, then `Ensembl ID`, `Gene name`, 12 sample columns) and
  `metadata/sample_info.xlsx` (`Sample ID`, `Tissue`, `OA grade`, `Donor`, `Age`). 6 donors x
  cartilage/meniscus; 20 genes 4x up and 12 down in meniscus, 10 genes 3x up with high OA grade.
- `sparc/`: Discover 427's layout: `primary/pool-<id>/<id>_Genes_ReadCount.txt` (htseq-style)
  with `samples.xlsx` (two animals per pool) and `subjects.xlsx` (sex, group, species).
- `genesets/`: tiny stand-ins for the image's MSigDB files (set `GENESETS_DIR` to it).

## Real data

Discover 427 (`workflow-bulk-rnaseq/scripts/fetch-sample-data.sh`): mount the dataset folder and
a sheet, e.g.

```sh
docker run --rm --platform linux/amd64 --memory 3g --network none --read-only --tmpfs /tmp \
  --env-file dev.env -e DESIGN="~ sex + condition" -e CONTRASTS="condition:T:N,sex:M:F" \
  -e SAMPLE_SHEET=design.csv \
  -v /path/to/data/tg-tumour:/data/input/dataset:ro \
  -v /path/to/design.csv:/data/input/design.csv:ro \
  -v "$PWD/data/output:/data/output" bulk-de:dev
```

`design.csv` has `library,sex,condition` rows such as `FN1,F,N` (derived from the SPARC sheets;
the male pools are `N1-N3`, `T1`, `T2`, `MT`). Without it, the SPARC sheets are used directly:
`-e DESIGN="~ subject_experimental_group"`.
