# Third-party software and data

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following third-party
components, installed from PyPI at the pinned versions in `requirements.txt`, and the MSigDB gene
sets below; their licences apply to the image. Nothing in the image is GPL.

## Software

| Component | What | Licence | Source |
|---|---|---|---|
| PyDESeq2 0.5.4 | DESeq2's size factors, dispersion and fold-change estimation, Wald tests, Cook's and independent filtering, LFC shrinkage, VST, in Python. | MIT | https://github.com/owkin/PyDESeq2 |
| GSEApy 1.3.2 | Preranked GSEA (Rust core). | BSD-3-Clause | https://github.com/zqfang/GSEApy |
| formulaic 1.2.2, formulaic-contrasts 1.0.0 | Design matrices for PyDESeq2. | MIT / BSD-3-Clause | https://github.com/matthewwardrop/formulaic |
| anndata 0.13.4 (+ h5py, zarr, numcodecs) | PyDESeq2's data container. | BSD-3-Clause | https://github.com/scverse/anndata |
| NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.0, scikit-learn 1.9.1 | Arrays, statistics, tables. | BSD-3-Clause | numpy.org, scipy.org, pandas.pydata.org, scikit-learn.org |
| Matplotlib 3.11.2 | PNG figures. | Matplotlib licence (PSF-based) | https://matplotlib.org |
| openpyxl 3.1.5 (+ et_xmlfile 2.0.0) | Reading XLSX count tables and sample sheets. | MIT | https://foss.heptapod.net/openpyxl/openpyxl |
| defusedxml 0.7.1 | Safe XML parsing for openpyxl (no entity expansion). | PSF-2.0 | https://github.com/tiran/defusedxml |
| requests (+ certifi, charset-normalizer, idna, urllib3) | Dependency of GSEApy; never used at run time (the container has no network). | Apache-2.0 / MPL-2.0 (certifi) / MIT | PyPI |
| annotated-types, array-api-compat, cloudpickle, contourpy, cycler, donfig, fonttools, google-crc32c, interface_meta, joblib, kiwisolver, legacy-api-wrap, msgspec, narwhals, natsort, packaging, pillow, pydantic (+ pydantic-core, pydantic-settings, typing-inspection), pyparsing, python-dateutil, python-dotenv, PyYAML, scverse-misc, session-info, session-info2, six, stdlib-list, threadpoolctl, typing_extensions, tzdata, wrapt | Dependencies of the above (versions in `requirements.txt`). | MIT / BSD-2/3-Clause / MIT-CMU (Pillow) / Apache-2.0 (google-crc32c, tzdata) / Apache-2.0 or BSD-2-Clause (packaging, python-dateutil) / MPL-2.0 (legacy-api-wrap, session-info2; unmodified) / PSF-2.0 (typing_extensions) | PyPI |
| awslambdaric 4.2.3 (+ simplejson 4.2.0, snapshot-restore-py 1.0.0) | AWS Lambda runtime interface client. | Apache-2.0 (simplejson: MIT / AFL-2.1) | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| python:3.12.16-slim-bookworm | Base image (Debian 12 + CPython 3.12). | PSF (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

## Data baked into the image (`/opt/genesets`, built by `genesets/build.py`)

All from MSigDB release 2026.2 (October 2026), downloaded at build time and checked against the
sha256 pinned in `genesets/build.py`. **Licence: CC BY 4.0** — the MSigDB licence terms
(https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp, read 2026-10-31) put MSigDB v2022.1
and later under the Creative Commons Attribution 4.0 International License, except KEGG
(CC BY-SA 4.0 / KEGG terms) and BioCarta sets, which are **not** included here. Attribution:
Liberzon et al. 2015 (Hallmark), Castanza et al. 2023 (mouse MSigDB), Subramanian et al. 2005.

| File | MSigDB source | Contents |
|---|---|---|
| `hallmark.human.gmt.gz` | `h.all.v2026.2.Hs.symbols.gmt` | 50 Hallmark sets, human symbols |
| `hallmark.mouse.gmt.gz` | `mh.all.v2026.2.Mm.symbols.gmt` | 50 mouse Hallmark sets |
| `go-bp.human.gmt.gz` | `c5.go.bp.v2026.2.Hs.symbols.gmt` | GO biological process sets. Derived from the Gene Ontology (CC BY 4.0, Gene Ontology Consortium; GO-basic released 2026-07-26 with NCBI gene2go annotations of 2026-09-30, per the MSigDB 2026.2 release notes) |
| `go-bp.mouse.gmt.gz` | `m5.go.bp.v2026.2.Mm.symbols.gmt` | the same for mouse |
| `map.human.tsv.gz` | `Human_Ensembl_Gene_ID_MSigDB.v2026.2.Hs.chip`, `Human_Gene_Symbol_with_Remapping_MSigDB.v2026.2.Hs.chip` | Ensembl id / symbol / alias -> MSigDB human symbol |
| `map.mouse.tsv.gz` | `Mouse_Ensembl_Gene_ID_MSigDB.v2026.2.Mm.chip`, `Mouse_Gene_Symbol_Remapping_MSigDB.v2026.2.Mm.chip` | the same for mouse (Ensembl 116) |
| `map.rat.tsv.gz` | `Rat_Ensembl_Gene_ID_Mouse_Orthologs_MSigDB.v2026.2.Mm.chip`, `Rat_Gene_Symbol_Remapping_Mouse_Orthologs_MSigDB.v2026.2.Mm.chip` | rat Ensembl id / symbol -> mouse ortholog symbol (Alliance of Genome Resources orthology 9.1.0) |

## Not bundled

- **DESeq2** (R/Bioconductor, LGPL-3.0) and **apeglm** (GPL-2) are the methods PyDESeq2
  reimplements; no R code is in the image.
- **Demo data** (not in the image or the repository): Discover dataset 427, Tumanov AV, Shein SA,
  Lai Z, Chen Y, Ruparel SB, doi:10.26275/gaap-zlpj, CC BY 4.0 (see `CITATION.cff`).
