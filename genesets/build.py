"""
Build-time only (Dockerfile): download the MSigDB 2026.2 gene sets and chip files,
check each against its pinned sha256, and write the compact files processor/gsea.py
reads into OUT (default /opt/genesets):

  hallmark.{human,mouse}.gmt.gz   MSigDB H / MH
  go-bp.{human,mouse}.gmt.gz      MSigDB C5 GO:BP / M5 GO:BP
  map.{human,mouse,rat}.tsv.gz    lower-case Ensembl gene id or symbol -> MSigDB symbol
                                  (rat -> mouse ortholog symbols)
  genesets.json                   release, sources, sha256, licence

MSigDB is CC BY 4.0 (https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp, checked
2026-10-31; the KEGG and BioCarta sets with other terms are not used).

    python genesets/build.py /opt/genesets
"""
import gzip
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

BASE = "https://data.broadinstitute.org/gsea-msigdb/msigdb"
RELEASE = "2026.2"
FILES = {
    "h.all.Hs": (f"release/{RELEASE}.Hs/h.all.v{RELEASE}.Hs.symbols.gmt",
                 "56bba1300555fb6302844fc8cdca63406eb50fa1a2472ec86468747efb3fb89b"),
    "c5.go.bp.Hs": (f"release/{RELEASE}.Hs/c5.go.bp.v{RELEASE}.Hs.symbols.gmt",
                    "c52dedc2e3aa722406a2a0829b29e11d22deb962714a90b763a0f8cc45c17c9e"),
    "mh.all.Mm": (f"release/{RELEASE}.Mm/mh.all.v{RELEASE}.Mm.symbols.gmt",
                  "a3068abd9d5cba5c4415109dbe5233de27e5e51b027ad7d1e8451251ce67491f"),
    "m5.go.bp.Mm": (f"release/{RELEASE}.Mm/m5.go.bp.v{RELEASE}.Mm.symbols.gmt",
                    "9f81390893bed8e0ad838c1c0cac8f27ff71583c99b73cc9277ab31b94fea73e"),
    "human.ensembl": (f"annotations/human/Human_Ensembl_Gene_ID_MSigDB.v{RELEASE}.Hs.chip",
                      "0ad1246a70df668afe0f752d8ee2299fba6d849b561b3c5c5d7f47a75cfbe17b"),
    "human.symbol": (f"annotations/human/Human_Gene_Symbol_with_Remapping_MSigDB.v{RELEASE}.Hs.chip",
                     "c9ca95d4e06eade8345eeb5f7963ad7733c328429751ffc708ee9185a692159e"),
    "mouse.ensembl": (f"annotations/mouse/Mouse_Ensembl_Gene_ID_MSigDB.v{RELEASE}.Mm.chip",
                      "e723090866bef01dde9458efaabcc4690de54b77156f0834d6300b3ea1d4a2ea"),
    "mouse.symbol": (f"annotations/mouse/Mouse_Gene_Symbol_Remapping_MSigDB.v{RELEASE}.Mm.chip",
                     "e0af57a83991f58225d8072b9b0d2191f6705f407ab5c7bd01bbec4ce4acab90"),
    "rat.ensembl": (f"annotations/mouse/Rat_Ensembl_Gene_ID_Mouse_Orthologs_MSigDB.v{RELEASE}.Mm.chip",
                    "bf45071d2704fd3f811646777e0b135644e102ce1b0903c19031d2450599be77"),
    "rat.symbol": (f"annotations/mouse/Rat_Gene_Symbol_Remapping_Mouse_Orthologs_MSigDB.v{RELEASE}.Mm.chip",
                   "40b903fd8d355f7d566b8f457415259763eb555ffccec6988b6e891e5607916e"),
}
GMT = {"hallmark.human": "h.all.Hs", "go-bp.human": "c5.go.bp.Hs", "hallmark.mouse": "mh.all.Mm",
       "go-bp.mouse": "m5.go.bp.Mm"}


def fetch(name: str) -> bytes:
    path, sha = FILES[name]
    with urllib.request.urlopen(f"{BASE}/{path}", timeout=300) as r:
        data = r.read()
    got = hashlib.sha256(data).hexdigest()
    if got != sha:
        sys.exit(f"{path}: sha256 {got} != pinned {sha}")
    return data


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    raw = {name: fetch(name) for name in FILES}
    for target, name in GMT.items():
        lines = [ln for ln in raw[name].decode().splitlines() if ln.strip()]
        (out / f"{target}.gmt.gz").write_bytes(gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0))
    for org in ("human", "mouse", "rat"):
        mapping: dict[str, str] = {}
        for kind in ("ensembl", "symbol"):
            for ln in raw[f"{org}.{kind}"].decode().splitlines()[1:]:
                parts = ln.split("\t")
                if len(parts) >= 2 and parts[0].strip() and parts[1].strip():
                    mapping.setdefault(parts[0].strip().lower(), parts[1].strip())
        text = "".join(f"{k}\t{v}\n" for k, v in sorted(mapping.items()))
        (out / f"map.{org}.tsv.gz").write_bytes(gzip.compress(text.encode(), mtime=0))
    (out / "genesets.json").write_text(json.dumps({
        "source": f"MSigDB {RELEASE} (Hs, Mm)", "release": RELEASE, "licence": "CC BY 4.0",
        "go": "GO-basic obo released 2026-07-26 with NCBI gene2go annotations of 2026-09-30 (MSigDB release notes)",
        "orthology": "rat -> mouse: MSigDB chip files, Alliance of Genome Resources orthology 9.1.0",
        "files": {name: {"url": f"{BASE}/{p}", "sha256": s} for name, (p, s) in FILES.items()},
    }, indent=2) + "\n")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "/opt/genesets"))
