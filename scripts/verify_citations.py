#!/usr/bin/env python3
"""Fetch each DOI from Crossref (papers) or DataCite (datasets) and PubMed (PMID via esearch on the DOI,
then esummary), and write processor/references.json with full author lists. Run by hand when the list
changes; the output is committed (the app needs no network at run time).

    python scripts/verify_citations.py processor/references.json
"""
import html
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date

DOIS = {
    "pydeseq2": "10.1093/bioinformatics/btad547",
    "deseq2": "10.1186/s13059-014-0550-8",
    "apeglm": "10.1093/bioinformatics/bty895",
    "vst": "10.1186/gb-2010-11-10-r106",
    "bh": "10.1111/j.2517-6161.1995.tb02031.x",
    "gsea": "10.1073/pnas.0506580102",
    "gseapy": "10.1093/bioinformatics/btac757",
    "msigdb_hallmark": "10.1016/j.cels.2015.12.004",
    "msigdb_mouse": "10.1038/s41592-023-02014-7",
    "go": "10.1038/75556",
    "go_2023": "10.1093/genetics/iyad031",
    "love_workflow": "10.12688/f1000research.7035.1",
    "conesa": "10.1186/s13059-016-0881-8",
    "tg_paper": "10.1038/s41598-023-40380-6",
    "tg_data": "10.26275/gaap-zlpj",
    "oa_data": "10.26275/robc-8zol",
    "oa_paper": "10.1016/j.joca.2025.08.013",
}
JOURNAL_SHORT = {"Bioinformatics": "Bioinformatics", "Genome Biology": "Genome Biol",
                 "Proceedings of the National Academy of Sciences": "Proc Natl Acad Sci U S A",
                 "Cell Systems": "Cell Syst", "Nature Methods": "Nat Methods", "Nature Genetics": "Nat Genet",
                 "GENETICS": "Genetics", "Genetics": "Genetics", "F1000Research": "F1000Res",
                 "Scientific Reports": "Sci Rep", "Osteoarthritis and Cartilage": "Osteoarthritis Cartilage",
                 "Journal of the Royal Statistical Society Series B: Statistical Methodology": "J R Stat Soc Series B"}


def get(url, tries=5):
    req = urllib.request.Request(url, headers={"User-Agent": "pennsieve-bulk-de citation check (mailto:support@pennsieve.io)"})
    for k in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as exc:
            if exc.code != 429 or k == tries - 1:
                raise
            time.sleep(2 ** k)


def short(authors):
    names = [f"{a['family']} {''.join(p[0] for p in a.get('given', '').replace('-', ' ').replace('.', ' ').split())}".strip()
             for a in authors]
    return ", ".join(names) if len(names) <= 6 else ", ".join(names[:3]) + ", et al"


def pubmed(doi):
    q = urllib.parse.quote(f"{doi}[doi]")
    ids = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={q}&retmode=json")["esearchresult"]["idlist"]
    time.sleep(0.4)
    if len(ids) != 1:
        return None
    s = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={ids[0]}&retmode=json")["result"][ids[0]]
    time.sleep(0.4)
    return {"pmid": ids[0], "title": s["title"], "n_authors": len(s["authors"]), "first": s["authors"][0]["name"],
            "year": s["pubdate"][:4], "source": s["source"], "volume": s.get("volume"), "pages": s.get("pages")}


def crossref(doi):
    m = get(f"https://api.crossref.org/works/{urllib.parse.quote(doi)}")["message"]
    authors = []
    for a in m.get("author", []):     # Crossref lists some consortium papers' authors twice
        a = {"family": a.get("family") or a.get("name"), "given": a.get("given", "")}
        if a not in authors:
            authors.append(a)
    journal = html.unescape((m.get("container-title") or [""])[0])
    year = (m.get("published-print") or m.get("published") or m.get("issued"))["date-parts"][0][0]
    return {"doi": doi, "source": "crossref", "title": html.unescape(" ".join(m["title"][0].split())), "authors": authors,
            "journal": journal, "journal_short": JOURNAL_SHORT.get(journal, journal), "year": year,
            "volume": m.get("volume", ""), "issue": m.get("issue", ""),
            "pages": m.get("page") or m.get("article-number", ""), "type": "article"}


def datacite(doi):
    a = get(f"https://api.datacite.org/dois/{urllib.parse.quote(doi)}")["data"]["attributes"]
    authors = [{"family": c.get("familyName") or c["name"], "given": c.get("givenName", ""),
                "orcid": next((i["nameIdentifier"] for i in c.get("nameIdentifiers", [])
                               if "orcid" in i.get("nameIdentifierScheme", "").lower()), "")} for c in a["creators"]]
    return {"doi": doi, "source": "datacite", "title": a["titles"][0]["title"], "authors": authors,
            "journal": a.get("publisher") if isinstance(a.get("publisher"), str) else a["publisher"].get("name"),
            "year": a["publicationYear"], "volume": "", "issue": "", "pages": "", "type": "data",
            "licence": ", ".join(r.get("rights", "") for r in a.get("rightsList", []))}


out = {}
for key, doi in DOIS.items():
    r = datacite(doi) if doi.startswith("10.26275") else crossref(doi)
    r["journal_short"] = r.get("journal_short") or r["journal"]
    r["short_authors"] = short(r["authors"])
    pm = pubmed(doi) if r["type"] == "article" else None
    r["pmid"] = pm["pmid"] if pm else ""
    r["checked"] = date.today().isoformat()
    out[key] = r
    print(f"{key}: {r['short_authors']}. {r['title']}. {r['journal_short']} {r['year']};{r['volume']}({r['issue']}):{r['pages']}. "
          f"doi:{doi} n_authors={len(r['authors'])}")
    if pm:
        print(f"   PubMed {pm['pmid']}: {pm['title']} | {pm['first']} +{pm['n_authors'] - 1} | {pm['year']} {pm['source']} {pm['volume']}:{pm['pages']}")
    elif r["type"] == "article":
        print("   not in PubMed")
    time.sleep(0.3)
with open(sys.argv[1], "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
    fh.write("\n")
