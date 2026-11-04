"""Plan-phase grounding: raw inventory of amino-acid pathway genes per Prochlorococcus genome.

No scoring, no completeness calls. Pulls, per genome:
  1. genes annotated to each KEGG 'Amino acid metabolism' pathway (subcategory 09105), with each gene's KO(s)
  2. genes annotated to GO-BP '<amino acid> biosynthetic process' / '... catabolic process' (20 amino acids)
Genomes: every Prochlorococcus genome strain in the KG.
Adapted from analyses/2026-11-03-alteromonas_aa_pathways/grounding/scripts/pull_inventory.py.
"""
from pathlib import Path
import re
import pandas as pd
import multiomics_explorer as m

OUT = Path(__file__).resolve().parents[1] / "data"
calls, skipped = [], []

def log(tool, **kw):
    calls.append(dict(tool=tool, **kw))

# --- genomes ---------------------------------------------------------------
org = m.list_organisms(organism_names=["Prochlorococcus"], limit=None)
genomes = sorted(r["organism_name"] for r in org["results"] if r["organism_type"] == "genome_strain")
pd.DataFrame([{k: r.get(k) for k in ("organism_name", "strain", "clade", "gene_count", "experiment_count", "omics_types")}
              for r in org["results"]]).to_csv(OUT / "genomes.csv", index=False, encoding="utf-8")
log("list_organisms", query="Prochlorococcus", total_matching=org["total_matching"], kept=len(genomes))

# --- KEGG amino-acid pathways ---------------------------------------------
sub = m.ontology_term_details(term_ids=["kegg.subcategory:09105"])
pathways = {c["id"]: c["name"] for c in sub["results"][0]["children"]}
log("ontology_term_details", query="kegg.subcategory:09105", total_matching=len(pathways))

kegg_rows, ko_rows = [], []
for g in list(genomes):
    try:
        r = m.genes_by_ontology(ontology="kegg", organism=g, term_ids=list(pathways),
                                min_gene_set_size=1, max_gene_set_size=None, limit=None, verbose=True)
    except Exception as e:  # e.g. organism name matching several strains
        skipped.append(dict(organism=g, error=str(e)[:200])); genomes.remove(g); continue
    assert r["returned"] == r["total_matching"], (g, r["returned"], r["total_matching"])
    log("genes_by_ontology[kegg pathways]", organism=g, total_matching=r["total_matching"], returned=r["returned"])
    for x in r["results"]:
        kegg_rows.append(dict(organism=g, pathway_id=x["term_id"], pathway=x["term_name"],
                              locus_tag=x["locus_tag"], gene_name=x.get("gene_name"), product=x.get("product")))
    tags = sorted({x["locus_tag"] for x in r["results"]})
    if tags:
        t = m.gene_ontology_terms(locus_tags=tags, organism=g, ontology="kegg", level=3, limit=None)
        assert t["returned"] == t["total_matching"], (g, "KO")
        log("gene_ontology_terms[kegg level3]", organism=g, total_matching=t["total_matching"], returned=t["returned"])
        for x in t["results"]:
            ko_rows.append(dict(organism=g, locus_tag=x["locus_tag"], ko_id=x.get("term_id"), ko_name=x.get("term_name")))

kegg = pd.DataFrame(kegg_rows); ko = pd.DataFrame(ko_rows)
kegg.to_csv(OUT / "kegg_pathway_genes.csv", index=False, encoding="utf-8")
ko.to_csv(OUT / "kegg_gene_ko.csv", index=False, encoding="utf-8")

# --- GO-BP per amino acid -------------------------------------------------
AA = ["alanine", "arginine", "asparagine", "aspartate", "cysteine", "glutamate", "glutamine",
      "glycine", "histidine", "isoleucine", "leucine", "lysine", "methionine", "phenylalanine",
      "proline", "serine", "threonine", "tryptophan", "tyrosine", "valine"]
go_terms = []
for aa in AA:
    for kind in ("biosynthetic", "catabolic"):
        s = m.search_ontology(search_text=f"{aa} AND {kind}", ontology="go_bp", limit=50)
        pat = re.compile(rf"^(L-)?{aa} {kind} process$", re.I)
        hits = [x for x in s["results"] if pat.match(x["name"])]
        log("search_ontology[go_bp]", query=f"{aa} {kind}", total_matching=s["total_matching"], kept=len(hits),
            kept_names=";".join(h["name"] for h in hits))
        go_terms += [dict(amino_acid=aa, direction=kind, go_id=h["id"], go_name=h["name"]) for h in hits]
go_terms = pd.DataFrame(go_terms)
go_terms.to_csv(OUT / "go_terms_used.csv", index=False, encoding="utf-8")

go_rows = []
for g in genomes:
    r = m.genes_by_ontology(ontology="go_bp", organism=g, term_ids=go_terms.go_id.tolist(),
                            min_gene_set_size=1, max_gene_set_size=None, limit=None)
    assert r["returned"] == r["total_matching"], (g, "GO")
    log("genes_by_ontology[go_bp aa terms]", organism=g, total_matching=r["total_matching"], returned=r["returned"])
    go_rows += [dict(organism=g, go_id=x["term_id"], go_name=x["term_name"], locus_tag=x["locus_tag"],
                     gene_name=x.get("gene_name"), product=x.get("product")) for x in r["results"]]
go = pd.DataFrame(go_rows).merge(go_terms[["go_id", "amino_acid", "direction"]], on="go_id", how="left")
go.to_csv(OUT / "go_aa_genes.csv", index=False, encoding="utf-8")
pd.DataFrame(calls).to_csv(OUT / "kg_calls.csv", index=False, encoding="utf-8")
pd.DataFrame(skipped, columns=["organism", "error"]).to_csv(OUT / "skipped_genomes.csv", index=False, encoding="utf-8")

# --- count matrices (raw counts only) -------------------------------------
kc = kegg.groupby(["pathway", "organism"]).locus_tag.nunique().unstack(fill_value=0).reindex(columns=genomes, fill_value=0)
kc.to_csv(OUT / "counts_kegg_pathway_by_genome.csv", encoding="utf-8")
ko_n = kegg.merge(ko, on=["organism", "locus_tag"]).groupby(["pathway", "organism"]).ko_id.nunique().unstack(fill_value=0).reindex(columns=genomes, fill_value=0)
ko_n.to_csv(OUT / "counts_kegg_distinct_ko_by_genome.csv", encoding="utf-8")
go["row"] = go.amino_acid + " " + go.direction
gc = go.groupby(["row", "organism"]).locus_tag.nunique().unstack(fill_value=0)
gc = gc.reindex(index=[f"{a} {k}" for a in AA for k in ("biosynthetic", "catabolic")], columns=genomes, fill_value=0)
gc.to_csv(OUT / "counts_go_aa_by_genome.csv", encoding="utf-8")
print("skipped:", skipped)
print(pd.DataFrame(calls).groupby("tool").size())
print(kc.T.describe().T[["min", "50%", "max"]].to_string())
print(gc.T.describe().T[["min", "50%", "max"]].to_string())
