"""Which Prochlorococcus genomes carry each KO used on the pathway map.

1. ko_names.csv      : KG name of every KO in pathway_steps.py (checks the KO assignments; not_found = KO absent from KG)
2. ko_presence.csv   : organism x KO -> locus tags. 44 genomes via genes_by_ontology; PAC1 (cannot be addressed by
                       organism= because the name also matches EQPAC1) via run_cypher on Gene_has_kegg_ko. The Cypher
                       query is run for all 45 genomes and must agree with genes_by_ontology on the 44 (asserted).
3. gene_presence.csv / route_presence.csv : strain counts per map gene and per route (all genes present).
"""
from pathlib import Path
import sys
import pandas as pd
import multiomics_explorer as m

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pathway_steps import ROUTES, all_kos, all_groups  # noqa: E402

OUT = HERE.parent / "data"
kos = all_kos()
calls = []

# 1. KO names
names = {}
for i in range(0, len(kos), 40):  # batched: the tool pages at 50 terms
    t = m.ontology_term_details(term_ids=kos[i:i + 40])
    assert t["returned"] == t["total_matching"], ("ko names truncated", t["returned"], t["total_matching"])
    names.update({r["term_id"]: r["name"] for r in t["results"]})
pd.DataFrame([dict(ko_id=k, kg_name=names.get(k)) for k in kos]).to_csv(OUT / "ko_names.csv", index=False, encoding="utf-8")
calls.append(dict(tool="ontology_term_details[map KOs]", n=len(kos), found=len(names), not_found=len(kos) - len(names)))

# 2a. API, 44 genomes
gen = pd.read_csv(OUT / "genomes.csv", encoding="utf-8")
genomes = sorted(gen.organism_name)
api_rows = []
for g in genomes:
    try:
        r = m.genes_by_ontology(ontology="kegg", organism=g, term_ids=kos, min_gene_set_size=1, max_gene_set_size=None, limit=None)
    except Exception as e:
        calls.append(dict(tool="genes_by_ontology[map KOs]", organism=g, error=str(e)[:120])); continue
    assert r["returned"] == r["total_matching"]
    calls.append(dict(tool="genes_by_ontology[map KOs]", organism=g, total_matching=r["total_matching"]))
    api_rows += [dict(organism=g, ko_id=x["term_id"], locus_tag=x["locus_tag"]) for x in r["results"]]
api = pd.DataFrame(api_rows)

# 2b. Cypher, all genomes (for PAC1, and as a cross-check)
ko_list = ", ".join(f"'{k}'" for k in kos)
org_list = ", ".join(f"'{g}'" for g in genomes)
q = (f"MATCH (g:Gene)-[:Gene_has_kegg_ko]->(k:KeggTerm) WHERE g.organism_name IN [{org_list}] AND k.id IN [{ko_list}] "
     "RETURN g.organism_name AS organism, k.id AS ko_id, g.locus_tag AS locus_tag")
c = m.run_cypher(query=q, limit=None)
cy = pd.DataFrame(c["results"])
calls.append(dict(tool="run_cypher[Gene_has_kegg_ko]", total_matching=len(cy), truncated=c["truncated"]))
assert not c["truncated"]

done = set(api.organism)
key = lambda d: set(map(tuple, d[["organism", "ko_id", "locus_tag"]].values))
a, b = key(api), key(cy[cy.organism.isin(done)])
assert a == b, (len(a - b), len(b - a), list(a - b)[:5], list(b - a)[:5])
calls.append(dict(tool="cross-check api vs cypher", organisms=len(done), rows=len(a), agree=True))
pres = pd.concat([api, cy[~cy.organism.isin(done)]]).drop_duplicates()
# ortholog-group stand-ins (genes whose KO is mis-assigned in the KG; see pathway_steps.py)
for og in all_groups():
    h = m.genes_by_homolog_group(group_ids=[og], limit=None)
    assert h["returned"] == h["total_matching"]
    rows = [dict(organism=x["organism_name"], ko_id=og, locus_tag=x["locus_tag"]) for x in h["results"] if x["organism_name"] in genomes]
    calls.append(dict(tool="genes_by_homolog_group", query=og, total_matching=h["total_matching"], kept_prochlorococcus=len(rows)))
    pres = pd.concat([pres, pd.DataFrame(rows)])
pres.to_csv(OUT / "ko_presence.csv", index=False, encoding="utf-8")
N = pres.organism.nunique()
assert N == len(genomes), (N, len(genomes))

# 3. gene / route presence
has = pres.groupby("organism").ko_id.apply(set).to_dict()
grows, rrows = [], []
for r in ROUTES:
    complete = set(genomes)
    for label, alts in r["genes"]:
        orgs = {o for o in genomes if has.get(o, set()) & set(alts)}
        complete &= orgs
        per_ko = {k: sum(k in has.get(o, set()) for o in genomes) for k in alts}
        grows.append(dict(route=r["id"], kind=r["kind"], src=r["src"], dst=r["dst"], gene=label, n_strains=len(orgs), of=N,
                          kos=";".join(alts), per_ko=";".join(f"{k.split(':')[1]}={v}" for k, v in per_ko.items()),
                          kg_names=" | ".join(f"{k.split(':')[1]}: {names.get(k, 'ortholog group' if not k.startswith('kegg.') else 'NOT IN KG')}" for k in alts),
                          missing_in=";".join(sorted(set(genomes) - orgs)) if 0 < len(orgs) < N else ""))
    rrows.append(dict(route=r["id"], kind=r["kind"], src=r["src"], dst=r["dst"], n_complete=len(complete), of=N,
                      incomplete_in=";".join(sorted(set(genomes) - complete)) if 0 < len(complete) < N else ""))
gp, rp = pd.DataFrame(grows), pd.DataFrame(rrows)
gp.to_csv(OUT / "map_gene_presence.csv", index=False, encoding="utf-8")
rp.to_csv(OUT / "map_route_presence.csv", index=False, encoding="utf-8")
pd.DataFrame(calls).to_csv(OUT / "kg_calls_map.csv", index=False, encoding="utf-8")

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 140); pd.set_option("display.max_rows", 500)
print(pd.DataFrame(calls).tail(4).to_string())
print(gp[["route", "gene", "n_strains", "per_ko", "kg_names"]].to_string())
