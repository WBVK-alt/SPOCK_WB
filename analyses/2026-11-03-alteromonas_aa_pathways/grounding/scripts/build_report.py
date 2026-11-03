"""Render the Plan-phase grounding inventory as one self-contained HTML page.

Reads only grounding/data/*.csv (written by pull_inventory.py). Raw counts, no scoring.
Writes grounding/reports/inventory.html.
"""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA, REP = ROOT / "data", ROOT / "reports"
REP.mkdir(exist_ok=True)


def short(o: str) -> str:
    return (o.replace("Alteromonas macleodii ", "").replace("Alteromonas mediterranea ", "")
             .replace("Prochlorococcus ", ""))


kegg = pd.read_csv(DATA / "kegg_pathway_genes.csv", encoding="utf-8")
ko = pd.read_csv(DATA / "kegg_gene_ko.csv", encoding="utf-8")
go = pd.read_csv(DATA / "go_aa_genes.csv", encoding="utf-8")
calls = pd.read_csv(DATA / "kg_calls.csv", encoding="utf-8")
c_kegg = pd.read_csv(DATA / "counts_kegg_pathway_by_genome.csv", index_col=0, encoding="utf-8")
c_ko = pd.read_csv(DATA / "counts_kegg_distinct_ko_by_genome.csv", index_col=0, encoding="utf-8")
c_go = pd.read_csv(DATA / "counts_go_aa_by_genome.csv", index_col=0, encoding="utf-8")
genomes = list(c_kegg.columns)
groups = ["Prochlorococcus" if g.startswith("Prochlorococcus") else "Alteromonas" for g in genomes]

ko_by_gene = ko.groupby(["organism", "locus_tag"]).ko_id.apply(lambda s: ";".join(sorted(set(s.str.replace("kegg.orthology:", ""))))).to_dict()


def label(r):
    n = r.gene_name if isinstance(r.gene_name, str) and r.gene_name else ""
    p = r["product"] if isinstance(r["product"], str) else ""
    return f"{n} ({r.locus_tag})" + (f" — {p}" if p else "") if n else f"{r.locus_tag}" + (f" — {p}" if p else "")


def cell_genes(df, rowcol, rowval_fn):
    out = {}
    for (rv, o), sub in df.groupby([rowcol, "organism"]):
        sub = sub.drop_duplicates("locus_tag").sort_values("locus_tag")
        out[f"{rowval_fn(rv)}|{o}"] = [label(r) + (f" [KO {ko_by_gene.get((o, r.locus_tag), '')}]" if rowcol == "pathway" else "")
                                       for _, r in sub.iterrows()]
    return out


go["row"] = go.amino_acid + " " + go.direction
panels = [
    dict(id="go", title="GO biological process: genes per amino acid, made vs broken down",
         note="Genes annotated to the GO term '<amino acid> biosynthetic process' or '… catabolic process' "
              "(exact term-name match; terms used are listed at the bottom). Rows with no matching GO term are marked.",
         unit="genes", rows=list(c_go.index), m=c_go.values.tolist(), genes=cell_genes(go, "row", str)),
    dict(id="kegg", title="KEGG amino-acid metabolism pathways: genes per genome",
         note="Genes annotated to each of the 14 KEGG pathways under 'Amino acid metabolism'. A gene can sit in several pathways.",
         unit="genes", rows=list(c_kegg.index), m=c_kegg.values.tolist(), genes=cell_genes(kegg, "pathway", str)),
    dict(id="ko", title="KEGG amino-acid metabolism pathways: distinct gene families (KOs) per genome",
         note="Distinct KEGG orthology (KO) families among those genes — paralogs collapse to one family.",
         unit="KO families", rows=list(c_ko.index), m=c_ko.values.tolist(), genes=cell_genes(kegg, "pathway", str)),
]

go_terms = pd.read_csv(DATA / "go_terms_used.csv", encoding="utf-8")
nomatch = calls[(calls.tool == "search_ontology[go_bp]") & (calls.kept == 0)]["query"].tolist()
payload = dict(genomes=genomes, short=[short(g) for g in genomes], groups=groups, panels=panels,
               nomatch=[q.replace(" biosynthetic", " biosynthetic").strip() for q in nomatch])

calls_html = calls.fillna("").to_html(index=False, classes="tbl", border=0)
terms_html = go_terms.to_html(index=False, classes="tbl", border=0)

HTML = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>Amino-acid pathway genes — raw inventory (Plan-phase grounding)</title>
<style>
:root{--surface:#fcfcfb;--ink:#1f1f1d;--ink2:#52514e;--muted:#8a8984;--rule:#e4e3df;--zero:#ffffff;}
@media (prefers-color-scheme:dark){:root{--surface:#1a1a19;--ink:#ecebe8;--ink2:#b9b8b3;--muted:#85847f;--rule:#383835;--zero:#1a1a19;}}
body{background:var(--surface);color:var(--ink);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif;margin:24px 32px;max-width:1300px}
h1{font-size:20px;margin:0 0 4px} h2{font-size:16px;margin:32px 0 2px} p.note{color:var(--ink2);margin:2px 0 10px;max-width:900px}
.meta{color:var(--muted);font-size:12px}
svg text{font-family:inherit}
.legend{display:flex;align-items:center;gap:10px;font-size:12px;color:var(--ink2);margin:4px 0 6px}
.legend .sw{display:inline-block;width:18px;height:12px;border-radius:2px;vertical-align:middle}
#tip{position:fixed;pointer-events:none;background:var(--surface);color:var(--ink);border:1px solid var(--rule);border-radius:6px;
 padding:8px 10px;font-size:12px;max-width:520px;box-shadow:0 4px 14px rgba(0,0,0,.15);display:none;z-index:9}
#tip b{font-weight:600} #tip ul{margin:4px 0 0 16px;padding:0} #tip li{margin:0}
.tbl{border-collapse:collapse;font-size:12px;margin:6px 0} .tbl th,.tbl td{border-bottom:1px solid var(--rule);padding:3px 8px;text-align:left}
details{margin:18px 0} summary{cursor:pointer;color:var(--ink2)}
.pin{font-size:12px;color:var(--ink2);white-space:pre-wrap;border-left:3px solid var(--rule);padding-left:10px;margin-top:6px}
</style></head><body>
<h1>Amino-acid pathway genes in Alteromonas and its Prochlorococcus partners — raw inventory</h1>
<div class="meta">Plan-phase grounding · counts straight from the KG annotation, no completeness scoring ·
9 Alteromonas genomes (AltDE excluded — the tool cannot address it apart from AltDE1) + 3 Prochlorococcus coculture partners ·
hover a cell for its genes; click to pin the list below the chart</div>
<div id="root"></div>
<details><summary>GO terms used (exact name match)</summary>__TERMS__</details>
<details><summary>Every KG call behind this page</summary>__CALLS__</details>
<div id="tip"></div>
<script>
const D = __DATA__;
const PURPLES=["#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d"];
function col(v,max){ if(v===0) return null; const t=Math.min(1,v/max); return PURPLES[Math.min(PURPLES.length-1,Math.floor(t*(PURPLES.length-1e-9)))]; }
function ink(v,max){ const t=v/max; return t>0.45?"#ffffff":"#1f1f1d"; }
const tip=document.getElementById("tip");
const root=document.getElementById("root");
D.panels.forEach(P=>{
  const h=document.createElement("h2"); h.textContent=P.title; root.appendChild(h);
  const n=document.createElement("p"); n.className="note"; n.textContent=P.note; root.appendChild(n);
  const max=Math.max(...P.m.flat());
  const lg=document.createElement("div"); lg.className="legend";
  lg.innerHTML=`<span>${P.unit}:</span><span class="sw" style="background:var(--zero);outline:1.5px dashed var(--muted);outline-offset:-1.5px"></span> 0`+
    `<span class="sw" style="background:${PURPLES[0]}"></span>1`+
    `<span style="display:inline-block;width:120px;height:12px;border-radius:2px;background:linear-gradient(90deg,${PURPLES.join(",")})"></span>${max}`;
  root.appendChild(lg);
  const cw=58,ch=22,left=P.id==="go"?190:330,top=58,gap=2,sep=14;
  const nAlt=D.groups.filter(g=>g==="Alteromonas").length;
  const x=j=>left+j*cw+(j>=nAlt?sep:0);
  const W=x(D.genomes.length)+10, H=top+P.rows.length*ch+10;
  const NS="http://www.w3.org/2000/svg";
  const svg=document.createElementNS(NS,"svg"); svg.setAttribute("width",W); svg.setAttribute("height",H); svg.setAttribute("role","img");
  svg.setAttribute("aria-label",P.title);
  const T=(x0,y0,s,a={})=>{const t=document.createElementNS(NS,"text");t.setAttribute("x",x0);t.setAttribute("y",y0);t.textContent=s;for(const k in a)t.setAttribute(k,a[k]);svg.appendChild(t);return t;};
  T(x(0),14,"Alteromonas",{"font-size":12,"font-weight":600,fill:"var(--ink2)"});
  T(x(nAlt),14,"Prochlorococcus",{"font-size":12,"font-weight":600,fill:"var(--ink2)"});
  D.short.forEach((s,j)=>T(x(j)+cw/2,top-8,s,{"font-size":11,"text-anchor":"middle",fill:"var(--ink2)"}));
  const ln=document.createElementNS(NS,"line");ln.setAttribute("x1",x(nAlt)-sep/2);ln.setAttribute("x2",x(nAlt)-sep/2);
  ln.setAttribute("y1",20);ln.setAttribute("y2",H-6);ln.setAttribute("stroke","var(--muted)");ln.setAttribute("stroke-width",1);svg.appendChild(ln);
  P.rows.forEach((r,i)=>{
    const y=top+i*ch;
    const nm = P.id==="go" && D.nomatch.includes(r);
    T(left-8,y+ch/2+4,r+(nm?"  (no GO term)":""),{"font-size":12,"text-anchor":"end",fill:nm?"var(--muted)":"var(--ink)"});
    D.genomes.forEach((g,j)=>{
      const v=P.m[i][j], c=col(v,max);
      const rc=document.createElementNS(NS,"rect");
      rc.setAttribute("x",x(j)+gap/2);rc.setAttribute("y",y+gap/2);rc.setAttribute("width",cw-gap);rc.setAttribute("height",ch-gap);rc.setAttribute("rx",3);
      if(c){rc.setAttribute("fill",c);} else {rc.setAttribute("fill","var(--zero)");rc.setAttribute("stroke","var(--muted)");rc.setAttribute("stroke-dasharray","3 2");rc.setAttribute("stroke-width",1);}
      rc.style.cursor="pointer";
      const genes=P.genes[r+"|"+g]||[];
      const show=e=>{tip.style.display="block";
        tip.innerHTML=`<b>${D.short[j]}</b> · ${r}<br>${v} ${P.unit}`+(genes.length?`<ul>${genes.slice(0,25).map(s=>`<li>${s}</li>`).join("")}</ul>${genes.length>25?`<i>+${genes.length-25} more — click to list all</i>`:""}`:"");
        const px=Math.min(e.clientX+14,window.innerWidth-540);tip.style.left=px+"px";tip.style.top=(e.clientY+14)+"px";};
      rc.addEventListener("mousemove",show); rc.addEventListener("mouseleave",()=>tip.style.display="none");
      rc.addEventListener("click",()=>{pin.textContent=`${D.genomes[j]} · ${r} · ${v} ${P.unit}\n`+genes.join("\n");});
      svg.appendChild(rc);
      T(x(j)+cw/2,y+ch/2+4,String(v),{"font-size":11,"text-anchor":"middle",fill:c?ink(v,max):"var(--muted)","pointer-events":"none"});
    });
  });
  root.appendChild(svg);
  const pin=document.createElement("div"); pin.className="pin"; pin.textContent="(click a cell to list its genes here)"; root.appendChild(pin);
});
</script></body></html>"""

out = HTML.replace("__DATA__", json.dumps(payload)).replace("__CALLS__", calls_html).replace("__TERMS__", terms_html)
(REP / "inventory.html").write_text(out, encoding="utf-8")
print("wrote", REP / "inventory.html", len(out), "bytes")
