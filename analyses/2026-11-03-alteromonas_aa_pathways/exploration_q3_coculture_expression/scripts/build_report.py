"""Q3 infographic: Alteromonas amino-acid genes in coculture with Prochlorococcus vs axenic.

Reads only exploration_q3_coculture_expression/data/*.csv. Writes reports/q3_infographic.html.
"""
from pathlib import Path
import json
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
D, R = HERE / "data", HERE / "reports"
R.mkdir(exist_ok=True)

rows = pd.read_csv(D / "de_rows.csv", encoding="utf-8")
memb = pd.read_csv(D / "gene_set_membership.csv", encoding="utf-8")
summ = pd.read_csv(D / "set_by_timepoint_summary.csv", encoding="utf-8")
gw = pd.read_csv(D / "genome_wide_status.csv", encoding="utf-8").sort_values("timepoint_order")
ora = pd.read_csv(D / "ora_kegg_aa_pathways.csv", encoding="utf-8")
days = gw.timepoint.tolist()

# --- tiles ---------------------------------------------------------------
n_genes = rows.locus_tag.nunique()
tiles = []
for _, g in gw.iterrows():
    r = rows[rows.timepoint == g.timepoint]
    up, dn = int((r.expression_status == "significant_up").sum()), int((r.expression_status == "significant_down").sum())
    tiles.append(dict(day=g.timepoint, up=up, down=dn, n=int(r.locus_tag.nunique()),
                      pct=round(100 * (up + dn) / r.locus_tag.nunique(), 1),
                      gw_pct=round(100 * (g.significant_up + g.significant_down) / g.genes_measured, 1),
                      gw_up=int(g.significant_up), gw_down=int(g.significant_down), gw_n=int(g.genes_measured)))
ora_n = len(ora); ora_sig = int((ora.p_adjust < 0.05).sum()); ora_min = ora.sort_values("p_adjust").iloc[0]

# --- KEGG pathway diverging bars -----------------------------------------
kp = summ[summ.set_kind == "KEGG pathway"]
pw_order = kp.groupby("set_name").n_genes.first().sort_values(ascending=False).index.tolist()
bars = [dict(name=p, n=int(kp[kp.set_name == p].n_genes.iloc[0]),
             days=[dict(day=d, up=int(x.n_up), down=int(x.n_down), up_genes=x.up_genes if isinstance(x.up_genes, str) else "",
                        down_genes=x.down_genes if isinstance(x.down_genes, str) else "")
                   for d in days for _, x in kp[(kp.set_name == p) & (kp.timepoint == d)].iterrows()])
        for p in pw_order]

# --- gene heatmap: genes significant at >=1 day --------------------------
names = rows.groupby("locus_tag").gene_name.first()
go = memb[memb.set_kind == "GO process"]
kg = memb[memb.set_kind == "KEGG pathway"]
def group_of(t):
    g = go[go.locus_tag == t]
    if len(g):
        return "; ".join(sorted(g.set_name.unique()))
    return "KEGG only: " + "; ".join(sorted(kg[kg.locus_tag == t].set_name.unique()))
sig_tags = sorted(rows[rows.expression_status != "not_significant"].locus_tag.unique())
genes = []
for t in sig_tags:
    r = rows[rows.locus_tag == t].set_index("timepoint")
    nm = names.get(t)
    genes.append(dict(tag=t, label=(f"{nm} ({t})" if isinstance(nm, str) else t), group=group_of(t),
                      product=str(memb[memb.locus_tag == t]["product"].dropna().iloc[0]) if memb[memb.locus_tag == t]["product"].notna().any() else "",
                      sets="; ".join(sorted(memb[memb.locus_tag == t].set_name.unique())),
                      cells=[dict(day=d, lfc=round(float(r.loc[d, "log2fc"]), 2), padj=float(r.loc[d, "padj"]),
                                  status=r.loc[d, "expression_status"]) for d in days]))
genes.sort(key=lambda g: (g["group"], g["label"]))

# --- full table ------------------------------------------------------------
wide = rows.pivot_table(index="locus_tag", columns="timepoint", values="log2fc").reindex(columns=days)
st = rows.pivot_table(index="locus_tag", columns="timepoint", values="expression_status", aggfunc="first").reindex(columns=days)
table = []
for t in wide.index:
    nm = names.get(t)
    table.append(dict(tag=t, name=nm if isinstance(nm, str) else "", sets="; ".join(sorted(memb[memb.locus_tag == t].set_name.unique())),
                      lfc=[round(float(wide.loc[t, d]), 2) for d in days], st=[st.loc[t, d] for d in days]))

pairs = [
    dict(pair="HOT1A3 + MED4", has=True, exp="Coculture vs axenic HOT1A3, RNA-seq, days 11 / 18 / 31, low-N medium, continuous light; both arms sampled the same day",
         src="Weissberg, Aharonovich & Sher 2025 (bioRxiv)"),
    dict(pair="MIT1002 + NATL2A", has=False, exp="KG holds: 24 h and 48 h vs 12 h of coculture (UP genes only published); extended darkness vs diel; dark-tolerant vs parental coculture. No axenic MIT1002 arm.",
         src="Biller 2016, 2018; Coe 2024"),
    dict(pair="EZ55 + MIT9312", has=False, exp="KG holds: MIT9312 partner vs the mean of three partners (not vs axenic); elevated vs ambient CO2 within coculture. The paper's general coculture-vs-axenic table (Table S5) is not among the KG experiments.",
         src="Barreto Filho 2022; Hennon 2018"),
]
payload = dict(days=days, tiles=tiles, n_genes=n_genes, bars=bars, genes=genes, table=table, pairs=pairs,
               ora=dict(n=ora_n, sig=ora_sig, min_term=ora_min.term_name, min_day=ora_min.timepoint,
                        min_dir=ora_min.direction, min_padj=round(float(ora_min.p_adjust), 3)))

HTML = r"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Q3 — Alteromonas amino-acid genes in coculture with Prochlorococcus</title>
<style>
:root{--surface:#fcfcfb;--card:#f5f4f1;--ink:#1f1f1d;--ink2:#52514e;--muted:#8a8984;--rule:#e4e3df;--up:#eb6834;--down:#2a78d6;--mid:#f0efec;--gap:#c9c8c3}
@media (prefers-color-scheme:dark){:root{--surface:#1a1a19;--card:#232322;--ink:#ecebe8;--ink2:#b9b8b3;--muted:#85847f;--rule:#383835;--up:#d95926;--down:#3987e5;--mid:#383835;--gap:#55544f}}
body{background:var(--surface);color:var(--ink);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif;margin:24px 32px;max-width:1200px}
h1{font-size:21px;margin:0 0 4px} h2{font-size:16px;margin:34px 0 4px} .sub{color:var(--ink2);margin:0 0 8px;max-width:920px}
.meta{color:var(--muted);font-size:12px}
.tiles{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:18px 0}
.tile{background:var(--card);border-radius:10px;padding:14px 16px}
.tile .k{font-size:12px;color:var(--ink2)} .tile .v{font-size:30px;font-weight:650;margin:2px 0} .tile .s{font-size:12px;color:var(--ink2)}
.ud{display:flex;gap:10px;font-size:13px;margin-top:4px}.ud span{display:inline-flex;align-items:center;gap:4px}
.dot{width:10px;height:10px;border-radius:2px;display:inline-block}
.pairs{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}
.pair{border-radius:10px;padding:12px 14px;border:1.5px solid var(--rule)} .pair.gap{border-style:dashed;color:var(--ink2)}
.pair b{font-size:15px} .badge{display:inline-block;font-size:11px;font-weight:600;border-radius:4px;padding:1px 6px;margin-left:6px}
.badge.ok{background:var(--ink);color:var(--surface)} .badge.no{border:1px dashed var(--muted);color:var(--ink2)}
.legend{display:flex;gap:14px;align-items:center;font-size:12px;color:var(--ink2);margin:6px 0}
#tip{position:fixed;pointer-events:none;background:var(--surface);color:var(--ink);border:1px solid var(--rule);border-radius:6px;padding:8px 10px;
 font-size:12px;max-width:460px;box-shadow:0 4px 14px rgba(0,0,0,.15);display:none;z-index:9}
table.t{border-collapse:collapse;font-size:12px;width:100%} table.t th,table.t td{border-bottom:1px solid var(--rule);padding:3px 6px;text-align:left}
table.t th{cursor:pointer;color:var(--ink2);position:sticky;top:0;background:var(--surface)} td.num{text-align:right;font-variant-numeric:tabular-nums}
.s-up{color:var(--up);font-weight:600}.s-down{color:var(--down);font-weight:600}
input#q{font:inherit;padding:4px 8px;border:1px solid var(--rule);border-radius:6px;background:var(--surface);color:var(--ink);width:280px}
.tw{max-height:520px;overflow:auto;border:1px solid var(--rule);border-radius:8px}
svg text{font-family:inherit}
</style></head><body>
<h1>Are Alteromonas amino-acid genes switched on and changed when grown with <i>Prochlorococcus</i>?</h1>
<div class="meta">Exploration before a locked proposal · source: multiomics KG 0.1.0-alpha.9 · significance = each paper's own call as stored in the KG · no new statistics except the package's pathway over-representation test · gene sets = HOT1A3 genes in the 14 KEGG amino-acid pathways or GO '&lt;amino acid&gt; biosynthetic/catabolic process'</div>

<h2>1 · Which partner pairs can answer this?</h2>
<p class="sub">Q3 needs Alteromonas measured in coculture <i>and</i> alone. Only one of the three pairs has that contrast in the KG.</p>
<div class="pairs" id="pairs"></div>

<h2>2 · HOT1A3 + MED4: the answer in numbers</h2>
<div class="tiles" id="tiles"></div>

<h2>3 · Per KEGG pathway: genes up and down in coculture</h2>
<p class="sub">Bars left of the line = genes significantly lower in coculture; right = higher. Number in brackets = genes in the pathway; a gene can sit in several pathways. Hover for genes.</p>
<div class="legend"><span><span class="dot" style="background:var(--down)"></span>lower in coculture</span><span><span class="dot" style="background:var(--up)"></span>higher in coculture</span></div>
<div id="bars"></div>

<h2>4 · Every gene called significant on at least one day</h2>
<p class="sub">Colour = log2 fold change, coculture vs axenic, on all three days (including days it was not significant). A ring marks the days the paper called it significant. Grouped by GO amino-acid process (or KEGG pathway where no GO term).</p>
<div class="legend" id="hmleg"></div>
<div id="heat"></div>

<h2>5 · All amino-acid genes</h2>
<p class="sub"><input id="q" placeholder="filter: gene, locus tag, pathway, amino acid"></p>
<div class="tw"><table class="t" id="tbl"></table></div>
<div id="tip"></div>
<script>
const D=__DATA__;
const tip=document.getElementById("tip");
function showTip(e,h){tip.style.display="block";tip.innerHTML=h;tip.style.left=Math.min(e.clientX+14,innerWidth-480)+"px";tip.style.top=(e.clientY+14)+"px";}
function hideTip(){tip.style.display="none";}
const css=v=>getComputedStyle(document.documentElement).getPropertyValue(v).trim();

// 1 pairs
document.getElementById("pairs").innerHTML=D.pairs.map(p=>`<div class="pair ${p.has?"":"gap"}"><b>${p.pair}</b>
<span class="badge ${p.has?"ok":"no"}">${p.has?"coculture vs axenic":"no coculture-vs-axenic data"}</span>
<div style="font-size:12px;margin-top:6px">${p.exp}</div><div class="meta" style="margin-top:4px">${p.src}</div></div>`).join("");

// 2 tiles
const T=document.getElementById("tiles");
T.innerHTML=`<div class="tile"><div class="k">Amino-acid genes measured</div><div class="v">${D.n_genes}</div>
<div class="s">of ${D.n_genes} HOT1A3 amino-acid genes — all detected on all 3 days</div></div>`+
D.tiles.map(t=>`<div class="tile"><div class="k">${t.day}: changed in coculture</div><div class="v">${t.up+t.down}<span style="font-size:15px;font-weight:500;color:var(--ink2)"> / ${t.n} · ${t.pct}%</span></div>
<div class="ud"><span><span class="dot" style="background:var(--up)"></span>${t.up} up</span><span><span class="dot" style="background:var(--down)"></span>${t.down} down</span></div>
<div class="s">whole genome that day: ${t.gw_pct}% (${t.gw_up} up, ${t.gw_down} down of ${t.gw_n})</div></div>`).join("");
const o=document.createElement("p");o.className="sub";
o.innerHTML=`Pathway over-representation test (package <code>pathway_enrichment</code>, 14 pathways × 3 days × up/down = ${D.ora.n} tests, background = genes measured): <b>${D.ora.sig} of ${D.ora.n}</b> pass adjusted p &lt; 0.05. Lowest: ${D.ora.min_term}, ${D.ora.min_dir}, ${D.ora.min_day}, adjusted p = ${D.ora.min_padj}.`;
T.after(o);

// 3 diverging bars, small multiples per day
(function(){
 const NS="http://www.w3.org/2000/svg", rowH=24, left=360, colW=230, gapC=26, top=30;
 const maxN=Math.max(...D.bars.map(b=>b.n)); const half=(colW-10)/2; const sc=half/Math.max(8,...D.bars.flatMap(b=>b.days.map(d=>Math.max(d.up,d.down))));
 const W=left+D.days.length*(colW+gapC), H=top+D.bars.length*rowH+8;
 const svg=document.createElementNS(NS,"svg");svg.setAttribute("width",W);svg.setAttribute("height",H);
 const el=(n,a)=>{const e=document.createElementNS(NS,n);for(const k in a)e.setAttribute(k,a[k]);svg.appendChild(e);return e;};
 const tx=(x,y,s,a={})=>{const t=el("text",Object.assign({x,y},a));t.textContent=s;return t;};
 D.days.forEach((d,j)=>{const cx=left+j*(colW+gapC)+colW/2; tx(cx,16,d,{"text-anchor":"middle","font-size":12,"font-weight":600,fill:css("--ink2")});
   el("line",{x1:cx,x2:cx,y1:top-6,y2:H-4,stroke:css("--muted"),"stroke-width":1});});
 D.bars.forEach((b,i)=>{const y=top+i*rowH;
   tx(left-10,y+rowH/2+4,`${b.name} (${b.n})`,{"text-anchor":"end","font-size":12,fill:css("--ink")});
   b.days.forEach((d,j)=>{const cx=left+j*(colW+gapC)+colW/2, bh=12, by=y+(rowH-bh)/2;
     if(d.down){el("rect",{x:cx-1-d.down*sc,y:by,width:d.down*sc,height:bh,rx:3,fill:css("--down")});tx(cx-5-d.down*sc,by+10,d.down,{"text-anchor":"end","font-size":11,fill:css("--ink2")});}
     if(d.up){el("rect",{x:cx+1,y:by,width:d.up*sc,height:bh,rx:3,fill:css("--up")});tx(cx+5+d.up*sc,by+10,d.up,{"font-size":11,fill:css("--ink2")});}
     if(!d.up&&!d.down) tx(cx,by+10,"·",{"text-anchor":"middle","font-size":12,fill:css("--muted")});
     const hit=el("rect",{x:cx-half,y:y,width:2*half,height:rowH,fill:"transparent"});
     hit.addEventListener("mousemove",e=>showTip(e,`<b>${b.name}</b> · ${d.day}<br>${d.up} up: ${d.up_genes||"–"}<br>${d.down} down: ${d.down_genes||"–"}<br><span style="color:var(--muted)">${b.n} genes in pathway</span>`));
     hit.addEventListener("mouseleave",hideTip);});});
 document.getElementById("bars").appendChild(svg);
})();

// 4 heatmap
(function(){
 const lim=Math.ceil(Math.max(...D.genes.flatMap(g=>g.cells.map(c=>Math.abs(c.lfc))))*2)/2;
 const UP=[0xeb,0x68,0x34],DN=[0x2a,0x78,0xd6],MID=[0xf0,0xef,0xec];
 const mix=(a,b,t)=>`rgb(${a.map((v,i)=>Math.round(v+(b[i]-v)*t)).join(",")})`;
 const color=v=>v>=0?mix(MID,UP,Math.min(1,v/lim)):mix(MID,DN,Math.min(1,-v/lim));
 document.getElementById("hmleg").innerHTML=`<span>← lower in coculture</span><span style="display:inline-block;width:220px;height:12px;border-radius:2px;background:linear-gradient(90deg,${color(-lim)},${color(0)},${color(lim)})"></span><span>higher in coculture →</span><span>log2FC, coculture vs axenic, −${lim} … +${lim}</span><span>◯ = significant that day</span>`;
 const NS="http://www.w3.org/2000/svg", cw=64, ch=20, left=190, gl=330, top=28;
 const groups=[...new Set(D.genes.map(g=>g.group))];
 let y=top; const pos=[];
 groups.forEach(gr=>{D.genes.filter(g=>g.group===gr).forEach(g=>{pos.push([g,y]);y+=ch;}); y+=8;});
 const W=gl+left+D.days.length*cw+20, H=y+4;
 const svg=document.createElementNS(NS,"svg");svg.setAttribute("width",W);svg.setAttribute("height",H);
 const el=(n,a)=>{const e=document.createElementNS(NS,n);for(const k in a)e.setAttribute(k,a[k]);svg.appendChild(e);return e;};
 const tx=(x,y,s,a={})=>{const t=el("text",Object.assign({x,y},a));t.textContent=s;return t;};
 const x0=gl+left;
 D.days.forEach((d,j)=>tx(x0+j*cw+cw/2,16,d,{"text-anchor":"middle","font-size":12,"font-weight":600,fill:css("--ink2")}));
 let last=null;
 pos.forEach(([g,yy])=>{
   if(g.group!==last){tx(8,yy+14,g.group.length>52?g.group.slice(0,50)+"…":g.group,{"font-size":11,fill:css("--ink2")});
     el("line",{x1:8,x2:W-8,y1:yy-4,y2:yy-4,stroke:css("--rule"),"stroke-width":1});last=g.group;}
   tx(x0-8,yy+14,g.label,{"text-anchor":"end","font-size":12,fill:css("--ink")});
   g.cells.forEach((c,j)=>{const x=x0+j*cw;
     const r=el("rect",{x:x+1,y:yy+1,width:cw-2,height:ch-2,rx:3,fill:color(c.lfc)});
     if(c.status!=="not_significant") el("rect",{x:x+3,y:yy+3,width:cw-6,height:ch-6,rx:3,fill:"none",stroke:"#1f1f1d","stroke-width":1.5});
     tx(x+cw/2,yy+14,(c.lfc>0?"+":"")+c.lfc.toFixed(1),{"text-anchor":"middle","font-size":11,fill:Math.abs(c.lfc)/lim>0.6?"#ffffff":"#1f1f1d","pointer-events":"none"});
     r.addEventListener("mousemove",e=>showTip(e,`<b>${g.label}</b> · ${c.day}<br>${g.product}<br>log2FC ${c.lfc} · adj. p ${c.padj.toExponential(2)} · ${c.status.replace("_"," ")}<br><span style="color:var(--muted)">${g.sets}</span>`));
     r.addEventListener("mouseleave",hideTip);});});
 document.getElementById("heat").appendChild(svg);
})();

// 5 table
(function(){
 const tb=document.getElementById("tbl"); let key=null,asc=true;
 const cls=s=>s==="significant_up"?"s-up":s==="significant_down"?"s-down":"";
 function render(){const q=document.getElementById("q").value.toLowerCase();
  let rs=D.table.filter(r=>!q||(r.tag+" "+r.name+" "+r.sets).toLowerCase().includes(q));
  if(key!==null) rs=rs.slice().sort((a,b)=>{const va=key<0?(key===-1?a.tag:a.name):a.lfc[key],vb=key<0?(key===-1?b.tag:b.name):b.lfc[key];return (va>vb?1:va<vb?-1:0)*(asc?1:-1);});
  tb.innerHTML=`<tr><th data-k="-1">locus tag</th><th data-k="-2">gene</th>${D.days.map((d,j)=>`<th data-k="${j}">log2FC ${d}</th>`).join("")}<th>gene sets</th></tr>`+
   rs.map(r=>`<tr><td>${r.tag}</td><td>${r.name}</td>${r.lfc.map((v,j)=>`<td class="num ${cls(r.st[j])}">${v.toFixed(2)}${r.st[j]==="significant_up"?" ▲":r.st[j]==="significant_down"?" ▼":""}</td>`).join("")}<td style="color:var(--ink2)">${r.sets}</td></tr>`).join("");
  tb.querySelectorAll("th[data-k]").forEach(th=>th.onclick=()=>{const k=+th.dataset.k;asc=key===k?!asc:true;key=k;render();});}
 document.getElementById("q").oninput=render; render();
})();
</script></body></html>"""
(R / "q3_infographic.html").write_text(HTML.replace("__DATA__", json.dumps(payload)), encoding="utf-8")
print("wrote", R / "q3_infographic.html")
