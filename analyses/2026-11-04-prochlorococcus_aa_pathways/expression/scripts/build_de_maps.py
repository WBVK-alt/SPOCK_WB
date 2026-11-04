"""Pathway maps coloured by differential expression, one per condition set.

Same layout as grounding/scripts/build_pathway_map.py. Under each arrow, every gene gets a strip of small squares, one
per experiment of the condition set (numbered in the key): orange = significantly up at any timepoint, blue = down,
split = both, grey = measured not significant, white = not in the experiment's table, dot = strain lacks the gene.
Reads expression/data/de_overview_cells.csv + experiments.csv. Writes expression/reports/de_maps.html and
expression/figures/de_map_<set>.svg.
"""
from pathlib import Path
import html
import sys
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
D, R, F = HERE / "data", HERE / "reports", HERE / "figures"
sys.path.insert(0, str(HERE.parent / "grounding" / "scripts"))
from pathway_steps import NODES, DEG_NODES, ROUTES, CX, CY  # noqa: E402

P = "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_"
SETS = [
    ("nitrogen", "Nitrogen depletion / limitation", [
        (P + "rnaseq_axenic", "MED4 low-N starvation, RNA, axenic"),
        (P + "rnaseq_coculture", "MED4 low-N starvation, RNA, coculture"),
        (P + "proteomics_axenic", "MED4 low-N starvation, protein, axenic"),
        (P + "proteomics_coculture", "MED4 low-N starvation, protein, coculture"),
        ("10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq", "MED4 N-depleted vs N-replete, RNA ◆"),
        ("10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray", "MED4 N deprivation, microarray"),
        ("10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray", "MIT9313 N deprivation, microarray")]),
    ("lightdark", "Extended darkness / extended light", [
        ("10.64898/2026.04.15.718746_extended_darkness_med4", "MED4 extended darkness (vs diel light)"),
        ("10.1128/mSystems.00040-18_darkness_extended_darkness_natl2a_rnaseq_axenic", "NATL2A extended darkness, axenic ◆"),
        ("10.1128/mSystems.00040-18_darkness_extended_darkness_natl2a_rnaseq_coculture", "NATL2A extended darkness, coculture ◆"),
        ("10.1371/journal.pone.0165375_light_stress_constant_dark_med4_rnaseq_dark", "MED4 constant dark vs constant light ◆"),
        ("10.64898/2026.04.15.718746_extended_light_med4", "MED4 extended LIGHT (vs diel dark)")]),
    ("mitomycin", "Mitomycin C (DNA damage)", [
        ("10.1016/j.cell.2022.12.006_mmc_dna_damage_mit0604", "MIT0604 mitomycin C 2 h vs untreated")]),
]

cells = pd.read_csv(D / "de_overview_cells.csv", encoding="utf-8")
st = {(c.gene, c.experiment_id): (c.status, c.detail if isinstance(c.detail, str) else "") for c in cells.itertuples()}
exp_ids = set(pd.read_csv(D / "experiments.csv", encoding="utf-8").experiment_id)
for _, _, ex in SETS:
    for e, _ in ex:
        assert e in exp_ids, e

INK, INK2, MUTED, RULE, SURF, CARD = "#1f1f1d", "#52514e", "#8a8984", "#d6d5d0", "#fcfcfb", "#f1f0ec"
UP, DN, NS = "#eb6834", "#2a78d6", "#d9d8d3"
W, H = 1780, 1290
esc = html.escape


def box(label, kind):
    return 7.2 * len(label) + (22 if kind == "aa" else 18), (28 if kind == "aa" else 24)


def node_xy(n):
    return NODES[n][0] if n in NODES else DEG_NODES[n]


def trim(p, q, name, pad=5):
    w, h = box(name, NODES[name][1]) if name in NODES else (7.0 * len(name) + 8, 18)
    dx, dy = q[0] - p[0], q[1] - p[1]
    if dx == 0 and dy == 0:
        return p
    t = min((w / 2 + pad) / abs(dx) if dx else 1e9, (h / 2 + pad) / abs(dy) if dy else 1e9, 1)
    return (p[0] + dx * t, p[1] + dy * t)


def draw(key, title, exps):
    ne = len(exps)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="system-ui,Segoe UI,sans-serif">',
         f'<rect width="{W}" height="{H}" fill="{SURF}"/><defs>']
    add = o.append
    for nm, col in (("a-c", MUTED), ("a-b", INK), ("a-d", INK2)):
        add(f'<marker id="{key}{nm}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{col}"/></marker>')
    add('</defs>')
    add(f'<text x="30" y="34" font-size="20" font-weight="700" fill="{INK}">{esc(title)}</text>')
    add(f'<circle cx="{CX}" cy="{CY}" r="185" fill="none" stroke="#e4e3df" stroke-width="18" stroke-opacity=".45"/>')
    add(f'<text x="{CX}" y="{CY}" text-anchor="middle" font-size="16" font-weight="700" fill="{INK2}">TCA cycle</text>')
    mids = []
    for r in ROUTES:
        p, q, c = node_xy(r["src"]), node_xy(r["dst"]), r.get("via")
        a, b = trim(p, c or q, r["src"]), trim(q, c or p, r["dst"], pad=7)
        col, mk = {"central": (MUTED, "a-c"), "bio": (INK, "a-b"), "deg": (INK2, "a-d")}[r["kind"]]
        dash = ' stroke-dasharray="6 5"' if r["kind"] == "deg" else ""
        d = f"M{a[0]:.0f},{a[1]:.0f} Q{c[0]},{c[1]} {b[0]:.0f},{b[1]:.0f}" if c else f"M{a[0]:.0f},{a[1]:.0f} L{b[0]:.0f},{b[1]:.0f}"
        add(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="1.6"{dash} stroke-opacity=".7" marker-end="url(#{key}{mk})"/>')
        mids.append((r, (0.25 * a[0] + 0.5 * c[0] + 0.25 * b[0], 0.25 * a[1] + 0.5 * c[1] + 0.25 * b[1]) if c else ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)))
    sq, gapx = 7, 6
    for r, (mx, my) in mids:
        genes = [g for g, _ in r["genes"]]
        gw = [max(6.2 * len(g), ne * sq) for g in genes]
        lines, cur, cw_ = [], [], 0
        for g, w in zip(genes, gw):
            if cur and cw_ + w > 150:
                lines.append(cur); cur, cw_ = [], 0
            cur.append((g, w)); cw_ += w + gapx
        lines.append(cur)
        lh = 22
        y0 = my - (len(lines) * lh) / 2
        tot = len(lines) * lh + 4
        add(f'<rect x="{mx - 80}" y="{y0 - 3}" width="160" height="{tot}" rx="5" fill="{SURF}" fill-opacity=".8"/>')
        for li, ln in enumerate(lines):
            tw = sum(w for _, w in ln) + gapx * (len(ln) - 1)
            x = mx - tw / 2
            yy = y0 + li * lh
            for g, w in ln:
                cx = x + w / 2
                add(f'<text x="{cx}" y="{yy + 9}" text-anchor="middle" font-size="10" font-weight="600" fill="{INK}">{esc(g)}</text>')
                sx = cx - ne * sq / 2
                for k, (eid, elab) in enumerate(exps):
                    s, det = st.get((g, eid), ("no_data", ""))
                    tip = esc(f"{g} · [{k + 1}] {elab}: {s}" + (f" — {det}" if det else ""))
                    X, Y = sx + k * sq, yy + 12
                    if s == "up":
                        add(f'<rect class="c" data-t="{tip}" x="{X}" y="{Y}" width="{sq - 1}" height="{sq}" fill="{UP}"/>')
                    elif s == "down":
                        add(f'<rect class="c" data-t="{tip}" x="{X}" y="{Y}" width="{sq - 1}" height="{sq}" fill="{DN}"/>')
                    elif s == "both":
                        add(f'<g class="c" data-t="{tip}"><path d="M{X},{Y} h{sq - 1} L{X},{Y + sq} z" fill="{UP}"/><path d="M{X + sq - 1},{Y} v{sq} h{-(sq - 1)} z" fill="{DN}"/></g>')
                    elif s == "ns":
                        add(f'<rect class="c" data-t="{tip}" x="{X}" y="{Y}" width="{sq - 1}" height="{sq}" fill="{NS}"/>')
                    elif s == "not_measured":
                        add(f'<rect class="c" data-t="{tip}" x="{X + .5}" y="{Y + .5}" width="{sq - 2}" height="{sq - 1}" fill="#fff" stroke="{MUTED}" stroke-width=".7"/>')
                    else:  # no_gene / no_data: strain lacks the gene
                        add(f'<circle class="c" data-t="{tip}" cx="{X + sq / 2}" cy="{Y + sq / 2}" r="1.3" fill="{MUTED}"/>')
                x += w + gapx
    for name, ((x, y), kind) in NODES.items():
        w, h = box(name, kind)
        if kind == "aa":
            add(f'<rect x="{x - w / 2}" y="{y - h / 2}" width="{w}" height="{h}" rx="{h / 2}" fill="{INK}"/><text x="{x}" y="{y + 5}" text-anchor="middle" font-size="14" font-weight="700" fill="#fff">{esc(name)}</text>')
        else:
            add(f'<rect x="{x - w / 2}" y="{y - h / 2}" width="{w}" height="{h}" rx="6" fill="#fff" stroke="{MUTED if kind == "central" else RULE}" stroke-width="1.3"/>'
                f'<text x="{x}" y="{y + 4.5}" text-anchor="middle" font-size="12.5" font-weight="{600 if kind == "central" else 400}" fill="{INK if kind == "central" else INK2}">{esc(name)}</text>')
    for name, (x, y) in DEG_NODES.items():
        add(f'<text x="{x}" y="{y + 4}" text-anchor="middle" font-size="11.5" font-style="italic" fill="{INK2}">{esc(name)}</text>')
    # key
    lx, ly = 30, 1085
    add(f'<text x="{lx}" y="{ly}" font-size="13" font-weight="700" fill="{INK}">Squares under each gene, left to right = experiments:</text>')
    for k, (_, elab) in enumerate(exps):
        add(f'<text x="{lx}" y="{ly + 18 + 15 * k}" font-size="11.5" fill="{INK2}">[{k + 1}] {esc(elab)}</text>')
    ly2 = ly + 30 + 15 * ne
    for i, (txt, kind) in enumerate([("up", UP), ("down", DN), ("up + down (timepoints)", "both"), ("not significant", NS),
                                     ("not in table", "nm"), ("strain lacks gene", "dot")]):
        x = lx + i * 120 if i < 3 else lx + (i - 3) * 120
        y = ly2 + (0 if i < 3 else 18)
        if kind == "both":
            add(f'<path d="M{x},{y} h10 L{x},{y + 10} z" fill="{UP}"/><path d="M{x + 10},{y} v10 h-10 z" fill="{DN}"/>')
        elif kind == "nm":
            add(f'<rect x="{x}" y="{y}" width="10" height="10" fill="#fff" stroke="{MUTED}"/>')
        elif kind == "dot":
            add(f'<circle cx="{x + 5}" cy="{y + 5}" r="1.6" fill="{MUTED}"/>')
        else:
            add(f'<rect x="{x}" y="{y}" width="10" height="10" fill="{kind}"/>')
        add(f'<text x="{x + 15}" y="{y + 9}" font-size="11" fill="{INK2}">{txt}</text>')
    add(f'<text x="{lx}" y="{ly2 + 44}" font-size="11" fill="{INK2}">Significance = each paper\'s own call; any timepoint, any paralog. ◆ = partial source table (only significant genes listed).</text>')
    add('</svg>')
    return "\n".join(o)


parts = []
for key, title, exps in SETS:
    svg = draw(key, title, exps)
    (F / f"de_map_{key}.svg").write_text(svg, encoding="utf-8")
    parts.append((title, svg))
    print(key, {s: sum(1 for r in ROUTES for g, _ in r["genes"] for e, _ in exps if st.get((g, e), ("no_data",))[0] == s)
                for s in ("up", "down", "both", "ns", "not_measured")})
tabs = "".join(f'<button onclick="show({i})" id="b{i}">{html.escape(t)}</button>' for i, (t, _) in enumerate(parts))
body = "".join(f'<div class="fig" id="f{i}" style="display:{"block" if i == 0 else "none"}">{s}</div>' for i, (_, s) in enumerate(parts))
page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>Prochlorococcus amino-acid pathways — expression by condition</title>
<style>body{{background:{SURF};margin:16px 24px;font:14px system-ui,Segoe UI,sans-serif}}.fig{{overflow:auto;border:1px solid {RULE};border-radius:10px}}
button{{font:inherit;padding:6px 12px;margin:0 6px 10px 0;border:1px solid {RULE};border-radius:6px;background:#fff;cursor:pointer}}button.on{{background:{INK};color:#fff}}
#tip{{position:fixed;pointer-events:none;background:#fff;border:1px solid {RULE};border-radius:6px;padding:6px 9px;font-size:12px;max-width:560px;box-shadow:0 4px 14px rgba(0,0,0,.15);display:none}}
.c:hover{{outline:1.5px solid {INK}}}</style></head><body>
<h1 style="font-size:20px;margin:0 0 8px">How amino-acid pathway genes respond, mapped onto the pathways</h1>
<div>{tabs}</div>{body}<div id="tip"></div>
<script>function show(i){{document.querySelectorAll(".fig").forEach((f,j)=>f.style.display=j===i?"block":"none");document.querySelectorAll("button").forEach((b,j)=>b.className=j===i?"on":"")}}show(0);
const t=document.getElementById("tip");document.querySelectorAll(".c").forEach(r=>{{
r.addEventListener("mousemove",e=>{{t.style.display="block";t.textContent=r.dataset.t;t.style.left=Math.min(e.clientX+12,innerWidth-580)+"px";t.style.top=Math.min(e.clientY+12,innerHeight-120)+"px"}});
r.addEventListener("mouseleave",()=>t.style.display="none")}});</script></body></html>"""
(R / "de_maps.html").write_text(page, encoding="utf-8")
print("wrote", R / "de_maps.html")
