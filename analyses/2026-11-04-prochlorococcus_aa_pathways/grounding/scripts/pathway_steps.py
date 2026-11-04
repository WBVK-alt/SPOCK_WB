"""Step definitions for the amino-acid pathway map.

The pathway topology (which reaction links which metabolite, which enzyme catalyses it, and the KEGG KO
for that enzyme) is textbook biochemistry [interpretation], not KG data. The KG supplies only which
genomes carry each KO (pull_ko_presence.py). KO names are checked against the KG (ko_names.csv) so a
wrong KO assignment shows up as a name mismatch.

ROUTES: one arrow on the map. kind = central | bio | deg.
  genes: list of (label, [KO alternatives]) — the gene counts as present in a strain if ANY listed KO is.
  A route counts as complete in a strain when ALL its genes are present.
NODES: metabolite label -> (x, y) on the map canvas; kind = central | aa | other.
"""
import math

K = lambda *ids: [f"kegg.orthology:{i}" for i in ids]
# PSTA (3-phosphoserine transaminase, MED4 PMM0035) carries KO K00436 (hoxH) in the KG, so it is counted by its
# Cyanorak ortholog group instead (see gaps_and_friction.md).
OG_PSTA = "cyanorak:CK_00000338"

# --- TCA ring geometry -------------------------------------------------------
CX, CY, RR = 800, 990, 185
_ring = ["Citrate", "Isocitrate", "2-Oxoglutarate", "Succinyl-CoA", "Succinate", "Fumarate", "Malate", "Oxaloacetate"]
TCA = {n: (round(CX + RR * math.sin(math.radians(22.5 + 45 * i))), round(CY - RR * math.cos(math.radians(22.5 + 45 * i))))
       for i, n in enumerate(_ring)}

NODES = {
    # central carbon
    "CO₂ (Calvin cycle)": ((800, 70), "central"),
    "Sugar phosphates": ((800, 170), "central"),
    "Ribose-5-P": ((1060, 120), "central"),
    "Erythrose-4-P": ((1060, 220), "central"),
    "PRPP": ((1290, 120), "other"),
    "3-Phosphoglycerate": ((800, 290), "central"),
    "PEP": ((800, 450), "central"),
    "Pyruvate": ((800, 600), "central"),
    "Acetyl-CoA": ((800, 760), "central"),
    **{n: (xy, "central") for n, xy in TCA.items()},
    # Ser family
    "Ser": ((560, 290), "aa"), "Gly": ((330, 250), "aa"), "Cys": ((330, 340), "aa"),
    # aromatic family
    "Chorismate": ((1180, 450), "other"),
    "Phe": ((1450, 400), "aa"), "Tyr": ((1450, 470), "aa"), "Trp": ((1450, 540), "aa"),
    "His": ((1520, 120), "aa"),
    # pyruvate family
    "Ala": ((1060, 600), "aa"),
    "2-Oxoisovalerate": ((1060, 690), "other"),
    "Val": ((1290, 650), "aa"), "Leu": ((1290, 730), "aa"),
    # Glu family
    "Glu": ((1180, 1080), "aa"), "Gln": ((1420, 1010), "aa"), "Pro": ((1420, 1090), "aa"), "Arg": ((1420, 1170), "aa"),
    # Asp family
    "Asp": ((470, 860), "aa"), "Asn": ((250, 930), "aa"),
    "Aspartate semialdehyde": ((470, 740), "other"), "Lys": ((250, 740), "aa"),
    "Homoserine": ((470, 610), "other"), "Thr": ((250, 610), "aa"), "Met": ((470, 480), "aa"),
    "Ile": ((250, 480), "aa"),
}

# degradation end-products: small grey labels, positioned next to the amino acid
DEG_NODES = {
    "CO₂ + NH₃ + CH₂-THF": (110, 250),
    "→ pyruvate (from Ser)": (600, 395),
    "→ Gly + acetaldehyde": (95, 560),
    "→ fumarate (from Asp)": (470, 980),
    "→ Asp (from Asn)": (95, 990),
    "→ Glu (from Gln)": (1600, 1010),
    "→ putrescine (polyamines)": (1650, 1125),
    "→ ornithine + urea": (1650, 1225),
    "→ urocanate": (1700, 60),
    "→ indole": (1650, 560),
    "→ homogentisate": (1650, 470),
    "→ branched-chain acyl-CoA": (1580, 770),
    "→ cadaverine": (95, 790),
    "→ methanethiol + 2-oxobutanoate": (660, 470),
    "→ Glu (from Pro)": (1620, 1070),
    "→ Tyr (from Phe)": (1650, 400),
    "→ pyruvate (from Ala)": (1060, 530),
}

ROUTES = [
    # ---------------- central carbon -----------------------------------------------------------
    dict(id="rubisco", src="CO₂ (Calvin cycle)", dst="3-Phosphoglycerate", kind="central", via=(640, 180),
         genes=[("rbcL", K("K01601")), ("rbcS", K("K01602"))]),
    dict(id="tkt", src="Sugar phosphates", dst="Ribose-5-P", kind="central", genes=[("tktA", K("K00615"))]),
    dict(id="tkt2", src="Sugar phosphates", dst="Erythrose-4-P", kind="central", genes=[("tktA", K("K00615"))]),
    dict(id="prs", src="Ribose-5-P", dst="PRPP", kind="central", genes=[("prs", K("K00948"))]),
    dict(id="pgk", src="3-Phosphoglycerate", dst="Sugar phosphates", kind="central",
         genes=[("pgk", K("K00927")), ("gap2", K("K00134", "K00150"))]),
    dict(id="lowglyc", src="3-Phosphoglycerate", dst="PEP", kind="central",
         genes=[("gpm", K("K15633", "K01834")), ("eno", K("K01689"))]),
    dict(id="pyk", src="PEP", dst="Pyruvate", kind="central", genes=[("pyk", K("K00873"))]),
    dict(id="pdh", src="Pyruvate", dst="Acetyl-CoA", kind="central",
         genes=[("pdhA", K("K00161")), ("pdhB", K("K00162")), ("aceF", K("K00627")), ("lpd", K("K00382"))]),
    dict(id="gltA", src="Acetyl-CoA", dst="Citrate", kind="central", genes=[("gltA", K("K01647"))]),
    dict(id="acn", src="Citrate", dst="Isocitrate", kind="central", genes=[("acnB", K("K01682", "K01681"))]),
    dict(id="icd", src="Isocitrate", dst="2-Oxoglutarate", kind="central", genes=[("icd", K("K00031"))]),
    dict(id="ogdh", src="2-Oxoglutarate", dst="Succinyl-CoA", kind="central",
         genes=[("sucA", K("K00164")), ("sucB", K("K00658"))]),
    dict(id="scs", src="Succinyl-CoA", dst="Succinate", kind="central", genes=[("sucC", K("K01903")), ("sucD", K("K01902"))]),
    dict(id="sdh", src="Succinate", dst="Fumarate", kind="central", genes=[("sdhA", K("K00239")), ("sdhB", K("K00240"))]),
    dict(id="fum", src="Fumarate", dst="Malate", kind="central", genes=[("fumC", K("K01679", "K01676"))]),
    dict(id="mdh", src="Malate", dst="Oxaloacetate", kind="central", genes=[("mdh/mqo", K("K00024", "K00116"))]),
    dict(id="ppc", src="PEP", dst="Oxaloacetate", kind="central", via=(640, 600), genes=[("ppc", K("K01595"))]),
    # ---------------- biosynthesis -------------------------------------------------------------
    dict(id="serine", src="3-Phosphoglycerate", dst="Ser", kind="bio",
         genes=[("serA", K("K00058")), ("serC/PSTA", K("K00831") + [OG_PSTA]), ("serB/psp", K("K01079", "K22305"))]),
    dict(id="glyA", src="Ser", dst="Gly", kind="bio", genes=[("glyA", K("K00600"))]),
    dict(id="cys", src="Ser", dst="Cys", kind="bio", genes=[("cysE", K("K00640")), ("cysK", K("K01738", "K12339"))]),
    dict(id="shikimate", src="Erythrose-4-P", dst="Chorismate", kind="bio", note="+ PEP",
         genes=[("aroG", K("K01626", "K03856")), ("aroB", K("K01735")), ("aroQ", K("K03786", "K03785")),
                ("aroE", K("K00014")), ("aroK", K("K00891")), ("aroA", K("K00800")), ("aroC", K("K01736"))]),
    dict(id="phe", src="Chorismate", dst="Phe", kind="bio",
         genes=[("CM", K("K06208", "K04516", "K14170", "K04092", "K14187", "K04093")), ("pheA", K("K04518", "K05359", "K14170"))]),
    dict(id="tyr", src="Chorismate", dst="Tyr", kind="bio",
         genes=[("CM", K("K06208", "K04516", "K14170", "K04092", "K14187", "K04093")), ("tyrA", K("K04517", "K15226", "K00220", "K14187"))]),
    dict(id="trp", src="Chorismate", dst="Trp", kind="bio",
         genes=[("trpE", K("K01657")), ("trpG", K("K01658")), ("trpD", K("K00766")), ("trpC", K("K01609", "K13498")),
                ("trpF", K("K01817")), ("trpA", K("K01695")), ("trpB", K("K01696"))]),
    dict(id="his", src="PRPP", dst="His", kind="bio",
         genes=[("hisG", K("K00765")), ("hisIE", K("K01496", "K01523", "K11755")), ("hisA", K("K01814")),
                ("hisH", K("K02501")), ("hisF", K("K02500")), ("hisB", K("K01693", "K01089")), ("hisC", K("K00817")),
                ("hisN", K("K05602", "K04486", "K18649")), ("hisD", K("K00013"))], note="no annotated histidinol-phosphatase"),
    dict(id="ala", src="Pyruvate", dst="Ala", kind="bio", genes=[("ald/alaA/AGXT", K("K00259", "K14260", "K14261", "K00830"))]),
    dict(id="kiv", src="Pyruvate", dst="2-Oxoisovalerate", kind="bio", via=(930, 690),
         genes=[("ilvB", K("K01652")), ("ilvN", K("K01653")), ("ilvC", K("K00053")), ("ilvD", K("K01687"))]),
    dict(id="val", src="2-Oxoisovalerate", dst="Val", kind="bio", genes=[("ilvE", K("K00826"))]),
    dict(id="leu", src="2-Oxoisovalerate", dst="Leu", kind="bio",
         genes=[("leuA", K("K01649")), ("leuC", K("K01703")), ("leuD", K("K01704")), ("leuB", K("K00052")), ("ilvE", K("K00826"))]),
    dict(id="gogat", src="2-Oxoglutarate", dst="Glu", kind="bio", note="+ Gln (GS-GOGAT)",
         genes=[("glsF/gltB", K("K00284", "K00265"))]),
    dict(id="gln", src="Glu", dst="Gln", kind="bio", genes=[("glnA", K("K01915", "K01949"))]),
    dict(id="pro", src="Glu", dst="Pro", kind="bio", genes=[("proB", K("K00931")), ("proA", K("K00147")), ("proC", K("K00286"))]),
    dict(id="arg", src="Glu", dst="Arg", kind="bio",
         genes=[("argJ", K("K00620")), ("argB", K("K00930")), ("argC", K("K00145")), ("argD", K("K00818", "K00821")),
                ("carA", K("K01956")), ("carB", K("K01955")), ("argF", K("K00611")), ("argG", K("K01940")), ("argH", K("K01755"))]),
    dict(id="asp", src="Oxaloacetate", dst="Asp", kind="bio", genes=[("aspC", K("K00812", "K00813", "K11358"))]),
    dict(id="asn", src="Asp", dst="Asn", kind="bio", genes=[("asnB", K("K01953"))]),
    dict(id="asa", src="Asp", dst="Aspartate semialdehyde", kind="bio",
         genes=[("lysC", K("K00928", "K12524", "K12525", "K12526")), ("asd", K("K00133"))]),
    dict(id="lys", src="Aspartate semialdehyde", dst="Lys", kind="bio",
         genes=[("dapA", K("K01714")), ("dapB", K("K00215")), ("dapL", K("K10206")), ("dapF", K("K01778")), ("lysA", K("K01586"))]),
    dict(id="hom", src="Aspartate semialdehyde", dst="Homoserine", kind="bio", genes=[("hom", K("K00003", "K12524", "K12525"))]),
    dict(id="thr", src="Homoserine", dst="Thr", kind="bio", genes=[("thrB", K("K00872", "K02204")), ("thrC", K("K01733"))]),
    dict(id="met", src="Homoserine", dst="Met", kind="bio",
         genes=[("metX/metA", K("K00641", "K00651")), ("metY/metB", K("K01740", "K01739")), ("metE/metH", K("K00549", "K00548"))]),
    dict(id="ile", src="Thr", dst="Ile", kind="bio", note="+ pyruvate",
         genes=[("ilvA", K("K01754", "K17989")), ("ilvB", K("K01652")), ("ilvN", K("K01653")), ("ilvC", K("K00053")),
                ("ilvD", K("K01687")), ("ilvE", K("K00826"))]),
    # ---------------- degradation / exits ------------------------------------------------------
    dict(id="gcv", src="Gly", dst="CO₂ + NH₃ + CH₂-THF", kind="deg",
         genes=[("gcvT", K("K00605")), ("gcvP", K("K00281", "K00282")), ("gcvH", K("K02437"))]),
    dict(id="sdaA", src="Ser", dst="→ pyruvate (from Ser)", kind="deg", genes=[("sdaA", K("K01752"))]),
    dict(id="dadA", src="Ala", dst="→ pyruvate (from Ala)", kind="deg", genes=[("dadA", K("K00285"))]),
    dict(id="thrdeg", src="Thr", dst="→ Gly + acetaldehyde", kind="deg", genes=[("ltaE/tdh", K("K01620", "K00060"))]),
    dict(id="aspA", src="Asp", dst="→ fumarate (from Asp)", kind="deg", genes=[("aspA", K("K01744"))]),
    dict(id="ans", src="Asn", dst="→ Asp (from Asn)", kind="deg", genes=[("ansA/iaaA", K("K01424", "K13051"))]),
    dict(id="glsA", src="Gln", dst="→ Glu (from Gln)", kind="deg", genes=[("glsA", K("K01425"))]),
    dict(id="gdh", src="Glu", dst="2-Oxoglutarate", kind="deg", via=(1080, 1190), genes=[("gdhA", K("K00262", "K00261", "K15371"))]),
    dict(id="putA", src="Pro", dst="→ Glu (from Pro)", kind="deg", genes=[("putA", K("K13821", "K00318"))]),
    dict(id="spe", src="Arg", dst="→ putrescine (polyamines)", kind="deg", genes=[("speA", K("K01585")), ("speB", K("K01480"))]),
    dict(id="rocF", src="Arg", dst="→ ornithine + urea", kind="deg", genes=[("rocF", K("K01476"))]),
    dict(id="hutH", src="His", dst="→ urocanate", kind="deg", genes=[("hutH", K("K01745"))]),
    dict(id="tnaA", src="Trp", dst="→ indole", kind="deg", genes=[("tnaA", K("K01667"))]),
    dict(id="phhA", src="Phe", dst="→ Tyr (from Phe)", kind="deg", genes=[("phhA", K("K00500"))]),
    dict(id="hpd", src="Tyr", dst="→ homogentisate", kind="deg", genes=[("hpd", K("K00457"))]),
    dict(id="bkd", src="Leu", dst="→ branched-chain acyl-CoA", kind="deg",
         genes=[("bkdA1", K("K00166")), ("bkdA2", K("K00167")), ("bkdB", K("K09699"))]),
    dict(id="ldc", src="Lys", dst="→ cadaverine", kind="deg", genes=[("ldcC", K("K01582"))]),
    dict(id="mdeA", src="Met", dst="→ methanethiol + 2-oxobutanoate", kind="deg", genes=[("mdeA", K("K01761"))]),
]


def all_kos():
    return sorted({ko for r in ROUTES for _, alts in r["genes"] for ko in alts if ko.startswith("kegg.")})


def all_groups():
    return sorted({g for r in ROUTES for _, alts in r["genes"] for g in alts if not g.startswith("kegg.")})
