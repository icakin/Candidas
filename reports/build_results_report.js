// Results report — Candida thermal energetics. Builds a .docx with Figs 1 & 2.
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  ImageRun, PageBreak
} = require("docx");

const FIGDIR = "/sessions/sweet-ecstatic-davinci/mnt/Candidas/figures";
const OUT = "/sessions/sweet-ecstatic-davinci/mnt/Candidas/reports/Results_report.docx";

// ---- helpers ----------------------------------------------------------------
const body = (runs, opts = {}) =>
  new Paragraph({ spacing: { after: 160, line: 276 }, alignment: AlignmentType.JUSTIFIED, children: runs, ...opts });
const t  = (s, o = {}) => new TextRun({ text: s, size: 22, font: "Calibri", ...o });   // 11 pt
const it = (s) => t(s, { italics: true });
const b  = (s) => t(s, { bold: true });

const H1 = (s) => new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 240, after: 120 },
  children: [new TextRun({ text: s, size: 28, bold: true, font: "Calibri", color: "1a1a1a" })] });
const H2 = (s) => new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 200, after: 80 },
  children: [new TextRun({ text: s, size: 24, bold: true, font: "Calibri", color: "333333" })] });

function figure(file, widthIn) {
  const data = fs.readFileSync(path.join(FIGDIR, file));
  // native px from the file is 600 dpi; scale to widthIn inches (96 dpi doc units)
  const wPx = widthIn * 96;
  // read PNG dims
  const buf = data;
  const w = buf.readUInt32BE(16), h = buf.readUInt32BE(20);
  const hPx = wPx * (h / w);
  return new Paragraph({
    alignment: AlignmentType.CENTER, spacing: { before: 120, after: 60 },
    children: [new ImageRun({ type: "png", data, transformation: { width: Math.round(wPx), height: Math.round(hPx) } })]
  });
}

const legend = (title, runs) => new Paragraph({
  spacing: { after: 200, line: 240 }, alignment: AlignmentType.JUSTIFIED,
  children: [new TextRun({ text: title + " ", bold: true, size: 20, font: "Calibri" }), ...runs]
});
const lt = (s, o = {}) => new TextRun({ text: s, size: 20, font: "Calibri", ...o });   // 10 pt legend
const li = (s) => lt(s, { italics: true });

// Species name: genus (roman) + epithet (italic). Only the epithet is italic,
// e.g. "C. auris" -> C. roman, auris italic. Returns an array of runs.
const sp = (genus, epithet, tail = " ") => [t(genus + " "), it(epithet + tail)];
const spL = (genus, epithet, tail = " ") => [lt(genus + " "), li(epithet + tail)]; // legend size

// ---- content ----------------------------------------------------------------
const children = [];

children.push(new Paragraph({
  spacing: { after: 60 },
  children: [new TextRun({ text: "The metabolic cost of the febrile host across the priority Candida",
    bold: true, size: 32, font: "Calibri" })] }));
children.push(new Paragraph({ spacing: { after: 240 },
  children: [new TextRun({ text: "Results", italics: true, size: 24, color: "666666", font: "Calibri" })] }));

children.push(new Paragraph({ spacing: { after: 200 }, alignment: AlignmentType.JUSTIFIED, children: [
  t("Dissolved-oxygen respirometry was used to fit thermal performance curves for four "),
  ...sp("Candida", "auris"), t("clades (I-IV; three isolates each) and "),
  ...sp("C.", "parapsilosis"),
  t("(three isolates), across 22-44 °C (12 temperatures). Growth was modelled with a Sharpe-Schoolfield function and per-cell respiration with a Boltzmann-Arrhenius function, both fitted hierarchically (isolates nested within clades) in a Bayesian framework. Unless stated, values are posterior medians with 95% credible intervals (CrI).")
]}));

children.push(H2("Methods (summary)"));
children.push(body([
  t("O"),
  t("2", { subScript: true }),
  t(" concentration was recorded continuously with PreSens SensorDish readers across the temperature series (22-44 °C, 12 temperatures; one plate per group × temperature, three isolates × five replicate wells). Specific growth rate and per-cell respiration rate were recovered simultaneously from each O"),
  t("2", { subScript: true }),
  t(" time series using the mechanistic oxygen-dynamics model of Cakin "),
  it("et al. "),
  t("(2026), which couples exponential population growth with biomass-proportional oxygen consumption, O"),
  t("2", { subScript: true }),
  t("(t) = O"),
  t("2,0", { subScript: true }),
  t(" + (K/r)(1 − e"),
  t("rt", { superScript: true }),
  t("), where r is the specific growth rate and K the volumetric O"),
  t("2", { subScript: true }),
  t(" consumption rate at the window start. Per-cell respiration, carbon-use efficiency and derived quantities were computed as in that study.")
]));

children.push(body([
  b("Units. "),
  t("The oxygen-dynamics model recovers the specific growth rate r (h"),
  t("−1", { superScript: true }),
  t(") from the exponential term. A per-cell carbon assimilation flux was then formed as G = r q, where q is the per-cell carbon quota (a temperature-independent constant per taxon); the per-cell respiratory carbon flux R was obtained from the volumetric O"),
  t("2", { subScript: true }),
  t(" consumption rate divided by the biomass integral, converted to carbon. G and R are therefore both per-cell carbon fluxes (fg C cell"),
  t("−1", { superScript: true }),
  t(" h"),
  t("−1", { superScript: true }),
  t("), and carbon-use efficiency CUE = G/(G + R) is the dimensionless ratio of these two matched quantities. Because q is independent of temperature, the growth activation energy E"),
  t("G", { subScript: true }),
  t(" equals the activation energy of r (q shifts the intercept, not the Boltzmann slope) and is directly comparable to the respiration activation energy E"),
  t("R", { subScript: true }),
  t(". Growth thermal performance curves (Fig. 1a) and E"),
  t("G", { subScript: true }),
  t(" thus derive from a rate constant; the r-to-carbon conversion enters only where CUE and downstream ratios are formed.")
]));

// ---- Section 1 ----
children.push(H1("Growth and respiration decouple, placing the carbon-economy optimum below body temperature"));

children.push(body([
  t("Growth rate rose with temperature to a clade-level optimum near 33-36 °C and then declined, following the expected unimodal thermal performance curve (Fig. 1a). Per-cell respiration showed no such turnover: it increased monotonically across the entire measured range, and leave-one-out cross-validation favoured the monotonic Arrhenius model over a unimodal Sharpe-Schoolfield alternative (Fig. 1b). Growth and respiration are therefore governed by different thermal sensitivities.")
]));

children.push(body([
  t("This difference is quantified by the activation energies (Fig. 1c). In every taxon, the activation energy of growth exceeded that of respiration: E"),
  t("G", { subScript: true }),
  t(" ranged 0.81-1.14 eV and E"),
  t("R", { subScript: true }),
  t(" 0.30-0.52 eV, with the growth-respiration difference credibly greater than zero in all five taxa (0.39-0.73 eV; 95% CrI excluding zero throughout). Growth is thus roughly twice as temperature-sensitive as respiration.")
]));

children.push(body([
  t("Because respiration continues to rise while growth saturates and falls, carbon-use efficiency (CUE = G/(G+R)) is itself a unimodal function of temperature that peaks early and then collapses (Fig. 1d). The temperature of this economic optimum, T"),
  t("opt", { subScript: true }),
  t("(CUE) lay below human body temperature in every taxon: 31.3-34.1 °C, with posterior probability P(T"),
  t("opt", { subScript: true }),
  t("(CUE) < 37 °C) ≥ 0.999 for all five (Fig. 1e). The economic optimum sat 1.0-2.6 °C below the growth optimum and 2.9-5.7 °C below 37 °C. The mammalian host temperature therefore lies beyond the point at which warming still improves the carbon economy of any species tested.")
]));

children.push(figure("FIG1_decoupling.png", 6.4));
children.push(legend("Figure 1. Growth and respiration decouple with temperature, placing every carbon-economy optimum below body temperature.", [
  lt("("), lt("a", { bold: true }), lt(") Specific growth rate "), li("r"),
  lt(" (h⁻¹) and ("), lt("b", { bold: true }),
  lt(") per-cell respiration carbon flux (fg C cell⁻¹ h⁻¹) versus temperature (points, isolate-level observations; lines, clade-level posterior medians; bands, 95% CrI). Growth turns over; respiration does not (Arrhenius favoured over Sharpe-Schoolfield by leave-one-out cross-validation). ("),
  lt("c", { bold: true }),
  lt(") Activation energies of growth (circles) and respiration (triangles); dashed line, the 0.65 eV metabolic-theory value. Growth is credibly more temperature-sensitive than respiration in all five taxa. ("),
  lt("d", { bold: true }),
  lt(") Carbon-use efficiency, CUE = G/(G+R), collapses with warming and peaks early; filled points mark each curve’s peak. ("),
  lt("e", { bold: true }),
  lt(") Thermal optimum of CUE (filled points, 95% CrI; open points, growth optimum; small points, isolates), with the posterior probability that it lies below 37 °C. Shaded band, 37 to 40 °C (body to fever).")
]));

children.push(new Paragraph({ children: [new PageBreak()] }));

// ---- Section 2 ----
children.push(H1("A febrile host imposes a quantifiable, taxon-specific carbon tax"));

children.push(body([
  t("Because the economic optimum lies below 37 °C, every taxon operates at a carbon penalty in the host. We express this as the relative respiratory cost, i.e. respiration per unit growth relative to each taxon’s own cheapest temperature (Fig. 2a). At 40 °C (fever), this cost reached 1.34× for "),
  ...sp("C.", "auris"), t("Clade IV but 3.22× for "), ...sp("C.", "parapsilosis", ""), t(".")
]));

children.push(body([
  t("The three-degree step from body temperature (37 °C) to fever (40 °C) captures the clinical penalty most directly (Fig. 2c). Over this interval "),
  ...sp("C.", "parapsilosis"),
  t("retained only 57% of its 37 °C growth rate (0.57, 95% CrI 0.49-0.65) while its respiratory cost per unit growth rose 2.00× (1.76-2.32). The "),
  ...sp("C.", "auris"),
  t("clades were far less affected: Clade IV retained 89% of growth (0.89, 0.86-0.93) at a cost increase of only 1.25× (1.19-1.31), with Clades I-III intermediate (growth retained 0.67-0.84; cost 1.37×-1.72×). Fever thus suppresses growth without suppressing respiration, and it does so most severely in "),
  ...sp("C.", "parapsilosis", ""), t(".")
]));

children.push(body([
  b("A warmer environmental optimum predicts a lower cost of fever. "),
  t("Across the five taxa, the environmental growth optimum was strongly and negatively associated with the metabolic cost of fever: the more thermotolerant a taxon, the less a fever cost it (Fig. 2b). Computed as a posterior distribution over Spearman’s ρ (propagating parameter uncertainty rather than treating the five taxon medians as fixed), the correlation was ρ = −0.90 (95% CrI −1.00 to −0.30; P(ρ < 0) = 0.996). "),
  ...sp("C.", "parapsilosis"),
  t("occupied the cool-optimum, expensive-fever extreme and "),
  ...sp("C.", "auris"),
  t("Clade IV the warm-optimum, cheap-fever extreme. This links adaptation to a warm environment with tolerance of the febrile host as facets of a single thermal trait, consistent with the hypothesis that environmental thermal selection underlies fungal thermotolerance to the mammalian host. We note that T"),
  t("opt", { subScript: true }),
  t(" and fever cost are partly geometrically coupled through the growth curve; however, respiration (E"),
  t("R", { subScript: true }),
  t(") is an independent parameter that could have broken the relationship and did not. The correlation rests on only five taxa; it is therefore reported as suggestive and will be tested directly with additional species.")
]));

children.push(body([
  b("Clade-level differences within "), ...sp("C.", "auris", ""), b(". "),
  t("Pairwise comparison of the fever cost resolved 8 of 10 taxon contrasts (95% CrI excluding a ratio of 1; Fig. 2d). All four "),
  ...sp("C.", "auris"),
  t("clades credibly paid less at fever than "),
  ...sp("C.", "parapsilosis"),
  t("(cost ratios 0.42-0.70), and Clade IV was credibly the cheapest of all clades. Only two contrasts, Clade I versus Clade II (0.85, 0.68-1.05) and Clade I versus Clade III (1.18, 0.99-1.41), could not be separated at three isolates per clade.")
]));

children.push(figure("FIG2_the_bill.png", 6.4));
children.push(legend("Figure 2. The metabolic cost of a fever, and its relationship to environmental thermotolerance.", [
  lt("("), lt("a", { bold: true }),
  lt(") Relative respiratory cost (respiration per unit growth, relative to each taxon’s own cheapest temperature) versus temperature; filled points mark each minimum, which coincides with T"),
  new TextRun({ text: "opt", subScript: true, size: 20, font: "Calibri" }),
  lt("(CUE) in Fig. 1e. Dashed lines, body (37 °C) and fever (40 °C). ("),
  lt("b", { bold: true }),
  lt(") Environmental growth optimum versus the cost of a fever ([R/G] at 40 °C ÷ at 37 °C), five taxa (points, 95% CrI); the warmer a taxon’s optimum, the lower its fever cost (posterior Spearman ρ = −0.90). ("),
  lt("c", { bold: true }),
  lt(") The 37→40 °C bill: growth retained (growth at 40 °C ÷ at 37 °C, left) and the accompanying increase in respiratory cost per unit growth (right); points, posterior medians; bars, 95% CrI; 1 = no change. ("),
  lt("d", { bold: true }),
  lt(") Pairwise ratios of fever cost (I-IV, "), ...spL("C.", "auris"),
  lt("clades). Red, 95% CrI excludes 1; grey, unresolved.")
]));

// ---- References ----
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(H1("Reference"));
children.push(new Paragraph({
  spacing: { after: 160, line: 264 }, alignment: AlignmentType.JUSTIFIED,
  children: [
    t("Cakin I, Millington R, Pawar S, Buckling A, Smirnoff N, Padfield D, Duffy J, Yvon-Durocher G. "),
    t("A novel method to simultaneously estimate bacterial respiration and growth from oxygen dynamics. ",
      { italics: true }),
    t("ISME Communications "), t("2026;6(1):ycag024. "),
    t("https://doi.org/10.1093/ismeco/ycag024", { color: "0563C1" })
  ]
}));

// ---- build ----
const doc = new Document({
  creator: "Candida thermal energetics",
  styles: { default: { document: { run: { font: "Calibri", size: 22 } } } },
  sections: [{
    properties: { page: { margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    children
  }]
});

Packer.toBuffer(doc).then((buf) => { fs.writeFileSync(OUT, buf); console.log("wrote", OUT); });
