// Render the clinically advanced drug-target radial tree to a standalone SVG.
//
//   npm install --prefix <staged-drugs-directory> d3@7 jsdom@26
//   node 05_render_radial_tree.mjs
//
// Uses the same Tree() renderer the browser/Observable notebook uses
// (src/tree_radial.js), driven headlessly through jsdom so the figure can be
// regenerated from the command line without a browser.
//
// Reads   table/63_druggable_genes_advanced_55genes_193drugs.json
//         src/tree_radial.js
// Writes  figure/CRC_scTWAS_CellType_Gene_Drug_Indication.svg
import fs from "node:fs";
import { createRequire } from "node:module";
import path from "node:path";
import vm from "node:vm";
import { fileURLToPath, pathToFileURL } from "node:url";

const ROOT = process.env.SCTWAS_DRUGS_ROOT || process.cwd();
const dependencies = process.env.SCTWAS_NODE_DEPS || ROOT;
const require = createRequire(path.resolve(dependencies, "package.json"));
const d3 = await import(pathToFileURL(require.resolve("d3")).href);
const { JSDOM } = require("jsdom");
const JSON_PATH = `${ROOT}/table/63_druggable_genes_advanced_55genes_193drugs.json`;
const TREE_PATH = fileURLToPath(new URL("./tree_radial.js", import.meta.url));
const OUT_SVG = `${ROOT}/figure/CRC_scTWAS_CellType_Gene_Drug_Indication.svg`;

// --- browser-like globals for tree_radial.js --------------------------------
const dom = new JSDOM("<!DOCTYPE html><body></body>");
globalThis.document = dom.window.document;
globalThis.d3 = d3;
vm.runInThisContext(fs.readFileSync(TREE_PATH, "utf8"), { filename: TREE_PATH });

// --- render ------------------------------------------------------------------
const data = JSON.parse(fs.readFileSync(JSON_PATH, "utf8"));

// Two rings only: Colorectal cancer -> gene -> drug. The indication layer is not
// drawn; instead a drug that carries at least one cancer indication is flagged so
// its label prints bold. The JSON keeps the indications, so nothing is lost and
// the sourcedata step can still list them.
let nBold = 0;
for (const gene of data.children) {
  for (const drug of gene.children) {
    drug.hasCancerIndication = drug.children.length > 0;
    nBold += drug.hasCancerIndication ? 1 : 0;
    delete drug.children;
  }
}

// Figure-only short labels, carried in the JSON so the source data can report the
// same substitution (meta.display_overrides).
const SHORT = (data.meta && data.meta.display_overrides) || {};

const svgNode = Tree(data, {
  label: d => SHORT[d.name] || d.name,
  title: (d, n) => n.ancestors().reverse().map(a => a.data.name).join("  >  "),
  width: 3880,
  height: 3880,
  margin: 640,   // the one 70-char drug name wraps to a 40-char line at 28px = 616px
  // Bigger type needs arc length, and the cheapest source of it is the gap the
  // layout puts between one gene's drugs and the next: at d3's default (siblings
  // 1, cross-parent 2) the 54 inter-gene gaps eat ~40% of the circle. Cutting the
  // cross-parent term to 1.4 returns most of that to the labels and lifts the font
  // ceiling from 28/26 to 33/30 at this radius, at the cost of slightly less
  // obvious visual grouping per gene.
  separation: (a, b) => (a.parent == b.parent ? 1 : 1.4) / a.depth,
  // radius 1300, rings at 715 (gene) and 1300 (drug). Dropping the indication ring
  // makes the 220 drugs the leaves, so the layout is re-derived: at this radius the
  // minimum arc gap is 32.7px on the gene ring (28px type needs 32) and 29.7px on
  // the drug ring (24px type needs 28). Anything below radius 1300 collides.
  // gene ring at 0.57 (741px), not 0.55: at 32px type the labels need 36.8px of
  // arc and 0.55 gave 36.0, so 14 pairs were fractionally tight.
  ringFractions: [0, 0.57, 1.0],
  drugDepth: 2,          // must be explicit now: root.height - 1 would be the genes
  // Cancer-indication drugs: bold black. The rest: regular, mid grey. Weight on
  // its own was too weak to read at this density (and was being dropped entirely
  // in the PDF), so the two cues are applied together.
  boldFlag: "hasCancerIndication",
  labelFill: "#000000",
  labelFillUnflagged: "#6e6e6e",
  geneFontSize: 32,
  drugFontSize: 28,
  r: 11,                // node circles: 11 root / 9.5 gene / 8 drug
  maxLength: 40,        // wraps the one 70-char drug name onto two lines
  rootTitle: "Colorectal cancer",
  rootTitleFontSize: 38,
  legend: true,
  legendFontSize: 34,
  legendPosition: "bottom-right",
  legendNotes: [
    { text: "Black: \u22651 cancer indication", weight: "bold", fill: "#000000" },
    { text: "Grey: no cancer indication", weight: "normal", fill: "#6e6e6e" },
  ],
  pdfDownload: false,   // return the bare <svg> node, no download button
});

const xml = new dom.window.XMLSerializer().serializeToString(svgNode);
fs.mkdirSync(`${ROOT}/figure`, { recursive: true });
fs.writeFileSync(OUT_SVG, `<?xml version="1.0" encoding="UTF-8"?>\n${xml}\n`, "utf8");

// --- report what was drawn ---------------------------------------------------
const root = d3.hierarchy(data);
const drugDepth = 2;
const drugs = root.descendants().filter(d => d.depth === drugDepth);
const green = drugs.filter(d => String(d.data.value).toLowerCase() === "true").length;

console.log(
  `genes ${root.descendants().filter(d => d.depth === drugDepth - 1).length} | ` +
  `drug nodes ${drugs.length} (green ${green}, grey ${drugs.length - green}) | ` +
  `bold (has cancer indication) ${nBold} | ` +
  `svg ${(fs.statSync(OUT_SVG).size / 1024).toFixed(0)} kB`
);
console.log(OUT_SVG);
