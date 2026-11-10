// tree_radial.js -- radial drug-indication tree.
//
// Copy of tree.js with four changes, so that it renders BOTH hierarchy shapes:
//   5 levels  Colorectal cancer -> Cell type -> Gene -> Drug -> Indication
//   4 levels  Colorectal cancer -> Gene -> Drug -> Indication   (no cell-type ring)
//
//   1. The drug ring is found from the tree height (drugDepth = root.height - 1)
//      instead of being hard-coded to depth 3. In tree.js the colour map and the
//      "hide the circle" rule fired on absolute depths, so with a 4-level JSON the
//      green/grey map landed on the indication ring and every drug came out grey.
//   2. Circles are hidden on the indication ring (drugDepth + 1), not on depth 4.
//   3. Font sizes are keyed to the role (gene / drug / indication), not to absolute
//      depth, so labels keep the same sizes in both shapes.
//   4. Optional legend for the approval colours (legend = true).
//
// Drug node colour comes from data.value: "true" = approved (green),
// anything else = clinical-trial stage (grey).
function Tree(data, {
  path,
  id = Array.isArray(data) ? d => d.id : null,
  parentId = Array.isArray(data) ? d => d.parentId : null,
  children,
  tree = d3.tree,
  separation = tree === d3.tree
    ? (a, b) => (a.parent == b.parent ? 1 : 2) / a.depth
    : (a, b) => a.parent == b.parent ? 1 : 2,
  sort,
  label,
  title,
  link,
  linkTarget = "_blank",
  width = 1600,
  height = 1600,
  margin = 200,
  marginTop = margin,
  marginRight = margin,
  marginBottom = margin,
  marginLeft = margin,
  radius = Math.min(
    width - marginLeft - marginRight,
    height - marginTop - marginBottom
  ) / 2,
  r = 8,
  padding = 1,
  fill = "#808080",
  fillOpacity,
  stroke = "#555",
  strokeWidth = 1.5,
  strokeOpacity = 0.4,
  strokeLinejoin,
  strokeLinecap,
  halo = "#fff",
  haloWidth = 3,

  // Ring that carries the drugs. Default: the ring just inside the leaves.
  drugDepth = null,

  // Label sizes, per role.
  geneFontSize = 19,
  drugFontSize = 16,
  indicationFontSize = 13,

  // Leaf labels longer than this are wrapped onto two lines.
  maxLength = 90,

  // Name of a boolean property on the node data; where it is true the label is
  // drawn bold. Used to mark drugs that carry a cancer indication once the
  // indication ring itself is dropped. null = every label normal weight.
  // Weight alone is easy to miss at this density, so nodes CARRYING the property
  // are also split by colour: true -> labelFill, false -> labelFillUnflagged.
  // Nodes without the property (root, genes) always use labelFill.
  boldFlag = null,
  labelFill = "#000",
  labelFillUnflagged = "#000",

  // Put "Helvetica" first if the SVG will be converted to PDF by svglib: given
  // "Arial" it resolves a system Arial TTF and ignores font-weight entirely, so
  // every label lands in the PDF as regular. With Helvetica leading, the base-14
  // Helvetica / Helvetica-Bold pair is used and the weights survive. Renderers on
  // machines without Helvetica fall through to Arial, which is metrically the same.
  fontFamily = "Helvetica, Arial, sans-serif",

  // Radial position of each ring as a fraction of `radius`, indexed by depth
  // ([root, gene, drug, indication]). d3.tree spaces rings evenly, which wastes
  // the middle of the plot and leaves the inner rings too tight for large type;
  // pushing the gene and drug rings outward buys arc length for their labels at
  // the same canvas size. null = even spacing (d3 default).
  ringFractions = null,

  // Horizontal title drawn just above the centre node. The root's own label is
  // placed radially like every other node, which leaves it rotated and buried in
  // the link bundle, so it is suppressed when this is set.
  rootTitle = null,
  rootTitleFontSize = 32,
  rootTitleOffset = 28,

  // Approval legend
  legend = true,
  legendApprovedText = "Approved drugs",
  legendClinicalText = "Drugs in clinical trials (phase II+)",
  legendFontSize = 20,
  legendPosition = "top-left",   // top-left | top-right | bottom-left | bottom-right
  // Extra text-only legend rows, e.g. to explain a label-weight encoding:
  // [{text: "Bold: ...", weight: "bold", fill: "#000"}, ...]
  legendNotes = [],

  // PDF export options
  pdfDownload = true,
  pdfFilename = "CRC_CellType_Gene_Drug_Indication_tree.pdf",
  pdfButtonText = "Download PDF",
} = {}) {

  // -----------------------------------------------------------------------
  // Color mappings
  // -----------------------------------------------------------------------
  // Drug node: colored by approval status (value = "true"/"false")
  const approveColorMap = {
    true:  "#70c17f",
    false: "#808080",
  };

  // -----------------------------------------------------------------------
  // Build hierarchy
  // -----------------------------------------------------------------------
  const root = path != null ? d3.stratify().path(path)(data)
      : id != null || parentId != null ? d3.stratify().id(id).parentId(parentId)(data)
      : d3.hierarchy(data, children);

  if (sort != null) root.sort(sort);

  // Drug ring: one level above the deepest level (the indications).
  const DRUG_DEPTH = drugDepth == null ? Math.max(root.height - 1, 1) : drugDepth;
  const IND_DEPTH = DRUG_DEPTH + 1;
  const GENE_DEPTH = DRUG_DEPTH - 1;

  const descendants = root.descendants();
  const L = label == null ? null : descendants.map(d => label(d.data, d));

  tree().size([2 * Math.PI, radius]).separation(separation)(root);

  // Re-place the rings radially (angles are untouched, so no label can collide
  // with a new neighbour; only the ring-to-ring distances change).
  if (ringFractions != null) {
    root.each(d => {
      const f = ringFractions[Math.min(d.depth, ringFractions.length - 1)];
      d.y = radius * f;
    });
  }

  // -----------------------------------------------------------------------
  // Tight square SVG around the radial tree.
  // -----------------------------------------------------------------------
  const plotWidth = 2 * radius + marginLeft + marginRight;
  const plotHeight = 2 * radius + marginTop + marginBottom;
  const viewX = -marginLeft - radius;
  const viewY = -marginTop - radius;

  const svg = d3.create("svg")
      .attr("viewBox", [viewX, viewY, plotWidth, plotHeight])
      .attr("width", plotWidth)
      .attr("height", plotHeight)
      .attr("style", "max-width: 100%; height: auto;")
      .attr("font-family", fontFamily)
      .attr("font-size", 8);

  // -----------------------------------------------------------------------
  // Links
  // -----------------------------------------------------------------------
  svg.append("g")
      .attr("fill", "none")
      .attr("stroke", stroke)
      .attr("stroke-opacity", strokeOpacity)
      .attr("stroke-linecap", strokeLinecap)
      .attr("stroke-linejoin", strokeLinejoin)
      .attr("stroke-width", strokeWidth)
    .selectAll("path")
    .data(root.links())
    .join("path")
      .attr("d", d3.linkRadial()
          .angle(d => d.x)
          .radius(d => d.y));

  // -----------------------------------------------------------------------
  // Nodes
  // -----------------------------------------------------------------------
  const node = svg.append("g")
    .selectAll("a")
    .data(root.descendants())
    .join("a")
      .attr("xlink:href", link == null ? null : d => link(d.data, d))
      .attr("target", link == null ? null : linkTarget)
      .attr("transform", d =>
        `rotate(${d.x * 180 / Math.PI - 90}) translate(${d.y},0)`
      );

  node.append("circle")
      .attr("fill", d => {
        if (d.depth === DRUG_DEPTH) {
          return approveColorMap[String(d.data.value).toLowerCase()] || fill;
        }
        return fill;
      })
      .attr("r", d => {
        if (d.depth >= IND_DEPTH) return 0;
        return Math.max(r - d.depth * 1.5, 1);
      });

  if (title != null) {
    node.append("title").text(d => title(d.data, d));
  }

  // -----------------------------------------------------------------------
  // Labels
  // -----------------------------------------------------------------------
  if (L) {
    node.append("text")
      .attr("transform", d => `rotate(${d.x >= Math.PI ? 180 : 0})`)
      .attr("dy", "0.32em")
      .attr("x", d => d.x < Math.PI === !d.children ? 6 : -6)
      .attr("text-anchor", d => d.x < Math.PI === !d.children ? "start" : "end")
      .attr("paint-order", "stroke")
      .attr("stroke", halo)
      .attr("stroke-width", haloWidth)
      .each(function(d, i) {
        const text = d3.select(this);
        const labelText = L[i];
        const isLeaf = !d.children;

        if (d.depth === 0 && rootTitle != null) {
          text.text(null);   // drawn horizontally below instead
        } else if (isLeaf && labelText && labelText.length > maxLength) {
          const splitIndex = labelText.lastIndexOf(" ", maxLength) || maxLength;
          const part1 = labelText.substring(0, splitIndex);
          const part2 = labelText.substring(splitIndex + 1);

          text.text(null)
            .append("tspan")
              .attr("x", text.attr("x"))
              .attr("dy", "-0.1em")
              .text(part1);

          text.append("tspan")
              .attr("x", text.attr("x"))
              .attr("dy", "1.1em")
              .text(part2);
        } else {
          text.text(labelText);
        }
      })
      // font-weight and fill are set as presentation ATTRIBUTES, not in style:
      // svglib (the SVG -> PDF step) reads attributes but ignores font-weight
      // declared inside a style attribute, which silently flattened every label
      // to regular Arial in the PDF.
      .attr("font-weight", d =>
        boldFlag != null && d.data[boldFlag] ? "bold" : "normal")
      .attr("fill", d => {
        if (boldFlag == null || !Object.prototype.hasOwnProperty.call(d.data, boldFlag)) {
          return labelFill;
        }
        return d.data[boldFlag] ? labelFill : labelFillUnflagged;
      })
      .style("font-size", d => {
        if (d.depth <= GENE_DEPTH) return `${geneFontSize}px`;   // cell type / gene
        if (d.depth === DRUG_DEPTH) return `${drugFontSize}px`;  // drug
        return `${indicationFontSize}px`;                        // indication
      });
  }

  // -----------------------------------------------------------------------
  // Horizontal title above the centre node
  // -----------------------------------------------------------------------
  if (rootTitle != null) {
    svg.append("text")
        .attr("x", 0)
        .attr("y", -rootTitleOffset)
        .attr("text-anchor", "middle")
        .attr("font-size", rootTitleFontSize)
        .attr("paint-order", "stroke")
        .attr("stroke", halo)
        .attr("stroke-width", haloWidth)
        .style("fill", "#000")
        .text(rootTitle);
  }

  // -----------------------------------------------------------------------
  // Legend for the drug-node colours
  // -----------------------------------------------------------------------
  if (legend) {
    const legendRows = [[approveColorMap.true, legendApprovedText],
                        [approveColorMap.false, legendClinicalText]];
    const rowStep = legendFontSize * 1.6;
    const pad = 80;
    // Width estimate for the block, so a right-hand corner stays inside the canvas.
    const allText = legendRows.map(([, t]) => t).concat(legendNotes.map(n => n.text));
    const blockWidth = Math.max(...allText.map(t => t.length))
                       * 0.55 * legendFontSize + 2 * legendFontSize;
    const nRows = legendRows.length + legendNotes.length;

    const lx = legendPosition.endsWith("right")
      ? viewX + plotWidth - pad - blockWidth
      : viewX + pad;
    const ly = legendPosition.startsWith("bottom")
      ? viewY + plotHeight - pad - (nRows - 1) * rowStep
      : viewY + pad + legendFontSize;

    const lg = svg.append("g")
        .attr("transform", `translate(${lx},${ly})`)
        .attr("font-size", legendFontSize);

    legendRows.forEach(([colour, text], i) => {
      const row = lg.append("g")
          .attr("transform", `translate(0,${i * rowStep})`);
      row.append("circle")
          .attr("r", Math.max(r - DRUG_DEPTH * 1.5, 1))
          .attr("fill", colour);
      row.append("text")
          .attr("x", 16)
          .attr("dy", "0.32em")
          .attr("fill", "#000")
          .text(text);
    });

    legendNotes.forEach((note, i) => {
      lg.append("text")
          .attr("x", 16)
          .attr("y", (legendRows.length + i) * rowStep)
          .attr("dy", "0.32em")
          .attr("font-weight", note.weight || "normal")
          .attr("fill", note.fill || "#000")
          .text(note.text);
    });
  }

  // -----------------------------------------------------------------------
  // Return SVG only if PDF export is disabled.
  // -----------------------------------------------------------------------
  if (!pdfDownload) return svg.node();

  // -----------------------------------------------------------------------
  // Observable/browser wrapper + automatic PDF download button.
  // -----------------------------------------------------------------------
  const wrapper = document.createElement("div");
  wrapper.style.display = "inline-flex";
  wrapper.style.flexDirection = "column";
  wrapper.style.alignItems = "flex-start";
  wrapper.style.gap = "10px";

  const svgNode = svg.node();
  wrapper.appendChild(svgNode);

  const button = document.createElement("button");
  button.textContent = pdfButtonText;
  button.style.padding = "7px 12px";
  button.style.fontSize = "14px";
  button.style.cursor = "pointer";

  button.addEventListener("click", async () => {
    const oldText = button.textContent;
    button.disabled = true;
    button.textContent = "Generating PDF...";

    try {
      const [{ jsPDF }, svg2pdfModule] = await Promise.all([
        import("https://cdn.jsdelivr.net/npm/jspdf@2.5.1/+esm"),
        import("https://cdn.jsdelivr.net/npm/svg2pdf.js@2.2.4/+esm")
      ]);

      const svg2pdf = svg2pdfModule.svg2pdf || svg2pdfModule.default;
      if (typeof svg2pdf !== "function") {
        throw new Error("Could not load svg2pdf.js");
      }

      const exportSvg = svgNode.cloneNode(true);
      exportSvg.setAttribute("width", plotWidth);
      exportSvg.setAttribute("height", plotHeight);
      exportSvg.setAttribute("viewBox", `${viewX} ${viewY} ${plotWidth} ${plotHeight}`);
      exportSvg.style.maxWidth = "none";
      exportSvg.style.height = "auto";

      const pdf = new jsPDF({
        orientation: plotWidth >= plotHeight ? "landscape" : "portrait",
        unit: "pt",
        format: [plotWidth, plotHeight],
        compress: true
      });

      await svg2pdf(exportSvg, pdf, {
        x: 0,
        y: 0,
        width: plotWidth,
        height: plotHeight
      });

      pdf.save(pdfFilename);
      button.textContent = "PDF downloaded";
      setTimeout(() => {
        button.textContent = oldText;
      }, 1500);

    } catch (error) {
      console.error(error);
      button.textContent = "PDF export failed";
      alert(`PDF export failed: ${error.message}`);
      setTimeout(() => {
        button.textContent = oldText;
      }, 2500);
    } finally {
      button.disabled = false;
    }
  });

  wrapper.appendChild(button);
  return wrapper;
}

if (typeof module !== "undefined" && module.exports) {
  module.exports = Tree;
}
