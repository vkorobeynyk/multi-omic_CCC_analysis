# Builds evidence_graph.html from evidence.yaml: one network per LR pair, with a
# selector at the top to switch between them.
#
# Edit evidence.yaml (the single document), then:  Rscript build_graph.R
# Needs: yaml, jsonlite

library(yaml)
library(jsonlite)

ev <- read_yaml("evidence.yaml")

# written into the project folder so the report links (e.g. "11_Proteomics_LR_AgeTrends.html") resolve
out_html <- file.path("..", "evidence_graph.html")

effect_color <- c(increase = "#d55e00", decrease = "#0072b2", no_change = "#999999",
                  not_detected = "#cccccc")
col_for <- function(effect) if (!is.null(effect) && effect %in% names(effect_color)) effect_color[[effect]] else "#333333"
or_empty <- function(x) if (is.null(x)) "" else as.character(x)

build_pair <- function(p) {
  nodes <- list(); edges <- list()
  add_node <- function(id, label, shape, bg, border, title, size = 14) {
    if (!id %in% vapply(nodes, `[[`, "", "id"))
      nodes[[length(nodes) + 1]] <<- list(id = id, label = label, shape = shape, title = title,
                                           color = list(background = bg, border = border),
                                           font = list(size = size))
  }
  add_edge <- function(from, to, label, color, report, title, arrows = "to", dashes = FALSE, width = 2) {
    edges[[length(edges) + 1]] <<- list(id = length(edges) + 1, from = from, to = to, label = label,
                                        color = list(color = color), report = report, title = title,
                                        arrows = arrows, dashes = dashes, width = width)
  }

  # ligand - receptor
  add_node(p$ligand,   p$ligand,   "ellipse", "#dbe9f6", "#2a78d6",
           paste0("Ligand", if (!is.null(p$sender)) paste0(" (", p$sender, ")")), size = 22)
  add_node(p$receptor, p$receptor, "ellipse", "#dbe9f6", "#2a78d6",
           paste0("Receptor", if (!is.null(p$receiver)) paste0(" (", p$receiver, ")")), size = 22)
  add_edge(p$ligand, p$receptor, "LR pair", "#333333", or_empty(p$report),
           paste0(p$name, if (!is.null(p$notes)) paste0("\n", p$notes)), arrows = "", width = 5)

  for (e in p$evidence) {
    src_id <- paste0("src:", e$source)
    add_node(src_id, e$source, "box", "#f2f2f2", "#555555", e$source)
    add_edge(src_id, e$gene, or_empty(e$label), col_for(e$effect), or_empty(e$report),
             paste0(e$source, " -> ", e$gene, "\n", or_empty(e$label),
                    if (!is.null(e$notes)) paste0("\n", e$notes),
                    if (!is.null(e$report)) paste0("\nclick to open: ", e$report)))
  }

  for (s in p$stainings) {
    st_id <- paste0("st:", s$name)
    add_node(st_id, s$name, "diamond", "#fde3cf", "#d55e00",
             paste0(s$name, " (", or_empty(s$status), ")", if (!is.null(s$notes)) paste0("\n", s$notes)))
    planned <- identical(s$status, "planned")
    for (g in s$genes) {
      add_edge(st_id, g, if (planned) "planned" else or_empty(s$label), col_for(s$effect),
               or_empty(s$report),
               paste0(s$name, " -> ", g, if (!is.null(s$notes)) paste0("\n", s$notes)),
               dashes = planned)
    }
  }
  list(name = p$name, nodes = nodes, edges = edges)
}

networks <- lapply(ev$pairs, build_pair)
stopifnot(length(networks) > 0)

html <- sprintf('<!doctype html>
<html><head><meta charset="utf-8"><title>LR evidence graphs</title>
<script src="https://unpkg.com/vis-network@9.1.9/standalone/umd/vis-network.min.js"></script>
<style>
  body { font-family: sans-serif; margin: 0; }
  #bar { padding: 8px 14px; font-size: 13px; color: #444; border-bottom: 1px solid #ddd; }
  #bar select { font-size: 14px; margin-right: 16px; }
  #graph { width: 100vw; height: 90vh; }
  .k { display: inline-block; width: 22px; height: 3px; vertical-align: middle; margin: 0 4px 0 12px; }
</style></head><body>
<div id="bar">
  <b>LR pair:</b> <select id="pair"></select>
  built %s &nbsp;|
  <span class="k" style="background:#d55e00"></span>increase with age
  <span class="k" style="background:#0072b2"></span>decrease
  <span class="k" style="background:#999"></span>no change
  <span class="k" style="background:#ccc"></span>not detected
  <span class="k" style="background:#333;height:5px"></span>LR pair
  &nbsp; dashed = planned &nbsp; (click an arrow to open its report)
</div>
<div id="graph"></div>
<script>
const networks = %s;
const sel = document.getElementById("pair");
networks.forEach((n, i) => sel.add(new Option(n.name, i)));
let net = null;
function show(i) {
  const n = networks[i];
  const edges = new vis.DataSet(n.edges);
  if (net) net.destroy();
  net = new vis.Network(document.getElementById("graph"), { nodes: new vis.DataSet(n.nodes), edges }, {
    physics: { solver: "forceAtlas2Based", forceAtlas2Based: { springLength: 180, gravitationalConstant: -80 },
               stabilization: { iterations: 300 } },
    edges: { font: { size: 12, align: "middle", background: "white" }, smooth: { type: "dynamic" } },
    interaction: { hover: true, tooltipDelay: 150 }
  });
  net.on("click", p => {
    if (p.edges.length && !p.nodes.length) {
      const e = edges.get(p.edges[0]);
      if (e.report) window.open(e.report, "_blank");
    }
  });
}
sel.onchange = () => show(sel.value);
show(0);
</script></body></html>',
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  toJSON(networks, auto_unbox = TRUE))

writeLines(html, out_html)
cat("Wrote", normalizePath(out_html), "with", length(networks), "LR pair network(s):",
    paste(vapply(networks, `[[`, "", "name"), collapse = ", "), "\n")
