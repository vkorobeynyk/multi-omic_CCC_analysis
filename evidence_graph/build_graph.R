# Builds evidence_graph.html from evidence.yaml: one network per LR pair (dropdown at the top).
#   ligand -- [interaction] -- receptor, plus one arrow per evidence entry from its technology
#   to the interaction (target: pair) or to a single gene (target: <gene>).
# Click a node or an arrow: the side panel shows technology, samples, method, description and
# links to the reports.
#
# Usage:  Rscript build_graph.R   (run from this folder; writes evidence_graph.html here)   Needs: yaml, jsonlite

library(yaml)
library(jsonlite)

# make sure non-ASCII text (en dashes, µ, ’) survives even if R runs in a non-UTF-8 locale
if (!l10n_info()$`UTF-8`) for (loc in c("C.UTF-8", "en_US.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break

ev <- yaml.load(paste(readLines("evidence.yaml", encoding = "UTF-8", warn = FALSE), collapse = "\n"))
out_html   <- "evidence_graph.html"   # written next to this script
report_dir <- ".."                    # where the HTML reports live, relative to out_html (the project folder)

effect_color <- c(increase = "#d55e00", decrease = "#0072b2", no_change = "#999999",
                  not_detected = "#cccccc", weak = "#c9a227")
modality_style <- list(transcriptomics = list(bg = "#e3f2e1", border = "#2e7d32"),
                       proteomics      = list(bg = "#efe6f8", border = "#6a3fa0"),
                       staining        = list(bg = "#fde3cf", border = "#d55e00"))
col_for  <- function(x) if (!is.null(x) && x %in% names(effect_color)) effect_color[[x]] else "#333333"
or_empty <- function(x) if (is.null(x)) "" else paste(x, collapse = " ")
as_list  <- function(x) if (is.null(x)) list() else as.list(x)
# report paths in evidence.yaml are relative to the project folder; make them relative to out_html
report_links <- function(x) lapply(as_list(x), function(r)
  if (grepl("^(/|[a-z]+://)", r)) r else file.path(report_dir, r))

build_pair <- function(p) {
  nodes <- list(); edges <- list()
  add_node <- function(id, label, shape, bg, border, info, size = 14) {
    if (!id %in% vapply(nodes, `[[`, "", "id"))
      nodes[[length(nodes) + 1]] <<- list(id = id, label = label, shape = shape, info = info,
                                           color = list(background = bg, border = border),
                                           font = list(size = size))
  }
  add_edge <- function(from, to, label, color, info, arrows = "to", dashes = FALSE, width = 2) {
    edges[[length(edges) + 1]] <<- list(id = length(edges) + 1, from = from, to = to, label = label,
                                        color = list(color = color), info = info,
                                        arrows = arrows, dashes = dashes, width = width,
                                        # wrapped labels are ~22 characters per line; give each line room
                                        length = if (nchar(label) == 0) 160 else 220 + 25 * ceiling(nchar(label) / 22))
  }
  pair_id <- paste0("pair:", p$name)
  target_id <- function(t) if (is.null(t) || t == "pair") pair_id else t

  add_node(p$ligand, p$ligand, "ellipse", "#dbe9f6", "#2a78d6",
           list(title = p$ligand, fields = list(role = "ligand", sender = or_empty(p$sender))), size = 22)
  add_node(p$receptor, p$receptor, "ellipse", "#dbe9f6", "#2a78d6",
           list(title = p$receptor, fields = list(role = "receptor", receiver = or_empty(p$receiver))), size = 22)
  add_node(pair_id, "interaction", "dot", "#333333", "#333333",
           list(title = p$name, fields = list(sender = or_empty(p$sender), receiver = or_empty(p$receiver),
                                              verdict = or_empty(p$verdict), comments = or_empty(p$comments)),
                literature = lapply(as_list(p$literature), function(l)
                  list(citation = or_empty(l$citation), link = or_empty(l$link), summary = or_empty(l$summary)))),
           size = 12)
  add_edge(p$ligand, pair_id, "", "#333333", NULL, arrows = "", width = 5)
  add_edge(pair_id, p$receptor, "", "#333333", NULL, arrows = "to", width = 5)

  for (e in p$evidence) {
    src_id <- paste0("tech:", e$technology)
    tech <- ev$technologies[[e$technology]]
    # modality: from `technologies:` in evidence.yaml if given, else proteomics if the name says so
    modality <- if (!is.null(tech$modality)) tech$modality else
      if (grepl("proteom", e$technology, ignore.case = TRUE)) "proteomics" else "transcriptomics"
    mstyle <- modality_style[[modality]]
    if (is.null(mstyle)) mstyle <- list(bg = "#f2f2f2", border = "#555555")
    add_node(src_id, e$technology, "box", mstyle$bg, mstyle$border,
             list(title = e$technology, fields = list(citation = or_empty(tech$citation)),
                  links = if (!is.null(tech$paper)) list(list(label = "Paper", url = tech$paper)) else list()))
    info <- list(title = paste0(e$technology, " -> ", if (is.null(e$target) || e$target == "pair") p$name else e$target),
                 fields = list(effect = or_empty(e$effect), samples = or_empty(e$samples),
                               method = or_empty(e$method), description = or_empty(e$description)),
                 reports = report_links(e$reports))
    add_edge(src_id, target_id(e$target), or_empty(e$label), col_for(e$effect), info)
  }

  for (s in p$stainings) {
    st_id <- paste0("st:", s$name)
    planned <- identical(s$status, "planned")
    info <- list(title = s$name, fields = list(status = or_empty(s$status), samples = or_empty(s$samples),
                                               description = or_empty(s$description)),
                 reports = report_links(s$reports))
    add_node(st_id, s$name, "diamond", "#fde3cf", "#d55e00", info)
    add_edge(st_id, target_id(s$target), if (planned) "planned" else or_empty(s$label),
             col_for(s$effect), info, dashes = planned)
  }
  list(name = p$name, nodes = nodes, edges = edges)
}

networks <- lapply(ev$pairs, build_pair)
stopifnot(length(networks) > 0)

html <- sprintf('<!doctype html>
<html><head><meta charset="utf-8"><title>LR evidence graphs</title>
<script src="https://unpkg.com/vis-network@9.1.9/standalone/umd/vis-network.min.js"></script>
<style>
  body { font-family: sans-serif; margin: 0; display: grid; grid-template-rows: auto 1fr; height: 100vh; }
  #bar { padding: 8px 14px; font-size: 13px; color: #444; border-bottom: 1px solid #ddd; }
  #bar select { font-size: 14px; margin-right: 14px; }
  #main { display: grid; grid-template-columns: 1fr 380px; min-height: 0; }
  #graph { min-height: 0; }
  #panel { border-left: 1px solid #ddd; padding: 14px 16px; overflow-y: auto; font-size: 13px; line-height: 1.45; color: #222; }
  #panel h3 { margin: 0 0 10px; font-size: 15px; }
  #panel .f { margin-bottom: 8px; } #panel .f b { color: #555; text-transform: capitalize; }
  .sw { display: inline-block; width: 16px; height: 11px; border: 1.5px solid; border-radius: 2px; vertical-align: middle; margin: 0 4px 0 10px; }
  .lit { border-left: 3px solid #2a78d6; background: #f6f8fb; padding: 6px 8px; margin: 6px 0 10px 0; }
  .lit .cit { color: #333; }
  .lit .sum { margin-top: 6px; padding-top: 6px; border-top: 1px dashed #c9d3e0; color: #222; }
  .lit .sum b { color: #2a78d6; }
  .k { display: inline-block; width: 22px; height: 3px; vertical-align: middle; margin: 0 4px 0 10px; }
</style></head><body>
<div id="bar">
  <b>LR pair:</b> <select id="pair"></select> built %s |
  <span class="k" style="background:#d55e00"></span>increase
  <span class="k" style="background:#0072b2"></span>decrease
  <span class="k" style="background:#c9a227"></span>weak / unclear
  <span class="k" style="background:#999"></span>no change
  <span class="k" style="background:#ccc"></span>not detected
  &nbsp; dashed = planned &nbsp;|
  <span class="sw" style="background:#e3f2e1;border-color:#2e7d32"></span>transcriptomics
  <span class="sw" style="background:#efe6f8;border-color:#6a3fa0"></span>proteomics
</div>
<div id="main"><div id="graph"></div><div id="panel"><i>Click a node or an arrow for details.</i></div></div>
<script>
const networks = %s;
const sel = document.getElementById("pair"), panel = document.getElementById("panel");
networks.forEach((n, i) => sel.add(new Option(n.name, i)));
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;","\\"":"&quot;"}[c]));
function showInfo(info) {
  if (!info) { panel.innerHTML = "<i>No details for this element.</i>"; return; }
  let h = "<h3>" + esc(info.title) + "</h3>";
  for (const [k, v] of Object.entries(info.fields || {})) if (v) h += "<div class=f><b>" + esc(k) + ":</b> " + esc(v) + "</div>";
  for (const l of [].concat(info.links || [])) h += "<div class=f><b>" + esc(l.label) + ":</b> <a href=\\"" + esc(l.url) + "\\" target=_blank>" + esc(l.url) + "</a></div>";
  const lit = [].concat(info.literature || []);
  if (lit.length) h += "<div class=f><b>Literature:</b></div>" + lit.map(l =>
    "<div class=lit><div class=cit>" + esc(l.citation) +
    (l.link ? " <a href=\\"" + esc(l.link) + "\\" target=_blank>[link]</a>" : "") + "</div>" +
    (l.summary ? "<div class=sum><b>Summary:</b> " + esc(l.summary) + "</div>" : "") + "</div>").join("");
  const reps = [].concat(info.reports || []);
  if (reps.length) h += "<div class=f><b>Reports:</b><br>" + reps.map(r => "<a href=\\"" + esc(r) + "\\" target=_blank>" + esc(r) + "</a>").join("<br>") + "</div>";
  panel.innerHTML = h;
}
let net = null;
function show(i) {
  const n = networks[i];
  const nodes = new vis.DataSet(n.nodes), edges = new vis.DataSet(n.edges);
  if (net) net.destroy();
  net = new vis.Network(document.getElementById("graph"), { nodes, edges }, {
    physics: { solver: "forceAtlas2Based",
               forceAtlas2Based: { springLength: 220, gravitationalConstant: -250, avoidOverlap: 0.6 },
               stabilization: { iterations: 500 } },
    edges: { font: { size: 12, align: "horizontal", background: "white", strokeWidth: 0 },
             widthConstraint: { maximum: 160 }, smooth: { type: "continuous" } },
    layout: { randomSeed: 1 },
    nodes: { size: 14 },
    interaction: { hover: true }
  });
  net.on("click", p => {
    if (p.nodes.length) showInfo(nodes.get(p.nodes[0]).info);
    else if (p.edges.length) showInfo(edges.get(p.edges[0]).info);
  });
  const pairNode = n.nodes.find(x => x.id.startsWith("pair:"));
  if (pairNode && pairNode.info.fields.verdict) showInfo(pairNode.info);
  else panel.innerHTML = "<i>Click a node or an arrow for details.</i>";
}
sel.onchange = () => show(sel.value);
show(0);
</script></body></html>',
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  toJSON(networks, auto_unbox = TRUE, null = "null"))

writeLines(enc2utf8(html), out_html, useBytes = TRUE)
cat("Wrote", normalizePath(out_html), "with", length(networks), "LR pair network(s):",
    paste(vapply(networks, `[[`, "", "name"), collapse = ", "), "\n")
