const data = window.PROJECT_DATA || {};
const summary = data.summary || {};
const sigCount = summary["Adjusted P-value < 0.05"] ?? 1;
document.getElementById("sig-count").textContent = sigCount.toLocaleString();
document.getElementById("sig-probes").textContent = sigCount.toLocaleString();

const probes = data.significantProbes || [];
if (probes.length) {
  const p = probes[0];
  document.getElementById("probe-id").textContent = p.Probe_ID ?? "—";
  document.getElementById("probe-logfc").textContent = Number(p.logFC).toFixed(4);
  const padj = document.getElementById("probe-pad j");
  if (padj) padj.textContent = Number(p["adj.P.Val"]).toExponential(2);
}
const pathways = data.topPathways || [];
if (pathways.length) {
  const max = Math.max(...pathways.map(x => Number(x.NES) || 0));
  document.getElementById("pathway-list").innerHTML = pathways.slice(0, 8).map(p => {
    const width = Math.max(3, (Number(p.NES) / max) * 100);
    return `<div class="pathway-row"><span>${p.Description}</span><div class="bar-track"><i style="width:${width}%"></i></div><b>${Number(p.NES).toFixed(2)}</b></div>`;
  }).join("");
}
const toggle = document.querySelector(".menu-toggle");
const nav = document.getElementById("nav");
toggle.addEventListener("click", () => {
  const opened = nav.classList.toggle("open");
  toggle.setAttribute("aria-expanded", String(opened));
});
nav.querySelectorAll("a").forEach(a => a.addEventListener("click", () => nav.classList.remove("open")));
