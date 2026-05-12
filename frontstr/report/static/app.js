/* FRONTStr report — minimal vanilla interactivity.
 * Profile table: sort + free-text filter + status filter.
 * Per-locus details: keep state on open/close so anchor navigation works.
 *
 * No frameworks. No fetches. Works on a file:// URL.
 */
(() => {
  "use strict";

  /* ---------- profile table: sort + filter ---------- */
  const profile = document.querySelector("table.profile");
  const searchInput = document.querySelector("#profile-search");
  const statusSelect = document.querySelector("#profile-status");
  const countPill = document.querySelector("#profile-count");

  if (profile) {
    const tbody = profile.tBodies[0];
    const allRows = Array.from(tbody.rows);

    const visibleCount = () => allRows.filter((r) => r.style.display !== "none").length;
    const updatePill = () => {
      if (countPill) countPill.textContent = `${visibleCount()} / ${allRows.length} markers`;
    };

    const applyFilters = () => {
      const q = (searchInput?.value || "").trim().toLowerCase();
      const status = statusSelect?.value || "";
      for (const tr of allRows) {
        const text = tr.dataset.search || "";
        const chip = tr.dataset.chip || "";
        const matchesText = !q || text.includes(q);
        const matchesStatus = !status || chip === status;
        tr.style.display = matchesText && matchesStatus ? "" : "none";
      }
      updatePill();
    };
    searchInput?.addEventListener("input", applyFilters);
    statusSelect?.addEventListener("change", applyFilters);

    const headers = profile.tHead.querySelectorAll("th[data-sortable]");
    headers.forEach((th, idx) => {
      th.addEventListener("click", () => {
        const dir = th.dataset.sortDir === "asc" ? "desc" : "asc";
        headers.forEach((h) => h.removeAttribute("data-sort-dir"));
        th.dataset.sortDir = dir;
        const type = th.dataset.sortable;
        const sorted = [...allRows].sort((a, b) => {
          let av = a.cells[idx]?.dataset.sortValue ?? a.cells[idx]?.textContent?.trim() ?? "";
          let bv = b.cells[idx]?.dataset.sortValue ?? b.cells[idx]?.textContent?.trim() ?? "";
          if (type === "num") {
            av = parseFloat(av) || -Infinity;
            bv = parseFloat(bv) || -Infinity;
          }
          if (av < bv) return dir === "asc" ? -1 : 1;
          if (av > bv) return dir === "asc" ? 1 : -1;
          return 0;
        });
        for (const r of sorted) tbody.appendChild(r);
      });
    });

    updatePill();
  }

  /* ---------- locus jump-anchor smooth scroll ---------- */
  document.querySelectorAll("a[href^='#locus-']").forEach((a) => {
    a.addEventListener("click", (e) => {
      const id = a.getAttribute("href").slice(1);
      const tgt = document.getElementById(id);
      if (tgt && tgt.tagName === "DETAILS") {
        e.preventDefault();
        tgt.open = true;
        tgt.scrollIntoView({ behavior: "smooth", block: "start" });
      }
    });
  });

  /* ---------- expand/collapse all loci ---------- */
  const expandAll = document.querySelector("#expand-all");
  const collapseAll = document.querySelector("#collapse-all");
  const allLoci = document.querySelectorAll("details.locus");
  expandAll?.addEventListener("click", () => allLoci.forEach((d) => (d.open = true)));
  collapseAll?.addEventListener("click", () => allLoci.forEach((d) => (d.open = false)));

  /* ---------- "Copy hash" buttons in audit page ---------- */
  document.querySelectorAll("[data-copy]").forEach((btn) => {
    btn.addEventListener("click", async () => {
      try {
        await navigator.clipboard.writeText(btn.dataset.copy);
        const t = btn.textContent;
        btn.textContent = "copied";
        setTimeout(() => (btn.textContent = t), 1200);
      } catch (_) {}
    });
  });
})();
