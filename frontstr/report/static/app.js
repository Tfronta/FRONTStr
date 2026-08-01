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
  const flaggedSelect = document.querySelector("#profile-flagged");
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
      const flagged = flaggedSelect?.value || "";
      for (const tr of allRows) {
        const text = tr.dataset.search || "";
        const chip = tr.dataset.chip || "";
        const matchesText = !q || text.includes(q);
        const matchesStatus = !status || chip === status;
        // The search box already matches flag codes via data-search, so typing
        // "allele_imbalance" filters too; this select is the no-typing path.
        const matchesFlagged = !flagged || (tr.dataset.flagged || "0") === flagged;
        tr.style.display = matchesText && matchesStatus && matchesFlagged ? "" : "none";
      }
      updatePill();
    };
    searchInput?.addEventListener("input", applyFilters);
    statusSelect?.addEventListener("change", applyFilters);
    flaggedSelect?.addEventListener("change", applyFilters);

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

  /* ---------- any table marked data-sortable-table ----------
   * The block above binds to the first table.profile on the page, which is the
   * CE profile. Tables added later need sorting without inheriting that one's
   * filters, so they opt in by attribute. */
  document.querySelectorAll("table[data-sortable-table]").forEach((table) => {
    const tbody = table.tBodies[0];
    if (!tbody || !table.tHead) return;
    const rows = Array.from(tbody.rows);
    const headers = table.tHead.querySelectorAll("th[data-sortable]");
    headers.forEach((th, idx) => {
      th.addEventListener("click", () => {
        const dir = th.dataset.sortDir === "asc" ? "desc" : "asc";
        headers.forEach((h) => h.removeAttribute("data-sort-dir"));
        th.dataset.sortDir = dir;
        const numeric = th.dataset.sortable === "num";
        const sorted = [...rows].sort((a, b) => {
          let av = a.cells[idx]?.dataset.sortValue ?? a.cells[idx]?.textContent?.trim() ?? "";
          let bv = b.cells[idx]?.dataset.sortValue ?? b.cells[idx]?.textContent?.trim() ?? "";
          if (numeric) {
            av = parseFloat(av);
            bv = parseFloat(bv);
            if (Number.isNaN(av)) av = -Infinity;
            if (Number.isNaN(bv)) bv = -Infinity;
          }
          if (av < bv) return dir === "asc" ? -1 : 1;
          if (av > bv) return dir === "asc" ? 1 : -1;
          return 0;
        });
        for (const r of sorted) tbody.appendChild(r);
      });
    });
  });

  /* ---------- sequencing table: filter + sort + sync from CE table ---------- */
  const seqTable = document.querySelector("table.seqtable");
  const seqSearch = document.querySelector("#seq-search");
  const seqClear = document.querySelector("#seq-clear");
  const seqCount = document.querySelector("#seq-count");

  if (seqTable) {
    const sbody = seqTable.tBodies[0];
    const seqRows = Array.from(sbody.rows);

    const updateSeqPill = () => {
      if (seqCount) {
        const shown = seqRows.filter((r) => r.style.display !== "none").length;
        seqCount.textContent = `${shown} / ${seqRows.length} alleles`;
      }
    };
    const applySeqFilter = () => {
      const q = (seqSearch?.value || "").trim().toLowerCase();
      for (const tr of seqRows) {
        const text = tr.dataset.search || "";
        tr.style.display = !q || text.includes(q) ? "" : "none";
      }
      updateSeqPill();
    };
    seqSearch?.addEventListener("input", applySeqFilter);
    seqClear?.addEventListener("click", () => {
      if (seqSearch) seqSearch.value = "";
      applySeqFilter();
      seqSearch?.focus();
    });

    const seqHeaders = seqTable.tHead.querySelectorAll("th[data-sortable]");
    seqHeaders.forEach((th, idx) => {
      th.addEventListener("click", () => {
        const dir = th.dataset.sortDir === "asc" ? "desc" : "asc";
        seqHeaders.forEach((h) => h.removeAttribute("data-sort-dir"));
        th.dataset.sortDir = dir;
        const type = th.dataset.sortable;
        const sorted = [...seqRows].sort((a, b) => {
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
        for (const r of sorted) sbody.appendChild(r);
      });
    });

    updateSeqPill();

    /* click a marker in the CE table → filter this table to that marker */
    document.querySelectorAll("a.marker-link[data-marker]").forEach((a) => {
      a.addEventListener("click", (e) => {
        e.preventDefault();
        if (seqSearch) {
          seqSearch.value = a.dataset.marker || "";
          applySeqFilter();
        }
        document.getElementById("sequences")?.scrollIntoView({ behavior: "smooth", block: "start" });
      });
    });
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

  /* ---------- NGS panel: sync table rows ↔ stacked chart segments ---------- */
  document.querySelectorAll(".ngs-panel").forEach((panel) => {
    const rows = panel.querySelectorAll("tr.ngs-row[data-row-id]");
    const segments = panel.querySelectorAll("rect.ngs-segment[data-row-id]");

    const clear = () => {
      panel.querySelectorAll(".is-selected").forEach((el) => el.classList.remove("is-selected"));
    };

    const selectRowId = (id) => {
      clear();
      rows.forEach((tr) => {
        if (tr.dataset.rowId === id) tr.classList.add("is-selected");
      });
      segments.forEach((seg) => {
        if (seg.dataset.rowId === id) seg.classList.add("is-selected");
      });
    };

    segments.forEach((seg) => {
      seg.addEventListener("click", () => selectRowId(seg.dataset.rowId || ""));
      seg.addEventListener("keydown", (e) => {
        if (e.key === "Enter" || e.key === " ") {
          e.preventDefault();
          selectRowId(seg.dataset.rowId || "");
        }
      });
    });

    rows.forEach((tr) => {
      tr.addEventListener("click", () => selectRowId(tr.dataset.rowId || ""));
      tr.addEventListener("keydown", (e) => {
        if (e.key === "Enter" || e.key === " ") {
          e.preventDefault();
          selectRowId(tr.dataset.rowId || "");
        }
      });
    });
  });
})();
