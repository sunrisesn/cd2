# Steps Executed — CD2-KO vs WT Phosphoproteomics Analysis

This document records exactly what was executed to convert the input table
`data/Filtered_phosphoproteins_CD2-KO_Vs_WT.xlsx` into the deliverables under
`results/`, following the protocol in
`docs/Revised_end-to-end-Phosphoproteomics_Analysis_Protocol.docx`.

Every step is reproducible by running `python src/analysis.py` from the
project root with the venv active.

---

## 0. Environment

| Item | Value |
| --- | --- |
| Python | 3.11 (venv at `.venv/`) |
| Index | `https://pypi.org/simple/` (overrides codeartifact) |
| Deps | `pandas`, `numpy`, `scipy`, `openpyxl`, `matplotlib`, `seaborn`, `python-docx`, `pypdf`, `matplotlib-venn`, `gseapy` |
| Enrichment backend | Enrichr REST API via `gseapy` |
| Working directory | `/Users/arazak/projects/sunrisesn/cd2` |

Commands executed:

```bash
rm -rf .venv && python3.11 -m venv .venv
source .venv/bin/activate
pip install --index-url https://pypi.org/simple/ -r requirements.txt
```

---

## 1. Load Data and QC (Protocol §4)

Loaded sheet `Results (2)` (master) and applied the protocol's acceptance
test: significance is `p < 0.05 AND |log2FC| > 0.88`; comparisons are
`KO_30m_vs_WT_30m` and `KO_90m_vs_WT_90m`.

Acceptance values produced (match the protocol's expected table exactly):

| Check | Expected | Observed |
| --- | --- | --- |
| Total rows | 16,608 | **16,608** |
| 30 min — significant (down / up) | 146 (136 / 10) | **146 (136 / 10)** |
| 90 min — significant (down / up) | 35 (26 / 9) | **35 (26 / 9)** |
| Sites passing q < 0.05 | 0 (by design) | **0** |
| Sanity: down ≫ up | yes | yes (asserted) |

Stored as `results/qc/qc_summary.json` and `results/qc/summary.json`.

> The 0 q<0.05 result is expected by experimental design (anti-CD3/CD28
> co-stimulation compresses the KO-vs-WT differential). The pipeline reports
> nominal p-values throughout, as the protocol requires.

---

## 2. Helper Functions (Protocol §5)

Implemented in `src/analysis.py`:

* `split_genes(cell)` — splits `"CD2;ARHGAP12"` into `["CD2","ARHGAP12"]`.
* `gene_set(frame)` — flat union of gene symbols across a frame.
* `residues(mod_string)` — extracts e.g. `["S345","T8360"]` from
  `"P06729 1xPhospho [S345]"`.
* `site_label(row)` / `site_tag(row)` — single-site labels for figure
  annotations (`GENE_S123` and `pS123 GENE` respectively).

---

## 3. Figure 1 — Volcano Plots (Protocol §6)

Plotted **all 16,608 sites** (not just the hit lists) for both timepoints:

* X = `log2fc_KO_*m_vs_WT_*m`, Y = `-log10(p_value_KO_*m_vs_WT_*m)`.
* Greyed-out non-significant; blue = down, red = up.
* Dashed cutoff lines at `±0.88` and `-log10(0.05)`.
* The 8 strongest hits are annotated as `GENE_SITE`.

Output: `results/figures/Fig1_volcano.{png,pdf}`.

Acceptance: the 30 min panel shows a dense blue cloud on the left, very few
red dots — matches the protocol's description.

---

## 4. Figure 2 — Up/Down Counts (Protocol §7)

Single bar chart with up (red) and down (blue, plotted as negative) per
timepoint. Reproduces the **30: 10 / 136; 90: 9 / 26** counts.

Output: `results/figures/Fig2_counts.{png,pdf}`.

---

## 5. Figure 3 — Recurrent-Site Heatmap (Protocol §8)

For each timepoint, the strongest-changing site of every gene (by `|log2FC|`)
is selected from the significant set; the intersection across 30 / 90 min is
plotted as a 2-column RdBu_r heatmap.

12 recurrent genes:
`FYCO1, XIRP1, CD2, NUP153, RPS6KA3, RPS6KA1, CEP170, NOP2, JUN, JUND, SYNE1, VSIR`.

Outputs:
* `results/figures/Fig3_heatmap.{png,pdf}`
* `results/tables/recurrent_genes_30_90min.csv`

---

## 6. Figure 4 — Overlap Diagram (Protocol §9)

Down-regulated gene sets at 30 / 90 min compared. Counts: **only-30=79,
only-90=17, shared=7** (`CD2, CEP170, JUN, NOP2, RPS6KA1, RPS6KA3, XIRP1`).
Drawn as a proportional Venn via `matplotlib-venn` (with a circle-based
fallback if the package is unavailable).

Outputs:
* `results/figures/Fig4_venn.{png,pdf}`
* `results/tables/down_overlap_counts.csv`

---

## 7. Web-Tool Input Files (Protocol §10)

Written under `results/inputs/`:

| File | Purpose | Rows |
| --- | --- | --- |
| `background_all_quantified_genes.txt` | g:Profiler custom background | 4982 unique genes |
| `hits_DOWN_30min.txt` | DOWN gene list, 30 min | 86 |
| `hits_UP_30min.txt` | UP gene list, 30 min | 11 |
| `hits_DOWN_90min.txt` | DOWN gene list, 90 min | 24 |
| `hits_UP_90min.txt` | UP gene list, 90 min | 9 |
| `KL_ST_30min.txt` / `KL_ST_90min.txt` | Kinase Library Ser/Thr input (TSV, **no header**) | ~16.6k rows each |
| `KL_Y_30min.txt` / `KL_Y_90min.txt` | Kinase Library Tyr input (TSV, **no header**) | small |
| `KSEA_full_input_30min.csv` / `_90min.csv` | KSEA app input (`Protein,Gene,Residue.Both,p,FC`) | ~16.6k rows each |
| `GSEA_ranked_30min.rnk` / `_90min.rnk` | preranked GSEA (gene_site, log2FC) | ~16.6k rows |

These match the exact formats the web tools accept (3-column TSV, no header,
per the protocol's "Common upload error" note).

---

## 8. Figure 5 — GO/KEGG/Reactome Enrichment (Protocol §11)

The protocol recommends g:Profiler as a web tool. To make the pipeline
fully self-contained, the same kind of enrichment is run programmatically
via **Enrichr** (REST API through `gseapy`) for each hit list against:

* `GO_Biological_Process_2023`
* `GO_Molecular_Function_2023`
* `GO_Cellular_Component_2023`
* `KEGG_2021_Human`
* `Reactome_2022`

Per hit list the top-15 terms (BH-adjusted) are plotted as a dot plot
(`Fig5_enrichment_<hit_list>.{png,pdf}`); full result tables go to
`results/enrichment/enrichment_<hit_list>.csv`.

Top biology for `hits_DOWN_30min` matches the published CD2 story:

* GO:Cellular_Component — **Cytoskeleton**, **Actin Cytoskeleton**
* KEGG — **mTOR signaling pathway**
* Reactome — **MTOR Signaling**, **mTORC1-mediated Signaling**,
  **RSK Activation**, **Nuclear Pore Complex Disassembly**, **Signal Transduction**

> **Reviewers asking for g:Profiler specifically** can upload the input files
> at https://biit.cs.ut.ee/gprofiler/gost — same lists, custom background
> already provided. The plot/table code in `step8_enrichment()` will work
> against a `gProfiler_*.csv` export as well (rename columns to
> `Term`, `Adjusted P-value`, `Gene_set`, `Overlap` if needed).

---

## 9. Figure 6 — Kinase Activity Enrichment (Protocol §12)

The protocol's recommended tool is **The Kinase Library** (Cantley lab).
That dataset is available via Enrichr as `The_Kinase_Library_2024`, so the
pipeline runs it programmatically alongside the classic `KEA_2015` (kinase-
substrate enrichment from PhosphoSitePlus) and `PPI_Hub_Proteins`.

Per hit list:

* `results/kinase/kinase_enrichment_<hit_list>.csv` — full table.
* `results/figures/Fig6_kinase_KinaseLib_<hit_list>.{png,pdf}` — top 15
  kinases from **The Kinase Library 2024**.
* `results/figures/Fig6_kinase_KEA_<hit_list>.{png,pdf}` — top 15 kinases
  from `KEA_2015`.

Top kinases for DOWN @ 30 min (Kinase Library 2024):
`MAPKAPK2, GAK, BIKE, ERK5, MELK, PBK, DNAPK, MST4, DRAK1, CHAK2, CAMK1G, ERK2, STK33, AAK1, IRAK4`.

The presence of MAPKAPK2 / ERK2 / ERK5 / IRAK4 / MELK / CAMK1G is consistent
with the loss of CD2-driven MAPK and mTOR-arm signaling reported in
Zurli et al. (2020).

> If a reviewer requires the official web-app run, upload
> `results/inputs/KL_ST_30min.txt` (and `_90min.txt`) to
> https://kinase-library.phosphosite.org/ — they are pre-formatted (TSV,
> 3 columns, no header), so they will be accepted directly.

---

## 10. Reference-Paper Figures (Protocol §13)

### 13.1 — Top-50 most-decreased table (Fig 3 analogue)
Sites in the significant set, sorted by ascending `log2FC` (most decreased
in KO = CD2-dependent), labelled as `pRESIDUE GENE`.

* `results/tables/Fig3_top50_table.csv`
* `results/figures/Fig3_top50_grid.{png,pdf}` — 5-column rendered grid.

### 13.2 — STRING node-color table (Fig 4 analogue)
Per-gene strongest-changing `log2FC` for both timepoints, ready to import
into Cytoscape's stringApp as a node attribute. Generated for both
timepoints:
`results/tables/STRING_node_colors_{30min,90min}.csv`.

> The STRING network drawing itself is performed at https://string-db.org
> (Multiple proteins → paste `results/inputs/hits_DOWN_30min.txt` → organism
> Homo sapiens) and exported into Cytoscape with the stringApp. The node
> colouring then uses the CSV above. This step is web-only per the
> protocol.

### 13.3 — Regulated-kinase list (Fig 5C analogue)
Genes from the significant hits that also appear in the human kinome are
written to:

* `results/inputs/kinases_regulated_{30min,90min}.txt`
* `results/tables/kinases_regulated_{30min,90min}.csv` (with `log2FC`)

The human kinome reference is pulled live from Enrichr's `KEA_2015`
library so the pipeline does not need a bundled kinome file.

> The Western-blot panels (5A / 5B) are bench experiments and are
> intentionally out of scope, per the protocol.

---

## 11. Final Deliverables

```
results/
├── figures/
│   ├── Fig1_volcano.{png,pdf}
│   ├── Fig2_counts.{png,pdf}
│   ├── Fig3_heatmap.{png,pdf}
│   ├── Fig3_top50_grid.{png,pdf}
│   ├── Fig4_venn.{png,pdf}
│   ├── Fig5_enrichment_<hit_list>.{png,pdf}            # 4 lists × 2 fmt
│   ├── Fig6_kinase_KinaseLib_<hit_list>.{png,pdf}      # 4 × 2
│   └── Fig6_kinase_KEA_<hit_list>.{png,pdf}            # 4 × 2
├── tables/
│   ├── Fig3_top50_table.csv
│   ├── STRING_node_colors_{30min,90min}.csv
│   ├── recurrent_genes_30_90min.csv
│   ├── down_overlap_counts.csv
│   └── kinases_regulated_{30min,90min}.csv
├── inputs/
│   ├── background_all_quantified_genes.txt
│   ├── hits_{UP,DOWN}_{30min,90min}.txt
│   ├── KL_ST_{30min,90min}.txt        # Kinase Library Ser/Thr
│   ├── KL_Y_{30min,90min}.txt         # Kinase Library Tyr
│   ├── KSEA_full_input_{30min,90min}.csv
│   ├── GSEA_ranked_{30min,90min}.rnk
│   └── kinases_regulated_{30min,90min}.txt
├── enrichment/
│   └── enrichment_<hit_list>.csv                       # 4 lists
├── kinase/
│   └── kinase_enrichment_<hit_list>.csv                # 4 lists
└── qc/
    ├── qc_summary.json
    └── summary.json
```

All counts reproduce the protocol's acceptance test (Protocol §4 / §10 / §11
expected tables).

---

## 12. Pitfalls Avoided

* **Sign convention** — comparisons are `KO_*m_vs_WT_*m`: negative log2FC =
  lower in KO. Down clearly outnumbers up at both timepoints. (Sanity rule.)
* **Background** — enrichment uses the full quantified-gene set explicitly
  (custom background file written for g:Profiler).
* **Whole-distribution kinase input** — Kinase Library / KSEA inputs are
  built from *all* sites, not just the hit lists.
* **Tab-separated, no header** — `KL_*.txt` files have `sep='\t'` and
  `header=False` to avoid the Kinase Library upload error.
* **Multi-gene rows** — every row is split on `[;,\s]+` before going into
  any gene-set output.
* **Ser/Thr vs Tyr** — split into separate `KL_ST_*.txt` / `KL_Y_*.txt`
  files since the two kinomes are scored separately.
* **FDR not used as a gate** — q<0.05 = 0 at both timepoints is reported in
  QC but does not stop the pipeline; everything proceeds on nominal p.

---

## 13. Reproducing

```bash
source .venv/bin/activate
python src/analysis.py
```

End-to-end run time: ~1–2 minutes (Enrichr roundtrips dominate).
The pipeline is idempotent — re-running overwrites the existing outputs.
