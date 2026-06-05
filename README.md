# CD2-KO vs WT Phosphoproteomics — Downstream Analysis

End-to-end Python pipeline that turns the quantified phosphosite table from
the CD2 knockout vs wild-type regulatory-T-cell study into publication-ready
figures and result tables. Implements the protocol in
[`docs/Revised_end-to-end-Phosphoproteomics_Analysis_Protocol.docx`](docs/Revised_end-to-end-Phosphoproteomics_Analysis_Protocol.docx).

## Repository layout

```
.
├── data/                  # Input Excel (master + four hit sheets)
├── docs/
│   ├── Revised_end-to-end-Phosphoproteomics_Analysis_Protocol.docx
│   └── Steps_Executed.md  # Step-by-step record of what was run
├── sample-publication/    # Reference paper (Zurli et al., Sci. Signal. 2020)
├── src/
│   └── analysis.py        # Single-script pipeline
├── results/               # All generated outputs (figures, tables, inputs, QC)
├── requirements.txt
└── README.md
```

## Quick start

```bash
# 1. Create env (Python 3.10+; tested on 3.11)
python3.11 -m venv .venv
source .venv/bin/activate
pip install --index-url https://pypi.org/simple/ -r requirements.txt

# 2. Run the full pipeline
python src/analysis.py
```

Run time is ~1–2 minutes. Enrichr (used for GO/KEGG/Reactome and kinase
enrichment) is the only external call.

## What gets produced

| Deliverable | File |
| --- | --- |
| Volcano plots (30 + 90 min) | [results/figures/Fig1_volcano.png](results/figures/Fig1_volcano.png) |
| Up/down count chart | [results/figures/Fig2_counts.png](results/figures/Fig2_counts.png) |
| Recurrent-site heatmap | [results/figures/Fig3_heatmap.png](results/figures/Fig3_heatmap.png) |
| Top-50 most-decreased grid | [results/figures/Fig3_top50_grid.png](results/figures/Fig3_top50_grid.png) |
| Down-overlap Venn | [results/figures/Fig4_venn.png](results/figures/Fig4_venn.png) |
| GO/KEGG/Reactome dot plots | [results/figures/](results/figures/) (`Fig5_*`) |
| Kinase enrichment bar plots | [results/figures/](results/figures/) (`Fig6_*`) |
| Gene lists, background, kinase-library inputs | [results/inputs/](results/inputs/) |
| Enrichment tables (BH-adjusted) | [results/enrichment/](results/enrichment/) |
| Kinase enrichment tables | [results/kinase/](results/kinase/) |
| QC + final summary | [results/qc/](results/qc/) |

A guided walk-through of every step is in
[`docs/Steps_Executed.md`](docs/Steps_Executed.md).

## Key results

* **16,608** phosphosites quantified across 4,982 unique gene symbols.
* Significance threshold (protocol-defined): `p < 0.05 AND |log2FC| > 0.88`.

| Timepoint | Significant | Down in KO | Up in KO | q<0.05 |
| --- | --- | --- | --- | --- |
| 30 min | 146 | 136 | 10 | 0 |
| 90 min |  35 |  26 |  9 | 0 |

* **7 genes** are recurrently down at both timepoints:
  `CD2, CEP170, JUN, NOP2, RPS6KA1, RPS6KA3, XIRP1`.
* Top down-enriched pathways (30 min): **cytoskeleton / actin**, **mTOR signaling**,
  **mTORC1-mediated signaling**, **RSK activation**, **nuclear pore complex
  disassembly** — consistent with the Zurli et al. (2020) reference paper.
* Top kinases inferred by **The Kinase Library 2024** (down @ 30 min):
  `MAPKAPK2, GAK, BIKE, ERK5, MELK, PBK, DNAPK, MST4, DRAK1, CHAK2, CAMK1G,
  ERK2, STK33, AAK1, IRAK4`.

The 0 q<0.05 result is **expected by experimental design** (anti-CD3/CD28
co-stimulation compresses the KO-vs-WT differential); the protocol uses
nominal p throughout and reports hits as hypothesis-generating.

## Steps that still need a human / web tool

The protocol explicitly delegates these to interactive tools (database
freshness, hand-curated layouts). Pre-formatted input files for each are
already in [results/inputs/](results/inputs/):

| Step | Tool | Input file |
| --- | --- | --- |
| g:Profiler (if reviewer requires the official tool) | https://biit.cs.ut.ee/gprofiler/gost | `hits_*.txt` + `background_all_quantified_genes.txt` |
| Kinase Library web app (if reviewer requires it) | https://kinase-library.phosphosite.org/ | `KL_ST_*.txt`, `KL_Y_*.txt` |
| KSEA web app | https://casecpb.shinyapps.io/ksea/ | `KSEA_full_input_*.csv` |
| STRING network (Fig 4A/C analogue) | https://string-db.org → Cytoscape stringApp | `hits_DOWN_30min.txt`, `STRING_node_colors_30min.csv` |
| Kinome map (Fig 5C analogue) | STRING + Cytoscape | `kinases_regulated_*.txt` |
| Western blot panels (Fig 5A/B) | **Bench experiment** — out of scope | — |

## Reproducibility notes

* The pipeline is one file (`src/analysis.py`); each step is a function
  whose name matches the protocol section. Re-running overwrites results.
* Acceptance values (rows, hit counts) are asserted at load time; a sign
  flip or mis-mapped sheet will fail loudly rather than silently.
* Plots are saved at 300 dpi PNG **and** vector PDF (`pdf.fonttype=42`,
  editable in Illustrator).
* Enrichment is currently performed via Enrichr (consistent libraries,
  reproducible). The protocol's reference tool g:Profiler can be used
  instead — input files are already in the right format.

## License & citation

The reference paper used for figure analogues is:
Zurli V. *et al.* "Phosphoproteomics of CD2 signaling reveals AMPK-mediated
regulation of lytic granule release." *Sci. Signal.* (2020).
See [`sample-publication/`](sample-publication/).
