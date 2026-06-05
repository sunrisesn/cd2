"""
CD2-KO vs WT phosphoproteomics — end-to-end downstream analysis.

Implements the protocol in docs/Revised_end-to-end-Phosphoproteomics_Analysis_Protocol.docx.
Outputs are written under results/.
"""
from __future__ import annotations

import re
import sys
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "Filtered_phosphoproteins_CD2-KO_Vs_WT.xlsx"
RES = ROOT / "results"
FIG = RES / "figures"
TAB = RES / "tables"
INP = RES / "inputs"
QC = RES / "qc"
ENR = RES / "enrichment"
KIN = RES / "kinase"
for p in (FIG, TAB, INP, QC, ENR, KIN):
    p.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"pdf.fonttype": 42, "font.size": 9})

COMPARISONS = {
    "30min": "KO_30m_vs_WT_30m",
    "90min": "KO_90m_vs_WT_90m",
}
HIT_SHEETS = {
    "hits_UP_30min":   "30'_KO_Up(P<0.05,Log2FC>=0.88)",
    "hits_DOWN_30min": "30'_KODown(P<0.05,Log2FC<-0.99)",
    "hits_UP_90min":   "90'_KO_Up(P<0.05,Log2FC>=0.88)",
    "hits_DOWN_90min": "90'_KODown(P<0.05,Log2FC<-0.99)",
}
P_THR = 0.05
FC_THR = 0.88


def split_genes(cell) -> list[str]:
    if pd.isna(cell):
        return []
    return [g for g in re.split(r"[;,\s]+", str(cell)) if g and g != "nan"]


def gene_set(frame: pd.DataFrame, col: str = "Gene_Symbol") -> set[str]:
    out: set[str] = set()
    for cell in frame[col].dropna():
        out.update(split_genes(cell))
    return out


def residues(mod_string) -> list[str]:
    out = []
    for block in re.findall(r"\[([^\]]+)\]", str(mod_string)):
        for tok in re.split(r"[;,]", block):
            m = re.match(r"\s*([STY])(\d+)", tok.strip())
            if m:
                out.append(m.group(1) + m.group(2))
    return out


def site_label(row) -> str:
    g = split_genes(row["Gene_Symbol"])
    r = residues(row["Modifications_in_Master_Proteins"])
    return (g[0] + "_" + r[0]) if g and r else (g[0] if g else "")


def site_tag(row) -> str:
    g = split_genes(row["Gene_Symbol"])
    r = residues(row["Modifications_in_Master_Proteins"])
    g = g[0] if g else "?"
    r = r[0] if r else "?"
    return f"p{r} {g}"


def step1_load_and_qc() -> tuple[pd.DataFrame, dict]:
    print(f"[Step 1] Loading master sheet from {DATA.name}")
    df = pd.read_excel(DATA, sheet_name="Results (2)")

    qc = {
        "rows_total": int(len(df)),
        "unique_gene_symbols": int(df["Gene_Symbol"].nunique()),
        "by_timepoint": {},
    }
    for label, comp in COMPARISONS.items():
        fc = "log2fc_" + comp
        p = "p_value_" + comp
        q = "q_value_" + comp
        sub = df.dropna(subset=[fc, p])
        sig = (sub[p] < P_THR) & (sub[fc].abs() > FC_THR)
        qc["by_timepoint"][label] = {
            "quantified": int(len(sub)),
            "significant": int(sig.sum()),
            "down": int((sig & (sub[fc] < 0)).sum()),
            "up": int((sig & (sub[fc] > 0)).sum()),
            "q_under_0_05": int((sub[q] < 0.05).sum()),
        }

    print(json.dumps(qc, indent=2))
    (QC / "qc_summary.json").write_text(json.dumps(qc, indent=2))

    # Expected: 16608, 30min sig=146 (136 down, 10 up), 90min sig=35 (26 down, 9 up),
    # q<0.05 = 0 at both timepoints (by design).
    assert qc["rows_total"] == 16608, "Row count differs from protocol expectation"
    for label in COMPARISONS:
        b = qc["by_timepoint"][label]
        assert b["down"] > b["up"], f"{label}: down should outnumber up — sign may be flipped"
    print("[Step 1] QC checks passed.")
    return df, qc


def step3_volcano(df: pd.DataFrame) -> None:
    print("[Step 3] Volcano plots (Fig 1)")
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.3))
    for ax, (label, comp) in zip(axes, COMPARISONS.items()):
        fc = df["log2fc_" + comp]
        p = df["p_value_" + comp]
        nl = -np.log10(p)
        sig = (p < P_THR) & (fc.abs() > FC_THR)
        ax.scatter(fc[~sig], nl[~sig], s=4, c="#cccccc", alpha=0.4,
                   edgecolor="none", rasterized=True)
        ax.scatter(fc[sig & (fc < 0)], nl[sig & (fc < 0)], s=10, c="#2166ac",
                   label="Down in KO")
        ax.scatter(fc[sig & (fc > 0)], nl[sig & (fc > 0)], s=10, c="#b2182b",
                   label="Up in KO")
        ax.axvline(-FC_THR, ls="--", lw=0.5, c="k")
        ax.axvline(FC_THR, ls="--", lw=0.5, c="k")
        ax.axhline(-np.log10(P_THR), ls="--", lw=0.5, c="k")
        top = df[sig].reindex(
            df[sig]["log2fc_" + comp].abs().sort_values(ascending=False).index
        ).head(8)
        for _, r in top.iterrows():
            ax.annotate(site_label(r),
                        (r["log2fc_" + comp], -np.log10(r["p_value_" + comp])),
                        fontsize=5.5, xytext=(2, 1), textcoords="offset points")
        ax.set_title(f"CD2-KO vs WT, {label} (n={len(df)})", fontsize=9)
        ax.set_xlabel("log2 fold change (KO/WT)")
        ax.set_ylabel("-log10 p-value")
        ax.legend(frameon=False, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "Fig1_volcano.png", dpi=300)
    fig.savefig(FIG / "Fig1_volcano.pdf")
    plt.close(fig)


def step4_count_chart(df: pd.DataFrame) -> dict:
    print("[Step 4] Up/down count chart (Fig 2)")
    counts = {}
    for label, comp in COMPARISONS.items():
        fc = df["log2fc_" + comp]
        p = df["p_value_" + comp]
        sig = (p < P_THR) & (fc.abs() > FC_THR)
        counts[label] = (int((sig & (fc > 0)).sum()), int((sig & (fc < 0)).sum()))

    labels = list(counts)
    x = np.arange(len(labels))
    ups = [counts[l][0] for l in labels]
    dns = [counts[l][1] for l in labels]
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.bar(x, ups, color="#b2182b", label="Up in KO")
    ax.bar(x, [-d for d in dns], color="#2166ac", label="Down in KO")
    for i, (u, d) in enumerate(zip(ups, dns)):
        ax.text(i, u + 1, str(u), ha="center")
        ax.text(i, -d - 3, str(d), ha="center")
    ax.axhline(0, c="k", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Regulated phosphosites")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "Fig2_counts.png", dpi=300)
    fig.savefig(FIG / "Fig2_counts.pdf")
    plt.close(fig)
    return counts


def step5_heatmap(df: pd.DataFrame) -> list[str]:
    print("[Step 5] Recurrent-site heatmap (Fig 3)")

    def gene_to_fc(comp: str) -> dict[str, float]:
        fc = "log2fc_" + comp
        p = "p_value_" + comp
        sub = df[(df[p] < P_THR) & (df[fc].abs() > FC_THR)]
        best: dict[str, float] = {}
        for _, r in sub.iterrows():
            for g in split_genes(r["Gene_Symbol"]):
                if g not in best or abs(r[fc]) > abs(best[g]):
                    best[g] = r[fc]
        return best

    m30 = gene_to_fc(COMPARISONS["30min"])
    m90 = gene_to_fc(COMPARISONS["90min"])
    both = sorted(set(m30) & set(m90), key=lambda g: m30[g])
    mat = np.array([[m30[g], m90[g]] for g in both])

    if len(both) == 0:
        print("  no genes recurrent at both timepoints — skipping heatmap")
        return both

    fig, ax = plt.subplots(figsize=(3.6, max(3, 0.28 * len(both))))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["30 min", "90 min"])
    ax.set_yticks(range(len(both)))
    ax.set_yticklabels(both, fontsize=7)
    cb = fig.colorbar(im, ax=ax, fraction=0.06)
    cb.set_label("log2FC (KO/WT)")
    fig.tight_layout()
    fig.savefig(FIG / "Fig3_heatmap.png", dpi=300)
    fig.savefig(FIG / "Fig3_heatmap.pdf")
    plt.close(fig)

    pd.DataFrame({
        "gene": both,
        "log2FC_30min": [m30[g] for g in both],
        "log2FC_90min": [m90[g] for g in both],
    }).to_csv(TAB / "recurrent_genes_30_90min.csv", index=False)
    return both


def step6_overlap(df: pd.DataFrame) -> dict:
    print("[Step 6] Down-overlap diagram (Fig 4)")
    dn30 = gene_set(pd.read_excel(DATA, sheet_name=HIT_SHEETS["hits_DOWN_30min"]))
    dn90 = gene_set(pd.read_excel(DATA, sheet_name=HIT_SHEETS["hits_DOWN_90min"]))
    only30, only90, shared = len(dn30 - dn90), len(dn90 - dn30), len(dn30 & dn90)

    try:
        from matplotlib_venn import venn2
        fig, ax = plt.subplots(figsize=(4.5, 4))
        venn2(subsets=(only30, only90, shared),
              set_labels=("Down 30 min", "Down 90 min"),
              set_colors=("#2166ac", "#67a9cf"), alpha=0.55, ax=ax)
        fig.tight_layout()
    except Exception:
        fig, ax = plt.subplots(figsize=(4.2, 4))
        ax.add_patch(Circle((0.37, 0.5), 0.33, alpha=0.45, color="#2166ac"))
        ax.add_patch(Circle((0.63, 0.5), 0.33, alpha=0.45, color="#67a9cf"))
        ax.text(0.22, 0.5, only30, ha="center", va="center", fontsize=13)
        ax.text(0.78, 0.5, only90, ha="center", va="center", fontsize=13)
        ax.text(0.5, 0.5, shared, ha="center", va="center", fontsize=13, fontweight="bold")
        ax.text(0.22, 0.88, "Down 30 min", ha="center")
        ax.text(0.78, 0.88, "Down 90 min", ha="center")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        fig.tight_layout()

    fig.savefig(FIG / "Fig4_venn.png", dpi=300)
    fig.savefig(FIG / "Fig4_venn.pdf")
    plt.close(fig)

    pd.DataFrame({
        "category": ["only_30min", "only_90min", "shared"],
        "count": [only30, only90, shared],
    }).to_csv(TAB / "down_overlap_counts.csv", index=False)

    return {"only_30min": only30, "only_90min": only90, "shared": shared,
            "shared_genes": sorted(dn30 & dn90)}


def step7_export_inputs(df: pd.DataFrame) -> None:
    print("[Step 7] Exporting gene-list / kinase-library inputs")

    bg = sorted(gene_set(df))
    (INP / "background_all_quantified_genes.txt").write_text("\n".join(bg))

    for out_name, sheet in HIT_SHEETS.items():
        genes = sorted(gene_set(pd.read_excel(DATA, sheet_name=sheet)))
        (INP / f"{out_name}.txt").write_text("\n".join(genes))

    for label, comp in COMPARISONS.items():
        fc = "log2fc_" + comp
        rows: list[tuple[str, str, float]] = []
        for _, r in df.dropna(subset=[fc]).iterrows():
            genes = split_genes(r["Gene_Symbol"])
            res = residues(r["Modifications_in_Master_Proteins"])
            for g in genes:
                for site in res:
                    rows.append((g, site, float(r[fc])))
        out = pd.DataFrame(rows, columns=["gene", "site", "log2fc"])
        st = out[out["site"].str[0].isin(["S", "T"])]
        ty = out[out["site"].str[0] == "Y"]
        st.to_csv(INP / f"KL_ST_{label}.txt", sep="\t", header=False, index=False)
        ty.to_csv(INP / f"KL_Y_{label}.txt", sep="\t", header=False, index=False)

        ksea = []
        for _, r in df.dropna(subset=[fc, "p_value_" + comp]).iterrows():
            genes = split_genes(r["Gene_Symbol"])
            res = residues(r["Modifications_in_Master_Proteins"])
            acc = split_genes(r["Accession"])
            prot = acc[0] if acc else ""
            for g in genes:
                for site in res:
                    ksea.append((prot, g, site, float(r["p_value_" + comp]), float(r[fc])))
        pd.DataFrame(ksea, columns=["Protein", "Gene", "Residue.Both", "p", "FC"]).to_csv(
            INP / f"KSEA_full_input_{label}.csv", index=False
        )

        rnk = out.assign(label=out["gene"] + "_" + out["site"])[["label", "log2fc"]]
        rnk = rnk.groupby("label", as_index=False)["log2fc"].mean()
        rnk.to_csv(INP / f"GSEA_ranked_{label}.rnk", sep="\t", header=False, index=False)


def step13_1_top50(df: pd.DataFrame) -> None:
    print("[Step 13.1] Top-50 most-decreased table/grid")
    comp = COMPARISONS["30min"]
    fc = "log2fc_" + comp
    p = "p_value_" + comp
    sig = df[(df[p] < P_THR) & (df[fc].abs() > FC_THR)].copy()
    top = sig.sort_values(fc).head(50)
    top["label"] = top.apply(site_tag, axis=1)
    top[["label", fc, p]].to_csv(TAB / "Fig3_top50_table.csv", index=False)

    labels = list(top["label"])
    ncol = 5
    rows = [labels[i:i + ncol] for i in range(0, len(labels), ncol)]
    fig, ax = plt.subplots(figsize=(8.5, 0.32 * len(rows) + 0.6))
    ax.axis("off")
    tbl = ax.table(cellText=rows, cellLoc="left", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1, 1.3)
    ax.set_title("50 phosphosites most decreased in CD2-KO (30 min)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "Fig3_top50_grid.png", dpi=300)
    fig.savefig(FIG / "Fig3_top50_grid.pdf")
    plt.close(fig)


def step13_2_string_table(df: pd.DataFrame) -> None:
    print("[Step 13.2] STRING node-color table")
    for label, comp in COMPARISONS.items():
        fc = "log2fc_" + comp
        p = "p_value_" + comp
        sig = df[(df[p] < P_THR) & (df[fc].abs() > FC_THR)]
        best: dict[str, float] = {}
        for _, r in sig.iterrows():
            for g in split_genes(r["Gene_Symbol"]):
                if g not in best or abs(r[fc]) > abs(best[g]):
                    best[g] = float(r[fc])
        pd.DataFrame({"gene": list(best), "log2FC": list(best.values())}).to_csv(
            TAB / f"STRING_node_colors_{label}.csv", index=False
        )


def step8_enrichment() -> None:
    print("[Step 8] GO/KEGG enrichment via Enrichr (gseapy)")
    try:
        import gseapy as gp
    except ImportError:
        print("  gseapy not available — skipping; upload hit lists to g:Profiler instead.")
        return

    libraries = [
        "GO_Biological_Process_2023",
        "GO_Molecular_Function_2023",
        "GO_Cellular_Component_2023",
        "KEGG_2021_Human",
        "Reactome_2022",
    ]
    for hit_key in HIT_SHEETS:
        gene_file = INP / f"{hit_key}.txt"
        genes = [g for g in gene_file.read_text().splitlines() if g.strip()]
        if not genes:
            print(f"  {hit_key}: empty list — skipping")
            continue
        print(f"  {hit_key}: {len(genes)} genes -> Enrichr ({len(libraries)} libs)")
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=libraries,
                organism="human",
                outdir=None,
                no_plot=True,
            )
        except Exception as e:
            print(f"  Enrichr call failed for {hit_key}: {e}")
            continue
        res = enr.results
        if res is None or len(res) == 0:
            print(f"  {hit_key}: no results returned")
            continue
        res = res.sort_values("Adjusted P-value")
        res.to_csv(ENR / f"enrichment_{hit_key}.csv", index=False)

        plot_df = res.copy()
        plot_df = plot_df.sort_values("Adjusted P-value").head(15).iloc[::-1]
        if len(plot_df) == 0:
            continue
        plot_df["neglogP"] = -np.log10(plot_df["Adjusted P-value"].clip(lower=1e-300))
        plot_df["Overlap_size"] = plot_df["Overlap"].apply(
            lambda s: int(str(s).split("/")[0]) if "/" in str(s) else np.nan
        )

        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        sc = ax.scatter(plot_df["neglogP"], range(len(plot_df)),
                        s=plot_df["Overlap_size"].fillna(5) * 18,
                        c=plot_df["neglogP"], cmap="plasma",
                        edgecolor="k", linewidth=0.4)
        ax.set_yticks(range(len(plot_df)))
        terms = [f"[{src.split('_')[0]}] {t}" for src, t in
                 zip(plot_df["Gene_set"], plot_df["Term"])]
        terms = [t if len(t) < 75 else t[:72] + "..." for t in terms]
        ax.set_yticklabels(terms, fontsize=7.5)
        ax.set_xlabel("-log10 adjusted p-value (BH)")
        fig.colorbar(sc, label="-log10 FDR")
        ax.set_title(f"GO / KEGG / Reactome enrichment — {hit_key}", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG / f"Fig5_enrichment_{hit_key}.png", dpi=300)
        fig.savefig(FIG / f"Fig5_enrichment_{hit_key}.pdf")
        plt.close(fig)


def step9_kinase_enrichment() -> None:
    print("[Step 9] Kinase substrate enrichment via Enrichr (KEA)")
    try:
        import gseapy as gp
    except ImportError:
        print("  gseapy not available — skipping")
        return

    kinase_libs = ["KEA_2015", "The_Kinase_Library_2024", "PPI_Hub_Proteins"]
    for hit_key in ("hits_DOWN_30min", "hits_DOWN_90min",
                    "hits_UP_30min", "hits_UP_90min"):
        gene_file = INP / f"{hit_key}.txt"
        genes = [g for g in gene_file.read_text().splitlines() if g.strip()]
        if not genes:
            continue
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=kinase_libs,
                organism="human",
                outdir=None,
                no_plot=True,
            )
        except Exception as e:
            print(f"  KEA call failed for {hit_key}: {e}")
            continue
        res = enr.results
        if res is None or len(res) == 0:
            continue
        res = res.sort_values("Adjusted P-value")
        res.to_csv(KIN / f"kinase_enrichment_{hit_key}.csv", index=False)

        for lib_name, fig_tag in (("The_Kinase_Library_2024", "KinaseLib"),
                                  ("KEA_2015", "KEA")):
            sub = res[res["Gene_set"] == lib_name].head(15).iloc[::-1]
            if len(sub) == 0:
                continue
            sub = sub.copy()
            sub["neglogP"] = -np.log10(sub["Adjusted P-value"].clip(lower=1e-300))
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.barh(range(len(sub)), sub["neglogP"], color="#3690c0")
            ax.set_yticks(range(len(sub)))
            ax.set_yticklabels(sub["Term"], fontsize=8)
            ax.set_xlabel("-log10 adjusted p-value (BH)")
            ax.set_title(f"Kinase enrichment ({lib_name}) — {hit_key}", fontsize=10)
            fig.tight_layout()
            fig.savefig(FIG / f"Fig6_kinase_{fig_tag}_{hit_key}.png", dpi=300)
            fig.savefig(FIG / f"Fig6_kinase_{fig_tag}_{hit_key}.pdf")
            plt.close(fig)


def step13_3_regulated_kinases(df: pd.DataFrame) -> None:
    """Sec 13.3 — extract regulated kinase genes for Cytoscape kinome panel."""
    print("[Step 13.3] Regulated-kinase gene list (data-derived layer)")
    try:
        import gseapy as gp
    except ImportError:
        return

    try:
        kinome = gp.get_library("KEA_2015", organism="human")
    except Exception:
        kinome = {}
    kinome_genes = set(kinome.keys()) if kinome else set()

    for label, comp in COMPARISONS.items():
        fc = "log2fc_" + comp
        p = "p_value_" + comp
        sig = df[(df[p] < P_THR) & (df[fc].abs() > FC_THR)]
        best: dict[str, float] = {}
        for _, r in sig.iterrows():
            for g in split_genes(r["Gene_Symbol"]):
                if g not in best or abs(r[fc]) > abs(best[g]):
                    best[g] = float(r[fc])
        if kinome_genes:
            kinases = {g: v for g, v in best.items() if g in kinome_genes}
        else:
            kinases = {}
        (INP / f"kinases_regulated_{label}.txt").write_text("\n".join(sorted(kinases)))
        pd.DataFrame(
            {"gene": list(kinases), "log2FC": list(kinases.values())}
        ).to_csv(TAB / f"kinases_regulated_{label}.csv", index=False)


def write_final_summary(qc: dict, counts: dict, recurrent: list[str], overlap: dict) -> None:
    summary = {
        "qc": qc,
        "fig2_counts_up_down": {k: {"up": u, "down": d} for k, (u, d) in counts.items()},
        "fig3_recurrent_genes_both_timepoints": recurrent,
        "fig4_down_overlap": overlap,
    }
    (QC / "summary.json").write_text(json.dumps(summary, indent=2))
    print("\n=== Final summary ===")
    print(json.dumps(summary, indent=2))


def main() -> None:
    df, qc = step1_load_and_qc()
    step3_volcano(df)
    counts = step4_count_chart(df)
    recurrent = step5_heatmap(df)
    overlap = step6_overlap(df)
    step7_export_inputs(df)
    step13_1_top50(df)
    step13_2_string_table(df)
    step8_enrichment()
    step9_kinase_enrichment()
    step13_3_regulated_kinases(df)
    write_final_summary(qc, counts, recurrent, overlap)
    print("\nAll outputs written under results/.")


if __name__ == "__main__":
    sys.exit(main())
