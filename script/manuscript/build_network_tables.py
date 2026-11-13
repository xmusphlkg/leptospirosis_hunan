"""Publish tables from the sentinel-network evaluation without changing analyses."""
from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "output/hunan/current/network_redesign_2026-11-13"
TABLES = ROOT / "manuscript/tables"
DATA = ROOT / "data/processed/harmonized_2026-10-22"
SITE_NAMES = {
    "君山区": "Junshan", "涟源市": "Lianyuan", "醴陵市": "Liling",
    "浏阳市": "Liuyang", "宁乡县": "Ningxiang", "宁乡市": "Ningxiang",
    "双峰县": "Shuangfeng", "湘潭县": "Xiangtan", "沅江市": "Yuanjiang",
}
RULE_NAMES = {
    "Culture": "Kidney culture", "Capture": "Capture yield",
    "PastHuman": "Prior human serology", "HumanPooled": "Pooled human serology",
    "Random": "Uniform random expectation",
}


def count(value):
    return f"{round(value):,}" if np.isclose(value, round(value), atol=1e-9, rtol=0) else f"{value:,.1f}"


def markdown(headers, rows):
    return "\n".join([
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
        *["| " + " | ".join(map(str, row)) + " |" for row in rows],
    ])


def main():
    TABLES.mkdir(exist_ok=True)
    pooled = pd.read_csv(RESULTS / "pooled_metrics.csv")
    folds = pd.read_csv(RESULTS / "fold_metrics.csv")
    sites = pd.read_csv(RESULTS / "site_summary.csv").sort_values("human_fraction", ascending=False)
    table1_rows = []
    for row in sites.itertuples():
        table1_rows.append([
            row.site,
            f"{count(row.human_pos)}/{count(row.human_n)} ({100*row.human_fraction:.2f})",
            f"{count(row.culture_pos)}/{count(row.culture_n)} ({100*row.culture_fraction:.2f})",
            f"{count(row.caught)}/{count(row.traps)} ({100*row.capture_yield:.2f})",
        ])
    total = sites[["human_pos", "human_n", "culture_pos", "culture_n", "caught", "traps"]].sum()
    table1_rows.append([
        "Total",
        f"{count(total.human_pos)}/{count(total.human_n)} ({100*total.human_pos/total.human_n:.2f})",
        f"{count(total.culture_pos)}/{count(total.culture_n)} ({100*total.culture_pos/total.culture_n:.2f})",
        f"{count(total.caught)}/{count(total.traps)} ({100*total.caught/total.traps:.2f})",
    ])
    table1 = markdown([
        "Site", "Human positive / tested (%)", "Culture positive / effective (%)",
        "Rodents caught / trap-nights (%)",
    ], table1_rows)
    (TABLES / "table1.md").write_text(table1 + "\n")
    sites.to_csv(TABLES / "table1_network_site_summary.csv", index=False)

    rows = []
    for design, label, rules in [
        ("loyo87", "LOYO, 2015–2025", ["Culture", "Capture", "Random"]),
        ("forward_expanding", "Forward, 2018–2025", ["Culture", "Capture", "PastHuman", "Random"]),
    ]:
        for rule in rules:
            row = pooled.loc[(pooled.design == design) & (pooled.rule == rule) & (pooled.k == 3)].iloc[0]
            rows.append([
                label, RULE_NAMES[rule],
                f"{count(row.positive_selected)}/{count(row.positive_all)}",
                f"{100*row.coverage:.1f}", f"{100*row.testing_share:.1f}", f"{row.enrichment:.2f}",
            ])
    posterior = pd.read_csv(TABLES / "table_s13_posterior_decision_utility.csv")
    for strategy in ["Culture", "Capture"]:
        row = posterior.loc[posterior.strategy == strategy].iloc[0]
        rows.append(["Modelled, 2015–2025", RULE_NAMES[strategy], "—", f"{100*row['mean']:.1f}", "—", "—"])
    table2 = markdown([
        "Evaluation", "Selection rule", "Positives covered", "Coverage (%)", "Testing share (%)", "Enrichment",
    ], rows)
    (TABLES / "table2_network_evaluation.md").write_text(table2 + "\n")
    pooled[pooled.k == 3].to_csv(TABLES / "table2_observed_evaluation.csv", index=False)

    ledger = pd.read_csv(DATA / "panel_all_site_years.csv")
    ledger["site"] = ledger.site.replace(SITE_NAMES)
    columns = [
        "site", "year", "positive", "tested", "leptospira_culture_positive", "effective_cultures",
        "rodents_caught", "trap_setting_count", "eligible_human", "eligible_culture", "eligible_capture",
    ]
    ledger[columns].sort_values(["year", "site"]).to_csv(TABLES / "table_s27_complete_human_ledger.csv", index=False)
    for source, destination in {
        "pooled_metrics.csv": "table_s28_network_sizes.csv",
        "fold_metrics.csv": "table_s29_network_folds.csv",
        "selection_weights.csv": "table_s30_selection_weights.csv",
        "period_metrics.csv": "table_s31_network_periods.csv",
        "publication_overlap_counts.csv": "table_s32_publication_overlap.csv",
    }.items():
        shutil.copy2(RESULTS / source, TABLES / destination)
    forward3 = folds[(folds.design == "forward_expanding") & (folds.k == 3)]
    forward3.to_csv(TABLES / "table_s33_forward_three_sites.csv", index=False)

    wide = []
    for design, label in [
        ("loyo87", "Leave-one-year-out"), ("forward_expanding", "Forward, cumulative"),
        ("forward_rolling3", "Forward, recent three years"), ("pooled_fixed87", "Full-period fixed"),
    ]:
        for k in range(1, 8):
            subset = pooled[(pooled.design == design) & (pooled.k == k)].set_index("rule")
            human_rule = "HumanPooled" if design == "pooled_fixed87" else "PastHuman"
            values = [f"{100*subset.loc[r, 'coverage']:.1f}" if r in subset.index else "—"
                      for r in ["Culture", "Capture", human_rule, "Random"]]
            wide.append([label, k, *values])
    (TABLES / "table_s28_network_sizes.md").write_text(markdown([
        "Evaluation", "Sites selected", "Culture (%)", "Capture (%)", "Human rule (%)", "Random (%)",
    ], wide) + "\n")
    annual = []
    for year in sorted(forward3.evaluation_year.unique()):
        sub = forward3[forward3.evaluation_year == year].set_index("rule")
        def fractions(rule):
            row = sub.loc[rule]
            return f"{count(row.positive_selected)}/{count(row.positive_all)}"
        annual.append([year, fractions("Culture"), fractions("Capture"), fractions("PastHuman")])
    (TABLES / "table_s33_forward_three_sites.md").write_text(markdown([
        "Evaluation year", "Culture", "Capture", "Prior human serology",
    ], annual) + "\n")
    copied = [TABLES / f"table_s{i}_{suffix}.csv" for i, suffix in [
        (27, "complete_human_ledger"), (28, "network_sizes"), (29, "network_folds"),
        (30, "selection_weights"), (31, "network_periods"), (32, "publication_overlap"),
        (33, "forward_three_sites"),
    ]]
    manifest = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in copied}
    (TABLES / "network_table_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Wrote Table 1, Table 2, and Tables S27–S33 from the network evaluation.")


if __name__ == "__main__":
    main()
