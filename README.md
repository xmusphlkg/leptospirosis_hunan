# Leptospirosis sentinel-site selection in Hunan

Public, reproducible materials for the current manuscript, **Selecting leptospirosis sentinel sites using rodent surveillance: a retrospective evaluation in Hunan, China, 2015–2025**. The release evaluates how rodent kidney-culture positivity and capture yield select sentinel sites and how their representation of recorded human seropositive results changes with network size.

## Public contents

- `data/analysis_panel.csv`: 87 aggregated site–year records used for the pooled indicator, LOYO and model analyses.
- `data/complete_human_ledger.csv`: the complete 88-row human ledger with animal-record availability for the forward evaluation.
- `manuscript/main.md`: current manuscript source.
- `manuscript/supplementary_priority_methods.md`: current supplementary methods and Tables S1–S33.
- `manuscript/tables/`: current table data, including network-size, forward-evaluation and publication-reconciliation tables.
- `manuscript/figures/`: current main and supplementary figures plus figure source data.
- `results/network_redesign_2026-11-13/`: current network-evaluation outputs, protocols and independent validation.
- `results/priority_87_2026-10-22/`: model-support outputs retained for the current manuscript; binary NetCDF posterior draws are omitted.
- `script/hunan/`: analysis, verification and figure-generation scripts for the current network evaluation.
- `script/manuscript/`: table, figure-source and validation utilities.

All shared records are aggregated by site and year. The repository does not contain source surveillance workbooks, individual-level records, correction workbooks or binary posterior draws. The source ledgers are held by the Hunan Center for Disease Control and Prevention.

## Reproduce the current network evaluation

Install dependencies from `requirements.txt` in an isolated environment, then run:

```bash
python script/hunan/analyze_sentinel_network.py
python script/hunan/verify_sentinel_network.py
python script/manuscript/build_network_tables.py
```

The saved outputs in `results/network_redesign_2026-11-13/` reproduce the reported selection weights, coverage curves, forward evaluation and independent validation. Figure scripts require the R packages documented in their headers.
