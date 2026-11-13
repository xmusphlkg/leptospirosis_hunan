# Sentinel-network redesign: 2026-11-13

Run `.venv/bin/python script/hunan/analyze_sentinel_network.py` from the project root.

The new descriptive evaluations are recomputed from harmonised counts. The original 87 jointly observed site-years support LOYO and fixed historical rankings. Forward evaluation scores all 64 human site-years in 2018–2025, including Yuanjiang 2018, and uses only earlier valid records for each indicator.

`fold_metrics.csv` contains per-year counts and metrics; `selection_weights.csv` gives every site's exact tie-averaged inclusion weight. `pooled_metrics.csv` pools numerators and denominators, while `period_metrics.csv` partitions evaluation years without reselecting within periods. `HumanPooled` is the in-sample human-rate ranking; `PastHuman` is the forward human-history comparator. `Random` is the expectation over eligible portfolios, not one random draw. Undefined annual coverage is stored as a missing value when no human positives occur.

`pooled_fixed87` uses evaluation_year=0 to denote the full 2015–2025 record. Its selection is in-sample. `loyo87` trains on all other years and is historical representation. `forward_expanding` and `forward_rolling3` obey calendar order. All methods evaluate k=1–7.

`exact_permutation_k3.csv` and its JSON summary reproduce the primary historical culture-versus-capture comparison with 5,040 complete-trajectory permutations. No multiple-k search for statistical significance is performed.

`posterior_support.csv`, `posterior_site_rates.csv` and `posterior_paired_draws.csv` are reconstructed from the four saved primary model fits under `priority_87_2026-10-22`; no new posterior sampling is performed. The model endpoint remains the 87-row panel and should not be presented as a model of forward outcomes.

`analysis_protocol.json` records the chosen design; `analysis_summary.json` gives main findings; `validation.json` records checks and source/output hashes.
