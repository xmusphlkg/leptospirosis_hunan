# Leptospirosis sentinel-site prioritisation in Hunan

Public, reproducible materials for the 2015–2025 analysis of rodent culture and capture yield for sentinel-site prioritisation. The release contains aggregated site–year counts, analysis and figure-generation code, posterior draws, manuscript tables, figures, and machine-readable results for the current 87 site–year panel.

## Contents

- `data/analysis_panel.csv`: 87 aggregated site–year records from eight sites. The ten fields are defined in `data/README.md`.
- `script/hunan/reanalyze_priority_87.py`: descriptive, model, and posterior-summary pipeline.
- `script/hunan/plot_priority_figures.R` and `priority_figure_data.R`: rebuild the figures from saved results.
- `results/priority_87_2026-10-22/`: current analysis summaries, diagnostics, and machine-readable results.
- `manuscript/tables/` and `manuscript/figures/`: supplementary tables, figure source data, and final figures.
- `manuscript/supplementary_priority_methods.md`: methods, estimands, diagnostics, and result interpretation.

The shared data are aggregated by site and year. The repository does not contain source workbooks, person-level records, correction evidence, archived analyses, or review drafts.

## Reproduce

Install Python dependencies from `requirements.txt` in an isolated environment. The full model refit uses PyMC, nutpie, and ArviZ. Binary NetCDF posterior draws are omitted; the public panel and code can regenerate them by refitting all models.

```bash
python script/hunan/reanalyze_priority_87.py --stage all
python script/hunan/prepare_satellite_map.py
Rscript script/hunan/plot_priority_figures.R manuscript/figures
```

The map-preparation script retrieves external mapping material when run; the imagery cache is excluded from this repository. The R figure script requires `ggplot2`, `patchwork`, `dplyr`, `tidyr`, `readr`, `scales`, `jsonlite`, `ncdf4`, and `digest`.
