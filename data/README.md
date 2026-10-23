# Released analysis data

`analysis_panel.csv` contains the 87 observed site–year records used in the current Hunan sentinel-priority analysis (eight sites, 2015–2025). Rows are aggregate counts, not individual-level observations. The source surveillance workbooks and row-level audit trail are not included.

| Field | Definition |
|---|---|
| `site` | English site label used in the analysis |
| `year` | Surveillance year |
| `positive` | Recorded human seropositive count |
| `tested` | Human tests conducted |
| `leptospira_culture_positive` | Recorded positive rodent cultures |
| `effective_cultures` | Effective culture denominator used in the primary analysis |
| `rodents_caught` | Recorded rodents caught |
| `trap_setting_count` | Recorded trap-setting effort |
| `cultures_total` | Total cultures in the site–year aggregate, retained for the arithmetic-denominator sensitivity analysis |
| `contaminated` | Contaminated cultures in the site–year aggregate, retained for the arithmetic-denominator sensitivity analysis |

The primary manuscript data table is `../manuscript/tables/table_s2_site_year_counts.csv`. The two additional culture fields support the denominator sensitivity analysis in the released code.
