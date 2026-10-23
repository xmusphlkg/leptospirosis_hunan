# Supplementary Methods and Results

Accompanies *Rodent kidney culture versus capture yield for leptospirosis sentinel-site prioritisation: historical coverage and decision uncertainty in Hunan, China*.

## Key estimands

The main text separates three questions. Numbers below are the primary 87 site–year results. Definitions and sensitivity analyses follow in S1–S6.

| Estimand | What varies | What is held fixed | Primary summary |
|---|---|---|---|
| Rank concordance | Site labels | Pooled human, culture, and capture rates | Human–culture Spearman ρ = 0.714 (exact P = 0.0576). Human–capture ρ = 0.429 (P = 0.299). Culture–capture ρ = 0.286 (P = 0.501) |
| Historical coverage | Human trajectories, within the observation pattern | Animal-selected three-site sets from the other years | Culture 69.4% (E = 1.80) versus capture 47.8% (E = 1.28). Enrichment P = 0.505. The net +63 positives are concentrated in 2015–2017 and at Xiangtan |
| Decision utility | Year-standardised human rates | The two trios already ranked from the pooled animal ledgers | Beta-binomial ΔU = 10.4 percentage points (95% CrI −6.7 to 27.1; Pr(ΔU > 0) = 88.5%). Identical top-three probability 0.278, against 0.750 under the binomial likelihood |

LOYO trains on years both before and after the held-out year, so the coverage estimand is historical representation, not a forecast. Utility does not re-select animal portfolios inside each posterior draw, and neither utility nor regret includes laboratory cost. The environmental domain is not an estimand in this table.

## S1. Panel construction and measurement

The primary panel contains every site–year from 2015 to 2025 in which human serology, rodent capture, and kidney culture were all recorded: 87 site–years at eight sites. Human, culture, and capture tables were joined by site and year. Numerator and denominator totals are 291/8,834 for humans, 317/6,824 for effective cultures, and 7,333/270,469 for capture. Sampling and laboratory procedures followed the National Leptospirosis Surveillance Protocol across sites and years (Supplementary Methods S8). Indicator numerators and denominators were taken from the recorded counts. Species composition from the trapping logs describes the Xiangtan–Yuanjiang contrast and is not a model covariate.

Historical summary ledgers were audited before modelling. A site–year entered the joint panel only when it had a site-allocated human count, a recorded effective-culture count, and a trapping effort. No missing label was imputed, and no outcome was filled in as zero. The inclusion decisions are as follows.

| Period | Decision | Reason |
|---|---|---|
| 2014 | Excluded | No animal stream |
| 2015 | Included | The culture row of 0 positives among 99 effective cultures is Ningxiang |
| 2016–2017 | Included | Matched human, culture, and capture records |
| 2018 | Seven sites included | Yuanjiang conducted no animal surveillance. Its human result of 1/100 stays in the all-recorded table and is not in the joint panel. It is not coded as zero effort or a negative culture |
| 2019 | Included | Matched records |
| 2020 | Included | The annual human benchmark is the site sum, 24/805 |
| 2021–2022 | Included | Matched records |
| 2023–2025 | Included | Positives are divided by recorded effective cultures in every year |

The 2018 annual human benchmark is 39/803, which includes Yuanjiang. The joint panel for 2018 is the other seven sites, 38/703. Those denominators are not interchangeable. The 2020 benchmark, 24/805, is entirely inside the joint panel. A balanced sensitivity analysis drops 2018 entirely (80 site–years).

Culture positivity is recorded positives divided by recorded effective cultures in every year, including 2023–2025. Source reports for those three years had printed percentages of 1.92%, 11.11%, and 4.91%, matching 2/104, 13/117, and 11/224 submitted cultures. The analysis uses the recorded effective-culture totals for the same positives: 2/606, 13/628, and 11/603. On four site–years in the 87-row panel, the recorded effective count differs from total cultures minus contaminated tubes. The primary analysis keeps the recorded effective count. An arithmetic sensitivity refits the culture model with total cultures minus contaminated tubes (S5). Specimens outside the effective-culture denominator are not treated as negative cultures. No individual-level record linkage was available.

## S2. Historical selection and exact comparison

In each leave-one-year-out fold, animal selection used ratios of counts pooled across the other observed years, which include years both before and after the held-out year. The procedure asks whether an animal ranking formed inside this series represents the omitted year. It does not forecast a later year, and it does not lock a rule before any subsequent animal data are seen. Scores were sorted in descending order. The selected sets are stored in Table S4. In 2018 only seven sites were eligible, because Yuanjiang was not monitored. For each fold, selected human positives divided by all-site positives gave C, selected tests divided by all-site tests gave S, and E = C/S. Primary summaries sum counts across the 11 folds before calculating C, S, and E. Averaging annual fractions is a different estimand and gives extra weight to years with one or a few positives. Pooling counts weights coverage towards years with more positives: 2015–2017 contributed 199/291 (68.4%). In those three folds, culture selection covered 146/199 human positives and capture selection covered 96/199. In 2023–2025, inside this same full-panel evaluation, each rule covered 5/11. These are arithmetic decompositions of the existing folds.

| Rule | Human positives covered | Tests covered | Pooled C | Pooled S | E |
|---|---:|---:|---:|---:|---:|
| Capture yield | 139/291 | 3,306/8,834 | 0.4777 | 0.3742 | 1.2764 |
| Culture positivity | 202/291 | 3,415/8,834 | 0.6942 | 0.3866 | 1.7957 |
| Culture positives / trap-night | 202/291 | 3,415/8,834 | 0.6942 | 0.3866 | 1.7957 |

The compound rule selected the same sites as culture in every fold, so it does not independently replicate the culture rule. Culture selection returned the set Lianyuan, Shuangfeng, and Xiangtan in every fold. That set did not change when any one year was omitted, which is why culture's pooled LOYO coverage equals the coverage of the single full-period culture ranking (202/291). The order inside the set was not identical in every fold: when 2020 was omitted, culture ranked Shuangfeng ahead of Lianyuan, and the order of sites below the trio also moved in some folds. Capture yield produced one ordering of all eight sites in each of the ten folds in which Yuanjiang was eligible: Yuanjiang, Liuyang, Shuangfeng, then Lianyuan. Capture selection therefore returned Liuyang, Shuangfeng, and Yuanjiang except in 2018, when Yuanjiang was unavailable and Lianyuan, fourth in every other fold, entered. That one-year substitution is why the fixed capture trio covers 135/291 (46.4%) of human positives, four fewer than the pooled LOYO capture total. Those four positives are Lianyuan's human positives in 2018. Uniform random selection of three of eight sites has expected positive-count coverage 3/8 = 37.5%. That expectation is not an enrichment null and is not a cost-matched comparator.

For the paired comparison, each permutation assigns complete human positive and tested trajectories to animal-site labels while preserving the observation pattern. Yuanjiang's 10-year trajectory is not exchanged with the seven complete 11-year trajectories, so the exact reference distribution contains 7! = 5,040 assignments rather than 8! = 40,320. Animal selection remains fixed. The same assignment is used in every year and for both rules. The two-sided test uses |Δ − median(Δ)|, including ties within 10⁻¹². The null median of ΔE is 0.173, not zero. The observed ΔE = 0.519 gives two-sided P = 0.505. The central 95% null interval is −0.680 to 0.950. The observed ΔC = 0.216 gives two-sided P = 0.471; its null median is 0.065 and its central 95% null interval is −0.258 to 0.375. These intervals are null reference ranges, not confidence intervals for the observed differences.

Exactness is conditional on exchangeability within the observation pattern. Preserving human trajectories retains within-site temporal pairing. It does not establish exchangeability of sites with different populations or exposure histories, and it does not let Yuanjiang's missing year be filled by another site. Animal-only selection prevents direct selection on the evaluation outcome. The same-site historical design still limits generalisation.

## S3. Model specification and posterior summaries

For each stream separately, an observed count y has a beta-binomial distribution with denominator n and shape parameters μφ and (1 − μ)φ. Its variance is nμ(1 − μ)(n + φ)/(1 + φ). Logit μ is an intercept plus site and year contrasts. The orthonormal contrast basis Q was obtained by QR decomposition of a matrix formed by stacking Iₙ₋₁ above a row of −1 values. Coefficients are independent normal with scale τ, and effect vectors Qθ sum to zero. Consequently, τ is the coefficient prior scale, and the marginal prior standard deviation of an individual site effect is τ√(7/8), not precisely τ. Site coefficients use scale 0.70 and year coefficients use scale 0.50 in the primary fits. Human and culture intercept priors are N(−3.4, 2²) and N(−3.0, 2²). The beta-binomial concentration prior is Gamma(shape 2, rate 0.1), with prior mean 20.

Year-standardised rates are the equal-weight mean of the inverse-logit predictions across all 11 calendar years, 2015–2025. For every site and every year, observed or not, the summand is the inverse logit of the intercept, the site effect, and the year effect. It is a conditional mean, not an observed proportion and not a new beta-binomial predictive draw. Yuanjiang has no 2018 count. Its 2018 term uses that same linear predictor: the intercept, Yuanjiang's site effect, and the 2018 year effect estimated from the seven monitored sites. Because the inverse logit is strictly increasing and the intercept plus the year effect is common to all sites within a year, the standardised rate has the same site order as the site effect in every draw. The average supplies probability-scale rates for utility, with each year weighted equally. For each draw, sites were ranked from the highest to the lowest standardised rate. Under the binomial model with the same predictor, the year-standardised human medians were 7.23% at Shuangfeng, 6.67% at Xiangtan, 3.61% at Lianyuan, and 3.11% at Liuyang, with marginal top-three probabilities 1.000, 1.000, 0.750, and 0.247. Table S11 reports the beta-binomial rates. Marginal top-three probability is the frequency of ranks 1–3. Independently fitted stream draws were paired in chain-then-draw order, with 6,000 retained draws for each primary model. The induced rank correlation and joint rank probabilities therefore use the product of the two posterior distributions, without an estimated cross-stream covariance. Capture yield remains an empirical comparator and has no sampling model of its own.

| Model | Posterior rank correlation, median (95% credible interval) | Probability of an identical top-three set | Pr(ρ > 0) |
|---|---|---:|---:|
| Beta-binomial, site coefficient scale 0.70 | 0.667 (0.143, 0.929) | 0.278 | 0.992 |
| Binomial, site coefficient scale 0.70 | 0.810 (0.667, 0.952) | 0.750 | 1.000 |

No Laplace approximation was used for these rankings. Prior and window checks of the human decision contrast are in S5 and S6.2. They do not replace this likelihood comparison.

| Site | Beta-binomial joint top-three probability |
|---|---:|
| Shuangfeng | 0.870 |
| Xiangtan | 0.691 |
| Lianyuan | 0.596 |
| Liuyang | 0.054 |
| Yuanjiang | 0.001 |
| Ningxiang | 0.0002 |
| Liling | 0 |
| Junshan | 0 |

Joint site probabilities need not sum to three. The probability that the complete top-three sets coincide, 0.278 under the beta-binomial model, is a different quantity. Two-signal non-dominance probabilities were not recomputed for this panel.

## S4. Numerical diagnostics and posterior predictive checks

Sampling used the No-U-Turn Sampler in PyMC (v5.28.5) with the nutpie backend, four chains, and target acceptance 0.95. Each primary model retained 1,500 draws per chain after 1,500 tuning draws (6,000 draws). Each window and prior-sensitivity model retained 1,000 draws per chain after 1,000 tuning draws (4,000 draws). Diagnostics were read in ArviZ (v0.23.4).

| Model / stream | Divergences | Maximum R-hat | Minimum bulk ESS | Minimum tail ESS | Draws |
|---|---:|---:|---:|---:|---:|
| Binomial / human | 0 | 1.002 | 4,679 | 4,000 | 6,000 |
| Binomial / culture | 0 | 1.003 | 2,417 | 3,543 | 6,000 |
| Beta-binomial / human | 0 | 1.002 | 4,516 | 4,091 | 6,000 |
| Beta-binomial / culture | 0 | 1.002 | 4,271 | 4,097 | 6,000 |

Across all 14 stored fits, there were no divergent transitions. The maximum R-hat was 1.004, the minimum bulk effective sample size was 2,415, and the minimum tail effective sample size was 2,280. Those two minima are the 2015–2022 human fit and the balanced human fit, respectively. R-hat is a diagnostic summary, not proof of convergence. Beta-binomial concentration medians were 15.4 (95% credible interval 8.4–26.8) for human counts and 6.7 (3.5–11.9) for culture counts.

| Model / stream | Observed zero cells | Replicated zero cells, central 95% | Observed total | Replicated total, central 95% | Cells inside pointwise 95% intervals | Pearson tail fraction |
|---|---:|---|---:|---|---:|---:|
| Binomial / human | 43 | 19–35 | 291 | 247–338 | 81/87 | 0.001 |
| Binomial / culture | 56 | 31–46 | 317 | 272–362 | 75/87 | 0 |
| Beta-binomial / human | 43 | 30–51 | 291 | 191–444 | 85/87 | 0.887 |
| Beta-binomial / culture | 56 | 44–63 | 317 | 157–489 | 87/87 | 0.825 |

The Pearson discrepancy uses the posterior mean μ̄ and the binomial-scale denominator nμ̄(1 − μ̄) for observed and replicated counts: Σ(y − nμ̄)² / [nμ̄(1 − μ̄) + 10⁻¹²]. For the beta-binomial model this is a descriptive discrepancy, not a variance-standardised goodness-of-fit test. Tail fractions compare replicated and observed discrepancies. They are not classical calibrated P values. Replicated counts were drawn from the saved posterior with NumPy `default_rng(20261522)`. Beta-binomial simulation first drew p from Beta(μφ, (1 − μ)φ), then y from Binomial(n, p). The wider beta-binomial intervals include the observed totals and nearly all site–year counts. This is in-sample compatibility. Pointwise coverage is not an external calibration estimate, and broad intervals can conceal model deficiencies. The primary model is the beta-binomial specification because the binomial checks do not reproduce the zero counts. The binomial fit is retained as the likelihood comparison.

![Figure S1](figures/figure_s1_predictive_checks.png)

**Figure S1. Posterior predictive distributions of totals and zero counts.** Rows are human serology and rodent culture. Columns are the replicated total number of positives and the number of zero-positive site–years. Curves are binomial (blue) and beta-binomial (orange) densities from the saved draws. Dashed lines mark the observed values. The checks use the 87 observed site–years. Pointwise predictive intervals for individual cells are in Table S10 and are not simultaneous.

![Figure S2](figures/figure_s2_period_prior_sensitivity.png)

**Figure S2. Period structure and sensitivity of the culture-minus-capture utility difference.** (A) Pooled human seroreactivity against recorded culture yield for 2015–2022 (63 site–years) and for 2023–2025 alone (24 site–years). (B) Posterior mean and 95% credible interval of ΔU under the primary beta-binomial model, the binomial likelihood, five one-at-a-time prior changes, the balanced 80-row panel, and the 2015–2022 window. The dashed line is zero.

## S5. Window and denominator checks

The primary analysis does not allocate extra human positives to sites. The 2018 and 2020 annual benchmarks were set to the sums of the recorded site counts, 39/803 and 24/805. The earlier enumeration of unallocated annual positives is not part of this panel. What changes the joint analysis is which observed site–years are included, and which culture denominator is used.

| Analysis | Rows | Human–culture Spearman ρ (exact P) | LOYO culture vs capture coverage | Enrichment P | Posterior mean ΔU, percentage points (95% CrI; probability) | Posterior-mean optimum |
|---|---:|---|---|---:|---|---|
| Primary, 2015–2025 | 87 | 0.714 (0.0576) | 202/291 vs 139/291 | 0.505 | 10.4 (−6.7 to 27.1; 88.5%) | Lianyuan, Shuangfeng, Xiangtan |
| Drop 2018 | 80 | 0.762 (0.0368) | 177/253 vs 121/253 | 0.338 | 7.9 (−10.3 to 26.0; 80.2%) | Lianyuan, Shuangfeng, Xiangtan |
| Years through 2022 | 63 | 0.903 (0.0054) | 197/280 vs 134/280 | 0.488 | 12.0 (−5.2 to 28.7; 91.0%) | Liuyang, Shuangfeng, Xiangtan |
| Exclude 2015 | 79 | 0.707 (0.0595) | 126/195 vs 80/195 | 0.246 | Descriptive only | Not refitted |
| 2023–2025 only | 24 | 0.347 (0.387) | 4/11 vs 0/11 | 0.187 | Descriptive only | Not refitted |
| Arithmetic effective cultures | 87 | 0.714 (0.0576) | Same LOYO coverage as the primary panel | 0.505 | Culture model refitted; see Table S17 | Not a new human optimum |

The 2023–2025 row is a separate evaluation trained only inside those three years. It is not the contribution of those years to the full-panel LOYO, which was 5/11 under both rules (net zero; Table S6). Dropping 2018 makes every remaining site observed in every remaining year, so that permutation uses all 40,320 trajectory assignments and both null medians are zero. The observed enrichment difference is 0.534 (two-sided P = 0.338; central 95% null interval −0.968 to 0.968) and the observed coverage difference is 22.1 percentage points (P = 0.314; central 95% null interval −37.9 to 37.9 percentage points). The 2015–2022 and primary permutations still hold Yuanjiang's trajectory fixed. Posterior rank correlations for the balanced, early, and arithmetic-denominator fits stayed positive, with Pr(ρ > 0) of 0.981, 0.996, and 0.993 (Table S17).

## S6. Descriptive portfolios, decision utility, and the design simulation

All combinations of k = 1 to 7 sites were enumerated from pooled observed counts. Table S7 records the 56 three-site sets. Table S9 records all 254 non-empty sets smaller than the full network. The mean observed positive coverage across sets at each k equals k/8. The best observed fixed set uses the same human counts to choose and to score a portfolio, so its coverage is an in-sample reference, not an estimate of future performance. Culture- and capture-ranked portfolios on the full period are likewise in-sample. The animal-only LOYO comparison in S2 is the historical cross-evaluation. No budget, site-cost, geographic-coverage, or feasibility constraint was applied.

Tables S5 and S6 decompose the LOYO human-positive difference. Each contribution is (I^culture − I^capture) × Y. Site totals and year totals each sum to +63, the difference between 202 and 139 positives. The spatial pieces are Xiangtan +79, Lianyuan +37, Liuyang −35, and Yuanjiang −18. This is an arithmetic decomposition, not another hypothesis test.

Pooled animal rates were untied at every k (Table S8). Culture ranking covered more human positives at k = 1 to 4 and matched the best observed trio at k = 3 (69.4%). Capture ranking covered more from k = 5 to 7. At k = 5 the culture set adds Junshan (83.5%) and the capture set adds Xiangtan (87.6%). At k = 6 the coverages are 89.7% and 91.1%. At k = 7 they are 93.1% and 97.9%. Figure 3B shows these single traces. There is no tie range to average.

### S6.1 Posterior decision utility and regret

For each human posterior draw, year-standardised site rates are weighted by observed testing totals: 1,006 tests at Yuanjiang and 1,100–1,210 at the other sites. The culture and capture sets stay at the portfolios ranked from the pooled animal ledgers. They are not re-selected inside the posterior loop. Staff time, consumables, transport, and molecular-testing costs are outside both utility and regret. For a fixed three-site set A, U(A) = Σ_{s in A} n_s p_s / Σ_s n_s p_s. Table S13 reports posterior utility, the probability of being the draw-wise best set, and expected regret relative to the maximum across all 56 sets. Under the primary beta-binomial model, the culture-selected set had mean utility 0.552 and expected regret 0.034. The capture-selected set had mean utility 0.448 and expected regret 0.139. Uniform random selection among the 56 sets has conditional expected utility 0.375 in every draw, expected regret 0.211, and probability 1/56 of being the draw-wise best set. That degenerate utility interval describes the conditional average, not the variability of one randomly realised portfolio.

Table S14 compares the two fixed portfolios within each human draw: mean ΔU = 0.104, 95% equal-tailed credible interval −0.067 to 0.271, and Pr(ΔU > 0 | D) = 0.885. Being the single best of 56 sets is a different event. Its probability was 0.334 for the culture trio and 0.007 for the capture trio. The draw-wise maximum changes and cannot be implemented as one policy. Expected regret satisfies R_capture − R_culture = E(U_culture − U_capture). The posterior-mean utility maximiser is the culture-selected set and was identified using all human data in the primary fit. No minimax, geographic-coverage, or core-plus-rotation optimisation was performed.

The three culture-coverage quantities have different targets. For culture they happen to share the same three sites. For capture they do not.

| Quantity | How the portfolio is chosen | What is evaluated | Culture | Capture |
|---|---|---|---:|---:|
| LOYO observed coverage | Animal data from the other observed years | Held-out human positives, pooled across 11 folds | 69.4% | 47.8% |
| Full-period fixed coverage | Pooled animal ranking | The same observed human positives | 69.4% | 46.4% |
| Posterior mean utility | Fixed pooled animal-ranked trio | Year-standardised latent human rates, weighted by observed tests | 55.2% | 44.8% |

Utility uses historical testing totals. It does not add future-year variability, equal sampling across sites, laboratory cost, or geographic representativeness.

### S6.2 Targeted decision-prior sensitivity

Five additional human beta-binomial fits changed one prior at a time. The baseline intercept remained N(−3.4, 2²). The alternatives were site coefficient scales 0.35 and 1.40, a year coefficient scale of 1.00, and concentration priors Gamma(2, rate 0.2) and Gamma(2, rate 0.05), with prior means 10 and 40 around the baseline mean of 20. Each used four chains, 1,000 tuning draws, 1,000 retained draws, and target acceptance 0.95. Fixed pooled animal-ranked sets do not depend on the human prior, so only the human model was refitted. Table S18 reports the paired differences. These one-at-a-time comparisons do not cover combined prior changes or future-year prediction.

| Specification | Mean ΔU, percentage points | 95% CrI | Pr(culture > capture) |
|---|---:|---|---:|
| Primary: site SD 0.70, year SD 0.50, concentration mean 20 | 10.4 | −6.7 to 27.1 | 0.885 |
| Site SD 0.35 | 6.7 | −6.8 to 19.9 | 0.844 |
| Site SD 1.40 | 12.3 | −6.7 to 30.0 | 0.898 |
| Year SD 1.00 | 12.7 | −3.3 to 27.8 | 0.943 |
| Concentration mean 10 | 9.9 | −7.5 to 26.8 | 0.869 |
| Concentration mean 40 | 10.7 | −6.8 to 27.6 | 0.893 |

Every interval includes zero, and every posterior-mean optimum in this table is Lianyuan, Shuangfeng, and Xiangtan. The calendar window is less stable than these priors: the 2015–2022 fit moves the mean-utility optimum to Liuyang, Shuangfeng, and Xiangtan (S5).

### S6.3 Decision-utility sensitivity to likelihood family

Table S15 reports the paired utility difference for every stored human fit. The likelihood comparison below uses the primary binomial and beta-binomial rows. Both fits have 6,000 draws and the same mean-structure priors. Paired differences are computed within each human posterior separately.

| Quantity | Beta-binomial | Binomial |
|---|---:|---:|
| Culture mean utility (95% CrI), % | 55.2 (43.7–66.2) | 67.9 (62.6–72.9) |
| Capture mean utility (95% CrI), % | 44.8 (33.8–55.9) | 45.0 (39.8–50.2) |
| Paired mean difference (95% CrI), percentage points | 10.4 (−6.7 to 27.1) | 22.9 (14.8–30.8) |
| Probability culture utility exceeds capture | 0.885 | 1.000 |
| Culture expected regret | 0.034 | 0.004 |
| Capture expected regret | 0.139 | 0.232 |
| Probability of being best among 56 sets | 0.334 | 0.757 |

Both posterior-mean maxima are Lianyuan, Shuangfeng, and Xiangtan. The likelihood changes the size and uncertainty of the advantage. These calculations condition on the observed animal-selected sets and do not propagate repeated animal sampling or reselection. The LOYO permutation and this comparison answer different questions.

### S6.4 Conditional design simulation

A supporting simulation fixed the observed 2015–2025 site denominators, used logistic-normal site effects with standard deviation 0.70, and set the latent human–culture site-effect correlation to 0, 0.4, or 0.8. Each scenario had 1,500 binomial replicates. Each replicate was tested with 999 random site-label permutations. The generator was NumPy `default_rng(20261022)`. Rejection proportions at nominal two-sided 0.05 were 5.7%, 11.6%, and 48.5% (Table S1). These are properties of the specified binomial simulation. They do not estimate power under the beta-binomial model, calibrate the posterior ranks, or verify the data-generating mechanism.

## S7. Reproducibility

Scripts, the released aggregated site–year counts, posterior summaries, and the outputs that reproduce this analysis are available in the public [leptospirosis_hunan repository](https://github.com/xmusphlkg/leptospirosis_hunan). Binary NetCDF posterior draws are not included; the analysis code can regenerate them from the released panel.

The current analysis is a retrospective evaluation of the established sentinel network. It fixes the 87-row panel, the pooled Spearman comparison, animal-only three-site historical evaluation, the pattern-preserving trajectory permutation, and beta-binomial uncertainty as the primary presentation. Binomial summaries are the likelihood comparison.

| Calculation | Seed / sampling specification |
|---|---|
| Primary human beta-binomial | seed 20293088; 4 × 1,500 retained draws; nutpie |
| Primary culture beta-binomial | seed 20358260; 4 × 1,500 retained draws; nutpie |
| Primary human binomial | seed 20264042; 4 × 1,500 retained draws; nutpie |
| Primary culture binomial | seed 20288853; 4 × 1,500 retained draws; nutpie |
| Human prior refits, in the order site SD 0.35, site SD 1.40, year SD 1.00, concentration means 10 and 40 | seeds 20297592, 20285636, 20293241, 20261647, 20280162; 4 × 1,000 retained draws |
| Balanced 80-row human and culture fits | seeds 20316860 and 20286958 |
| 2015–2022 human and culture fits | seeds 20321969 and 20356439 |
| Arithmetic effective-culture fit | seed 20263082; culture model only |
| Posterior predictive replicates in S4 | `default_rng(20261522)`, saved draws only |
| Design simulation in S6.4 | `default_rng(20261022)`; 1,500 replicates; 999 permutations |
| Exact permutations | Exhaustive enumeration; no random-number generator |

Recreate the descriptive analyses, refit all 14 models, and rebuild posterior summaries from the released panel:

```bash
python script/hunan/reanalyze_priority_87.py --stage all
```

Redraw Figures 1–4, S1, and S2 from the saved outputs:

```bash
python script/hunan/prepare_satellite_map.py
Rscript script/hunan/plot_priority_figures.R manuscript/figures
```

Figures are 320 dpi PNG files and vector PDFs drawn with ggplot2, patchwork, and `theme_bw()`. Figures 1 and 4 have six and five panels. Figure 2 has five panels. Figure 3 has three panels. Figure S1 has four predictive panels, and Figure S2 has two. Panel letters are the only in-panel annotations. Site order follows pooled human seroreactivity. Wilson intervals in Figure 1 are descriptive binomial intervals. Model uncertainty is in Figure 4. With eight sites, the posterior Spearman correlation has discrete support, so Figure 4E shows that probability mass without density smoothing.

The earlier scripts `script/manuscript/build_priority_manuscript_assets.py` and `script/hunan/decision_extension.py` rebuild the five-year archive. They are not the entry points for this 87-row manuscript.

## Supplementary table files

- [Table S1: conditional design simulation at the 2015–2025 site denominators](tables/table_s1_design_simulation.csv).
- [Table S2: 87 site–year count records](tables/table_s2_site_year_counts.csv).
- [Table S3: leave-one-site-out rank influence](tables/table_s3_leave_one_site_out.csv).
- [Table S4: annual historical selections, coverage, and enrichment](tables/table_s4_historical_evaluation.csv).
- [Table S5: site contributions to the held-out coverage difference](tables/table_s5_site_contributions.csv).
- [Table S6: year contributions to the held-out coverage difference](tables/table_s6_year_contributions.csv).
- [Table S7: all 56 observed three-site portfolios](tables/table_s7_all_portfolios.csv).
- [Table S8: fixed-portfolio coverage frontier by network size](tables/table_s8_portfolio_frontier.csv).
- [Table S9: all fixed portfolios for k = 1–7](tables/table_s9_all_portfolios_k1_k7.csv).
- [Table S10: beta-binomial pointwise predictive count intervals](tables/table_s10_predictive_cells.csv).
- [Table S11: beta-binomial marginal site rates and top-three probabilities](tables/table_s11_beta_binomial_site_rates.csv).
- [Table S12: beta-binomial joint top-three probabilities](tables/table_s12_joint_priorities.csv).
- [Table S13: posterior decision utility and regret](tables/table_s13_posterior_decision_utility.csv).
- [Table S14: paired posterior utility comparison](tables/table_s14_paired_posterior_utility.csv).
- [Table S15: decision-utility sensitivity across stored human fits](tables/table_s15_decision_model_sensitivity.csv).
- [Table S16: sampling diagnostics for the 14 stored fits](tables/table_s16_nuts_laplace.csv). The filename retains an earlier suffix; the file records NUTS diagnostics and does not compare a Laplace approximation.
- [Table S17: posterior rank correlation under the balanced, early, and arithmetic-denominator fits](tables/table_s17_binomial_prior_sensitivity.csv). The filename retains an earlier suffix.
- [Table S18: human decision utility under the primary fit and the prior refits](tables/table_s18_decision_prior_sensitivity.csv).

## S8. Published surveillance-method framework

Provincial sentinel serosurveys were conducted each year in the pre-epidemic season (April–May). At each site, stratified cluster random sampling collected 20 sera from healthy residents in each of five age groups (10–19, 20–29, 30–39, 40–49, and ≥50 years). Persons with a history of leptospirosis or leptospirosis vaccination were excluded. MAT was performed against 15 major serovars, and a reciprocal titre ≥1:50 was classified as seropositive; the seropositivity rate was the number of seropositive persons divided by the number tested. These procedures implement the national trial surveillance protocol, which specifies pre-season collection of 100 healthy-person sera per sentinel site across five age bands and MAT using 15 representative strains, with titres read at the 50% agglutination endpoint. Its animal-culture appendix describes aseptic kidney inoculation into Korthof medium, incubation at 28°C, microscopy every 5–7 days and observation for 60 days before a negative determination; serogrouping is a subsequent laboratory step. The historical site–year ledgers used in this evaluation recorded only aggregate positive and tested counts, without serovar-level disaggregation. Trapping effort in the rodent density surveys is recorded as trap-nights (standardised night-trap deployments). Source: [national surveillance protocol](https://www.chinacdc.cn/jkyj/crb2/yl/gdlxtb/jswj_gdlxtb/202409/W020240909361196300710.pdf), PDF pages 6, 16 and 20 (printed pages 150, 160 and 164).
