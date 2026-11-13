"""Rebuild sentinel-site coverage analyses under the 2026-11-13 design.

The 87 jointly observed site-years support the historical LOYO comparison.
Forward evaluation scores all human records in 2018--2025 and never uses the
evaluation year's animal availability to determine candidate eligibility.
Existing posterior draws are re-summarised; no model is refitted.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import rankdata


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/processed/harmonized_2026-10-22"
LEGACY = ROOT / "output/hunan/current/priority_87_2026-10-22"
DEFAULT_OUT = ROOT / "output/hunan/current/network_redesign_2026-11-13"
LABELS = {
    "双峰县": "Shuangfeng", "君山区": "Junshan", "宁乡县": "Ningxiang",
    "沅江市": "Yuanjiang", "浏阳市": "Liuyang", "涟源市": "Lianyuan",
    "湘潭县": "Xiangtan", "醴陵市": "Liling",
}
SITES = sorted(LABELS.values())
RULES = {
    "Culture": ("leptospira_culture_positive", "effective_cultures", "eligible_culture"),
    "Capture": ("rodents_caught", "trap_setting_count", "eligible_capture"),
    "PastHuman": ("positive", "tested", "eligible_human"),
    "HumanPooled": ("positive", "tested", "eligible_human"),
}
COUNTS = ["positive", "tested", "leptospira_culture_positive", "effective_cultures",
          "rodents_caught", "trap_setting_count"]
METRIC_COUNTS = ["positive_selected", "positive_all", "tested_selected", "tested_all"]
CONFIG = {
    "protocol_date": "2026-11-13", "k_values": list(range(1, 8)),
    "historical_design": "loyo87", "forward_initialisation": [2015, 2016, 2017],
    "forward_evaluation_years": list(range(2018, 2026)),
    "forward_updates": ["expanding", "rolling3"],
    "ties": "equal weight over every boundary-tied admissible portfolio",
    "random": "equal weight over every k-site portfolio among eligible sites",
    "pooling": "sum selected and total counts across years before taking ratios",
    "zero_positive_year": "coverage and enrichment undefined; retain numerator and denominator",
    "forward_eligibility": "valid training records only, independent of current animal availability",
    "primary_permutation": "loyo87 k=3, human trajectories exchanged only within observation-pattern groups",
    "models": "re-summarise saved 87-site-year posterior draws without new sampling",
}


def dump(value, path):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False,
                                    default=lambda x: x.item() if isinstance(x, np.generic) else str(x)) + "\n")


def save(frame, out, name):
    frame.to_csv(out / name, index=False)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_panel(name):
    frame = pd.read_csv(DATA / name)
    frame["site"] = frame.site.map(LABELS)
    assert frame.site.notna().all() and not frame.duplicated(["site", "year"]).any()
    return frame.sort_values(["year", "site"]).reset_index(drop=True)


def ratios(ps, pa, ns, na):
    coverage = float(ps / pa) if pa else np.nan
    share = float(ns / na) if na else np.nan
    enrichment = float(coverage / share) if share > 0 and pa else np.nan
    return {"coverage": coverage, "testing_share": share, "enrichment": enrichment}


def ranked_portfolios(training, rule, candidates, k):
    """Use exact rational rates so equivalent fractions remain genuine ties."""
    if rule == "Random":
        choices = list(itertools.combinations(sorted(candidates), k))
        weights = {site: k / len(candidates) for site in candidates}
        return choices, weights, None
    num, den, valid = RULES[rule]
    grouped = training.loc[training[valid]].groupby("site")[[num, den]].sum()
    candidates = sorted(site for site in candidates if site in grouped.index and grouped.at[site, den] > 0)
    assert len(candidates) >= k
    scores = {site: Fraction(int(grouped.at[site, num]), int(grouped.at[site, den])) for site in candidates}
    threshold = sorted(scores.values(), reverse=True)[k - 1]
    required = [site for site in candidates if scores[site] > threshold]
    boundary = [site for site in candidates if scores[site] == threshold]
    choices = [tuple(sorted(required + list(extra)))
               for extra in itertools.combinations(boundary, k - len(required))]
    weights = {site: sum(site in portfolio for portfolio in choices) / len(choices) for site in candidates}
    return choices, weights, grouped


def one_fold(design, year, training, outcome, candidates, rules, outcome_panel):
    fold_rows, weight_rows = [], []
    outcome = outcome.set_index("site")
    years = sorted(int(y) for y in training.year.unique())
    for rule in rules:
        for k in range(1, 8):
            choices, weights, grouped = ranked_portfolios(training, rule, candidates, k)
            assert np.isclose(sum(weights.values()), k)
            assert set(weights).issubset(outcome.index)
            ps = sum(weights[site] * outcome.at[site, "positive"] for site in weights)
            ns = sum(weights[site] * outcome.at[site, "tested"] for site in weights)
            pa, na = int(outcome.positive.sum()), int(outcome.tested.sum())
            fold_rows.append({
                "design": design, "rule": rule, "k": k, "evaluation_year": year,
                "training_start": min(years), "training_end": max(years),
                "training_years": ";".join(map(str, years)), "outcome_panel": outcome_panel,
                "positive_selected": float(ps), "positive_all": pa,
                "tested_selected": float(ns), "tested_all": na,
                **ratios(ps, pa, ns, na), "n_tied_sets": len(choices),
                "n_candidates": len(weights),
                "selected_sites": ";".join(sorted(site for site, weight in weights.items() if weight > 0)),
                "admissible_sets": "|".join(";".join(portfolio) for portfolio in choices),
            })
            for site in SITES:
                num_value = den_value = rate = np.nan
                if grouped is not None and site in grouped.index:
                    num, den, _ = RULES[rule]
                    num_value, den_value = int(grouped.at[site, num]), int(grouped.at[site, den])
                    rate = num_value / den_value if den_value > 0 else np.nan
                weight_rows.append({
                    "design": design, "rule": rule, "k": k, "evaluation_year": year,
                    "site": site, "eligible": site in weights, "weight": weights.get(site, 0.0),
                    "training_numerator": num_value, "training_denominator": den_value,
                    "training_rate": rate,
                })
    return fold_rows, weight_rows


def evaluate(joint, full):
    fold_rows, weight_rows = [], []
    for year in range(2015, 2026):
        training = joint[joint.year != year]
        outcome = joint[joint.year == year]
        rows, weights = one_fold("loyo87", year, training, outcome, sorted(outcome.site),
                                 ["Culture", "Capture", "Random"], "joint87")
        fold_rows.extend(rows)
        weight_rows.extend(weights)
    for design in ["forward_expanding", "forward_rolling3"]:
        for year in range(2018, 2026):
            start = 2015 if design == "forward_expanding" else year - 3
            training = full[(full.year >= start) & (full.year < year)]
            # Candidate status is determined before looking at outcome data.
            candidates = sorted(training.site.unique())
            assert candidates == SITES and training.year.max() < year
            for _, _, valid in RULES.values():
                assert set(training.loc[training[valid], "site"]) == set(candidates)
            outcome = full.loc[(full.year == year) & full.eligible_human]
            assert set(outcome.site) == set(candidates)
            rows, weights = one_fold(design, year, training, outcome, candidates,
                                     ["Culture", "Capture", "PastHuman", "Random"], "human88")
            fold_rows.extend(rows)
            weight_rows.extend(weights)
    # Fixed full-period rankings are an in-sample descriptive reference.
    outcome = joint.groupby("site")[["positive", "tested"]].sum().reset_index()
    rows, weights = one_fold("pooled_fixed87", 0, joint, outcome, SITES,
                             ["Culture", "Capture", "HumanPooled", "Random"], "joint87")
    fold_rows.extend(rows)
    weight_rows.extend(weights)
    return pd.DataFrame(fold_rows), pd.DataFrame(weight_rows)


def pool_folds(folds):
    rows = []
    for (design, rule, k), part in folds.groupby(["design", "rule", "k"], sort=False):
        totals = part[METRIC_COUNTS].sum().to_dict()
        years = (2015, 2025) if design == "pooled_fixed87" else (int(part.evaluation_year.min()), int(part.evaluation_year.max()))
        rows.append({"design": design, "rule": rule, "k": int(k),
                     "year_start": years[0], "year_end": years[1], "n_folds": len(part),
                     "outcome_panel": part.outcome_panel.iloc[0], **totals,
                     **ratios(*(totals[col] for col in METRIC_COUNTS))})
    return pd.DataFrame(rows)


def site_summary(joint):
    result = joint.groupby("site").agg(
        human_pos=("positive", "sum"), human_n=("tested", "sum"),
        culture_pos=("leptospira_culture_positive", "sum"), culture_n=("effective_cultures", "sum"),
        caught=("rodents_caught", "sum"), traps=("trap_setting_count", "sum"),
        observed_years=("year", "nunique")).reset_index()
    for stream, num, den in [("human", "human_pos", "human_n"), ("culture", "culture_pos", "culture_n"),
                             ("capture", "caught", "traps")]:
        p, n = result[num] / result[den], result[den]
        result[f"{stream}_fraction"] = p
        result[f"{stream}_rank"] = rankdata(-p, method="average")
        z = 1.959963984540054
        middle = (p + z * z / (2 * n)) / (1 + z * z / n)
        half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
        result[f"{stream}_low"], result[f"{stream}_high"] = middle - half, middle + half
    result["capture_yield"] = result.capture_fraction
    result["trap_effort"] = result.traps
    return result


def exact_rank_correlation(x, y):
    xr, yr = rankdata(x), rankdata(y)
    xc, yc = xr - xr.mean(), yr - yr.mean()
    scale = np.linalg.norm(xc) * np.linalg.norm(yc)
    observed = float(xc @ yc / scale)
    values = np.asarray(list(itertools.permutations(yc))) @ xc / scale
    return {"rho": observed, "p_two_sided": float(np.mean(abs(values) >= abs(observed) - 1e-12)),
            "n_permutations": len(values)}


def trajectory_permutation(joint, weights, pooled):
    years = list(range(2015, 2026))
    y = joint.pivot(index="year", columns="site", values="positive").reindex(index=years, columns=SITES)
    n = joint.pivot(index="year", columns="site", values="tested").reindex(index=years, columns=SITES)
    patterns = {}
    for column, site in enumerate(SITES):
        patterns.setdefault(tuple(y[site].notna()), []).append(column)
    groups = list(patterns.values())
    mappings = []
    for permutations in itertools.product(*(itertools.permutations(group) for group in groups)):
        mapping = np.arange(len(SITES))
        for group, permutation in zip(groups, permutations):
            mapping[group] = permutation
        mappings.append(mapping)
    mappings = np.asarray(mappings)
    totals = {}
    yv, nv = y.fillna(0).to_numpy(), n.fillna(0).to_numpy()
    for rule in ["Culture", "Capture"]:
        w = weights[(weights.design == "loyo87") & (weights.rule == rule) & (weights.k == 3)]
        w = w.pivot(index="evaluation_year", columns="site", values="weight").reindex(index=years, columns=SITES).to_numpy()
        selected_y = np.einsum("yps,ys->p", yv[:, mappings], w)
        selected_n = np.einsum("yps,ys->p", nv[:, mappings], w)
        totals[rule] = (selected_y, selected_n)
    delta_coverage = (totals["Culture"][0] - totals["Capture"][0]) / yv.sum()
    delta_enrichment = (totals["Culture"][0] / totals["Culture"][1] -
                        totals["Capture"][0] / totals["Capture"][1]) * nv.sum() / yv.sum()
    primary = pooled[(pooled.design == "loyo87") & (pooled.k == 3)].set_index("rule")
    result = {"design": "loyo87", "k": 3, "n_permutations": len(mappings),
              "groups": [[SITES[column] for column in group] for group in groups]}
    for name, values in [("coverage", delta_coverage), ("enrichment", delta_enrichment)]:
        observed = float(primary.at["Culture", name] - primary.at["Capture", name])
        center = float(np.median(values))
        result[name] = {"observed_difference": observed, "null_median": center,
                        "null_interval": np.quantile(values, [.025, .975]).tolist(),
                        "p_two_sided": float(np.mean(abs(values - center) >= abs(observed - center) - 1e-12))}
    return pd.DataFrame({"delta_enrichment": delta_enrichment, "delta_coverage": delta_coverage}), result


def periods_and_contributions(folds, weights, joint, full):
    periods = []
    for name, start, end in [("2015-2017", 2015, 2017), ("2018-2025", 2018, 2025), ("2023-2025", 2023, 2025)]:
        subset = folds[folds.evaluation_year.between(start, end)]
        for row in pool_folds(subset).to_dict("records"):
            periods.append({"period": name, **row})
    contributions = []
    for design in ["loyo87", "forward_expanding", "forward_rolling3"]:
        frame = joint if design == "loyo87" else full[full.year >= 2018]
        w = weights[(weights.design == design) & (weights.k == 3) & weights.rule.isin(["Culture", "Capture"])]
        w = w.pivot(index=["evaluation_year", "site"], columns="rule", values="weight")
        for row in frame.itertuples():
            contribution = float((w.at[(row.year, row.site), "Culture"] - w.at[(row.year, row.site), "Capture"]) * row.positive)
            contributions.append({"design": design, "year": int(row.year), "site": row.site,
                                  "human_positive": int(row.positive), "contribution": contribution})
    return pd.DataFrame(periods), pd.DataFrame(contributions)


def posterior_support(joint, sites, out):
    """Reconstruct standardised probabilities and fixed-trio utilities from NC."""
    import xarray as xr
    rows, site_rows, hashes, rate_cache = [], [], {}, {}
    for stream in ["human", "culture"]:
        for family in ["binomial", "betabinomial"]:
            label = f"{stream}_{family}"
            metadata_path, nc_path = LEGACY / f"{label}_fit.json", LEGACY / f"{label}.nc"
            metadata = json.loads(metadata_path.read_text())
            model_sites, years = metadata["sites"], metadata["years"]
            assert model_sites == SITES and years == list(range(2015, 2026))
            num = "positive" if stream == "human" else "leptospira_culture_positive"
            den = "tested" if stream == "human" else "effective_cultures"
            observed = joint.sort_values(["year", "site"])
            with xr.open_dataset(nc_path, group="observed_data") as data:
                assert np.array_equal(data.y.values, observed[num].to_numpy())
            with xr.open_dataset(nc_path, group="posterior") as post:
                sb = np.linalg.qr(np.vstack([np.eye(7), -np.ones((1, 7))]))[0][:, :7]
                yb = np.linalg.qr(np.vstack([np.eye(10), -np.ones((1, 10))]))[0][:, :10]
                site_effects = post.site_theta.values.reshape(-1, 7) @ sb.T
                year_effects = post.year_theta.values.reshape(-1, 10) @ yb.T
                rates = expit(post.intercept.values.reshape(-1, 1, 1) + site_effects[:, :, None] + year_effects[:, None, :]).mean(axis=2)
            # Check that the saved fit was actually conditioned on these counts.
            count_text = observed[["site", "year", num, den]].to_csv(index=False)
            assert hashlib.sha256(count_text.encode()).hexdigest() == metadata["counts_sha256"]
            rate_cache[label] = rates
            hashes.update({str(metadata_path.relative_to(ROOT)): digest(metadata_path), str(nc_path.relative_to(ROOT)): digest(nc_path)})
            ranks = np.argsort(np.argsort(-rates, axis=1), axis=1) + 1
            for index, site in enumerate(SITES):
                values = rates[:, index]
                site_rows.append({"stream": stream, "family": family, "site": site,
                                  "mean": values.mean(), "median": np.median(values),
                                  "low": np.quantile(values, .025), "high": np.quantile(values, .975),
                                  "top3": np.mean(ranks[:, index] <= 3)})
    trios = {"Culture": ["Lianyuan", "Shuangfeng", "Xiangtan"], "Capture": ["Liuyang", "Shuangfeng", "Yuanjiang"]}
    testing = sites.set_index("site").loc[SITES, "human_n"].to_numpy()
    draw_rows = []
    for family in ["binomial", "betabinomial"]:
        rates = rate_cache[f"human_{family}"]
        weighted = rates * testing
        utility = {}
        for rule, selected in trios.items():
            values = weighted[:, [SITES.index(site) for site in selected]].sum(axis=1) / weighted.sum(axis=1)
            utility[rule] = values
            rows.append({"family": family, "estimand": "fixed_portfolio_expected_coverage", "rule": rule,
                         "k": 3, "draws": len(values), "mean": values.mean(), "median": np.median(values),
                         "low": np.quantile(values, .025), "high": np.quantile(values, .975),
                         "probability_gt_zero": np.nan, "probability_gt_10pp": np.nan})
        difference = utility["Culture"] - utility["Capture"]
        rows.append({"family": family, "estimand": "paired_expected_coverage_difference", "rule": "Culture-Capture",
                     "k": 3, "draws": len(difference), "mean": difference.mean(), "median": np.median(difference),
                     "low": np.quantile(difference, .025), "high": np.quantile(difference, .975),
                     "probability_gt_zero": np.mean(difference > 0), "probability_gt_10pp": np.mean(difference > .1)})
        for index, value in enumerate(difference):
            draw_rows.append({"family": family, "draw": index, "delta_coverage": value})
        if family == "betabinomial":
            existing = pd.read_csv(LEGACY / "paired_utility_draws.csv").delta_utility
            assert np.allclose(difference, existing, atol=1e-12, rtol=0)
    save(pd.DataFrame(rows), out, "posterior_support.csv")
    save(pd.DataFrame(site_rows), out, "posterior_site_rates.csv")
    save(pd.DataFrame(draw_rows), out, "posterior_paired_draws.csv")
    return rows, hashes


def validate_contract():
    # Exact equal-rate fractions must not be resolved using row order.
    toy = pd.DataFrame({"site": ["a", "b", "c", "d"], "positive": [1, 2, 3, 0],
                        "tested": [2, 4, 4, 2], "eligible_human": [True] * 4})
    sets, weights, _ = ranked_portfolios(toy, "PastHuman", list("abcd"), 2)
    assert sets == [("a", "c"), ("b", "c")]
    assert weights == {"a": .5, "b": .5, "c": 1.0, "d": 0.0}
    reverse = ranked_portfolios(toy.iloc[::-1], "PastHuman", list("dcba"), 2)
    assert reverse[0] == sets and reverse[1] == weights
    sets, weights, _ = ranked_portfolios(toy, "Random", list("abcd"), 2)
    assert len(sets) == math.comb(4, 2) and all(value == .5 for value in weights.values())
    assert math.isnan(ratios(0, 0, 100, 400)["coverage"])
    assert ratios(0, 3, 100, 400)["coverage"] == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    validate_contract()
    joint, full = load_panel("panel_complete_site_years.csv"), load_panel("panel_all_site_years.csv")
    assert len(joint) == 87 and len(full) == 88
    assert tuple(joint[COUNTS].sum().astype(int)) == (291, 8834, 317, 6824, 7333, 270469)
    assert tuple(full[["positive", "tested"]].sum().astype(int)) == (292, 8934)
    pd.testing.assert_frame_equal(joint[COUNTS + ["site", "year"]].astype({col: float for col in COUNTS}),
                                  full.loc[full.complete_site_year, COUNTS + ["site", "year"]].reset_index(drop=True).astype({col: float for col in COUNTS}))
    folds, weights = evaluate(joint, full)
    pooled, sites = pool_folds(folds), site_summary(joint)
    save(folds, out, "fold_metrics.csv")
    save(weights, out, "selection_weights.csv")
    save(pooled, out, "pooled_metrics.csv")
    save(sites, out, "site_summary.csv")
    # Write the core files first so the independent figure and audit steps can run.
    print("Core coverage CSVs written to", out, flush=True)
    period_metrics, contributions = periods_and_contributions(folds, weights, joint, full)
    save(period_metrics, out, "period_metrics.csv")
    save(contributions, out, "site_year_contributions.csv")
    permutations, permutation_summary = trajectory_permutation(joint, weights, pooled)
    save(permutations, out, "exact_permutation_k3.csv")
    dump(permutation_summary, out / "exact_permutation_summary.json")
    assert np.isclose(permutation_summary["enrichment"]["p_two_sided"], .5047619047619047, atol=1e-14)
    old = pd.read_csv(LEGACY / "primary87_loyo_by_year.csv")
    primary = folds[(folds.design == "loyo87") & (folds.k == 3)]
    for rule in ["Culture", "Capture"]:
        now = primary[primary.rule == rule].sort_values("evaluation_year")
        previous = old[old.rule == rule].sort_values("held_out_year")
        assert np.allclose(now[METRIC_COUNTS].to_numpy(), previous[METRIC_COUNTS].to_numpy())
    # Verify the previous compound indicator is not independent confirmation.
    compound = old[old.rule == "Compound"].sort_values("held_out_year")
    culture = primary[primary.rule == "Culture"].sort_values("evaluation_year")
    assert np.array_equal(compound.selected_sites, culture.selected_sites)
    assert np.allclose(compound[METRIC_COUNTS].to_numpy(), culture[METRIC_COUNTS].to_numpy())
    assert len(folds) == 707 and len(weights) == 5656 and len(pooled) == 105
    assert np.allclose(weights.groupby(["design", "rule", "k", "evaluation_year"]).weight.sum().to_numpy(),
                       weights.groupby(["design", "rule", "k", "evaluation_year"]).weight.sum().index.get_level_values("k"))
    forward = folds[folds.design.str.startswith("forward")]
    assert (forward.training_end < forward.evaluation_year).all()
    assert (forward.n_candidates == 8).all()
    assert (weights.loc[(weights.design.str.startswith("forward")) & (weights.evaluation_year == 2018) &
                        (weights.site == "Yuanjiang"), "eligible"]).all()
    for design in ["loyo87", "forward_expanding", "forward_rolling3", "pooled_fixed87"]:
        random = pooled[(pooled.design == design) & (pooled.rule == "Random")]
        if design != "loyo87":
            assert np.allclose(random.coverage, random.k / 8)
            assert np.allclose(random.testing_share, random.k / 8)
            assert np.allclose(random.enrichment, 1)
    by_year = folds[(folds.design == "forward_expanding") & (folds.k == 3)].pivot(index="evaluation_year", columns="rule", values="selected_sites")
    assert (by_year.Culture == by_year.PastHuman).all()
    posterior_rows, posterior_hashes = posterior_support(joint, sites, out)
    correlations = {
        "human_culture": exact_rank_correlation(sites.human_fraction, sites.culture_fraction),
        "human_capture": exact_rank_correlation(sites.human_fraction, sites.capture_fraction),
        "culture_capture": exact_rank_correlation(sites.culture_fraction, sites.capture_fraction),
    }
    comparison = pooled[pooled.rule.isin(["Culture", "Capture"])].pivot(index=["design", "k"], columns="rule", values="coverage")
    comparison["culture_minus_capture"] = comparison.Culture - comparison.Capture
    save(comparison.reset_index(), out, "coverage_contrasts.csv")
    summary = {
        "protocol": CONFIG,
        "panels": {"joint87": {"rows": 87, "human_positive": 291, "human_tested": 8834,
                                 "culture_positive": 317, "effective_cultures": 6824, "captures": 7333, "trap_nights": 270469},
                   "human88": {"rows": 88, "human_positive": 292, "human_tested": 8934},
                   "forward2018_2025": {"rows": 64, "human_positive": 93, "human_tested": 6525}},
        "k3_results": pooled[pooled.k == 3].to_dict("records"),
        "rank_correlations": correlations, "trajectory_permutation": permutation_summary,
        "k_by_coverage_direction": {design: {
            "culture_higher": [int(k) for k, row in sub.droplevel("design").iterrows() if row.culture_minus_capture > 1e-12],
            "capture_higher": [int(k) for k, row in sub.droplevel("design").iterrows() if row.culture_minus_capture < -1e-12],
        } for design, sub in comparison.groupby(level="design")},
        "forward_culture_past_human_same_k3_years": int((by_year.Culture == by_year.PastHuman).sum()),
        "nonrandom_tied_folds": int(((folds.rule != "Random") & (folds.n_tied_sets > 1)).sum()),
        "posterior_support": pd.DataFrame(posterior_rows).replace({np.nan: None}).to_dict("records"),
    }
    dump(CONFIG, out / "analysis_protocol.json")
    dump(summary, out / "analysis_summary.json")
    inputs = {str(path.relative_to(ROOT)): digest(path) for path in [DATA / "panel_complete_site_years.csv", DATA / "panel_all_site_years.csv", Path(__file__)]}
    inputs.update(posterior_hashes)
    dump({"status": "passed", "checks": ["source_panels_agree", "panel_totals", "legacy_k3_reproduced",
          "exact_5040_trajectory_permutations", "ties_use_exact_rational_rates", "tie_order_invariance",
          "random_portfolio_expectation", "zero_positive_year_contract", "forward_training_precedes_evaluation",
          "forward_yuanjiang_2018_eligible", "selection_weights_sum_to_k", "saved_posterior_counts_and_utilities"],
          "fold_rows": len(folds), "selection_weight_rows": len(weights), "pooled_rows": len(pooled),
          "input_sha256": inputs, "output_sha256": {path.name: digest(path) for path in sorted(out.glob("*.csv"))}}, out / "validation.json")
    (out / "README.md").write_text(
        "# Sentinel-network redesign: 2026-11-13\n\n"
        "Run `.venv/bin/python script/hunan/analyze_sentinel_network.py` from the project root.\n\n"
        "The new descriptive evaluations are recomputed from harmonised counts. The original 87 jointly observed site-years support LOYO and fixed historical rankings. Forward evaluation scores all 64 human site-years in 2018–2025, including Yuanjiang 2018, and uses only earlier valid records for each indicator.\n\n"
        "`fold_metrics.csv` contains per-year counts and metrics; `selection_weights.csv` gives every site's exact tie-averaged inclusion weight. `pooled_metrics.csv` pools numerators and denominators, while `period_metrics.csv` partitions evaluation years without reselecting within periods. `HumanPooled` is the in-sample human-rate ranking; `PastHuman` is the forward human-history comparator. `Random` is the expectation over eligible portfolios, not one random draw. Undefined annual coverage is stored as a missing value when no human positives occur.\n\n"
        "`pooled_fixed87` uses evaluation_year=0 to denote the full 2015–2025 record. Its selection is in-sample. `loyo87` trains on all other years and is historical representation. `forward_expanding` and `forward_rolling3` obey calendar order. All methods evaluate k=1–7.\n\n"
        "`exact_permutation_k3.csv` and its JSON summary reproduce the primary historical culture-versus-capture comparison with 5,040 complete-trajectory permutations. No multiple-k search for statistical significance is performed.\n\n"
        "`posterior_support.csv`, `posterior_site_rates.csv` and `posterior_paired_draws.csv` are reconstructed from the four saved primary model fits under `priority_87_2026-10-22`; no new posterior sampling is performed. The model endpoint remains the 87-row panel and should not be presented as a model of forward outcomes.\n\n"
        "`analysis_protocol.json` records the chosen design; `analysis_summary.json` gives main findings; `validation.json` records checks and source/output hashes.\n",
        encoding="utf-8")
    print(pooled[pooled.k == 3][["design", "rule"] + METRIC_COUNTS + ["coverage", "testing_share", "enrichment"]].to_string(index=False))
    print("All checks passed.")


if __name__ == "__main__":
    main()
