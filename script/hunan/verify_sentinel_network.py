"""Independently verify the 2026-11-13 sentinel-network analysis.

The oracle enumerates every possible k-site combination and maximises the sum
of exact rational rates. It does not import or call the analysis implementation.
It also reconstructs the k=3 trajectory permutation distribution directly.

Run from the repository root with .venv/bin/python. --data-only writes the
published-data reconciliation before the network outputs are available.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/processed/harmonized_2026-10-22"
OUT = ROOT / "output/hunan/current/network_redesign_2026-11-13"
OLD = ROOT / "output/hunan/current/priority_87_2026-10-22"
LABELS = {
    "双峰县": "Shuangfeng", "君山区": "Junshan", "宁乡县": "Ningxiang",
    "沅江市": "Yuanjiang", "浏阳市": "Liuyang", "涟源市": "Lianyuan",
    "湘潭县": "Xiangtan", "醴陵市": "Liling",
}
STREAMS = {
    "Culture": ("leptospira_culture_positive", "effective_cultures", "eligible_culture"),
    "Capture": ("rodents_caught", "trap_setting_count", "eligible_capture"),
    "PastHuman": ("positive", "tested", "eligible_human"),
    "HumanPooled": ("positive", "tested", "eligible_human"),
}
KEYS = ["design", "rule", "k", "evaluation_year"]
NUMERIC = ["positive_selected", "positive_all", "tested_selected", "tested_all",
           "coverage", "testing_share", "enrichment"]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    full = pd.read_csv(DATA / "panel_all_site_years.csv")
    full["site"] = full.site.replace(LABELS)
    complete = pd.read_csv(DATA / "panel_complete_site_years.csv")
    complete["site"] = complete.site.replace(LABELS)
    assert len(full) == 88 and len(complete) == 87
    assert not full.duplicated(["site", "year"]).any()
    assert tuple(full[["positive", "tested"]].sum()) == (292, 8934)
    assert tuple(complete[["positive", "tested"]].sum()) == (291, 8834)
    assert tuple(complete[["leptospira_culture_positive", "effective_cultures",
                           "rodents_caught", "trap_setting_count"]].sum()) == (
                               317, 6824, 7333, 270469)
    missing = full.loc[~full.complete_site_year]
    assert list(zip(missing.site, missing.year)) == [("Yuanjiang", 2018)]
    assert tuple(missing[["positive", "tested"]].iloc[0]) == (1, 100)
    pd.testing.assert_frame_equal(
        full.loc[full.complete_site_year].sort_values(["site", "year"]).reset_index(drop=True),
        complete.sort_values(["site", "year"]).reset_index(drop=True),
        check_dtype=False,
    )
    return full, complete


def reconcile_publication(full: pd.DataFrame, complete: pd.DataFrame) -> None:
    # Transcribed from Tables 2 and 3 of the publisher's three-page PDF.
    # The percentages printed alongside human totals are not recalculated here.
    published = pd.DataFrame({
        "year": [2018, 2019, 2020, 2021, 2022],
        "published_human_positive": [58, 15, 25, 3, 1],
        "published_human_tested": [803, 834, 805, 810, 820],
        "published_human_rate_percent": [4.86, 1.80, 2.98, 0.37, 0.12],
        "published_culture_positive": [29, 14, 26, 28, 10],
        "published_effective_cultures": [500, 599, 598, 596, 607],
        "published_rodents_caught": [529, 671, 635, 636, 615],
        "published_trap_setting_count": [23812, 25798, 20568, 30094, 24401],
    })
    columns = ["positive", "tested", "leptospira_culture_positive", "effective_cultures",
               "rodents_caught", "trap_setting_count"]
    for frame, prefix in [(full, "current_all"), (complete, "current_joint")]:
        totals = frame.loc[frame.year.between(2018, 2022)].groupby("year")[columns].sum()
        totals = totals.rename(columns={name: f"{prefix}_{name}" for name in columns})
        published = published.merge(totals, on="year", validate="one_to_one")
    for original, current in [
        ("published_human_tested", "current_all_tested"),
        ("published_culture_positive", "current_all_leptospira_culture_positive"),
        ("published_effective_cultures", "current_all_effective_cultures"),
        ("published_rodents_caught", "current_all_rodents_caught"),
        ("published_trap_setting_count", "current_all_trap_setting_count"),
    ]:
        assert np.array_equal(published[original], published[current]), original
    np.testing.assert_allclose(
        100 * published.current_all_positive / published.current_all_tested,
        published.published_human_rate_percent, atol=0.005, rtol=0,
    )
    assert published.published_human_positive.sum() == 102
    assert published.current_all_positive.sum() == 82
    assert published.current_joint_positive.sum() == 81
    published["published_minus_current_human_positive"] = (
        published.published_human_positive - published.current_all_positive
    )
    published.to_csv(OUT / "publication_overlap_counts.csv", index=False)
    pdf = ROOT / "manuscript/literature/zhang_2025_hunan_2018_2022.pdf"
    table = published[["year", "published_human_positive", "current_all_positive",
                       "current_all_tested", "current_joint_positive", "current_joint_tested"]]
    note = f"""# Published-data reconciliation

The 2025 Hunan report (DOI 10.3969/j.issn.1006-3110.2025.04.015)
was checked against its three-page publisher PDF, Tables 2 and 3. This document
records data provenance; it does not change the source workbooks or infer why
the originally printed human counts differ.

{table.to_markdown(index=False)}

For 2018–2022, the published rodent capture counts, trapping effort, effective
culture counts and culture positives agree with the present records in every
year (107/2,900 cultures; 3,086 captures/124,673 trap-nights).
The human testing total also agrees (4,072). The paper printed 102 human positives,
whereas the current site and age ledgers agree on 82. The paper's printed total
rates for 2018 (4.86%) and 2020 (2.98%) agree with 39/803 and 24/805, respectively,
rather than its printed positive counts of 58 and 25. The published overall
2.50% instead corresponds to 102/4,072. These inconsistencies are documented,
without attributing the difference to cross-reactivity or duplicated individuals.

The current 87-row joint panel omits Yuanjiang's 2018 human 1/100 because animal
surveillance was not conducted then. Its overlapping-period human total is
81/3,972. The forward design retains the human observation because prior animal
records are sufficient for selecting that site. It evaluates 93/6,525 human
results in 2018–2025.

The joint panel includes 39 site–years within the previous publication's period
and 48 outside that period. Source sharing is acknowledged; the new comparison
concerns site selection and coverage rather than reproducing the prior descriptive
epidemiology. The data-holder-approved counting choice and original backups are
documented in `data/corrections/2026-10-22/更正记录.md`.

Publisher PDF SHA-256: `{digest(pdf)}`.
Full panel SHA-256: `{digest(DATA / 'panel_all_site_years.csv')}`.
Joint panel SHA-256: `{digest(DATA / 'panel_complete_site_years.csv')}`.
"""
    (OUT / "data_reconciliation.md").write_text(note, encoding="utf-8")


def all_best_subsets(rates: dict[str, Fraction], k: int, random: bool = False):
    """Exhaustive optimisation; does not determine a kth-ranked threshold."""
    options = list(itertools.combinations(sorted(rates), k))
    assert options
    if not random:
        objectives = [sum((rates[site] for site in subset), start=Fraction(0))
                      for subset in options]
        optimum = max(objectives)
        options = [subset for subset, score in zip(options, objectives) if score == optimum]
    weights = {site: sum(site in subset for subset in options) / len(options)
               for site in rates}
    assert abs(sum(weights.values()) - k) < 1e-12
    return options, weights


def metrics(frame: pd.DataFrame, weights: dict[str, float]) -> dict[str, float]:
    selected_positive = sum(float(row.positive) * weights.get(row.site, 0)
                            for row in frame.itertuples())
    selected_tested = sum(float(row.tested) * weights.get(row.site, 0)
                          for row in frame.itertuples())
    all_positive, all_tested = frame[["positive", "tested"]].sum()
    coverage = selected_positive / all_positive if all_positive else np.nan
    share = selected_tested / all_tested
    return dict(positive_selected=selected_positive, positive_all=float(all_positive),
                tested_selected=selected_tested, tested_all=float(all_tested),
                coverage=coverage, testing_share=share, enrichment=coverage / share)


def oracle(full: pd.DataFrame, complete: pd.DataFrame):
    folds, weights = [], []
    designs = {
        "loyo87": (range(2015, 2026), ["Culture", "Capture", "Random"]),
        "forward_expanding": (range(2018, 2026), ["Culture", "Capture", "PastHuman", "Random"]),
        "forward_rolling3": (range(2018, 2026), ["Culture", "Capture", "PastHuman", "Random"]),
        "pooled_fixed87": ([0], ["Culture", "Capture", "HumanPooled", "Random"]),
    }
    all_sites = sorted(full.site.unique())
    for design, (years, rules) in designs.items():
        for year in years:
            if design == "loyo87":
                train = complete.loc[complete.year != year]
                target = complete.loc[complete.year == year]
            elif design == "pooled_fixed87":
                train = target = complete
            else:
                train = full.loc[full.year < year]
                if design == "forward_rolling3":
                    train = train.loc[train.year >= year - 3]
                target = full.loc[full.year == year]
                assert train.year.max() < year
            base_candidates = sorted(train.site.unique() if design.startswith("forward")
                                     else target.site.unique())
            for rule in rules:
                candidates = base_candidates.copy()
                if rule == "Random":
                    stats = {site: (np.nan, np.nan) for site in all_sites}
                    rates = {site: Fraction(0) for site in candidates}
                else:
                    numerator, denominator, valid = STREAMS[rule]
                    grouped = train.loc[train[valid]].groupby("site")[[numerator, denominator]].sum()
                    stats = {site: (int(grouped.loc[site, numerator]), int(grouped.loc[site, denominator]))
                             for site in grouped.index}
                    # Forward availability must be derived from past records.
                    if design.startswith("forward"):
                        candidates = sorted(site for site in all_sites if stats.get(site, (0, 0))[1] > 0)
                    rates = {site: Fraction(*stats[site]) for site in candidates}
                for k in range(1, 8):
                    subsets, selection = all_best_subsets(rates, k, random=rule == "Random")
                    key = dict(design=design, rule=rule, k=k, evaluation_year=year)
                    folds.append({**key, **metrics(target, selection),
                                  "n_tied_sets": len(subsets), "n_candidates": len(candidates),
                                  "training_start": int(train.year.min()),
                                  "training_end": int(train.year.max()),
                                  "training_years": ";".join(map(str, sorted(train.year.unique())))})
                    for site in all_sites:
                        n, d = stats.get(site, (0, 0))
                        weights.append({**key, "site": site, "eligible": site in candidates,
                                        "weight": selection.get(site, 0),
                                        "training_numerator": n, "training_denominator": d,
                                        "training_rate": n / d if d > 0 else np.nan})
    return pd.DataFrame(folds), pd.DataFrame(weights)


def compare_columns(actual: pd.DataFrame, expected: pd.DataFrame, keys: list[str],
                    columns: list[str], description: str) -> int:
    assert not actual.duplicated(keys).any(), (description, "duplicate keys")
    a = actual.set_index(keys).sort_index()
    e = expected.set_index(keys).sort_index()
    assert a.index.equals(e.index), (description, "key mismatch", a.shape, e.shape)
    for name in columns:
        np.testing.assert_allclose(a[name].to_numpy(dtype=float), e[name].to_numpy(dtype=float),
                                   rtol=1e-11, atol=1e-11, equal_nan=True,
                                   err_msg=f"{description}: {name}")
    return len(a) * len(columns)


def pooled(folds: pd.DataFrame) -> pd.DataFrame:
    result = folds.groupby(["design", "rule", "k"])[NUMERIC[:4]].sum().reset_index()
    result["coverage"] = result.positive_selected / result.positive_all
    result["testing_share"] = result.tested_selected / result.tested_all
    result["enrichment"] = result.coverage / result.testing_share
    return result


def permutation(complete: pd.DataFrame, expected_weights: pd.DataFrame) -> tuple[dict, np.ndarray]:
    sites = sorted(complete.site.unique())
    years = sorted(complete.year.unique())
    variable = [site for site in sites if site != "Yuanjiang"]
    outcome = {(row.year, row.site): (float(row.positive), float(row.tested))
               for row in complete.itertuples()}
    selectors = expected_weights.query('design == "loyo87" and k == 3 and rule != "Random"')
    selectors = selectors.set_index(["rule", "evaluation_year", "site"]).weight.to_dict()
    deltas = []
    for assignment in itertools.permutations(variable):
        mapping = dict(zip(variable, assignment))
        mapping["Yuanjiang"] = "Yuanjiang"
        totals = {rule: [0., 0.] for rule in ["Culture", "Capture"]}
        for year in years:
            for target in sites:
                source = mapping[target]
                if (year, source) not in outcome:
                    continue
                positives, tests = outcome[(year, source)]
                for rule in totals:
                    weight = selectors[(rule, year, target)]
                    totals[rule][0] += positives * weight
                    totals[rule][1] += tests * weight
        cp, cn = totals["Culture"]
        ap, an = totals["Capture"]
        deltas.append(((cp / cn - ap / an) * 8834 / 291, (cp - ap) / 291))
    deltas = np.asarray(deltas)
    reference = pooled(oracle_results_for_permutation(complete, expected_weights)).set_index("rule")
    observed = [reference.loc["Culture", "enrichment"] - reference.loc["Capture", "enrichment"],
                reference.loc["Culture", "coverage"] - reference.loc["Capture", "coverage"]]
    result = {"permutations": len(deltas)}
    for i, name in enumerate(["enrichment", "coverage"]):
        center = float(np.median(deltas[:, i]))
        p = float(np.mean(np.abs(deltas[:, i] - center) >= abs(observed[i] - center) - 1e-12))
        result[name] = dict(observed_difference=float(observed[i]), null_median=center,
                            null_interval=np.quantile(deltas[:, i], [.025, .975]).tolist(),
                            p_two_sided=p)
    assert result["permutations"] == 5040
    assert np.isclose(result["enrichment"]["p_two_sided"], 0.5047619047619047)
    assert np.isclose(result["coverage"]["p_two_sided"], 0.4714285714285714)
    return result, deltas


def verify_posterior_outputs() -> int:
    """Compare new summaries with the legacy outputs already rebuilt from NC draws."""
    new = pd.read_csv(OUT / "posterior_support.csv")
    original = pd.read_csv(OLD / "posterior_decisions.csv").set_index(["fit", "strategy"])
    contrasts = pd.read_csv(OLD / "posterior_decision_sensitivity.csv").set_index("fit")
    checked = 0
    assert len(new) == 6 and (new.draws == 6000).all() and (new.k == 3).all()
    for row in new.itertuples():
        fit = "human_" + row.family
        old = (contrasts.loc[fit] if row.estimand == "paired_expected_coverage_difference"
               else original.loc[(fit, row.rule)])
        for column in ["mean", "median", "low", "high"]:
            assert np.isclose(getattr(row, column), old[column], rtol=1e-11, atol=1e-12)
            checked += 1
        if row.estimand == "paired_expected_coverage_difference":
            assert np.isclose(row.probability_gt_zero, old.probability_culture_better)
            checked += 1
    new_sites = pd.read_csv(OUT / "posterior_site_rates.csv")
    old_sites = pd.read_csv(OLD / "posterior_site_rates.csv").rename(columns={"model": "family"})
    checked += compare_columns(new_sites, old_sites, ["stream", "family", "site"],
                               ["mean", "median", "low", "high", "top3"], "posterior site rates")
    return checked


def oracle_results_for_permutation(complete, weights):
    rows = []
    for (rule, year), frame in weights.query('design == "loyo87" and k == 3 and rule != "Random"').groupby(["rule", "evaluation_year"]):
        rows.append({"design": "loyo87", "rule": rule, "k": 3,
                     **metrics(complete.loc[complete.year == year], dict(zip(frame.site, frame.weight)))})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-only", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    full, complete = load_data()
    reconcile_publication(full, complete)
    if args.data_only:
        print("PASS: 87/88-row panels and published-period count reconciliation.")
        return
    expected, expected_weights = oracle(full, complete)
    actual = pd.read_csv(OUT / "fold_metrics.csv")
    actual_weights = pd.read_csv(OUT / "selection_weights.csv")
    actual_pooled = pd.read_csv(OUT / "pooled_metrics.csv")
    checked = compare_columns(actual, expected, KEYS, NUMERIC + ["n_tied_sets", "n_candidates", "training_start", "training_end"], "fold metrics")
    checked += compare_columns(actual_weights, expected_weights, KEYS + ["site"], ["eligible", "weight"], "selection weights")
    valid_a = actual_weights.loc[actual_weights.eligible & actual_weights.rule.ne("Random")]
    valid_e = expected_weights.loc[expected_weights.eligible & expected_weights.rule.ne("Random")]
    checked += compare_columns(valid_a, valid_e, KEYS + ["site"], ["training_numerator", "training_denominator", "training_rate"], "training rates")
    checked += compare_columns(actual_pooled, pooled(expected), ["design", "rule", "k"], NUMERIC, "pooled metrics")
    site_expected = complete.groupby("site").agg(
        human_pos=("positive", "sum"), human_n=("tested", "sum"),
        culture_pos=("leptospira_culture_positive", "sum"), culture_n=("effective_cultures", "sum"),
        caught=("rodents_caught", "sum"), traps=("trap_setting_count", "sum"),
        observed_years=("year", "nunique"),
    ).reset_index()
    for stream, n, d in [("human", "human_pos", "human_n"), ("culture", "culture_pos", "culture_n"),
                         ("capture", "caught", "traps")]:
        site_expected[stream + "_fraction"] = site_expected[n] / site_expected[d]
        site_expected[stream + "_rank"] = site_expected[stream + "_fraction"].rank(ascending=False)
    checked += compare_columns(pd.read_csv(OUT / "site_summary.csv"), site_expected,
                               ["site"], list(site_expected.columns[1:]), "site summaries")
    checked += verify_posterior_outputs()
    # Actual selection records must be independent of the evaluation-year animal gap.
    gap = actual_weights.query('design.str.startswith("forward") and evaluation_year == 2018 and site == "Yuanjiang"', engine="python")
    assert gap.eligible.all()
    forward = actual.loc[actual.design.str.startswith("forward")]
    assert (forward.training_end < forward.evaluation_year).all()
    for row in forward.itertuples():
        training_years = set(map(int, str(row.training_years).replace(",", ";").split(";")))
        assert max(training_years) < row.evaluation_year
        if row.design == "forward_rolling3":
            assert training_years == set(range(row.evaluation_year - 3, row.evaluation_year))
    all_selected = actual_weights.query('design == "loyo87" and evaluation_year == 2018 and k == 7')
    assert np.allclose(all_selected.loc[all_selected.eligible, "weight"], 1)
    assert np.allclose(all_selected.loc[~all_selected.eligible, "weight"], 0)
    assert (actual.loc[(actual.positive_all == 0), "coverage"].isna()).all()
    assert (actual.loc[(actual.positive_all == 0), "enrichment"].isna()).all()
    rng = actual.loc[actual.rule.eq("Random")]
    assert np.allclose(rng.coverage.dropna(), (rng.k / rng.n_candidates).loc[rng.coverage.notna()])
    assert np.allclose(rng.testing_share, rng.k / rng.n_candidates)
    # Meaningful edge case for the exhaustive oracle: equal fractions at a boundary.
    options, w = all_best_subsets({"a": Fraction(1, 2), "b": Fraction(1, 3), "c": Fraction(2, 6), "d": Fraction(0)}, 2)
    assert set(options) == {("a", "b"), ("a", "c")} and w == {"a": 1., "b": .5, "c": .5, "d": 0.}
    print(f"PASS: {len(expected)} folds, {len(expected_weights)} site weights, {checked} numerical checks.", flush=True)
    perm, distribution = permutation(complete, expected_weights)
    legacy = json.loads((OLD / "primary87_summary.json").read_text())["permutation"]
    current_permutation = json.loads((OUT / "exact_permutation_summary.json").read_text())
    assert current_permutation["n_permutations"] == perm["permutations"]
    current_distribution = pd.read_csv(OUT / "exact_permutation_k3.csv")
    for name in ["enrichment", "coverage"]:
        for field in ["observed_difference", "null_median", "null_interval", "p_two_sided"]:
            np.testing.assert_allclose(perm[name][field], legacy[name][field], rtol=1e-11, atol=1e-12)
            np.testing.assert_allclose(perm[name][field], current_permutation[name][field], rtol=1e-11, atol=1e-12)
    # Sort paired values to avoid depending on either implementation's enumeration order.
    current_pairs = current_distribution[["delta_enrichment", "delta_coverage"]].to_numpy()
    independent_pairs = distribution
    np.testing.assert_allclose(
        sorted(map(tuple, current_pairs)), sorted(map(tuple, independent_pairs)),
        rtol=1e-11, atol=1e-12,
    )
    report = dict(status="passed", implementation="independent exhaustive combination optimisation using rational rates",
                  folds_verified=len(expected), selection_weights_verified=len(expected_weights),
                  pooled_rows_verified=len(actual_pooled), numerical_checks=checked,
                  fractional_tie_folds=int(((actual.n_tied_sets > 1) & actual.rule.ne("Random")).sum()),
                  zero_positive_folds=int((actual.positive_all == 0).sum()),
                  forward_timing="all training years precede evaluation; Yuanjiang 2018 retained",
                  publication_overlap="rodent counts and human testing totals match; human numerator discrepancy documented",
                  posterior_support="new summaries agree with raw-NC-verified legacy outputs",
                  independent_trajectory_permutation=perm,
                  legacy_recomputation_report="manuscript/tables/supplementary_recompute_validation.json")
    paths = [Path(__file__), DATA / "panel_all_site_years.csv", DATA / "panel_complete_site_years.csv",
             OUT / "fold_metrics.csv", OUT / "selection_weights.csv", OUT / "pooled_metrics.csv",
             OUT / "publication_overlap_counts.csv", OUT / "site_summary.csv", OUT / "posterior_support.csv",
             OUT / "posterior_site_rates.csv", OUT / "exact_permutation_summary.json", OUT / "exact_permutation_k3.csv",
             ROOT / "manuscript/tables/supplementary_recompute_validation.json"]
    report["inputs_sha256"] = {str(p.relative_to(ROOT)): digest(p) for p in paths}
    (OUT / "independent_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
