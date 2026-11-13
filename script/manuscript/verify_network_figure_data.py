"""Independently verify ggplot2 figure provenance and plotted numerical data.

Checks the seven rendered figures, all labelled panels, copied inputs and source
tables. Derived coverage, membership, histogram, and posterior probabilities are
recalculated in Python from the analysis CSVs or saved posterior draws, without
calling the R plotting functions. No source data or analysis result is modified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from PIL import Image
from scipy.special import expit
import xarray as xr


ROOT = Path(__file__).resolve().parents[2]
NETWORK = ROOT / "output/hunan/current/network_redesign_2026-11-13"
LEGACY = ROOT / "output/hunan/current/priority_87_2026-10-22"
DEFAULT_SOURCE = ROOT / "manuscript/figures/source_data/network_ggplot_2026-11-13"
FIGURES = [
    "figure1_indicator_priorities", "figure2_network_coverage", "figure3_forward_evaluation",
    "figure_s1_predictive_checks", "figure_s2_period_prior_sensitivity",
    "figure_s3_network_sensitivity", "figure4_model_uncertainty",
]
MODEL_NAMES = {"binomial": "Binomial", "betabinomial": "Beta-binomial"}
METRICS = ["positive_selected", "positive_all", "tested_selected", "tested_all",
           "coverage", "testing_share", "enrichment"]
CHECKS = 0
SECTIONS: dict[str, int] = {}
USED: set[Path] = set()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> pd.DataFrame:
    USED.add(path)
    return pd.read_csv(path)


def numerical(actual, expected, label):
    global CHECKS
    a = np.asarray(actual, dtype=float)
    e = np.asarray(expected, dtype=float)
    assert a.shape == e.shape, (label, a.shape, e.shape)
    np.testing.assert_allclose(a, e, rtol=1e-10, atol=1e-11, equal_nan=True, err_msg=label)
    CHECKS += a.size


def frame(actual: pd.DataFrame, expected: pd.DataFrame, keys: list[str], columns: list[str], label: str):
    a = actual.sort_values(keys).reset_index(drop=True)
    e = expected.sort_values(keys).reset_index(drop=True)
    assert len(a) == len(e), (label, len(a), len(e))
    for col in dict.fromkeys(keys + columns):
        assert col in a and col in e, (label, col)
        if pd.api.types.is_numeric_dtype(e[col]) or pd.api.types.is_bool_dtype(e[col]):
            numerical(a[col], e[col], f"{label}: {col}")
        else:
            pd.testing.assert_series_equal(a[col].astype("string"), e[col].astype("string"),
                                           check_names=False, obj=f"{label}: {col}")


def wilson(y, n):
    y, n = np.asarray(y, float), np.asarray(n, float)
    z = 1.959963984540054
    p = y / n
    denominator = 1 + z*z/n
    center = (p + z*z/(2*n)) / denominator
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denominator
    return center-half, center+half


def verify_manifest(source: Path):
    import pymupdf

    path = source / "manifest.json"
    USED.add(path)
    manifest = json.loads(path.read_text())
    assert manifest["engine"] == "ggplot2" and manifest["composition"] == "patchwork"
    assert set(manifest["figures"]) == set(FIGURES)
    assert len(manifest["outputs_sha256"]) == 14
    assert manifest["script"] == "script/hunan/plot_sentinel_network.R"
    assert sha256(ROOT/manifest["script"]) == manifest["script_sha256"]
    for name, digest in manifest["scripts_sha256"].items():
        assert sha256(ROOT/name) == digest, f"Script hash differs: {name}"
    for name, digest in manifest["inputs_sha256"].items():
        assert sha256(ROOT/name) == digest, f"Input hash differs: {name}"
        assert sha256(source/manifest["input_copies"][name]) == digest, f"Input copy differs: {name}"
    for name, digest in manifest["panel_data_sha256"].items():
        assert sha256(source/name) == digest, f"Panel-data hash differs: {name}"
    for name, digest in manifest["outputs_sha256"].items():
        assert sha256(source.parent.parent/name) == digest, f"Figure hash differs: {name}"
    panels = manifest["panels"]
    assert len({(p["figure"], p["panel"]) for p in panels}) == len(panels)
    for stem, spec in manifest["figures"].items():
        tags = sorted(p["panel"] for p in panels if p["figure"] == stem)
        assert len(tags) >= 6 and tags == list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(tags)])
        assert tags == spec["panels"] and len(tags) == spec["panel_count"]
        assert spec["engine"] == "ggplot2" and spec["composition"] == "patchwork"
        for panel in [p for p in panels if p["figure"] == stem]:
            assert panel["engine"] == "ggplot2"
            assert panel["source"] in manifest["panel_data_sha256"]
            assert len(pd.read_csv(source/panel["source"])) == panel["rows"]
        image = Image.open(source.parent.parent/f"{stem}.png")
        assert abs(image.width - spec["width_inches"]*spec["dpi"]) <= 2
        assert abs(image.height - spec["height_inches"]*spec["dpi"]) <= 2
        pdf = pymupdf.open(source.parent.parent/f"{stem}.pdf")
        assert len(pdf) == 1
        words = {w[4] for w in pdf[0].get_text("words")}
        assert set(tags).issubset(words), (stem, "Rendered panel tags missing", set(tags)-words)
        assert abs(pdf[0].rect.width - spec["width_inches"]*72) < 1
        assert abs(pdf[0].rect.height - spec["height_inches"]*72) < 1
    return manifest


def verify_main_figures(source, sites, pooled, folds, weights):
    start = CHECKS
    panel = lambda stem, tag: read(source/f"{stem}_{tag}.csv")
    f1 = "figure1_indicator_priorities"
    sites = sites.copy()
    sites["member_group"] = np.select([
        (sites.culture_rank <= 3) & (sites.capture_rank <= 3),
        sites.culture_rank <= 3, sites.capture_rank <= 3,
    ], ["Shared", "Culture only", "Capture only"], default="Neither")
    rates = []
    for stream in ["human", "culture", "capture"]:
        data = sites[["site", "member_group"]].copy()
        data["stream"] = stream.title()
        for field in ["rank", "low", "high"]:
            data[field] = sites[f"{stream}_{field}"]
        data["rate"] = sites[f"{stream}_fraction"]
        data["selected"] = (data["rank"] <= 3).astype(int)
        rates.append(data)
    rates = pd.concat(rates, ignore_index=True)
    for tag in "AB":
        frame(panel(f1, tag), rates, ["site", "stream"],
              ["rate", "low", "high", "rank", "member_group"], f"Figure 1{tag}")
    frame(panel(f1, "C"), rates, ["site", "stream"], ["selected", "rank"], "Figure 1C")
    counts = read(ROOT/"manuscript/tables/table_s2_site_year_counts.csv")
    selected = weights.query('design=="loyo87" and k==3 and rule in ["Culture","Capture"]')
    selected = selected.pivot(index=["evaluation_year", "site"], columns="rule", values="weight").reset_index()
    contrib = counts.merge(selected, left_on=["year", "site"], right_on=["evaluation_year", "site"], validate="one_to_one")
    contrib["contribution"] = (contrib.Culture-contrib.Capture)*contrib.positive
    contribution = contrib.groupby("site", as_index=False).contribution.sum()
    numerical(contribution.contribution.sum(), 63, "LOYO contribution total")
    frame(panel(f1, "D"), contribution, ["site"], ["contribution"], "Figure 1D")
    groups = sites.groupby("member_group").agg(
        human_pos=("human_pos", "sum"), human_n=("human_n", "sum"),
        culture_pos=("culture_pos", "sum"), culture_n=("culture_n", "sum"),
        n_sites=("site", "size"), sites=("site", lambda x: ";".join(sorted(x))),
    ).reset_index()
    groups["rate"] = groups.human_pos/groups.human_n
    groups["low"], groups["high"] = wilson(groups.human_pos, groups.human_n)
    frame(panel(f1, "E"), groups, ["member_group"],
          ["sites", "n_sites", "human_pos", "human_n", "culture_pos", "culture_n", "rate", "low", "high"], "Figure 1E")
    composition = groups.melt(id_vars="member_group", value_vars=["human_pos", "culture_pos", "human_n"],
                              var_name="quantity", value_name="count")
    composition["total"] = composition.groupby("quantity")["count"].transform("sum")
    composition["share"] = composition["count"]/composition.total
    composition["quantity"] = composition.quantity.map({
        "human_pos": "Human\npositives", "culture_pos": "Culture\npositives", "human_n": "Human\ntests"})
    frame(panel(f1, "F"), composition, ["member_group", "quantity"], ["count", "total", "share"], "Figure 1F")
    SECTIONS["figure1_membership_and_contributions"] = CHECKS-start

    start = CHECKS
    f2 = "figure2_network_coverage"
    main = pooled.query('design=="loyo87" and rule in ["Culture","Capture","Random"]')
    for tag in "ABCD":
        frame(panel(f2, tag), main, ["rule", "k"], METRICS, f"Figure 2{tag}")
    increments = main.query('rule!="Random"').sort_values(["rule", "k"]).copy()
    increments["previous_coverage"] = increments.groupby("rule").coverage.shift().fillna(0)
    increments["added_coverage"] = increments.coverage-increments.previous_coverage
    increments["added_coverage_pp"] = increments.added_coverage*100
    frame(panel(f2, "E"), increments, ["rule", "k"],
          ["coverage", "previous_coverage", "added_coverage", "added_coverage_pp"], "Figure 2E")
    period_data = read(NETWORK/"period_metrics.csv").query('design=="loyo87" and rule in ["Culture","Capture"]')
    contrast = period_data.pivot(index=["period", "k", "positive_all"], columns="rule", values="coverage").reset_index()
    contrast["difference_pp"] = 100*(contrast.Culture-contrast.Capture)
    frame(panel(f2, "F"), contrast, ["period", "k"],
          ["positive_all", "Culture", "Capture", "difference_pp"], "Figure 2F")
    SECTIONS["figure2_size_and_period_derivations"] = CHECKS-start

    start = CHECKS
    f3 = "figure3_forward_evaluation"
    timeline = folds.query('design in ["forward_expanding","forward_rolling3"] and rule=="Culture" and k==3')
    frame(panel(f3, "A"), timeline, ["design", "evaluation_year"], ["training_start", "training_end"], "Figure 3A")
    assert (timeline.training_end < timeline.evaluation_year).all()
    selection = weights.query('design=="forward_expanding" and k==3 and rule in ["Culture","Capture","PastHuman"]')
    frame(panel(f3, "B"), selection, ["rule", "evaluation_year", "site"],
          ["eligible", "weight", "training_numerator", "training_denominator", "training_rate"], "Figure 3B")
    forward = folds.query('design=="forward_expanding" and k==3').copy()
    frame(panel(f3, "C"), forward.query('rule!="Random"'), ["rule", "evaluation_year"], METRICS, "Figure 3C")
    cumulative = forward.sort_values(["rule", "evaluation_year"])
    cumulative["cumulative_selected"] = cumulative.groupby("rule").positive_selected.cumsum()
    cumulative["cumulative_total"] = cumulative.groupby("rule").positive_all.cumsum()
    cumulative["cumulative_coverage"] = cumulative.cumulative_selected/cumulative.cumulative_total
    frame(panel(f3, "D"), cumulative, ["rule", "evaluation_year"],
          ["cumulative_selected", "cumulative_total", "cumulative_coverage"], "Figure 3D")
    human = pooled.query('design=="forward_expanding" and rule=="PastHuman"')[["k", "coverage"]].rename(columns={"coverage": "human_coverage"})
    relative = pooled.query('design=="forward_expanding" and rule in ["Culture","Capture"]').merge(human, on="k", validate="many_to_one")
    relative["difference_pp"] = 100*(relative.coverage-relative.human_coverage)
    frame(panel(f3, "E"), relative, ["rule", "k"], ["coverage", "human_coverage", "difference_pp"], "Figure 3E")
    update = pooled.query('design in ["forward_expanding","forward_rolling3"] and rule in ["Culture","Capture","PastHuman"]')
    update = update.pivot(index=["rule", "k"], columns="design", values="coverage").reset_index()
    update["difference_pp"] = 100*(update.forward_rolling3-update.forward_expanding)
    frame(panel(f3, "F"), update, ["rule", "k"],
          ["forward_expanding", "forward_rolling3", "difference_pp"], "Figure 3F")
    SECTIONS["figure3_chronological_derivations"] = CHECKS-start


def verify_support(source, pooled):
    panel = lambda stem, tag: read(source/f"{stem}_{tag}.csv")
    start = CHECKS
    counts = read(LEGACY/"analysis_panel.csv")
    s1 = "figure_s1_predictive_checks"
    for tag, stream, variable in [("A", "human", "total"), ("B", "human", "zeros"),
                                  ("D", "culture", "total"), ("E", "culture", "zeros")]:
        histogram = panel(s1, tag)
        observed_values = counts.positive if stream=="human" else counts.leptospira_culture_positive
        observed = observed_values.sum() if variable=="total" else (observed_values==0).sum()
        numerical(histogram.observed, np.repeat(observed, len(histogram)), f"Figure S1{tag} observed")
        for family, model in MODEL_NAMES.items():
            values = read(LEGACY/f"{stream}_{family}_predictive_totals.csv")[variable].to_numpy()
            h = histogram.loc[histogram.model==model].sort_values("lower")
            assert len(values)==6000
            numerical(h.lower.to_numpy()[1:], h.upper.to_numpy()[:-1], "Predictive histogram contiguous bins")
            edges = np.r_[h.lower.to_numpy()[0], h.upper.to_numpy()]
            # R hist uses (a,b] intervals, with the first left endpoint included.
            bin_index = np.searchsorted(edges, values, side="left")-1
            bin_index[values==edges[0]] = 0
            assert ((bin_index>=0)&(bin_index<len(h))).all()
            frequency = np.bincount(bin_index, minlength=len(h))
            numerical(h.draws, frequency, f"Figure S1{tag} {model} histogram counts")
            numerical(h.probability, frequency/len(values), f"Figure S1{tag} probabilities")
            numerical(h.density, frequency/len(values)/np.diff(edges), f"Figure S1{tag} densities")
    cells = read(LEGACY/"posterior_predictive_cells.csv")
    for tag, stream in [("C", "human"), ("F", "culture")]:
        expected = cells.query('stream==@stream and model=="betabinomial"').sort_values(["observed", "year", "site"]).copy()
        expected["ordered_cell"] = np.arange(1, len(expected)+1)
        expected["covered"] = (expected.observed>=expected.low)&(expected.observed<=expected.high)
        frame(panel(s1, tag), expected, ["site", "year"],
              ["observed", "low", "high", "ordered_cell", "covered"], f"Figure S1{tag}")
        truth = counts[["site", "year", "positive" if stream=="human" else "leptospira_culture_positive"]].copy()
        truth = truth.rename(columns={truth.columns[-1]: "observed"})
        frame(expected, truth, ["site", "year"], ["observed"], f"Figure S1{tag} original counts")
    SECTIONS["figure_s1_predictive_distributions"] = CHECKS-start

    start = CHECKS
    s2 = "figure_s2_period_prior_sensitivity"
    windows = []
    for name in ["primary87", "balanced80", "early63", "exclude2015", "recent24", "arithmetic_effective"]:
        d = read(LEGACY/f"{name}_loyo_pooled.csv").set_index("rule")
        windows.append({"window": name, "positive_all": d.loc["Culture", "positive_all"],
                        "culture_coverage": d.loc["Culture", "coverage"],
                        "capture_coverage": d.loc["Capture", "coverage"],
                        "difference_pp": 100*(d.loc["Culture", "coverage"]-d.loc["Capture", "coverage"])})
    frame(panel(s2, "A"), pd.DataFrame(windows), ["window"],
          ["positive_all", "culture_coverage", "capture_coverage", "difference_pp"], "Figure S2A")
    sensitivity = read(LEGACY/"posterior_decision_sensitivity.csv")
    prior_fits = ["human_betabinomial", "human_sensitivity_site035", "human_sensitivity_site140",
                  "human_sensitivity_year100", "human_sensitivity_phi_mean10", "human_sensitivity_phi_mean40"]
    prior = sensitivity.loc[sensitivity.fit.isin(prior_fits)].copy()
    prior[["mean", "low", "high"]] *= 100
    frame(panel(s2, "B"), prior, ["fit"], ["mean", "low", "high"], "Figure S2B")
    delta = read(NETWORK/"posterior_paired_draws.csv")
    delta["difference_pp"] = 100*delta.delta_coverage
    frame(panel(s2, "C"), delta, ["family", "draw"], ["delta_coverage", "difference_pp"], "Figure S2C")
    decisions = read(LEGACY/"posterior_decisions.csv")
    fits = ["human_betabinomial", "human_balanced80", "human_early63", "human_binomial"]
    decisions = decisions.loc[decisions.fit.isin(fits)].copy()
    decisions["rule"] = decisions.strategy.replace({"Uniform random": "Random"})
    expected = decisions.loc[decisions.rule.isin(["Culture", "Capture"])].copy()
    expected[["mean", "low", "high"]] *= 100
    frame(panel(s2, "D"), expected, ["fit", "rule"], ["mean", "low", "high"], "Figure S2D")
    regret = decisions.loc[decisions.rule.isin(["Culture", "Capture", "Random"])].copy()
    regret["regret_pp"] = 100*regret.expected_regret
    frame(panel(s2, "E"), regret, ["fit", "rule"], ["expected_regret", "regret_pp"], "Figure S2E")
    threshold = panel(s2, "F")
    for family, model in MODEL_NAMES.items():
        values = delta.loc[delta.family==family, "difference_pp"].to_numpy()
        sub = threshold.loc[threshold.model==model].sort_values("threshold_pp")
        numerical(sub.threshold_pp, np.arange(0, 30.001, .25), f"Figure S2F {model} threshold grid")
        probabilities = [(values>t).mean() for t in sub.threshold_pp]
        numerical(sub.probability, probabilities, f"Figure S2F {model} posterior probabilities")
        assert (np.diff(sub.probability)<=1e-12).all()
    SECTIONS["figure_s2_sensitivity_and_thresholds"] = CHECKS-start

    start = CHECKS
    s3 = "figure_s3_network_sensitivity"
    for tags, design in [("AB", "pooled_fixed87"), ("CD", "forward_expanding"), ("EF", "forward_rolling3")]:
        for tag in tags:
            frame(panel(s3, tag), pooled.query('design==@design'), ["rule", "k"], METRICS, f"Figure S3{tag}")
    SECTIONS["figure_s3_evaluation_designs"] = CHECKS-start


def posterior_rate_draws(stream, family):
    metadata_path = LEGACY/f"{stream}_{family}_fit.json"
    nc = LEGACY/f"{stream}_{family}.nc"
    USED.update([metadata_path, nc])
    meta = json.loads(metadata_path.read_text())
    sites = meta["sites"]
    years = meta["years"]
    ns, ny = len(sites), len(years)
    site_basis = np.linalg.qr(np.vstack([np.eye(ns-1), -np.ones((1, ns-1))]))[0][:, :ns-1]
    year_basis = np.linalg.qr(np.vstack([np.eye(ny-1), -np.ones((1, ny-1))]))[0][:, :ny-1]
    with xr.open_dataset(nc, group="posterior") as post:
        site_effect = post.site_theta.values.reshape(-1, ns-1)@site_basis.T
        year_effect = post.year_theta.values.reshape(-1, ny-1)@year_basis.T
        intercept = post.intercept.values.reshape(-1, 1, 1)
        rates = expit(intercept+site_effect[:, :, None]+year_effect[:, None, :]).mean(axis=2)
    assert rates.shape == (6000, 8)
    return sites, rates


def verify_rank_panels(source):
    start = CHECKS
    s4 = "figure4_model_uncertainty"
    panel = lambda tag: read(source/f"{s4}_{tag}.csv")
    rate_rows, rank_rows, inclusion_rows, mass_rows, agreements = [], [], [], [], []
    for family, model in MODEL_NAMES.items():
        draws, orders = {}, {}
        for stream in ["human", "culture"]:
            sites, rates = posterior_rate_draws(stream, family)
            draws[stream] = rates
            ranks = np.argsort(np.argsort(-rates, axis=1), axis=1)+1
            orders[stream] = ranks
            for j, site in enumerate(sites):
                rate_rows.append({"model": model, "stream": stream, "site": site,
                                  "mean": rates[:, j].mean(), "median": np.median(rates[:, j]),
                                  "low": np.quantile(rates[:, j], .025), "high": np.quantile(rates[:, j], .975)})
                for rank in range(1, 9):
                    rank_rows.append({"model": model, "stream": stream.title(), "site": site,
                                      "rank": rank, "probability": (ranks[:, j]==rank).mean()})
                inclusion_rows.append({"model": model, "stream": stream.title(), "site": site,
                                       "probability": (ranks[:, j]<=3).mean()})
        for j, site in enumerate(sites):
            inclusion_rows.append({"model": model, "stream": "Joint", "site": site,
                                   "probability": ((orders["human"][:, j]<=3)&(orders["culture"][:, j]<=3)).mean()})
        rho = ((orders["human"]-4.5)*(orders["culture"]-4.5)).sum(axis=1)/42
        values, frequencies = np.unique(np.round(rho, 12), return_counts=True)
        mass_rows.extend({"model": model, "rho": value, "draws": frequency, "probability": frequency/6000}
                         for value, frequency in zip(values, frequencies))
        same = ((orders["human"]<=3)==(orders["culture"]<=3)).all(axis=1)
        agreements.append({"model": model, "same_top3": same.mean(), "probability_positive": (rho>0).mean(),
                           "rho_mean": rho.mean(), "rho_median": np.median(rho),
                           "rho_low": np.quantile(rho, .025), "rho_high": np.quantile(rho, .975)})
    rate_data = pd.DataFrame(rate_rows)
    for tag, stream in [("A", "human"), ("B", "culture")]:
        frame(panel(tag), rate_data.query('stream==@stream'), ["model", "site"], ["median", "low", "high"], f"Figure S4{tag}")
    ranks = panel("C")
    frame(ranks, pd.DataFrame(rank_rows), ["model", "stream", "site", "rank"], ["probability"], "Figure S4C posterior ranks")
    numerical(ranks.groupby(["model", "stream", "site"]).probability.sum(), np.ones(32), "Rank-probability normalisation")
    frame(panel("D"), pd.DataFrame(inclusion_rows), ["model", "stream", "site"], ["probability"], "Figure S4D inclusion")
    frame(panel("E"), pd.DataFrame(mass_rows), ["model", "rho"], ["draws", "probability"], "Figure S4E discrete concordance")
    frame(panel("F"), pd.DataFrame(agreements), ["model"],
          ["same_top3", "probability_positive", "rho_mean", "rho_median", "rho_low", "rho_high"], "Figure S4F complete top-three agreement")
    SECTIONS["figure_s4_posterior_reconstruction"] = CHECKS-start


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args()
    source = args.source_dir.resolve()
    manifest = verify_manifest(source)
    sites = read(NETWORK/"site_summary.csv")
    pooled = read(NETWORK/"pooled_metrics.csv")
    folds = read(NETWORK/"fold_metrics.csv")
    weights = read(NETWORK/"selection_weights.csv")
    verify_main_figures(source, sites, pooled, folds, weights)
    verify_support(source, pooled)
    verify_rank_panels(source)
    report = {
        "status": "passed", "engine": "ggplot2", "composition": "patchwork",
        "figures_verified": len(manifest["figures"]), "outputs_verified": len(manifest["outputs_sha256"]),
        "panels_verified": len(manifest["panels"]),
        "minimum_panel_count": min(spec["panel_count"] for spec in manifest["figures"].values()),
        "rendered_pdf_tags": "all registered panel letters found in rendered vector PDFs",
        "numerical_checks": CHECKS, "numerical_checks_by_section": SECTIONS,
        "posterior_check": "reconstructed ranks, inclusion and complete top-three agreement from four saved NetCDF fits",
        "inputs_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in sorted(USED)},
        "verification_script_sha256": sha256(Path(__file__)),
    }
    (source/"verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k: v for k, v in report.items() if k!="inputs_sha256"}, indent=2))


if __name__ == "__main__":
    main()
