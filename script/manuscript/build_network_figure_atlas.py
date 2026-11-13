"""Combine the ggplot2 figure PDFs at native size, with a source-hash manifest."""
from pathlib import Path
import hashlib
import json

import pymupdf


ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "manuscript/figures"
SOURCE = FIGURES / "source_data/network_ggplot_2026-11-13"
ORDER = [
    ("figure1_indicator_priorities", "Figure 1 — Indicator choice and spatial priorities"),
    ("figure2_network_coverage", "Figure 2 — Network expansion and coverage gains"),
    ("figure3_forward_evaluation", "Figure 3 — Historical information, evaluation and updating"),
    ("figure_s1_predictive_checks", "Figure S1 — Sparse-count predictive assessment"),
    ("figure_s2_period_prior_sensitivity", "Figure S2 — Coverage sensitivity to analytical choices"),
    ("figure_s3_network_sensitivity", "Figure S3 — Alternative uses of surveillance history"),
    ("figure4_model_uncertainty", "Figure S4 — From site uncertainty to portfolio agreement"),
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    source_manifest = json.loads((SOURCE / "manifest.json").read_text())
    atlas = pymupdf.open()
    bookmarks, inputs = [], {}
    for stem, label in ORDER:
        path = FIGURES / f"{stem}.pdf"
        inputs[path.name] = digest(path)
        assert inputs[path.name] == source_manifest["outputs_sha256"][path.name]
        with pymupdf.open(path) as figure:
            assert len(figure) == 1
            bookmarks.append([1, label, len(atlas) + 1])
            atlas.insert_pdf(figure)
    atlas.set_toc(bookmarks)
    atlas.set_metadata({
        "title": "Sentinel-network figures: ggplot2 six-panel suite",
        "subject": "Three main and four supplementary figures; 42 panels",
        "creator": "ggplot2 / patchwork; atlas assembled with PyMuPDF",
    })
    output = FIGURES / "network_figure_atlas_2026-11-13.pdf"
    atlas.save(output, garbage=4, deflate=True)
    manifest = {
        "inputs_sha256": inputs, "atlas": output.name,
        "atlas_sha256": digest(output), "pages": len(atlas),
        "builder_sha256": digest(Path(__file__)),
    }
    (SOURCE / "atlas_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
