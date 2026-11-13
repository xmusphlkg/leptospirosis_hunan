# Network-size coverage and historical-window sensitivity.
# Called after network_plot_common.R by the main figure-suite driver.

build_figures_2_s3 <- function() {
  metrics <- read_network("pooled_metrics.csv")
  periods <- read_network("period_metrics.csv")
  stem <- "figure2_network_coverage"
  main <- metrics |>
    filter(design == "loyo87", rule %in% c("Culture", "Capture", "Random")) |>
    mutate(rule = factor(rule, c("Culture", "Capture", "Random")))
  stopifnot(nrow(main) == 21, all(table(main$rule) == 7))

  curve <- function(data, value, ylabel, tag, keep_legend = FALSE, strip = NULL) {
    legend_labels <- rule_labels
    if (!is.null(strip)) legend_labels[c("PastHuman", "HumanPooled")] <- "Human-serology reference"
    p <- ggplot(data, aes(k, .data[[value]], colour = rule, group = rule)) +
      annotate("rect", xmin = 2.80, xmax = 3.20, ymin = -Inf, ymax = Inf,
               fill = "#D9D9D9", alpha = .27) +
      geom_line(aes(linetype = rule), linewidth = .85) +
      geom_point(data = data |> filter(rule != "Random"), aes(shape = rule),
                 size = 2.6, stroke = .65) +
      scale_colour_manual(values = rule_colors, labels = legend_labels, name = NULL,
                          guide = if (keep_legend) "legend" else "none") +
      scale_shape_manual(values = rule_shapes, guide = "none") +
      scale_linetype_manual(values = c(Culture = "solid", Capture = "solid",
                                      PastHuman = "solid", HumanPooled = "solid",
                                      Random = "22"), guide = "none") +
      scale_x_continuous(breaks = 1:7, limits = c(.75, 7.25)) +
      labs(x = "Number of sites selected", y = ylabel)
    if (value %in% c("coverage", "testing_share")) {
      p <- p + scale_y_continuous(limits = c(0, 1), labels = percent)
    } else {
      upper <- max(2.5, ceiling(max(data[[value]], na.rm = TRUE) * 2) / 2)
      p <- p + scale_y_continuous(limits = c(0, upper)) +
        geom_hline(yintercept = 1, colour = "#969696", linewidth = .4, linetype = "dotted")
    }
    if (!is.null(strip)) {
      data$comparison <- strip
      # The single strip identifies the evaluation design without a plot title.
      p$data$comparison <- strip
      for (i in seq_along(p$layers)) {
        if (is.data.frame(p$layers[[i]]$data)) p$layers[[i]]$data$comparison <- strip
      }
      p <- p + facet_wrap(~comparison, nrow = 1)
    }
    if (keep_legend) {
      rules <- levels(droplevels(data$rule))
      legend_shapes <- unname(rule_shapes[rules])
      legend_shapes[rules == "Random"] <- NA
      p <- p + guides(colour = guide_legend(nrow = 1,
        override.aes = list(linewidth = .9, shape = legend_shapes,
                            linetype = ifelse(rules == "Random", "22", "solid"))))
    }
    panel(p, tag)
  }

  panel_data(main, stem, "A", "Held-out human-positive coverage as the LOYO network grows")
  a <- curve(main, "coverage", "Held-out human-positive coverage", "A", TRUE)
  panel_data(main, stem, "B", "Human testing share represented by the same LOYO portfolios")
  b <- curve(main, "testing_share", "Human testing share", "B")
  panel_data(main, stem, "C", "Coverage concentration relative to testing share across LOYO network sizes")
  c <- curve(main, "enrichment", "Enrichment (coverage / testing share)", "C")

  plane <- main |> mutate(reference_enrichment = enrichment)
  panel_data(plane, stem, "D", "Coverage-testing-share trade-off for every evaluated network size")
  reference <- tibble(enrichment = c(1, 1.5, 2)) |>
    mutate(xend = pmin(1, 1 / enrichment), yend = pmin(enrichment, 1))
  d <- ggplot(plane, aes(testing_share, coverage, colour = rule, group = rule)) +
    geom_segment(data = reference, aes(x = 0, y = 0, xend = xend, yend = yend),
                 inherit.aes = FALSE, colour = "#C4C4C4", linetype = "22", linewidth = .45) +
    geom_path(aes(linetype = rule), linewidth = .7) +
    geom_point(aes(size = k, shape = rule), stroke = .7) +
    scale_colour_manual(values = rule_colors, guide = "none") +
    scale_shape_manual(values = rule_shapes, guide = "none") +
    scale_linetype_manual(values = c(Culture = "solid", Capture = "solid", Random = "22"), guide = "none") +
    scale_size_continuous(name = "Sites selected", breaks = 1:7, limits = c(1, 7), range = c(1.5, 4.8)) +
    scale_x_continuous(limits = c(0, 1), labels = percent) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(x = "Human testing share", y = "Held-out human-positive coverage") +
    guides(size = guide_legend(nrow = 1, override.aes = list(colour = "#606060", shape = 16)))
  d <- panel(d, "D")

  increments <- main |>
    filter(rule != "Random") |>
    arrange(rule, k) |>
    group_by(rule) |>
    mutate(previous_coverage = lag(coverage, default = 0),
           added_coverage = coverage - previous_coverage,
           added_coverage_pp = 100 * added_coverage) |>
    ungroup()
  stopifnot(all(increments$added_coverage >= -1e-12))
  totals <- increments |> group_by(rule) |>
    summarise(total = sum(added_coverage), last_coverage = last(coverage), .groups = "drop")
  stopifnot(all(abs(totals$total - totals$last_coverage) < 1e-12))
  panel_data(increments, stem, "E", "Marginal coverage gained when one further site is included in the LOYO ranking")
  e <- ggplot(increments, aes(k, added_coverage_pp, fill = rule)) +
    geom_col(position = position_dodge(width = .72), width = .62, alpha = .95) +
    scale_fill_manual(values = rule_colors, guide = "none") +
    scale_x_continuous(breaks = 1:7, limits = c(.5, 7.5)) +
    scale_y_continuous(limits = c(0, 35)) +
    labs(x = "Number of sites after expansion", y = "Additional coverage (percentage points)")
  e <- panel(e, "E")

  period_difference <- periods |>
    filter(design == "loyo87", rule %in% c("Culture", "Capture")) |>
    select(period, k, rule, coverage, positive_all) |>
    pivot_wider(names_from = rule, values_from = coverage) |>
    mutate(difference_pp = 100 * (Culture - Capture),
           period_label = factor(paste0(gsub("-", "\u2013", period), "  (n = ", positive_all, ")"),
             levels = c("2023\u20132025  (n = 11)", "2018\u20132025  (n = 92)", "2015\u20132017  (n = 199)")),
           k_factor = factor(k, levels = 1:7))
  stopifnot(nrow(period_difference) == 21, !anyNA(period_difference$period_label))
  primary_early <- period_difference |> filter(period == "2015-2017", k == 3)
  stopifnot(abs(primary_early$difference_pp - 100 * 50 / 199) < 1e-10)
  panel_data(period_difference, stem, "F", "Period-specific coverage contrasts retaining the original LOYO selection rules")
  limit <- ceiling(max(abs(period_difference$difference_pp)) / 10) * 10
  f <- ggplot(period_difference, aes(k_factor, period_label, fill = difference_pp)) +
    geom_tile(colour = "#D5D5D5", linewidth = .35, width = .98, height = .88) +
    scale_fill_gradient2(low = rule_colors[["Capture"]], mid = "#FFFFFF", high = rule_colors[["Culture"]],
                         midpoint = 0, limits = c(-limit, limit),
                         name = "Culture \u2212 capture\n(percentage points)") +
    labs(x = "Number of sites selected", y = "Evaluation period (human positives)") +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(size = 9.5)) +
    guides(fill = guide_colourbar(title.position = "top", barwidth = grid::unit(44, "mm"),
                                  barheight = grid::unit(3, "mm")))
  f <- panel(f, "F")

  figure <- ((a | b | c) / (d | e | f)) +
    plot_layout(heights = c(1, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")
  save_network_figure(figure, stem, width = 18, height = 10.7)

  # Six panels, each with one comparison design and one coverage estimand.
  sensitivity_stem <- "figure_s3_network_sensitivity"
  designs <- c("pooled_fixed87", "forward_expanding", "forward_rolling3")
  strips <- c("Full-period fixed rankings", "Forward: cumulative history", "Forward: preceding three years")
  panels <- list()
  for (i in seq_along(designs)) {
    data <- metrics |>
      filter(design == designs[i]) |>
      mutate(rule = factor(rule, c("Culture", "Capture", "HumanPooled", "PastHuman", "Random")))
    stopifnot(nrow(data) == 28)
    human <- if (designs[i] == "pooled_fixed87") "HumanPooled" else "PastHuman"
    stopifnot(setequal(as.character(unique(data$rule)), c("Culture", "Capture", human, "Random")))
    tag1 <- LETTERS[2*i-1]
    tag2 <- LETTERS[2*i]
    panel_data(data, sensitivity_stem, tag1, paste(strips[i], "coverage across all network sizes"))
    panel_data(data, sensitivity_stem, tag2, paste(strips[i], "enrichment across all network sizes"))
    panels[[tag1]] <- curve(data, "coverage", "Human-positive coverage", tag1, i == 1, strips[i])
    panels[[tag2]] <- curve(data, "enrichment", "Enrichment (coverage / testing share)", tag2, FALSE, strips[i])
  }
  # The strip and caption distinguish in-sample and previous-history references.
  sensitivity <- ((panels$A | panels$B) / (panels$C | panels$D) / (panels$E | panels$F)) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save_network_figure(sensitivity, sensitivity_stem, width = 15.5, height = 13.2)
  invisible(list(figure2 = figure, figure_s3 = sensitivity))
}
