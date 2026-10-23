#!/usr/bin/env Rscript
# Rebuild the three manuscript figures from the current 87-row analysis, using ggplot2
# for every panel and patchwork for composition. Run from any working directory.
# Usage: Rscript script/hunan/plot_priority_figures.R [output_directory]
packages <- c("ggplot2", "patchwork", "dplyr", "tidyr", "readr",
              "scales", "jsonlite", "ncdf4", "digest")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
})
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), "../.."))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else file.path(root, "manuscript/figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "script/hunan/priority_figure_data.R"))
message("Reading recomputed 87-row counts and posterior summaries...")
d <- prepare_priority_data(root)

palette <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF",
             "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF")
ink <- "#222222"
stream_colors <- setNames(palette[c(1, 2, 4, 5)], c("Human", "Culture", "Capture", "Joint"))
rule_colors <- setNames(palette[c(2, 4, 5)], c("Culture", "Capture", "Compound"))
model_colors <- setNames(palette[c(2, 3)], c("Binomial", "Beta-binomial"))
percent <- scales::label_percent(accuracy = 1)
theme_set(theme_bw(base_size = 10, base_family = "DejaVu Sans") + theme(
  text = element_text(colour = ink),
  plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
  plot.tag = element_text(face = "bold", size = 16, colour = ink),
  plot.tag.position = "topleft", plot.title.position = "panel",
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#E8E8E8", linewidth = .3),
  panel.grid.major.y = element_blank(),
  panel.border = element_rect(colour = "#444444", fill = NA, linewidth = .5),
  axis.title = element_text(size = 9.5), axis.text = element_text(colour = ink, size = 9),
  legend.position = "bottom", legend.title = element_text(size = 9),
  legend.text = element_text(size = 8.5), legend.key.height = grid::unit(3.5, "mm"),
  legend.margin = margin(0, 0, 0, 0),
  strip.text = element_text(face = "bold", size = 9),
  strip.background = element_rect(fill = "white", colour = "#444444", linewidth = .5),
  panel.spacing = grid::unit(22, "pt"), plot.margin = margin(10, 12, 10, 10)
))
# Apply the same pretty-break contract to every quantitative axis and legend.
# Numeric encodings of site names/ranks/years remain ordinal, not continuous.
pretty_scale <- function(p) {
  built <- ggplot_build(p)
  for (axis in c("x", "y")) {
    sc <- p$scales$get_scales(axis)
    trained <- if (axis == "x") built$layout$panel_scales_x[[1]] else built$layout$panel_scales_y[[1]]
    if (is.null(sc) && !trained$is_discrete()) {
      p <- p + if (axis == "x") scale_x_continuous() else scale_y_continuous()
    }
  }
  for (i in seq_along(p$scales$scales)) {
    sc <- p$scales$scales[[i]]$clone()
    if (sc$is_discrete()) next
    ordinal <- is.character(sc$labels) ||
      (is.numeric(sc$breaks) && (identical(as.numeric(sc$breaks), as.numeric(1:8)) ||
                                all(sc$breaks >= 1900)))
    if (ordinal) next
    domain <- sc$limits
    sc$limits <- if (is.numeric(domain)) range(pretty(domain)) else function(x) range(pretty(x))
    sc$breaks <- function(x) pretty(x)
    sc$expand <- c(0, 0)
    p$scales$scales[[i]] <- sc
  }
  p
}
panel <- function(p, tag) pretty_scale(p) + labs(tag = tag, title = NULL, subtitle = NULL, caption = NULL)

# Figure 1: magnitudes, rank displacement, two associations, and a site fingerprint.
p1a <- ggplot(d$rates, aes(rate, site_label, colour = stream)) +
  geom_segment(aes(x = low, xend = high, yend = site_label), linewidth = 1) +
  geom_point(size = 2.5) + facet_wrap(~stream, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = stream_colors, guide = "none") +
  scale_x_continuous(labels = percent, expand = expansion(mult = c(.02, .12))) +
  labs(x = "Pooled rate", y = NULL)
p1a <- panel(p1a, "B")
span <- d$rates |> group_by(site_label) |> summarise(lo = min(rank), hi = max(rank), .groups = "drop")
p1b <- ggplot(d$rates, aes(rank, site_label)) +
  geom_segment(data = span, aes(x = lo, xend = hi, yend = site_label), colour = "#BDBDBD", linewidth = 1.2) +
  geom_point(aes(colour = stream, shape = stream), size = 3, alpha = .9) +
  scale_colour_manual(values = stream_colors, name = NULL) +
  scale_shape_manual(values = c(Human = 16, Culture = 17, Capture = 1), name = NULL) +
  scale_x_continuous(breaks = 1:8, limits = c(.65, 8.35)) +
  labs(x = "Rank", y = NULL)
p1b <- panel(p1b, "C")

scatter <- function(variable, tag, xlabel) {
  s <- d$site |> mutate(x = .data[[variable]],
                       human_top = rank(-human_fraction) <= 3,
                       animal_top = rank(-x) <= 3,
                       top3 = factor(case_when(human_top & animal_top ~ "Both",
                                               human_top | animal_top ~ "One", TRUE ~ "Neither"),
                                     c("Both", "One", "Neither")))
  p <- ggplot(s, aes(x, human_fraction)) +
    geom_point(aes(size = human_n, fill = top3), shape = 21, colour = "white", stroke = .6, show.legend = TRUE) +
    scale_fill_manual(values = setNames(palette[c(2, 3, 6)], c("Both", "One", "Neither")),
                      drop = FALSE, name = "Top 3 in the paired streams") +
    scale_size_continuous(range = c(3.5, 5.5), name = "Human tests", breaks = c(1000, 1100, 1200)) +
    scale_x_continuous(labels = percent, expand = expansion(mult = c(.12, .25))) +
    scale_y_continuous(labels = percent, breaks = seq(0, .08, .02), expand = expansion(mult = c(.12, .22))) +
    labs(x = xlabel, y = "Human seroreactivity") +
    guides(size = "none", fill = guide_legend(override.aes = list(size = 3)))
  panel(p, tag) + coord_cartesian(clip = "off")
}
p1c <- scatter("culture_fraction", "D", "Rodent culture positivity")
p1d <- scatter("capture_yield", "E", "Rodent capture yield")
# Quantitative rates: independent, explicitly labelled percentage scales.
rate_tile <- function(stream_name, upper, breaks, tag = NULL, show_sites = FALSE) {
  x <- d$rates |> filter(as.character(stream) == stream_name)
  p <- ggplot(x, aes(x = stream, y = site_label, fill = rate)) +
    geom_tile(colour = "white", linewidth = 1, width = .98, height = .98) +
    geom_text(aes(label = sprintf("%.2f", 100*rate),
                  colour = rate > upper*.55), size = 3) +
    scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "white"), guide = "none") +
    scale_fill_gradient(low = "white", high = stream_colors[[stream_name]],
                        limits = c(0, upper), breaks = breaks,
                        labels = function(x) sprintf("%g", 100*x),
                        name = paste0(stream_name, " (%)")) +
    guides(fill = guide_colourbar(barwidth = grid::unit(23, "mm"),
                                  barheight = grid::unit(2.5, "mm"),
                                  title.position = "top")) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size = 7), legend.title = element_text(size = 8),
          plot.margin = margin(10, 2, 10, 2))
  if (!show_sites) p <- p + theme(axis.text.y = element_blank())
  if (!is.null(tag)) p <- panel(p, tag) else p <- pretty_scale(p)
  p
}
p1e <- (rate_tile("Human", .08, c(0,.04,.08), "F", TRUE) |
        rate_tile("Culture", .16, c(0,.08,.16)) |
        rate_tile("Capture", .08, c(0,.04,.08))) + plot_layout(widths = c(1,1,1))
# Real satellite imagery with county-level reference locations (not sampling GPS).
map_dir <- file.path(root, "data/figure1_map_cache")
if (!file.exists(file.path(map_dir, "satellite.png")))
  stop("Run python script/hunan/prepare_satellite_map.py first")
boundary <- readr::read_csv(file.path(map_dir, "boundary.csv"), show_col_types = FALSE)
locations <- readr::read_csv(file.path(map_dir, "sites.csv"), show_col_types = FALSE)
labels <- tibble::tribble(
  ~site, ~label_lon, ~label_lat,
  "Junshan", 112.0, 29.7, "Yuanjiang", 111.7, 29.05,
  "Ningxiang", 111.5, 28.4, "Liuyang", 113.65, 28.6,
  "Xiangtan", 112.8, 27.05, "Liling", 113.8, 27.3,
  "Shuangfeng", 111.7, 26.85, "Lianyuan", 110.8, 27.65)
locations <- left_join(locations, labels, by = "site")
rgb_connection <- gzfile(file.path(map_dir, "satellite.rgb.gz"), "rb")
rgb_values <- matrix(as.integer(readBin(rgb_connection, "raw", n = 1800*1800*3)), nrow = 3)
close(rgb_connection)
satellite <- matrix(grDevices::rgb(rgb_values[1,], rgb_values[2,], rgb_values[3,], maxColorValue = 255),
                    nrow = 1800, byrow = TRUE)
p1map <- ggplot() +
  annotation_custom(grid::rasterGrob(as.raster(satellite),
                                    width = grid::unit(1, "npc"), height = grid::unit(1, "npc")),
                    xmin = 108, xmax = 115, ymin = 24, ymax = 31) +
  geom_path(data = boundary, aes(lon, lat, group = group), colour = "white", linewidth = .5) +
  geom_segment(data = locations, aes(lon, lat, xend = label_lon, yend = label_lat),
               colour = "white", linewidth = .35) +
  geom_point(data = locations, aes(lon, lat), shape = 21, fill = palette[1],
             colour = "white", size = 2.6, stroke = .6) +
  geom_label(data = locations, aes(label_lon, label_lat, label = site),
             size = 2.5, fill = scales::alpha("white", .9), label.padding = grid::unit(.12, "lines")) +
  annotate("text", x = 109, y = 29.65, label = "N ↑", colour = "white", size = 3.5) +
  coord_fixed(ratio = 1/cos(27.3*pi/180), xlim = range(pretty(c(108.5,114.5))), ylim = range(pretty(c(24.3,30.3))), expand = FALSE) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)") +
  annotate("text", x = 114.8, y = 24.15, hjust = 1, vjust = 0, size = 1.6, colour = "white",
           label = "Sources: Esri, Vantor, Earthstar Geographics, and the GIS User Community") +
  theme(panel.grid = element_blank())
p1map <- panel(p1map, "A")
# Wrap each row to keep its column sizing independent of the other row.
fig1_top <- (p1map | p1a | p1b) + plot_layout(widths = c(.95, 1.65, 1))
fig1_bottom <- (p1c | p1d | p1e) + plot_layout(widths = c(1, 1, 1.6))
fig1 <- wrap_elements(full = fig1_top) / wrap_elements(full = fig1_bottom) +
  plot_layout(heights = c(1, 1))

# Figure 2: fixed animal-only selection evaluated on held-out human counts.
p2a <- ggplot(d$membership, aes(factor(held_out_year), site_label, fill = status)) +
  geom_tile(colour = "white", linewidth = 1, width = .96, height = .96) +
  geom_point(data = filter(d$membership, selected), shape = 16, colour = "white", size = 1.7) +
  facet_wrap(~rule, nrow = 1) +
  scale_fill_manual(values = c(`Not selected` = "white", Selected = palette[2], `No monitoring` = "#BDBDBD"), name = NULL, drop = FALSE) +
  labs(x = "Held-out year", y = NULL) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=60,hjust=1))
p2a <- panel(p2a, "A")
annual <- d$folds |> filter(rule != "Compound")
p2b <- ggplot(annual, aes(y = factor(held_out_year))) +
  geom_segment(aes(x = testing_share_S, xend = positive_coverage_C, yend = factor(held_out_year)),
               colour = "#BDBDBD", linewidth = 1.3) +
  geom_point(aes(x = testing_share_S, shape = "Testing share S"), colour = "#777777", size = 2.7) +
  geom_point(aes(x = positive_coverage_C, shape = "Positive coverage C", colour = rule), size = 2.7) +
  facet_wrap(~rule, nrow = 1) +
  scale_colour_manual(values = rule_colors, guide = "none") +
  scale_shape_manual(values = c(`Testing share S` = 1, `Positive coverage C` = 16), name = NULL) +
  scale_x_continuous(labels = percent, limits = c(0, 1), breaks = c(0, .5, 1)) +
  labs(x = "Share of held-out human counts", y = NULL)
p2b <- panel(p2b, "B") + coord_cartesian(clip = "off")
year_labels <- d$folds |> distinct(held_out_year, positive_all_sites) |> arrange(held_out_year)
p2c <- ggplot(annual, aes(held_out_year, detection_enrichment_E_C_over_S, colour = rule)) +
  annotate("rect", xmin = 2020.6, xmax = 2025.4, ymin = -Inf, ymax = Inf, fill = scales::alpha(palette[6], .12)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#777777", linewidth = .45) +
  geom_line(linewidth = .8) + geom_point(size = 2.5) +
  scale_colour_manual(values = rule_colors, name = NULL) +
  scale_x_continuous(breaks = year_labels$held_out_year,
                     limits = c(2014.6, 2025.4)) +
  scale_y_continuous(limits = c(0, 3), breaks = 0:3) +
  labs(x = "Held-out year", y = "Detection enrichment E = C / S")
p2c <- panel(p2c, "C") + theme(axis.text.x=element_text(angle=60,hjust=1))

# Concentric symbols preserve identical culture/compound coordinates.
isolines <- tibble(E = factor(c("1", "1.5", "2"), c("1", "1.5", "2")), slope = c(1, 1.5, 2))
p2d <- ggplot(d$pooled |> arrange(desc(rule)), aes(S, C)) +
  geom_abline(data = isolines, aes(slope = slope, intercept = 0, linetype = E),
              colour = "#999999", linewidth = .5) +
  geom_point(aes(colour = rule, shape = rule, size = rule), stroke = 1) +
  scale_colour_manual(values = rule_colors, name = NULL) +
  scale_shape_manual(values = c(Culture = 16, Capture = 17, Compound = 1), name = NULL) +
  scale_size_manual(values = c(Culture = 3.3, Capture = 3.3, Compound = 6), guide = "none") +
  scale_linetype_manual(values = c(`1` = "dashed", `1.5` = "dotted", `2` = "dotdash"), name = "E") +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 3)),
         shape = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  scale_x_continuous(labels = percent, limits = c(.2, .6), breaks = seq(.2, .6, .1)) +
  scale_y_continuous(labels = percent, limits = c(.25, .85), breaks = seq(.3, .8, .1)) +
  labs(x = "Pooled testing share S", y = "Pooled positive coverage C") +
  theme(legend.box = "vertical")
p2d <- panel(p2d, "D")
ps <- d$perm_summary
# Split bin counts by exact tail membership, including bins crossing a cutoff.
# The P value always uses individual permutations, never binned counts.
hh <- hist(d$permutation$delta, breaks = 48, plot = FALSE)
hist_data <- d$permutation |>
  mutate(bin = findInterval(delta, hh$breaks, rightmost.closed = TRUE)) |>
  count(bin, tail, name = "count") |> mutate(x = hh$mids[bin])
p2e <- ggplot(hist_data, aes(x, count)) +
  annotate("rect", xmin = ps$interval[1], xmax = ps$interval[2], ymin = -Inf, ymax = Inf,
           fill = palette[6], alpha = .12) +
  geom_col(aes(fill = tail), width = diff(hh$breaks)[1]*.94) +
  geom_vline(xintercept = ps$median, colour = "#777777", linetype = "dotted", linewidth = .6) +
  geom_vline(xintercept = c(ps$median-abs(ps$observed-ps$median), ps$observed),
             colour = palette[1], linetype = c("dashed", "solid"), linewidth = .7) +
  scale_fill_manual(values = c(`FALSE` = palette[6], `TRUE` = palette[1]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .07))) +
  labs(x = "ΔE = Culture − Capture", y = "Permutations")
p2e <- panel(p2e, "E")
fig2 <- ((p2a | p2b) + plot_layout(widths = c(1.7, 1))) /
  ((p2c | p2d | p2e) + plot_layout(widths = c(1, 1, 1.05))) +
  plot_layout(heights = c(1.05, 1))

# Figure 4: uncertainty of rates, full rankings, individual and set-level decisions.
interval_panel <- function(stream, tag) {
  x <- d$intervals |> filter(.data$stream == .env$stream) |>
    mutate(y = as.numeric(site_label) + if_else(model == "Binomial", .13, -.13))
  p <- ggplot(x, aes(median, y, colour = model)) +
    geom_segment(aes(x = low, xend = high, yend = y), linewidth = .85) + geom_point(size = 2.3) +
    scale_colour_manual(values = model_colors, name = NULL) +
    scale_y_continuous(breaks = 1:8, labels = rev(d$display_order), limits = c(.55, 8.45)) +
    scale_x_continuous(labels = scales::label_percent(accuracy = .1), expand = expansion(mult = c(.02, .08))) +
    labs(x = "Equal-year standardised positivity", y = NULL)
  panel(p, tag)
}
p3a <- interval_panel("human", "A")
p3b <- interval_panel("culture", "B")
rank_data <- d$ranks |> mutate(stream = factor(stream, c("human", "culture"), c("Human", "Culture")))
p3c <- ggplot(rank_data, aes(rank, site_label, fill = probability)) +
  geom_tile(colour = "white", linewidth = .35) +
  facet_grid(model ~ stream) +
  scale_x_continuous(breaks = 1:8, expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = palette[2], limits = c(0, 1),
                      labels = percent, name = "Rank probability") +
  labs(x = "Posterior rank", y = NULL) + theme(panel.grid = element_blank())
p3c <- panel(p3c, "C") + guides(fill = guide_colourbar(barwidth = grid::unit(50, "mm")))
inc <- d$inclusion |> mutate(y = as.numeric(site_label) +
                              c(Human = .2, Culture = 0, Joint = -.2)[as.character(stream)])
p3d <- ggplot(inc, aes(probability, y, colour = stream)) +
  geom_point(aes(shape = stream), size = 2) + facet_wrap(~model, nrow = 1) +
  scale_colour_manual(values = stream_colors, name = NULL) +
  scale_shape_manual(values = c(Human = 16, Culture = 17, Joint = 18), name = NULL) +
  scale_y_continuous(breaks = 1:8, labels = rev(d$display_order), limits = c(.55, 8.45)) +
  scale_x_continuous(labels = percent, breaks = c(0, .5, 1), limits = c(0, 1)) +
  labs(x = "Pr(rank ≤ 3)", y = NULL)
p3d <- panel(p3d, "D") + coord_cartesian(clip = "off")

# With eight sites, Spearman rho is discrete. Plot empirical probability mass
# at its attainable values instead of inventing a smooth posterior density.
rho_mass <- d$rho |> mutate(rho = round(rho*42)/42) |>
  count(model, rho, name = "n") |> group_by(model) |> mutate(probability = n/sum(n)) |> ungroup()
rho_top <- max(rho_mass$probability)*1.12
p3e <- ggplot(rho_mass, aes(rho, probability, fill = model)) +
  annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf, fill = scales::alpha(palette[6], .12)) +
  geom_col(width = 1/42*.85) + facet_wrap(~model, nrow = 1) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = .4) +
  geom_segment(data = d$models, aes(x = low, xend = high, y = rho_top, yend = rho_top, colour = model),
               inherit.aes = FALSE, linewidth = 1) +
  geom_point(data = d$models, aes(x = median, y = rho_top, colour = model), inherit.aes = FALSE, size = 2.4) +
  scale_fill_manual(values = model_colors, guide = "none") +
  scale_colour_manual(values = model_colors, guide = "none") +
  scale_x_continuous(breaks = seq(-1, 1, .5), limits = c(-1, 1), oob = scales::squish) +
  scale_y_continuous(labels = percent, limits = c(0, rho_top*1.08)) +
  labs(x = "Posterior Spearman ρ", y = "Probability mass")
p3e <- panel(p3e, "E")
# The two complete-set overlap probabilities are reported in the Results;
# a separate two-point panel adds no information to the compact figure.
fig3_top <- (p3a | p3b | p3c) + plot_layout(widths = c(1, 1, 1.4))
fig3_bottom <- (p3d | p3e) + plot_layout(widths = c(1.3, 1))
fig3 <- wrap_elements(full = fig3_top) / wrap_elements(full = fig3_bottom) +
  plot_layout(heights = c(1.25, .85))

save_figure <- function(p, stem, width, height) {
  message("Rendering ", stem, "...")
  ggsave(file.path(out, paste0(stem, ".png")), p, width = width, height = height,
         dpi = 320, bg = "white", type = "cairo")
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width, height = height,
         bg = "white", device = grDevices::cairo_pdf)
}
save_figure(fig1, "figure1_site_signals", 18, 10)
save_figure(fig2, "figure2_historical_coverage", 15, 10)
save_figure(fig3, "figure4_model_uncertainty", 18, 10)

# Recomputed decision and predictive-check figures.
current <- file.path(root,'results/priority_87_2026-10-22')
read_current <- function(n) readr::read_csv(file.path(current,n),show_col_types=FALSE)
contributions <- read_current('primary87_contributions.csv') |> group_by(site) |>
  summarise(contribution=sum(contribution),.groups='drop') |> mutate(site=factor(site,levels=rev(d$display_order)))
f3a <- ggplot(contributions,aes(contribution,site)) + geom_vline(xintercept=0,colour='#777777') +
  geom_col(fill=palette[2],width=.65) + labs(x='Net human positives covered: Culture − Capture',y=NULL)
front <- read_current('pooled_rule_frontier.csv')
portfolios <- read_current('all_portfolios.csv') |> group_by(k) |> summarise(best=max(coverage),.groups='drop')
f3b <- ggplot(front,aes(k,coverage,colour=rule)) +
  geom_ribbon(aes(ymin=minimum,ymax=maximum,fill=rule),alpha=.12,colour=NA) +
  geom_line(linewidth=.8) + geom_point(size=2) +
  geom_line(data=portfolios,aes(k,best),inherit.aes=FALSE,colour=palette[5],linetype='dotted') +
  geom_line(data=portfolios,aes(k,k/8),inherit.aes=FALSE,colour=palette[6],linetype='dashed') +
  scale_colour_manual(values=rule_colors,name=NULL) + scale_fill_manual(values=rule_colors,guide='none') +
  scale_x_continuous(breaks=1:7) + scale_y_continuous(labels=percent,limits=c(0,1)) +
  labs(x='Portfolio size k',y='Observed positive coverage')
decisions <- read_current('posterior_decisions.csv') |> filter(fit=='human_betabinomial') |>
  mutate(strategy=factor(strategy,levels=c('Uniform random','Capture','Culture','Posterior optimum')))
f3c <- ggplot(decisions,aes(expected_regret,strategy,fill=strategy)) + geom_col(width=.65) +
  scale_fill_manual(values=c('Culture'=palette[2],'Capture'=palette[4],'Uniform random'=palette[6],'Posterior optimum'=palette[5]),guide='none') +
  scale_x_continuous(labels=percent) + labs(x='Posterior expected regret',y=NULL)
fig_decision <- (panel(f3a,'A') | panel(f3b,'B') | panel(f3c,'C')) + plot_layout(widths=c(1.2,1,1.2))
save_figure(fig_decision,'figure3_decision_analysis',18,5.6)
replicated <- bind_rows(lapply(c('human','culture'),function(s) bind_rows(lapply(c('binomial','betabinomial'),function(m)
  read_current(paste0(s,'_',m,'_predictive_totals.csv')) |> mutate(stream=s,model=m))))) |>
  pivot_longer(c(total,zeros),names_to='statistic',values_to='value') |>
  mutate(model=factor(model,c('binomial','betabinomial'),c('Binomial','Beta-binomial')),
         stream=factor(stream,c('human','culture'),c('Human','Culture')),
         statistic=factor(statistic,c('total','zeros'),c('Total positives','Zero-positive site–years')))
observed <- read_current('posterior_predictive_checks.csv') |> filter(fit %in% c('human_betabinomial','culture_betabinomial')) |>
  select(stream,total=observed_total,zeros=observed_zeros) |> pivot_longer(c(total,zeros),names_to='statistic',values_to='value') |>
  mutate(stream=factor(stream,c('human','culture'),c('Human','Culture')),
         statistic=factor(statistic,c('total','zeros'),c('Total positives','Zero-positive site–years')))
fs1 <- ggplot(replicated,aes(value,colour=model)) + geom_freqpoly(aes(y=after_stat(density)),bins=35,linewidth=.75) +
  geom_vline(data=observed,aes(xintercept=value),colour='#222222',linetype='dashed') +
  facet_wrap(stream~statistic,scales='free',ncol=2) + scale_colour_manual(values=model_colors,name=NULL) +
  labs(x='Replicated count',y='Density')
save_figure(pretty_scale(fs1),'figure_s1_predictive_checks',12,8)
era <- bind_rows(read_current('early63_site_summary.csv') |> mutate(era='2015–2022'),
                 read_current('recent24_site_summary.csv') |> mutate(era='2023–2025'))
fs2a <- ggplot(era,aes(culture_fraction,human_fraction,colour=site)) + geom_point(size=2.5) +
  facet_wrap(~era,scales='free') + scale_colour_manual(values=setNames(palette,sort(unique(era$site))),name=NULL) +
  scale_x_continuous(labels=percent) + scale_y_continuous(labels=scales::label_percent(accuracy=0.1)) + labs(x='Recorded culture yield',y='Human seroreactivity')
sensitivity <- read_current('posterior_decision_sensitivity.csv') |>
  mutate(label=recode(fit,human_betabinomial='Primary beta-binomial',human_binomial='Binomial',
    human_sensitivity_site035='Site SD 0.35',human_sensitivity_site140='Site SD 1.40',human_sensitivity_year100='Year SD 1.00',
    human_sensitivity_phi_mean10='Concentration prior mean 10',human_sensitivity_phi_mean40='Concentration prior mean 40',
    human_balanced80='Balanced 80 rows',human_early63='2015–2022: 63 rows'))
fs2b <- ggplot(sensitivity,aes(mean,reorder(label,mean))) + geom_vline(xintercept=0,colour='#777777',linetype='dashed') +
  geom_segment(aes(x=low,xend=high,yend=reorder(label,mean)),colour=palette[3],linewidth=.8) +
  geom_point(colour=palette[3],size=2.2) + scale_x_continuous(labels=percent) +
  labs(x='Posterior coverage utility: Culture − Capture',y=NULL)
fs2 <- (panel(fs2a,'A') / panel(fs2b,'B')) + plot_layout(heights=c(1,1.2))
save_figure(fs2,'figure_s2_period_prior_sensitivity',12,10)

source_dir <- file.path(out, "source_data")
dir.create(source_dir, showWarnings = FALSE)
tables <- c("rates", "correlations", "panel", "folds", "membership", "pooled",
            "permutation", "intervals", "ranks", "inclusion", "rho", "models")
for (name in tables) readr::write_csv(d[[name]], file.path(source_dir, paste0(name, ".csv")))
jsonlite::write_json(list(
  snapshot = "priority_87_2026-10-22", verified_inputs = d$sources,
  style = list(base_theme = "theme_bw", palette = palette, text_annotations = FALSE),
  permutation = d$perm_summary, posterior = d$models,
  notes = c("Current input files hashed in verified_inputs.",
            "Posterior summaries recomputed from current NUTS draws; rank probabilities sum to one.",
            "Posterior draws are paired by chain then draw, in the current analysis.",
            "All 5040 observation-pattern-preserving whole-trajectory permutations are enumerated; median-centred two-sided P.",
            "Wilson intervals are crude descriptive binomial intervals, not cluster-adjusted.",
            "No fitting or random posterior resampling is performed.")),
  file.path(source_dir, "build_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
writeLines(capture.output(sessionInfo()), file.path(source_dir, "sessionInfo.txt"))
message("Complete: ", normalizePath(out))
