# Six-panel supporting figures for the sentinel-network evaluation.
# Sourced by the shared ggplot2 driver after network_plot_common.R.
build_figures_support <- function() {
  model_name <- function(x) factor(x, c("binomial", "betabinomial"),
                                  c("Binomial", "Beta-binomial"))
  stream_name <- function(x) factor(x, c("human", "culture"), c("Human", "Culture"))
  counts <- read_legacy("analysis_panel.csv")
  cells <- read_legacy("posterior_predictive_cells.csv")
  posterior_checks <- read_legacy("posterior_predictive_checks.csv")
  stopifnot(nrow(counts)==87, sum(counts$positive)==291,
            sum(counts$leptospira_culture_positive)==317)
  replicate_counts <- bind_rows(lapply(c("human","culture"), function(s) {
    bind_rows(lapply(c("binomial","betabinomial"), function(m) {
      read_legacy(paste0(s,"_",m,"_predictive_totals.csv")) |>
        mutate(stream=s, model=model_name(m), draw=row_number())
    }))
  }))
  stopifnot(all((replicate_counts |> count(stream,model))$n==6000))

  # S1: Predictive calibration of sparse counts, then the actual site-year cells.
  s1 <- "figure_s1_predictive_checks"
  predictive_distribution <- function(s, variable, tag, xlabel) {
    d <- replicate_counts |> filter(stream==s)
    observed <- if(variable=="total") {
      if(s=="human") sum(counts$positive) else sum(counts$leptospira_culture_positive)
    } else {
      if(s=="human") sum(counts$positive==0) else sum(counts$leptospira_culture_positive==0)
    }
    d <- d |> mutate(value=.data[[variable]], observed=observed)
    breaks <- if(variable=="zeros") seq(-.5,max(d$value)+.5,by=1) else pretty(c(0,d$value),n=35)
    histograms <- bind_rows(lapply(levels(d$model),function(m) {
      h <- hist(d$value[d$model==m],breaks=breaks,plot=FALSE,include.lowest=TRUE)
      tibble(model=factor(m,levels=levels(d$model)),value=h$mids,
             lower=head(h$breaks,-1),upper=tail(h$breaks,-1),
             density=h$density,probability=h$counts/sum(h$counts),draws=h$counts,
             observed=observed)
    }))
    stopifnot(all(abs((histograms |> group_by(model) |> summarise(p=sum(probability),.groups="drop"))$p-1)<1e-10))
    panel_data(histograms,s1,tag,paste(s,variable,"posterior-predictive distribution"))
    p <- ggplot(histograms,aes(value,density,colour=model)) +
      geom_line(linewidth=.95) +
      geom_vline(xintercept=observed, colour=ink, linetype="dashed", linewidth=.65) +
      scale_colour_manual(values=model_colors,name=NULL) +
      scale_x_continuous(limits=c(0,max(histograms$upper))) +
      labs(x=xlabel,y="Posterior predictive density")
    panel(p,tag)
  }
  cell_intervals <- function(s, tag, ylabel) {
    d <- cells |> filter(stream==s,model=="betabinomial") |>
      arrange(observed,year,site) |> mutate(ordered_cell=row_number(),
                                        covered=observed>=low & observed<=high)
    check <- posterior_checks |> filter(fit==paste0(s,"_betabinomial"))
    stopifnot(nrow(d)==87, sum(d$covered)==check$cells_covered)
    d <- d |> mutate(interval_status=factor(if_else(covered,"Inside","Outside"),levels=c("Outside","Inside")))
    panel_data(d,s1,tag,paste(s,"observations with beta-binomial predictive intervals"))
    p <- ggplot(d,aes(ordered_cell,observed)) +
      geom_linerange(aes(ymin=low,ymax=high),colour=palette[3],alpha=.55,linewidth=.5) +
      geom_point(aes(shape=interval_status),colour=ink,size=1.6) +
      scale_shape_manual(values=c(Inside=16,Outside=4),breaks=c("Outside","Inside"),
                         limits=c("Outside","Inside"),drop=FALSE,
                         labels=c(Inside="Inside 95% interval",Outside="Outside 95% interval"),name=NULL) +
      labs(x="Site–years ordered by observed count",y=ylabel) +
      theme(legend.text=element_text(size=8.5))
    p <- panel(p,tag)
    if(s=="culture") p <- p + guides(shape="none")
    p
  }
  s1a <- predictive_distribution("human","total","A","Human positives across the network")
  s1b <- predictive_distribution("human","zeros","B","Human zero-positive site–years")
  s1c <- cell_intervals("human","C","Human positives per site–year")
  s1d <- predictive_distribution("culture","total","D","Culture positives across the network")
  s1e <- predictive_distribution("culture","zeros","E","Culture zero-positive site–years")
  s1f <- cell_intervals("culture","F","Culture positives per site–year")
  s1plot <- ((s1a|s1b|s1c)/(s1d|s1e|s1f)) + plot_layout(guides="collect") &
    theme(legend.position="bottom")
  save_network_figure(s1plot,s1,width=17,height=10.5)

  # S2: Sensitivity of the concrete portfolio comparison to distinct choices.
  s2 <- "figure_s2_period_prior_sensitivity"
  window_names <- c(primary87="All 87 site–years",balanced80="Balanced panel",
                    early63="2015–2022",exclude2015="Exclude 2015",
                    recent24="2023–2025",arithmetic_effective="Arithmetic culture\ndenominator")
  windows <- bind_rows(lapply(names(window_names),function(prefix) {
    d <- read_legacy(paste0(prefix,"_loyo_pooled.csv")) |>
      filter(rule %in% c("Culture","Capture"))
    stopifnot(nrow(d)==2,length(unique(d$positive_all))==1)
    tibble(window=prefix,label=window_names[[prefix]],
           positive_all=d$positive_all[1],
           culture_coverage=d$coverage[d$rule=="Culture"],
           capture_coverage=d$coverage[d$rule=="Capture"],
           difference_pp=100*(d$coverage[d$rule=="Culture"]-d$coverage[d$rule=="Capture"]))
  })) |> mutate(label=factor(label,levels=rev(unname(window_names))))
  panel_data(windows,s2,"A","Observed three-site LOYO coverage difference across record definitions")
  s2a <- ggplot(windows,aes(difference_pp,label)) +
    geom_segment(aes(x=0,xend=difference_pp,yend=label),colour=palette[6],linewidth=1) +
    geom_point(aes(size=positive_all),colour=palette[2]) +
    geom_vline(xintercept=0,colour="#888888",linewidth=.5) +
    scale_size_area(max_size=5,breaks=c(11,100,200,291),name="Human positives") +
    labs(x="Observed coverage difference (percentage points)",y=NULL)
  s2a <- panel(s2a,"A") + guides(size=guide_legend(nrow=1))

  sensitivity <- read_legacy("posterior_decision_sensitivity.csv")
  prior_labels <- c(human_betabinomial="Primary beta-binomial",
                    human_sensitivity_site035="Site contrast SD 0.35",
                    human_sensitivity_site140="Site contrast SD 1.40",
                    human_sensitivity_year100="Year contrast SD 1.00",
                    human_sensitivity_phi_mean10="Concentration mean 10",
                    human_sensitivity_phi_mean40="Concentration mean 40")
  priors <- sensitivity |> filter(fit %in% names(prior_labels)) |>
    mutate(label=factor(prior_labels[fit],levels=rev(unname(prior_labels))),
           across(c(mean,low,high),~100*.x),primary=fit=="human_betabinomial")
  stopifnot(nrow(priors)==6,all(priors$low<priors$mean),all(priors$high>priors$mean))
  panel_data(priors,s2,"B","Mean paired expected-coverage differences under prior specifications; 95% CrIs")
  s2b <- ggplot(priors,aes(mean,label)) +
    geom_vline(xintercept=0,linetype="dashed",colour="#888888") +
    geom_segment(aes(x=low,xend=high,yend=label),colour=palette[3],linewidth=.9) +
    geom_point(aes(shape=primary),colour=palette[3],size=3) +
    scale_shape_manual(values=c(`TRUE`=18,`FALSE`=16),guide="none") +
    labs(x="Expected coverage difference (percentage points)",y=NULL)
  s2b <- panel(s2b,"B")

  delta <- read_network("posterior_paired_draws.csv") |>
    mutate(model=model_name(family),difference_pp=100*delta_coverage)
  stopifnot(all((delta |> count(model))$n==6000))
  panel_data(delta,s2,"C","Paired posterior coverage differences by likelihood")
  s2c <- ggplot(delta,aes(difference_pp,colour=model)) +
    geom_freqpoly(aes(y=after_stat(density)),bins=45,linewidth=1) +
    geom_vline(xintercept=0,linetype="dashed",colour="#888888") +
    scale_colour_manual(values=model_colors,name=NULL) +
    labs(x="Expected coverage difference (percentage points)",y="Posterior density")
  s2c <- panel(s2c,"C") + guides(colour="none")

  scenario_labels <- c(human_betabinomial="Primary beta-binomial",human_balanced80="Balanced panel",
                       human_early63="2015–2022",human_binomial="Binomial likelihood")
  decisions <- read_legacy("posterior_decisions.csv") |>
    filter(fit %in% names(scenario_labels)) |>
    mutate(label=factor(scenario_labels[fit],levels=rev(unname(scenario_labels))),
           rule=recode(strategy,`Uniform random`="Random"))
  expected <- decisions |> filter(rule %in% c("Culture","Capture")) |>
    mutate(y=as.numeric(label)+if_else(rule=="Culture",.12,-.12),
           across(c(mean,low,high),~100*.x))
  panel_data(expected,s2,"D","Fixed-trio expected coverage across likelihood and evaluation-window fits")
  s2d <- ggplot(expected,aes(mean,y,colour=rule)) +
    geom_segment(aes(x=low,xend=high,yend=y),linewidth=.85) +
    geom_point(aes(shape=rule),size=2.8) +
    rule_colour_scale() + scale_shape_manual(values=rule_shapes,labels=rule_labels,name=NULL) +
    scale_y_continuous(breaks=1:4,labels=rev(unname(scenario_labels)),limits=c(.5,4.5)) +
    labs(x="Expected human-positive coverage (%)",y=NULL)
  s2d <- panel(s2d,"D") + guides(colour="none",shape="none")
  regret <- decisions |> filter(rule %in% c("Culture","Capture","Random")) |>
    mutate(y=as.numeric(label)+case_when(rule=="Culture"~.2,rule=="Capture"~0,TRUE~-.2),
           regret_pp=100*expected_regret)
  panel_data(regret,s2,"E","Expected regret relative to the draw-specific best three-site portfolio")
  s2e <- ggplot(regret,aes(regret_pp,y,colour=rule)) +
    geom_point(aes(shape=rule),size=3) +
    geom_segment(aes(x=0,xend=regret_pp,yend=y),alpha=.4,linewidth=.65) +
    rule_colour_scale() + scale_shape_manual(values=rule_shapes,labels=rule_labels,name=NULL) +
    scale_y_continuous(breaks=1:4,labels=rev(unname(scenario_labels)),limits=c(.5,4.5)) +
    labs(x="Expected coverage shortfall (percentage points)",y=NULL)
  s2e <- panel(s2e,"E")
  thresholds <- bind_rows(lapply(levels(delta$model),function(m) {
    v <- delta$difference_pp[delta$model==m]
    tibble(model=factor(m,levels=levels(delta$model)),threshold_pp=seq(0,30,.25),
           probability=vapply(seq(0,30,.25),function(t)mean(v>t),numeric(1)))
  }))
  check <- thresholds |> filter(model=="Beta-binomial",threshold_pp %in% c(0,10,15))
  stopifnot(isTRUE(all.equal(check$probability,c(.8851666666666667,.523,.2995),tolerance=.0001)))
  panel_data(thresholds,s2,"F","Exploratory posterior probability of coverage gain above each threshold")
  s2f <- ggplot(thresholds,aes(threshold_pp,probability,colour=model)) +
    geom_hline(yintercept=.5,colour="#aaaaaa",linetype="dotted") +
    geom_line(linewidth=1) + scale_colour_manual(values=model_colors,name=NULL) +
    scale_y_continuous(labels=percent,limits=c(0,1)) +
    labs(x="Coverage-gain threshold (percentage points)",y="Posterior probability of exceeding threshold")
  s2f <- panel(s2f,"F")
  s2plot <- ((s2a|s2b|s2c)/(s2d|s2e|s2f)) + plot_layout(guides="collect") &
    theme(legend.position="bottom")
  save_network_figure(s2plot,s2,width=18,height=11)

  # S4: Rate estimates -> rank uncertainty -> a complete shared priority set.
  s4 <- "figure4_model_uncertainty"
  sites <- read_legacy("posterior_site_rates.csv") |>
    mutate(model=model_name(model),site_label=site_factor(site))
  rate_intervals <- function(s,tag,xlabel) {
    d <- sites |> filter(stream==s) |>
      mutate(y=as.numeric(site_label)+if_else(model=="Binomial",.12,-.12))
    panel_data(d,s4,tag,paste(s,"year-standardised rates with 95% credible intervals"))
    p <- ggplot(d,aes(median,y,colour=model)) +
      geom_segment(aes(x=low,xend=high,yend=y),linewidth=.85) +
      geom_point(size=2.5) + scale_colour_manual(values=model_colors,name=NULL) +
      scale_y_continuous(breaks=1:8,labels=rev(display_order),limits=c(.5,8.5)) +
      scale_x_continuous(labels=percent) + labs(x=xlabel,y=NULL)
    panel(p,tag)
  }
  s4a <- rate_intervals("human","A","Year-standardised human seroreactivity")
  s4b <- rate_intervals("culture","B","Year-standardised culture positivity")
  ranks <- read_legacy("posterior_ranks.csv") |>
    mutate(model=model_name(model),stream=stream_name(stream),site_label=site_factor(site))
  stopifnot(all(abs((ranks |> group_by(model,stream,site) |> summarise(p=sum(probability),.groups="drop"))$p-1)<1e-10))
  panel_data(ranks,s4,"C","Complete posterior rank probabilities by stream and likelihood")
  s4c <- ggplot(ranks,aes(rank,site_label,fill=probability)) +
    geom_tile(colour="white",linewidth=.12) +
    facet_grid(stream~model) +
    scale_fill_gradient(low="white",high=palette[2],limits=c(0,1),labels=percent,name="Rank probability") +
    scale_x_continuous(breaks=1:8,expand=c(0,0)) +
    labs(x="Site rank",y=NULL) +
    theme(axis.text.y=element_text(size=8.5),panel.spacing=grid::unit(8,"pt"),
          strip.text=element_text(size=9),legend.key.width=grid::unit(18,"mm"))
  s4c <- panel(s4c,"C")

  marginal <- sites |> transmute(model,stream=tools::toTitleCase(stream),site,site_label,probability=top3)
  joint <- read_legacy("posterior_joint.csv") |>
    transmute(model=model_name(model),stream="Joint",site,site_label=site_factor(site),probability=joint_top3)
  inclusion <- bind_rows(marginal,joint) |>
    mutate(stream=factor(stream,c("Human","Culture","Joint")),
           y=as.numeric(site_label)+case_when(stream=="Human"~.18,stream=="Culture"~0,TRUE~-.18))
  panel_data(inclusion,s4,"D","Marginal and joint top-three inclusion probabilities")
  s4d <- ggplot(inclusion,aes(probability,y,colour=stream,shape=stream)) +
    geom_point(size=2.1) + facet_wrap(~model,nrow=1) +
    scale_colour_manual(values=stream_colors,name=NULL) +
    scale_shape_manual(values=c(Human=16,Culture=17,Joint=15),name=NULL) +
    scale_y_continuous(breaks=1:8,labels=rev(display_order),limits=c(.5,8.5)) +
    scale_x_continuous(labels=percent,limits=c(0,1)) +
    labs(x="Probability of top-three inclusion",y=NULL) +
    theme(axis.text.y=element_text(size=9),panel.spacing=grid::unit(9,"pt"))
  s4d <- panel(s4d,"D") +
    coord_cartesian(xlim=c(-.025,1.025),clip="off") +
    theme(panel.spacing=grid::unit(18,"pt"))
  inclusion_x_scale <- s4d$scales$get_scales("x")
  inclusion_x_scale$breaks <- pretty(c(0,1),n=2)
  rho <- read_legacy("posterior_rho_draws.csv") |>
    mutate(model=model_name(model),rho=round(rho,12)) |>
    count(model,rho,name="draws") |> group_by(model) |>
    mutate(probability=draws/sum(draws)) |> ungroup()
  stopifnot(all(abs((rho |> group_by(model) |> summarise(p=sum(probability),.groups="drop"))$p-1)<1e-10))
  panel_data(rho,s4,"E","Discrete posterior probability mass of cross-stream Spearman concordance")
  s4e <- ggplot(rho,aes(rho,probability,colour=model)) +
    geom_segment(aes(xend=rho,y=0,yend=probability),alpha=.45,linewidth=.45) +
    geom_point(size=1.8) + geom_vline(xintercept=0,colour="#aaaaaa",linetype="dashed") +
    scale_colour_manual(values=model_colors,name=NULL) +
    scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(labels=percent) +
    labs(x="Human–culture Spearman correlation",y="Posterior probability mass")
  s4e <- panel(s4e,"E") + guides(colour="none")
  agreements <- read_legacy("posterior_correlations.csv") |>
    mutate(model=model_name(model))
  stopifnot(abs(agreements$same_top3[agreements$model=="Binomial"]-.75)<1e-12,
            abs(agreements$same_top3[agreements$model=="Beta-binomial"]-.2778333333333333)<1e-12)
  panel_data(agreements,s4,"F","Posterior probability that the complete human and culture top-three sets agree")
  s4f <- ggplot(agreements,aes(model,same_top3,fill=model)) +
    geom_col(width=.55,colour="white") +
    scale_fill_manual(values=model_colors,guide="none") +
    scale_y_continuous(labels=percent,limits=c(0,1)) +
    labs(x="Count likelihood",y="Probability of identical top-three sets")
  s4f <- panel(s4f,"F")
  s4plot <- ((s4a|s4b|s4c)/(s4d|s4e|s4f)) + plot_layout(guides="collect") &
    theme(legend.position="bottom")
  save_network_figure(s4plot,s4,width=18,height=11)
  invisible(NULL)
}
