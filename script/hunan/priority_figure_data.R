# Current 87-row analysis. Numerical work is performed by reanalyze_priority_87.py.
prepare_priority_data <- function(root) {
  base <- file.path(root, 'results/priority_87_2026-10-22')
  used <- character()
  input <- function(name) { used <<- unique(c(used,name)); file.path(base,name) }
  read <- function(name) readr::read_csv(input(name),show_col_types=FALSE)
  panel <- read('analysis_panel.csv')
  stopifnot(nrow(panel)==87, sum(panel$positive)==291, sum(panel$tested)==8834,
            sum(panel$leptospira_culture_positive)==317, sum(panel$effective_cultures)==6824)
  site <- read('primary87_site_summary.csv') |> mutate(capture_yield=capture_fraction)
  display_order <- site$site[order(site$human_fraction,decreasing=TRUE)]
  site_factor <- function(x) factor(x,levels=rev(display_order))
  site$site_label <- site_factor(site$site)
  # tibble() masks the data-frame name once a column called site is created.
  site_rates <- site
  rates <- bind_rows(lapply(c('human','culture','capture'),function(s) tibble(
    site=site_rates$site,stream=tools::toTitleCase(s),rate=site_rates[[paste0(s,'_fraction')]],
    low=site_rates[[paste0(s,'_low')]],high=site_rates[[paste0(s,'_high')]]))) |>
    mutate(stream=factor(stream,c('Human','Culture','Capture')),site_label=site_factor(site)) |> group_by(stream) |> mutate(rank=rank(-rate)) |> ungroup()
  summary <- jsonlite::fromJSON(input('primary87_summary.json'))
  correlations <- bind_rows(lapply(c('culture','capture'),function(s)
    tibble(stream=s,rho=summary$correlations[[s]]$rho,p=summary$correlations[[s]]$p_two_sided)))
  folds <- read('primary87_loyo_by_year.csv') |> mutate(
    positive_coverage_C=coverage,testing_share_S=testing_share,
    detection_enrichment_E_C_over_S=enrichment,positive_all_sites=positive_all)
  membership <- read('primary87_loyo_selection.csv') |> mutate(
    held_out_year=year,selected=selection_probability>0,site_label=site_factor(site),
    status=factor(if_else(!available,'No monitoring',if_else(selected,'Selected','Not selected')),
                  levels=c('Not selected','Selected','No monitoring')))
  pooled <- read('primary87_loyo_pooled.csv') |> mutate(C=coverage,S=testing_share)
  pe <- summary$permutation$enrichment
  perm_summary <- list(observed=pe$observed_difference,median=pe$null_median,
                      interval=unlist(pe$null_interval),p=pe$p_two_sided,
                      permutations=summary$permutation$n_permutations)
  permutation <- read('primary87_permutations.csv') |> mutate(delta=delta_enrichment,
    tail=abs(delta-pe$null_median)>=abs(pe$observed_difference-pe$null_median)-1e-12)
  model_factor <- function(x) factor(x,c('binomial','betabinomial'),c('Binomial','Beta-binomial'))
  intervals <- read('posterior_site_rates.csv') |> mutate(model=model_factor(model),site_label=site_factor(site))
  ranks <- read('posterior_ranks.csv') |> mutate(model=model_factor(model),site_label=site_factor(site))
  inclusion <- intervals |> transmute(model,stream=tools::toTitleCase(stream),site,site_label,probability=top3)
  joint <- read('posterior_joint.csv') |> transmute(model=model_factor(model),stream='Joint',site,
                                                site_label=site_factor(site),probability=joint_top3)
  inclusion <- bind_rows(inclusion,joint) |> mutate(stream=factor(stream,c('Human','Culture','Joint')))
  rho <- read('posterior_rho_draws.csv') |> mutate(model=model_factor(model))
  models <- read('posterior_correlations.csv') |> mutate(model=model_factor(model),median=rho_median,low=rho_low,high=rho_high)
  stopifnot(nrow(permutation)==5040,all(abs((ranks |> group_by(model,stream,site) |>
                  summarise(total=sum(probability),.groups='drop'))$total-1)<1e-10))
  list(site=site,rates=rates,correlations=correlations,panel=panel,display_order=display_order,
       folds=folds,membership=membership,pooled=pooled,permutation=permutation,perm_summary=perm_summary,
       intervals=intervals,ranks=ranks,inclusion=inclusion,rho=rho,models=models,
       sources=tibble(path=used,sha256=vapply(file.path(base,used),digest::digest,'',algo='sha256',file=TRUE)))
}
