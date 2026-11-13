# Figure 1: indicator choice reshapes the selected network.
# Figure 3: historical information is carried forward and updated.
build_figures_13 <- function() {
  site <- read_network("site_summary.csv")
  selections <- read_network("selection_weights.csv")
  folds <- read_network("fold_metrics.csv")
  pooled <- read_network("pooled_metrics.csv")
  contributions <- read_network("site_year_contributions.csv")
  group_levels <- c("Shared", "Culture only", "Capture only", "Neither")
  group_colors <- setNames(c(palette[5], palette[2], palette[4], "#B6B6B6"), group_levels)
  group_labels <- c("Shared", "Culture only", "Capture only", "Neither")
  interval <- function(y,n) {
    z <- qnorm(.975); p <- y/n; div <- 1+z^2/n
    middle <- (p+z^2/(2*n))/div
    half <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))/div
    tibble(low=middle-half,high=middle+half)
  }
  stem1 <- "figure1_indicator_priorities"
  site <- site |> mutate(site_label=site_factor(site),
    member_group=factor(case_when(
      culture_rank<=3 & capture_rank<=3 ~ "Shared",
      culture_rank<=3 ~ "Culture only", capture_rank<=3 ~ "Capture only",
      TRUE ~ "Neither"),levels=group_levels))
  stopifnot(setequal(site$site[site$culture_rank<=3],site$site[site$human_rank<=3]),
    sum(site$culture_rank<=3 & site$capture_rank<=3)==1)

  site_rows <- site
  rates <- bind_rows(lapply(c("human","culture","capture"),function(s) tibble(
    site=site_rows$site,site_label=site_rows$site_label,stream=tools::toTitleCase(s),
    rate=site_rows[[paste0(s,"_fraction")]],low=site_rows[[paste0(s,"_low")]],
    high=site_rows[[paste0(s,"_high")]],rank=site_rows[[paste0(s,"_rank")]],
    member_group=site_rows$member_group))) |>
    mutate(stream=factor(stream,c("Human","Culture","Capture")))
  panel_data(rates,stem1,"A","Pooled indicator magnitudes and descriptive Wilson intervals")
  f1a <- ggplot(rates,aes(rate,site_label,colour=stream)) +
    geom_segment(aes(x=low,xend=high,yend=site_label),linewidth=.85) +
    geom_point(aes(shape=rank<=3),size=2.4,stroke=.9) +
    facet_wrap(~stream,nrow=1,scales="free_x") +
    scale_colour_manual(values=stream_colors,guide="none") +
    scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1),guide="none") +
    scale_x_continuous(labels=percent) + coord_cartesian(clip="off") +
    labs(x="Pooled proportion",y=NULL) +
    theme(panel.spacing=grid::unit(24,"pt"),axis.text.y=element_text(size=9.5))
  f1a <- panel(f1a,"A")
  # Facet endpoints need fewer ticks at journal-column width.
  a_scale <- f1a$scales$get_scales("x")
  a_scale$breaks <- function(x) pretty(x,n=3)

  panel_data(rates,stem1,"B","Site rank trajectories across human, culture and capture indicators")
  f1b <- ggplot(rates,aes(stream,rank,group=site,colour=member_group)) +
    geom_hline(yintercept=3.5,colour="#9A9A9A",linetype="dashed",linewidth=.45) +
    geom_line(linewidth=.9,alpha=.8) + geom_point(size=2.8) +
    scale_y_reverse(breaks=1:8,limits=c(8.35,.65)) +
    scale_colour_manual(values=group_colors,guide="none") +
    scale_x_discrete(expand=expansion(add=.22)) +
    labs(x=NULL,y="Site rank (1 = highest)") +
    theme(panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(colour="#ECECEC",linewidth=.25))
  f1b <- panel(f1b,"B")

  member <- rates |> mutate(selected=as.integer(rank<=3),
    tile=if_else(selected==1,as.character(stream),"Unselected"))
  panel_data(member,stem1,"C","Membership of pooled three-site portfolios for each indicator")
  f1c <- ggplot(member,aes(stream,site_label)) +
    geom_tile(aes(fill=tile),colour="white",linewidth=.8,width=.94,height=.92) +
    geom_point(data=member |> filter(selected==1),colour="white",size=1.7) +
    scale_fill_manual(values=c(stream_colors,Unselected="#F0F0F0"),guide="none") +
    labs(x="Pooled top-three membership",y=NULL) +
    theme(panel.grid=element_blank(),axis.ticks=element_blank())
  f1c <- panel(f1c,"C")

  contribution <- contributions |> filter(design=="loyo87") |>
    group_by(site) |> summarise(contribution=sum(contribution),.groups="drop") |>
    mutate(site_label=site_factor(site),direction=case_when(contribution>0~"Culture",
      contribution<0~"Capture",TRUE~"Neither"))
  stopifnot(sum(contribution$contribution)==63,
    contribution$contribution[contribution$site=="Xiangtan"]==79,
    contribution$contribution[contribution$site=="Lianyuan"]==37,
    contribution$contribution[contribution$site=="Liuyang"]==-35,
    contribution$contribution[contribution$site=="Yuanjiang"]==-18)
  panel_data(contribution,stem1,"D","Site contributions to the LOYO three-site positive-coverage difference")
  f1d <- ggplot(contribution,aes(contribution,site_label,colour=direction)) +
    geom_vline(xintercept=0,colour="#777777",linewidth=.5) +
    geom_segment(aes(x=0,xend=contribution,yend=site_label),linewidth=2.2) +
    geom_point(size=3) +
    scale_colour_manual(values=c(rule_colors,Neither="#B6B6B6"),guide="none") +
    labs(x="Additional human positives covered\n(culture minus capture, LOYO)",y=NULL)
  f1d <- panel(f1d,"D")

  grouped <- site |> group_by(member_group) |>
    summarise(sites=paste(sort(site),collapse=";"),n_sites=n(),
      human_pos=sum(human_pos),human_n=sum(human_n),
      culture_pos=sum(culture_pos),culture_n=sum(culture_n),.groups="drop") |>
    mutate(rate=human_pos/human_n,group_label=factor(member_group,levels=rev(group_levels)))
  grouped <- bind_cols(grouped,interval(grouped$human_pos,grouped$human_n))
  stopifnot(sum(grouped$human_pos)==291,sum(grouped$human_n)==8834,
    sum(grouped$culture_pos)==317,sum(grouped$culture_n)==6824)
  panel_data(grouped,stem1,"E","Human seroreactivity among four classes of pooled animal-portfolio membership")
  f1e <- ggplot(grouped,aes(rate,group_label,colour=member_group)) +
    geom_segment(aes(x=low,xend=high,yend=group_label),linewidth=1.1) +
    geom_point(aes(size=human_n),shape=16) +
    scale_colour_manual(values=group_colors,guide="none") +
    scale_size_continuous(range=c(3,6),breaks=c(1000,2000,3000),name="Human tests") +
    scale_x_continuous(labels=percent,limits=c(0,.10)) +
    labs(x="Human seroreactivity",y="Membership in pooled animal trios") +
    guides(size=guide_legend(override.aes=list(colour="#555555")))
  f1e <- panel(f1e,"E")

  composition <- grouped |> select(member_group,human_pos,culture_pos,human_n) |>
    pivot_longer(-member_group,names_to="quantity",values_to="count") |>
    group_by(quantity) |> mutate(share=count/sum(count),total=sum(count)) |> ungroup() |>
    mutate(quantity=factor(quantity,c("human_pos","culture_pos","human_n"),
      c("Human\npositives","Culture\npositives","Human\ntests")))
  stopifnot(all(abs((composition |> group_by(quantity) |> summarise(z=sum(share)))$z-1)<1e-12))
  panel_data(composition,stem1,"F","Distribution of positive results and human testing across animal-portfolio membership classes")
  f1f <- ggplot(composition,aes(quantity,share,fill=member_group)) +
    geom_col(width=.65,colour="white",linewidth=.5) +
    scale_fill_manual(values=group_colors,labels=group_labels,name=NULL,drop=FALSE) +
    scale_y_continuous(labels=percent,limits=c(0,1)) +
    labs(x=NULL,y="Share of network total") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme(panel.grid.major.x=element_blank())
  f1f <- panel(f1f,"F")
  fig1 <- ((f1a+f1b+f1c)+plot_layout(widths=c(1.45,1,1))) /
          ((f1d+f1e+f1f)+plot_layout(widths=c(1.45,1,1)))
  fig1 <- fig1 + plot_layout(heights=c(1,1),guides="keep")
  save_network_figure(fig1,stem1,width=18,height=10.8)

  stem3 <- "figure3_forward_evaluation"
  forward <- folds |> filter(design=="forward_expanding",k==3) |>
    mutate(rule=factor(rule,c("Culture","Capture","PastHuman","Random")))
  timeline <- folds |> filter(design %in% c("forward_expanding","forward_rolling3"),
      rule=="Culture",k==3) |>
    select(design,evaluation_year,training_start,training_end) |>
    mutate(evaluation_label=factor(evaluation_year),
      design_label=factor(design,c("forward_expanding","forward_rolling3"),
        c("Cumulative history","Recent three years")))
  stopifnot(nrow(timeline)==16,all(timeline$training_end<timeline$evaluation_year))
  panel_data(timeline,stem3,"A","Training windows and subsequent evaluation years for both updating approaches")
  f3a <- ggplot(timeline,aes(y=evaluation_label)) +
    geom_segment(aes(x=training_start,xend=training_end,yend=evaluation_label),
      linewidth=4,colour="#91A7B5",lineend="round") +
    geom_segment(aes(x=training_end,xend=evaluation_year,yend=evaluation_label),
      linewidth=.4,colour="#777777",linetype="dotted") +
    geom_point(aes(x=evaluation_year),shape=18,size=2.8,colour="#353535") +
    facet_wrap(~design_label,nrow=1) +
    scale_x_continuous(breaks=seq(2015,2025,2),labels=function(x)substr(x,3,4),
      limits=c(2014.5,2025.5),expand=c(0,0)) +
    labs(x="Training and evaluation year (20xx)",y="Evaluation year") +
    theme(panel.spacing=grid::unit(10,"pt"),panel.grid.major.x=element_line(colour="#EAEAEA",linewidth=.25))
  f3a <- panel(f3a,"A")

  forward_selection <- selections |> filter(design=="forward_expanding",k==3,
      rule %in% c("Culture","Capture","PastHuman")) |>
    mutate(site_label=site_factor(site),rule=factor(rule,c("Culture","Capture","PastHuman")),
      tile=if_else(weight>0,as.character(rule),"Unselected"))
  stopifnot(all(forward_selection$weight %in% c(0,1)),
    all((forward_selection |> group_by(rule,evaluation_year) |> summarise(w=sum(weight),.groups="drop"))$w==3))
  panel_data(forward_selection,stem3,"B","Three-site portfolio membership under cumulative-history selection")
  f3b <- ggplot(forward_selection,aes(factor(evaluation_year),site_label)) +
    geom_tile(aes(fill=tile),colour="white",linewidth=.6,width=.96,height=.94) +
    geom_point(data=forward_selection |> filter(weight>0),colour="white",size=.7) +
    facet_wrap(~rule,nrow=1,labeller=as_labeller(c(Culture="Culture",Capture="Capture",PastHuman="Prior human"))) +
    scale_fill_manual(values=c(rule_colors,Unselected="#F0F0F0"),guide="none") +
    scale_x_discrete(labels=function(x)substr(x,3,4)) +
    labs(x="Evaluation year (20xx)",y=NULL) +
    theme(panel.grid=element_blank(),axis.ticks=element_blank(),
      axis.text.x=element_text(size=8),axis.text.y=element_text(size=9),panel.spacing=grid::unit(9,"pt"))
  f3b <- panel(f3b,"B")

  annual <- forward |> filter(rule!="Random")
  panel_data(annual,stem3,"C","Annual positive-result coverage with point area reflecting evaluation-year positive count")
  f3c <- ggplot(annual,aes(evaluation_year,coverage,colour=rule,group=rule)) +
    geom_line(data=annual |> filter(rule!="PastHuman"),linewidth=.85) +
    geom_point(aes(shape=rule,size=positive_all),stroke=1) +
    scale_colour_manual(values=rule_colors,guide="none") +
    scale_shape_manual(values=c(Culture=16,Capture=15,PastHuman=1),guide="none") +
    scale_size_area(max_size=6.5,breaks=c(1,10,25,40),name="Human positives") +
    scale_x_continuous(breaks=2018:2025,labels=function(x)substr(x,3,4),
      limits=c(2017.7,2025.3),expand=c(0,0)) +
    scale_y_continuous(labels=percent,limits=c(0,1)) +
    coord_cartesian(clip="off") +
    labs(x="Evaluation year (20xx)",y="Annual positive coverage") +
    guides(size=guide_legend(override.aes=list(colour="#666666",shape=16)))
  f3c <- panel(f3c,"C")

  cumulative <- forward |> arrange(rule,evaluation_year) |> group_by(rule) |>
    mutate(cumulative_selected=cumsum(positive_selected),cumulative_total=cumsum(positive_all),
      cumulative_coverage=cumulative_selected/cumulative_total) |> ungroup()
  end <- cumulative |> filter(evaluation_year==2025) |> arrange(rule)
  stopifnot(end$cumulative_selected[end$rule=="Culture"]==56,
    end$cumulative_selected[end$rule=="Capture"]==40,
    all(end$cumulative_total==93))
  panel_data(cumulative,stem3,"D","Cumulative coverage as subsequent-year outcomes accumulate")
  f3d <- ggplot(cumulative,aes(evaluation_year,cumulative_coverage,colour=rule,group=rule)) +
    geom_line(data=cumulative |> filter(rule!="PastHuman"),aes(linetype=rule),linewidth=.95) +
    geom_point(data=cumulative |> filter(rule!="Random"),aes(shape=rule),size=2.6,stroke=.9) +
    scale_colour_manual(values=rule_colors,labels=rule_labels,name=NULL,
      limits=c("Culture","Capture","PastHuman","Random")) +
    scale_shape_manual(values=c(Culture=16,Capture=15,PastHuman=1),labels=rule_labels,name=NULL) +
    scale_linetype_manual(values=c(Culture="solid",Capture="solid",Random="dashed"),labels=rule_labels,name=NULL) +
    scale_x_continuous(breaks=2018:2025,labels=function(x)substr(x,3,4),
      limits=c(2017.7,2025.3),expand=c(0,0)) +
    scale_y_continuous(labels=percent,limits=c(0,1)) +
    labs(x="Evaluation through year (20xx)",y="Cumulative positive coverage") +
    guides(colour=guide_legend(nrow=2,byrow=TRUE,
      override.aes=list(shape=c(16,15,1,NA),linetype=c("solid","solid","blank","dashed"))),
      shape="none",linetype="none")
  f3d <- panel(f3d,"D")

  prior <- pooled |> filter(design=="forward_expanding",rule=="PastHuman") |>
    select(k,human_coverage=coverage)
  difference <- pooled |> filter(design=="forward_expanding",rule %in% c("Culture","Capture")) |>
    left_join(prior,by="k") |> mutate(difference_pp=100*(coverage-human_coverage))
  stopifnot(abs(difference$difference_pp[difference$rule=="Culture" & difference$k==3])<1e-12)
  panel_data(difference,stem3,"E","Animal-rule coverage relative to the prior-human-serology benchmark at every network size")
  f3e <- ggplot(difference,aes(k,difference_pp,colour=rule,shape=rule)) +
    geom_hline(yintercept=0,colour="#777777",linetype="dashed",linewidth=.5) +
    geom_line(linewidth=1) + geom_point(size=2.8) +
    scale_colour_manual(values=rule_colors,guide="none") +
    scale_shape_manual(values=rule_shapes,guide="none") +
    scale_x_continuous(breaks=1:7,limits=c(.8,7.2),expand=c(0,0)) +
    labs(x="Number of sites selected",y="Coverage minus prior human serology\n(percentage points)")
  f3e <- panel(f3e,"E")

  update_effect <- pooled |> filter(design %in% c("forward_expanding","forward_rolling3"),
      rule %in% c("Culture","Capture","PastHuman")) |>
    select(design,rule,k,coverage) |> pivot_wider(names_from=design,values_from=coverage) |>
    mutate(difference_pp=100*(forward_rolling3-forward_expanding))
  stopifnot(all(abs(update_effect$difference_pp[update_effect$k==3 & update_effect$rule!="PastHuman"])<1e-12),
    abs(update_effect$difference_pp[update_effect$k==3 & update_effect$rule=="PastHuman"]-200/93)<1e-10)
  panel_data(update_effect,stem3,"F","Change in coverage when rankings use recent-three-year rather than cumulative history")
  f3f <- ggplot(update_effect,aes(k,difference_pp,colour=rule,shape=rule)) +
    geom_hline(yintercept=0,colour="#777777",linetype="dashed",linewidth=.5) +
    geom_line(linewidth=1) + geom_point(size=2.8) +
    scale_colour_manual(values=rule_colors,guide="none") +
    scale_shape_manual(values=rule_shapes,guide="none") +
    scale_x_continuous(breaks=1:7,limits=c(.8,7.2),expand=c(0,0)) +
    labs(x="Number of sites selected",y="Recent minus cumulative coverage\n(percentage points)")
  f3f <- panel(f3f,"F")
  fig3 <- ((f3a+f3b+f3c)+plot_layout(widths=c(1.12,1.3,1))) /
          ((f3d+f3e+f3f)+plot_layout(widths=c(1.12,1.3,1)))
  fig3 <- fig3 + plot_layout(heights=c(1,1),guides="keep")
  save_network_figure(fig3,stem3,width=18,height=10.8)
  invisible(list(figure1=fig1,figure3=fig3))
}
