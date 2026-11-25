# Shared ggplot2 style and provenance for the sentinel-network figure suite.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
})
palette <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF",
             "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF")
rule_colors <- c(Culture=palette[2], Capture=palette[4], PastHuman=palette[1],
                 HumanPooled=palette[1], Random=palette[6])
stream_colors <- c(Human=palette[1], Culture=palette[2], Capture=palette[4], Joint=palette[5])
model_colors <- c(Binomial=palette[2], `Beta-binomial`=palette[3])
rule_labels <- c(Culture="Kidney culture", Capture="Capture yield",
                 PastHuman="Prior human serology", HumanPooled="Pooled human serology",
                 Random="Uniform random")
rule_shapes <- c(Culture=16, Capture=15, PastHuman=18, HumanPooled=18, Random=1)
ink <- "#252525"
theme_set(theme_bw(base_size=12, base_family="DejaVu Sans") + theme(
  text=element_text(colour=ink),
  plot.title=element_blank(), plot.subtitle=element_blank(), plot.caption=element_blank(),
  plot.tag=element_text(face="bold", size=18, colour=ink),
  plot.tag.position="topleft", plot.title.position="panel",
  panel.grid.minor=element_blank(),
  panel.grid.major=element_line(colour="#E8E8E8", linewidth=.28),
  panel.grid.major.y=element_blank(),
  panel.border=element_rect(colour="#555555", fill=NA, linewidth=.5),
  axis.title=element_text(size=11), axis.text=element_text(size=10, colour=ink),
  legend.position="bottom", legend.title=element_text(size=10),
  legend.text=element_text(size=9), legend.key.height=grid::unit(4,"mm"),
  legend.margin=margin(0,0,0,0),
  strip.text=element_text(face="bold",size=10),
  strip.background=element_rect(fill="white",colour="#666666",linewidth=.4),
  panel.spacing=grid::unit(16,"pt"), plot.margin=margin(12,12,12,12)
))
.figure_inputs <- character()
.panel_registry <- list()
.figure_outputs <- character()
.figure_specs <- list()

read_plot_csv <- function(path) {
  .figure_inputs <<- unique(c(.figure_inputs, normalizePath(path)))
  readr::read_csv(path, show_col_types=FALSE)
}
read_network <- function(name) read_plot_csv(file.path(network_dir,name))
read_legacy <- function(name) read_plot_csv(file.path(legacy_dir,name))
read_table <- function(name) read_plot_csv(file.path(root,"manuscript/tables",name))

panel_data <- function(data, figure, tag, description, suffix="") {
  name <- paste0(figure,"_",tag,if(nchar(suffix)) paste0("_",suffix) else "", ".csv")
  readr::write_csv(data, file.path(source_dir,name), na="")
  key <- paste(figure,tag,sep="_")
  .panel_registry[[key]] <<- list(figure=figure, panel=tag, description=description,
                                source=name, rows=nrow(data), engine="ggplot2")
  invisible(data)
}

pretty_scale <- function(p) {
  built <- ggplot_build(p)
  for(axis in c("x","y")) {
    sc <- p$scales$get_scales(axis)
    trained <- if(axis=="x") built$layout$panel_scales_x[[1]] else built$layout$panel_scales_y[[1]]
    if(is.null(sc) && !trained$is_discrete())
      p <- p + if(axis=="x") scale_x_continuous() else scale_y_continuous()
  }
  for(i in seq_along(p$scales$scales)) {
    sc <- p$scales$scales[[i]]$clone()
    if(sc$is_discrete()) next
    ordinal <- is.character(sc$labels) ||
      (is.numeric(sc$breaks) && (identical(as.numeric(sc$breaks),as.numeric(1:8)) ||
       identical(as.numeric(sc$breaks),as.numeric(1:7)) || all(sc$breaks>=1900)))
    if(ordinal) next
    domain <- sc$limits
    sc$limits <- if(is.numeric(domain)) range(pretty(domain)) else function(x) range(pretty(x))
    sc$breaks <- function(x) pretty(x)
    sc$expand <- c(0,0)
    p$scales$scales[[i]] <- sc
  }
  p
}
panel <- function(p,tag) pretty_scale(p) + labs(tag=tag,title=NULL,subtitle=NULL,caption=NULL)
percent <- scales::label_percent(accuracy=1)
site_factor <- function(x) factor(x,levels=rev(display_order))
rule_colour_scale <- function(...) scale_colour_manual(values=rule_colors,labels=rule_labels,name=NULL,...)
rule_fill_scale <- function(...) scale_fill_manual(values=rule_colors,labels=rule_labels,name=NULL,...)

save_network_figure <- function(plot,stem,width=17,height=10.5,main_figure=FALSE) {
  keys <- vapply(.panel_registry,function(x)x$figure,character(1))==stem
  tags <- unique(vapply(.panel_registry[keys],function(x)x$panel,character(1)))
  stopifnot(length(tags)>=6, all(LETTERS[seq_along(tags)] %in% tags))
  plot_tags <- function(p) {
    own <- p$labels$tag
    if(!is.character(own)) own <- character()
    if(inherits(p,"patchwork"))
      own <- c(own,unlist(lapply(p$patches$plots,plot_tags),use.names=FALSE))
    own
  }
  actual_tags <- unique(plot_tags(plot))
  stopifnot(setequal(actual_tags,tags))
  .figure_specs[[stem]] <<- list(panels=sort(tags), panel_count=length(tags),
                                width_inches=width,height_inches=height,dpi=320,
                                engine="ggplot2",composition="patchwork")
  for(ext in c("png","pdf",if(main_figure) "tif")) {
    target <- file.path(out,paste0(stem,".",ext))
    if(ext=="png") ggsave(target,plot,width=width,height=height,dpi=320,bg="white",limitsize=FALSE)
    else if(ext=="tif") ggsave(target,plot,width=width,height=height,device="tiff",
                              dpi=320,compression="lzw",bg="white",limitsize=FALSE)
    else ggsave(target,plot,width=width,height=height,device=grDevices::cairo_pdf,bg="white",limitsize=FALSE)
    .figure_outputs <<- c(.figure_outputs,normalizePath(target))
  }
}
