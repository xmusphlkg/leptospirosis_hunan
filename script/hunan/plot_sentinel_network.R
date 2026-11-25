#!/usr/bin/env Rscript
# Render all three main and four supplementary figures with ggplot2/patchwork.
# Usage: Rscript script/hunan/plot_sentinel_network.R [output_directory]
packages <- c("ggplot2","patchwork","dplyr","tidyr","readr","scales","jsonlite","digest")
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Missing R packages: ",paste(missing,collapse=", "))
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- normalizePath(file.path(dirname(script),"../.."))
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else file.path(root,"manuscript/figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out)
network_dir <- file.path(root,"results/network_redesign_2026-11-13")
legacy_dir <- file.path(root,"results/priority_87_2026-10-22")
source_dir <- file.path(out,"source_data/network_ggplot_2026-11-13")
dir.create(source_dir,recursive=TRUE,showWarnings=FALSE)
modules <- c("network_plot_common.R","network_figures_13.R",
             "network_figures_2_s3.R","network_figures_support.R")
for(name in modules) source(file.path(root,"script/hunan",name))
sites <- read_network("site_summary.csv")
display_order <- sites$site[order(sites$human_fraction,decreasing=TRUE)]
message("Rendering Figure 1 and Figure 3...")
build_figures_13()
message("Rendering Figure 2 and Figure S3...")
build_figures_2_s3()
message("Rendering Figures S1, S2 and S4...")
build_figures_support()

relative <- function(path) substring(normalizePath(path),nchar(root)+2)
hash <- function(path) digest::digest(file=path,algo="sha256")
input_names <- vapply(.figure_inputs,relative,character(1))
input_hashes <- setNames(lapply(.figure_inputs,hash),input_names)
input_copies <- setNames(as.list(file.path("inputs",input_names)),input_names)
for(i in seq_along(.figure_inputs)) {
  target <- file.path(source_dir,input_copies[[i]])
  dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  file.copy(.figure_inputs[i],target,overwrite=TRUE)
  stopifnot(identical(hash(target),input_hashes[[i]]))
}
panel_paths <- vapply(.panel_registry,function(x)x$source,character(1))
panel_hashes <- setNames(lapply(file.path(source_dir,panel_paths),hash),panel_paths)
script_paths <- c(file.path(root,"script/hunan/plot_sentinel_network.R"),
                  file.path(root,"script/hunan",modules))
output_hashes <- setNames(lapply(.figure_outputs,hash),basename(.figure_outputs))
stopifnot(length(.figure_specs)==7,length(.figure_outputs)==17,
          all(vapply(.figure_specs,function(x)x$panel_count,integer(1))>=6))
manifest <- list(
  engine="ggplot2",composition="patchwork",style="Hunan eight-colour theme_bw",
  script="script/hunan/plot_sentinel_network.R",script_sha256=hash(script_paths[1]),
  scripts_sha256=setNames(lapply(script_paths,hash),vapply(script_paths,relative,character(1))),
  inputs_sha256=input_hashes, input_copies=input_copies,
  panel_data_sha256=panel_hashes, outputs_sha256=output_hashes,
  figures=.figure_specs, panels=unname(.panel_registry)
)
jsonlite::write_json(manifest,file.path(source_dir,"manifest.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(source_dir,"sessionInfo.txt"))
message("Rendered ",length(.figure_specs)," figures; ",length(.panel_registry),
        " panels. Provenance: ",source_dir)
