#!/usr/bin/env Rscript
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", all_args[grepl("^--file=", all_args)])
root <- normalizePath(file.path(dirname(script_arg[[1]]), "..", ".."), mustWork = TRUE)
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
if (!requireNamespace("targets", quietly = TRUE)) {
  stop("Package `targets` is required. Add it to the local R package source, run `bash scripts/setup.sh r`, then retry.", call. = FALSE)
}
target <- NULL
if (length(args)) {
  if (length(args) != 2L || args[[1]] != "--target") stop("Usage: Rscript code/pipeline/run.R [--target TARGET]", call. = FALSE)
  target <- args[[2]]
}
targets::tar_make(names = target)
