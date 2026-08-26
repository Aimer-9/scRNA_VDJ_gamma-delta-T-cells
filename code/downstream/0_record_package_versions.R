# R version 4.5.2 (2025-10-31)
rm(list = ls())

project_dir <- "/path/to/project"
if (dir.exists(project_dir)) {
  setwd(project_dir)
}

code_dir <- file.path("code", "downstream")
output_path <- "renv.lock"
script_name <- "0_record_package_versions.R"

script_paths <- list.files(
  code_dir,
  pattern = "\\.[Rr]$",
  full.names = TRUE,
  recursive = TRUE
)
script_paths <- script_paths[basename(script_paths) != script_name]

strip_comments <- function(lines) {
  sub("#.*$", "", lines)
}

extract_function_packages <- function(lines, fun_name) {
  pattern <- paste0("\\b", fun_name, "\\s*\\(\\s*['\"]?([A-Za-z][A-Za-z0-9._]*)['\"]?")
  matches <- gregexpr(pattern, lines, perl = TRUE)
  found <- unlist(regmatches(lines, matches), use.names = FALSE)
  sub(pattern, "\\1", found, perl = TRUE)
}

extract_namespace_packages <- function(lines) {
  matches <- gregexpr("\\b([A-Za-z][A-Za-z0-9.]*)\\s*:::{0,1}\\s*[A-Za-z.][A-Za-z0-9._]*", lines, perl = TRUE)
  found <- unlist(regmatches(lines, matches), use.names = FALSE)
  sub("\\s*:::{0,1}.*$", "", found, perl = TRUE)
}

read_script_packages <- function(path) {
  lines <- strip_comments(readLines(path, warn = FALSE))
  packages <- unique(c(
    extract_function_packages(lines, "library"),
    extract_function_packages(lines, "require"),
    extract_namespace_packages(lines)
  ))
  packages <- packages[nzchar(packages)]
  setdiff(packages, c("FALSE", "TRUE", "NULL", "NA"))
}

json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  x <- gsub("\n", "\\\\n", x)
  x
}

json_string <- function(x) {
  paste0("\"", json_escape(x), "\"")
}

write_package_record <- function(package, installed_packages) {
  record <- installed_packages[package, , drop = FALSE]
  fields <- c(
    paste0("      \"Package\": ", json_string(package)),
    paste0("      \"Version\": ", json_string(record[, "Version"]))
  )
  repository <- if ("Repository" %in% colnames(record)) record[, "Repository"] else NA_character_
  if (!is.na(repository) && nzchar(repository)) {
    fields <- c(fields, paste0("      \"Repository\": ", json_string(repository)))
  }
  paste0(
    "    \"", json_escape(package), "\": {\n",
    paste(fields, collapse = ",\n"),
    "\n    }"
  )
}

used_packages <- unique(unlist(lapply(script_paths, read_script_packages), use.names = FALSE))
used_packages <- sort(used_packages[nzchar(used_packages)])

installed_packages <- utils::installed.packages()
installed_used_packages <- intersect(used_packages, rownames(installed_packages))
missing_packages <- setdiff(used_packages, rownames(installed_packages))

package_records <- vapply(
  installed_used_packages,
  write_package_record,
  character(1),
  installed_packages = installed_packages
)

lockfile <- c(
  "{",
  "  \"R\": {",
  paste0("    \"Version\": ", json_string(getRversion())),
  "  },",
  "  \"Packages\": {",
  paste(package_records, collapse = ",\n"),
  "  }",
  "}"
)

writeLines(lockfile, output_path)

message("Wrote renv-style package version lockfile: ", normalizePath(output_path, mustWork = FALSE))
if (length(missing_packages) > 0) {
  message("Packages used by scripts but not installed here: ", paste(missing_packages, collapse = ", "))
}
