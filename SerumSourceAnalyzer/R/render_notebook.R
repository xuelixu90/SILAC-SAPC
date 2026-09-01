# ============================================================
# Notebook rendering wrapper
# ============================================================

#' Locate the analysis notebook template
#'
#' Looks first at the installed package `inst/rmarkdown/` directory and
#' falls back to the source-tree equivalent (handy when developing with
#' devtools or loading the code without installing).
#'
#' @return path to `analysis_notebook.Rmd` (invisibly)
#' @keywords internal
find_notebook_rmd <- function() {
  candidates <- c(
    system.file("rmarkdown", "analysis_notebook.Rmd",
                package = "SerumSourceAnalyzer"),
    file.path(getwd(), "inst", "rmarkdown", "analysis_notebook.Rmd"),
    file.path(dirname(getwd()), "inst", "rmarkdown", "analysis_notebook.Rmd")
  )
  found <- candidates[file.exists(candidates) & nchar(candidates) > 0]
  if (length(found) == 0) {
    stop("Cannot find analysis_notebook.Rmd template. Looked in:\n  ",
         paste(candidates, collapse = "\n  "))
  }
  normalizePath(found[1], winslash = "/", mustWork = TRUE)
}

#' Render the execution notebook (HTML + Markdown)
#'
#' Runs the bundled `analysis_notebook.Rmd` against the real parameters of
#' the current analysis so the resulting files contain the full log: code,
#' console messages, intermediate tables, and every plot.
#'
#' The Markdown version is rendered so users can inspect the run without a
#' browser; the HTML version is a full self-contained notebook that can be
#' opened directly.
#'
#' @param input_file workbook path used by the run
#' @param output_dir destination directory for the generated files
#' @param sample_info sample mapping table
#' @param params a [serum_source_params()] object
#' @return a named list with two paths: `html` and `markdown`
#' @export
render_analysis_notebook <- function(input_file,
                                     output_dir,
                                     sample_info = build_sample_info(),
                                     params = serum_source_params()) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("Package 'rmarkdown' is not installed. Skipping notebook render. ",
            "Install with:\n  install.packages(c('rmarkdown','knitr'))")
    return(NULL)
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    warning("Package 'knitr' is not installed. Skipping notebook render.")
    return(NULL)
  }

  # RStudio bundles pandoc under resources/app/bin/quarto/bin/tools/ or
  # resources/app/bin/quarto/bin/.  Non-interactive Rscript sessions do
  # not automatically set RSTUDIO_PANDOC, so rmarkdown reports pandoc as
  # missing.  Probe common install locations and export RSTUDIO_PANDOC
  # to point at the pandoc binary directory before rendering.
  probe_pandoc <- function() {
    if (rmarkdown::pandoc_available()) return(TRUE)
    candidates <- c(
      "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/pandoc.exe",
      "C:/Program Files/RStudio/resources/app/bin/quarto/bin/pandoc.exe",
      "C:/Program Files (x86)/RStudio/resources/app/bin/quarto/bin/tools/pandoc.exe",
      "C:/Program Files/Pandoc/pandoc.exe",
      file.path(Sys.getenv("LOCALAPPDATA"), "Pandoc", "pandoc.exe", fsep = "/"),
      file.path(Sys.getenv("ProgramFiles"),     "Pandoc", "pandoc.exe", fsep = "/")
    )
    found <- candidates[file.exists(candidates)]
    if (length(found) == 0) return(FALSE)
    pandoc_dir <- dirname(normalizePath(found[1], winslash = "/", mustWork = TRUE))
    Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
    # Force rmarkdown to re-scan pandoc location
    tryCatch({
      unlink(Sys.getenv("RMARKDOWN_PANDOC"), recursive = TRUE)
    }, error = function(e) NULL)
    rmarkdown::pandoc_available()
  }
  have_pandoc <- probe_pandoc()
  if (!have_pandoc) {
    warning("pandoc >= 1.12.3 is required but not found. HTML notebook ",
            "cannot be rendered. Install RStudio or Pandoc (https://pandoc.org).")
    return(NULL)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rmd_src    <- find_notebook_rmd()
  rmd_copy   <- file.path(output_dir, "Analysis_Notebook.Rmd")
  file.copy(rmd_src, rmd_copy, overwrite = TRUE)

  # Detect the library the package was loaded from (if any).  This ensures
  # that the freshly rendered notebook loads the same build of the package
  # rather than whatever happens to be on .libPaths().
  lib_dir <- if (requireNamespace("SerumSourceAnalyzer", quietly = TRUE)) {
    dirname(system.file(package = "SerumSourceAnalyzer"))
  } else {
    .libPaths()[1]
  }

  params_list <- list(
    INPUT_FILE  = normalizePath(input_file, winslash = "/", mustWork = TRUE),
    OUTPUT_DIR  = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
    PARAMS      = params,
    SAMPLE_INFO = sample_info,
    LIB_DIR     = lib_dir
  )

  # (1) Render the html_notebook output -> *.nb.html
  message("Rendering HTML notebook ...")
  nb_html <- tryCatch({
    rmarkdown::render(
      input          = rmd_copy,
      output_format  = "html_notebook",
      output_file    = "Analysis_Notebook.nb.html",
      output_dir     = output_dir,
      params         = params_list,
      envir          = new.env(parent = globalenv()),
      quiet          = FALSE,
      knit_root_dir   = output_dir
    )
  }, error = function(e) {
    warning("html_notebook render failed: ", conditionMessage(e),
            "\nFalling back to html_document.")
    NULL
  })

  # (2) Render an html_document fallback and a plain md version.
  message("Rendering html_document (fallback + print-friendly) ...")
  reg_html <- tryCatch({
    rmarkdown::render(
      input          = rmd_copy,
      output_format  = "html_document",
      output_file    = "Analysis_Notebook.html",
      output_dir     = output_dir,
      params         = params_list,
      envir          = new.env(parent = globalenv()),
      quiet          = FALSE,
      knit_root_dir   = output_dir
    )
  }, error = function(e) {
    warning("html_document render failed: ", conditionMessage(e))
    NULL
  })

  message("Rendering Markdown notebook ...")
  md_file <- tryCatch({
    rmarkdown::render(
      input          = rmd_copy,
      output_format  = "md_document",
      output_file    = "Analysis_Notebook.md",
      output_dir     = output_dir,
      params         = params_list,
      envir          = new.env(parent = globalenv()),
      quiet          = FALSE,
      knit_root_dir   = output_dir
    )
  }, error = function(e) {
    warning("md_document render failed: ", conditionMessage(e))
    NULL
  })

  # Pick the "primary HTML" entry (prefer the notebook, fall back to the
  # regular html_document).
  html_out <- if (!is.null(nb_html) && file.exists(nb_html)) {
    nb_html
  } else if (!is.null(reg_html) && file.exists(reg_html)) {
    reg_html
  } else {
    NA_character_
  }

  out <- list(
    rmd      = rmd_copy,
    html     = if (is.na(html_out)) NULL else html_out,
    nb_html  = if (!is.null(nb_html) && file.exists(nb_html)) nb_html else NULL,
    html_doc = if (!is.null(reg_html) && file.exists(reg_html)) reg_html else NULL,
    markdown = if (!is.null(md_file)  && file.exists(md_file))  md_file  else NULL
  )
  class(out) <- "notebook_output"
  message("Notebook render finished.")
  message("  HTML (open in browser): ", if (!is.null(out$html)) out$html else "<failed>")
  message("  Markdown:              ", if (!is.null(out$markdown)) out$markdown else "<failed>")
  invisible(out)
}
