setwd("/home/ubuntu/vol_store/differential_digestion")


render_fun <- function(rmd) {
  rmarkdown::render(
    rmd,
    output_dir = paste0("analysis/", basename(rmd), "/"),
    clean = FALSE
  )
}

render_fun("scripts/initial_investigation.Rmd")
render_fun("scripts/digestion_coverage.Rmd")
