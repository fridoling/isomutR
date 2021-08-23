## code to prepare `gtf_scer` dataset goes here

gtf_scer <- system.file(
  "extdata",
  "Saccharomyces_cerevisiae.R64-1-1.102.gtf",
  package = "isomutR"
)
usethis::use_data(gtf_scer, overwrite = TRUE)
