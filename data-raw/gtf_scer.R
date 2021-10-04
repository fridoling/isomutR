## code to prepare `gtf_scer` dataset goes here

gtf_scer_path <- system.file(
  "extdata",
  "Saccharomyces_cerevisiae.R64-1-1.102.gtf",
  package = "isomutR"
)
gtf_scer <- rtracklayer::import(gtf_scer_path)
usethis::use_data(gtf_scer, overwrite = TRUE)
