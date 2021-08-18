## code to prepare `txdb_scer` dataset goes here
txdb_scer <- system.file(
  "extdata",
  "Saccharomyces_cerevisiae.R64-1-1.102.gtf",
  package = "isomutR"
)
usethis::use_data(txdb_scer, overwrite = TRUE)
