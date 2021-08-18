library(Rsamtools)
ref_file <- system.file(
  "extdata",
  "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",
  package = "isomutR"
)
ref_scer <- FaFile(ref_file)
usethis::use_data(ref_scer, overwrite = TRUE)
