## code to prepare `SGD_features` dataset goes here
SGD_features_file <- system.file(
  "extdata",
  "SGD_features.txt",
  package = "isomutR"
)
## load data and create GRanges object
SGD_features_df <- read.table(SGD_features_file, header = TRUE, sep = "\t")
SGD_features_df <- SGD_features_df[!SGD_features_df$chr %in% c("", "2-micron"),]
SGD_features_df <- SGD_features_df[!is.na(SGD_features_df$start),]
SGD_features_df <- SGD_features_df[!is.na(SGD_features_df$end),]
SGD_chr <- as.character(as.roman(SGD_features_df$chr))
SGD_start <- pmin(SGD_features_df$start, SGD_features_df$end)
SGD_end <- pmax(SGD_features_df$start, SGD_features_df$end)
SGD_features <- GenomicRanges::GRanges(SGD_chr,
                                       IRanges::IRanges(SGD_start, SGD_end),
                                       )
## add metadata
GenomicRanges::values(SGD_features) <- SGD_features_df[,1:7]

## filter for features to be removed
SGD_features <- SGD_features[
  grepl('centromere', SGD_features$FeatureType) |
  grepl('telomere', SGD_features$FeatureType) |
  grepl('telomeric',SGD_features$FeatureType) |
  grepl('LTR', SGD_features$FeatureType) |
  grepl('repeat', SGD_features$FeatureType)
]
usethis::use_data(SGD_features, overwrite = TRUE)
