#' Read and preprocess an isomut file
#'
#' `read_isomut` reads an isomut file, adds a `read` column, optionally
#' filters according to given criteria and adds extra columns.
#'
#' @param file the name of the file which the data are to be read from.
#' @param minReads integer indicating the minimum number of reads
#' (default = 0L).
#' @param minCoverage double indicating the minimum coverage
#' (default = 0).
#' @param minMutFreq double indicating the minimum mutation frequency
#' (default = 0).
#' @param removePatterns character vector providing patterns in
#' sample_names that are to be removed.
#' @param extraColumns named list to generate extra columns. The list names are the
#' names of the new columns. The list elements can either be a vector of length
#' one (each row gets the same value for this column), a named vector with names
#' being a subset of the `sample_name` column in the original isomut file, or a
#' function that accepts those names as arguments.
#'
#' @return a data.frame containing the filtered and preprocessed data.
#' @export
read_isomut <- function(file, minReads = 0L, minCoverage = 0, minMutFreq = 0,
                        removePatterns = NULL, extraColumns = NULL) {
  ## read data
  df <- read.table(file, header = TRUE)
  ## rename `sample_name` column
  df$file_name <- df$sample_name
  df$sample_name <- NULL
  ## remove samples that match certain patterns
  if(!is.null(remove_patterns)) {
    for(pattern in remove_patterns) {
      df <- df[!grepl(pattern, df$file_name),]
    }
  }
  ## add `read` column
  df[["read"]] <- with(df, round(coverage * mut_freq))
  ## remove rows with too few reads
  if(minReads > 0) {
    df <- df[df$read > minReads,]
  }
  ## remove rows with low coverage
  if(minCoverage > 0) {
    df <- df[df$coverage > minCoverage,]
  }
  ## remove rows with mutation frequency
  if(minMutFreq > 0) {
    df <- df[df$mut_freq > minMutFreq,]
  }
  ## replace cleanliness 42 with NA
  df$cleanliness[df$cleanliness==42] <- NA
  rownames(df) <- NULL
  ## add extra columns
  if(!is.null(extra_columns)) {
    stopifnot(is.list(extra_columns), !is.null(names(extra_columns)))
    for(col_name in names(extra_columns)) {
      col_value <- extra_columns[[col_name]]
      stopifnot(is.vector(col_value) || is.function(col_value))
      if(is.function(col_value)) {
        df[[col_name]] <- vapply(df$file_name, col_value, character(1))
      } else {
        if(length(col_value) == 1) {
          df[[col_name]] <- col_value
        } else {
          stopifnot(df$file_name %in% names(col_value))
          df[[col_name]] <- col_value[df$file_name]
        }
      }
    }
  }
  ## order
  df <- df[order(df$file_name, df$chr, df$pos),]
  return(df)
}
