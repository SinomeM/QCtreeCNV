
#' Extratc Raw Data in Specific Loci
#'
#' \code{extractRDS()} extract the raw LRR and BAF values in the loci from
#' the intensity files.
#'
#' @param loci lorem ipsum
#' @param samples_file lorem ipsum
#' @param rds_path lorem ipsum
#' @param overwrite lorem ipsum
#'
#' @export
#'
#' @import data.table


## TODO!
# 1. minimal documentation, in particular on the input files!
# 2. Add a progress bar?

# tar.gz are supported by default if R.utils is installed

extractRDS <- function(loci, samples_file, rds_path, overwrite=F) {

  message("Reading final Reports")

  for (s in unique(samples_file$sample_ID)) {

    if (overwrite == F & file.exists(paste0(rds_path, s, ".rds"))) next

    ff <- samples_file[sample_ID == s, file_path][1]

    if (length(ff) == 0) stop("Can't find sample ", s,
                              " in samples list file")

    tmp <- fread(ff, skip = "Sample ID")[,
             .(Chr, Position, `Log R Ratio`, `B Allele Freq`)]

    # select only the points within each locus
    dt <- data.table()
    for (l in 1:nrow(loci)) {
      cc <- loci$chr[l]
      if (cc %in% c(23,24)) {
        warning("Chr X and Y are not supported!")
        next
      }
      st <- loci$start[l]
      sp <- loci$end[l]
      dt <- rbind(dt, tmp[Chr == cc & between(Position, st, sp), ])
    }
    saveRDS(dt, paste0(rds_path, s, ".rds"))
    rm(tmp); rm(dt)
  }
message("Done!")
}
