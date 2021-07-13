
#' Extratc Raw Data in Specific Loci
#'
#' \code{extractRDS()} extract the raw LRR and BAF values in the loci from
#' the intensity files.
#'
#' @param loci lorem ipsum
#' @param samples_file lorem ipsum
#' @param rds_path lorem ipsum
#'
#' @export
#'
#' @import data.table


## TODO!
# 1. minimal documentation, in particular on the input files!
# 2. Add a progress bar?

# tar.gz are supported by default if R.utils is installed

extractRDS <- function(loci, samples_file, rds_path) {

  message("Reading final Reports")

  for (s in unique(samples_file$sample_ID)) {

    ff <- samples_file[sample_ID == s, file_path][1]

    if (length(ff) == 0) stop("Can't find sample ", s,
                              " in samples list file")

    tmp <- fread(ff, skip = "Sample ID")

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
      tmp1 <- tmp[Chr == cc & between(Position, st, sp), ]
      dt <- rbind(dt,tmp1)
    }

    dt <- dt[, .(Chr, Position, `Log R Ratio`, `B Allele Freq`)]

    saveRDS(dt, paste0(rds_path, s, ".rds"))
  }

message("Done!")

}
