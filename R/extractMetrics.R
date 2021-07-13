
#' Extratc QC metrics
#'
#' \code{extractMetrics()} extract the QC metrics from the raw data in RDS
#' format.
#'
#' @param loci lorem ipsum
#' @param cnvs lorem ipsum
#' @param pennQC lorem ipsum
#' @param int_rds_path lorem ipsum
#' @param tmp_rds_path lorem ipsum
#'
#' @export
#'
#' @import data.table


## TODO!
# 1. move here the code from the RMD file from GDK and standardize it
# 2. minimal documentation, in particular on the input files!

extractMetrics <- function(loci, cnvs, pennQC, int_rds_path, tmp_rds_path) {

  # initial checks
  # # TODO

  ids <- unique(cnvs$sample_ID)

  for (l in 1:nrow(loci)) {
    # info on locus
    loc <- getline_locus(loci[i])
    dt <- data.table(sample_ID = ids)
    message("Locus #", l, ": ", ll)

    for (s in ids) {
      tmp <-
        readRDS(file.path(int_rds_path,s,".rds"))[Chr == loc[2] &
                                                  between(Position, loc[3], loc[4]),]

      bafc <- nrow(tmp[between(`B Allele Freq`, 0.4, 0.6, incbounds=T), ]) /
                nrow(tmp)
      bafb <- nrow(tmp[between(`B Allele Freq`, 0.2, 0.4, incbounds=F) |
                       between(`B Allele Freq`, 0.6, 0.8, incbounds=F) |
                       `B Allele Freq` %in% c(0.2, 0.8), ]) / nrow(tmp)

      # sample QC measure from PennCNV
      qcline <- pennQC[sample_ID == s, ]
      dt[sample_ID == s, `:=` (BAFdrift = qcline$BAFdrift,
                               LRRSD = qcline$LRRSD,
                               GCWF = qcline$GCWF)]

      # locus measures (compute them regardless of the presence of a call)
      dt[sample_ID == s, `:=` (mLRRlocus = mean(tmp[, `Log R Ratio`], na.rm=T),
                         LRRSDlocus = sd(tmp[, `Log R Ratio`], na.rm=T),
                         BAFc = bafc, BAFb = bafb)]

      # check if this sample has a call in the locus
      put <- cnvs[Locus == ll & sample_ID == s,]
      if (nrow(put) == 0) {
        dt[sample_ID == s, `:=` (putCarrier= F, mLRRcall = NA_real_,
                                 centDistProp = NA_real_, overlapProp = NA_real_)]
      } else {
        # A sample can have only one call per locus
        if (nrow(put) > 1) stop(paste0("Sample ", s, " has more than one call in",
                                       " locus ", loc[1]))
        # info on putative call, if present
        putline <- getline_cnv(put)
        tmp1 <- tmp[between(Position, putline[6], putline[7]),]
        ov <- min(loc[4], putline[7]) - max(loc[3], putline[6]) +1
        if (ov < 0) stop("Something is wrong, overlap can't be negative)")
        dt[sample_ID == s, `:=` (putCarrier =T,
                                 mLRRcall = mean(tmp1[, `Log R Ratio`], na.rm=T),
                                 centDistProp = abs(loc[6] - putline[9]) / loc[5],
                                 overlapProp = putline[8] / loc[5])]
      }

      dt[,logr1 := log(abs(mLRRcall / mLRRlocus))]

      # rind all samples and loci together
      dtOUT <- rbind(dtOUT, dt)

    } # end samples loop
  } # end loci loop

  # sort columns
  dtOUT <- dtOUT[ , .(sample_ID, locus, putCarrier, LRRSD, BAFdrift, GCWF,
                      logr1, LRRSDlocus, BAFc, BAFb, centDistProp, overlapProp)]
  return(dtOUT)
}

