
#' Extratc QC metrics
#'
#' \code{extractMetrics()} extract the QC metrics from the raw data in RDS
#' format.
#'
#' @param loci lorem ipsum
#' @param cnvs lorem ipsum
#' @param pennQC lorem ipsum
#' @param int_rds_path lorem ipsum
#'
#' @export
#'
#' @import data.table


## TODO!
# 1. minimal documentation, in particular on the input files!

extractMetrics <- function(loci, cnvs, pennQC, int_rds_path) {

  # initial checks
  # # TODO

  ids <- unique(cnvs$sample_ID)

  for (l in 1:nrow(loci)) {
    # info on locus
    loc <- getline_locus(loci[l])
    dt <- data.table(sample_ID = ids)
    message("Locus #", l, ": ", loc[1])

    for (s in ids) {
      tmp <- readRDS(paste0(int_rds_path,"/",s,".rds"))[
               Chr == loc[2] & between(Position, as.integer(loc[3]), as.integer(loc[4])),]

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
        tmp1 <- tmp[between(Position, as.integer(putline[6]), as.integer(putline[7])),]
        ov <- min(as.integer(loc[4]), as.integer(putline[7])) - max(as.integer(loc[3]), as.integer(putline[6])) +1
        if (ov < 0) stop("Something is wrong, overlap can't be negative)")
        dt[sample_ID == s, `:=` (putCarrier =T,
                                 mLRRcall = mean(tmp1[, `Log R Ratio`], na.rm=T),
                                 centDistProp = abs(as.integer(loc[6]) - as.integer(putline[9])) / as.integer(loc[5]),
                                 overlapProp = as.integer(putline[8]) / as.integer(loc[5]))]
      }

      dt[,logr1 := log(abs(mLRRcall / mLRRlocus))]

      # rind all samples and loci together
      dtOUT <- rbind(dtOUT, dt)

    } # end samples loop
  } # end loci loop

  # sort columns
  dtOUT <- dtOUT[ , .(sample_ID, locus, putCarrier, LRRSD, BAFdrift, GCWF,
                      logr1, LRRSDlocus, BAFc, BAFb, centDistProp, overlapProp, mLRRlocus)]
  return(dtOUT)
}

