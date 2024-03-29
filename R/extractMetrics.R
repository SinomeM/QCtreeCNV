
#' Extratc QC metrics
#'
#' \code{extractMetrics()} extract the QC metrics from the raw data in RDS
#' format.
#'
#' @param loci lorem ipsum
#' @param cnvs lorem ipsum
#' @param pennQC lorem ipsum
#' @param samples_list lorem ipsum
#' @param snppos lorem ipsum
#'
#' @export
#'
#' @import data.table


## TODO!
# 1. minimal documentation, in particular on the input files!

extractMetrics <- function(loci, cnvs, pennQC, samples_list, snppos = NA) {

  # initial checks
  # # TODO

  ids <- unique(cnvs$sample_ID)
  dtOUT <- data.table()

  for (l in 1:nrow(loci)) {
    # info on locus
    loc <- getline_locus(loci[l])
    dt <- data.table(sample_ID = ids, locus = loc[1])
    message("Locus #", l, ": ", loc[1])

    for (s in ids) {

      f_path <- samples_list[sample_ID == s, file_path_tabix]
      if (length(f_path) > 1) {
            warning("sample ", s, " has more than one entry in samples_file")
            f_path <- f_path[1]
      }

      if (is.na(snppos)) tmp <- get_region_tabix(loc[2], loc[3], loc[4], f_path)
      else tmp <- get_region_tabix(loc[2], loc[3], loc[4], f_path, snppos)

      bafc <- nrow(tmp[between(BAF, 0.4, 0.6, incbounds=T), ]) /
                nrow(tmp)
      bafb <- nrow(tmp[between(BAF, 0.2, 0.4, incbounds=F) |
                       between(BAF, 0.6, 0.8, incbounds=F) |
                       BAF %in% c(0.2, 0.8), ]) / nrow(tmp)

      # sample QC measure from PennCNV
      qcline <- pennQC[sample_ID == s, ][1] # temporary fix!!!!
      dt[sample_ID == s, `:=` (BAFdrift = qcline$BAFdrift,
                               LRRSD = qcline$LRRSD,
                               GCWF = qcline$GCWF)]

      # locus measures (compute them regardless of the presence of a call)
      dt[sample_ID == s, `:=` (mLRRlocus = mean(tmp[, LRR], na.rm=T),
                         LRRSDlocus = sd(tmp[, LRR], na.rm=T),
                         BAFc = bafc, BAFb = bafb)]

      # check if this sample has a call in the locus
      put <- cnvs[locus == loc[1] & sample_ID == s,]
      if (nrow(put) == 0) {
        dt[sample_ID == s, `:=` (putCarrier= F, mLRRcall = NA_real_,
                                 centDistProp = NA_real_)]
      } else {
        # A sample can have only one call per locus
        if (nrow(put) > 1) stop(paste0("Sample ", s, " has more than one call in",
                                       " locus ", loc[1]))
        # info on putative call, if present
        putline <- getline_cnv(put)
        tmp1 <- tmp[between(position, as.integer(putline[5]), as.integer(putline[6])),]
        ov <- min(as.integer(loc[4]), as.integer(putline[6])) - max(as.integer(loc[3]), as.integer(putline[5])) +1
        if (ov < 0) stop("Something is wrong, overlap can't be negative)")
        dt[sample_ID == s, `:=` (putCarrier =T,
                                 mLRRcall = mean(tmp1[, LRR], na.rm=T),
                                 centDistProp = abs(as.integer(loc[6]) - as.integer(putline[8])) / as.integer(loc[5]))]
      }

      dt[,logr1 := log(abs(mLRRcall / mLRRlocus))]

    } # end samples loop
    # rind all samples and loci together
    dtOUT <- rbind(dtOUT, dt)
  } # end loci loop

  # sort columns
  dtOUT <- dtOUT[ , .(sample_ID, locus, putCarrier, LRRSD, BAFdrift, GCWF,
                      logr1, LRRSDlocus, BAFc, BAFb, centDistProp, mLRRlocus)]

  # add qs measures, i.e. merge the two tables
  setkeyv(cnvs, c("sample_ID", "locus"))
  setkeyv(dtOUT, c("sample_ID", "locus"))
  # join tables using data.table keys
  cnvsOUT <- dtOUT[cnvs]
  # add excl column to keep track of the process
  cnvsOUT[, excl := -1]

  return(cnvsOUT)
}

