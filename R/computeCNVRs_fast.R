#' Compute Copy Number Variable Regions (CNVRs)
#'
#' Faster CNVRs computation algorithm
#'
#' It is based on pamk clustering algorithm. It is a much simpler process that
#' takes advantage of the fact that this package is meant to deal only with groups
#' of calls overlapping a specific locus.
#'
#' @param cnvs lorem ipsum
#' @param chr_arms lorem ipsum
#' @param final_merge lorem ipsum
#'
#' @return a \code{list} of two elements. The first element is a \code{data.table}
#'   that contains the actual CNVR information, genomic location and frequency in
#'   the cohort. The second element is the \code{CNVresults}
#'
#' @export
#'
#' @import data.table


cnvr_fast <- function(put_cnvs) {

  cnvs_all <- data.table()
  cnvrs_all <- data.table()

  for (l in unique(put_cnvs$locus)) {
    cnvs <- put_cnvs[locus == l, ]
    cnvs[, `:=` (start = as.integer(start), end = as.integer(end),
                 length = as.integer(end - start + 1),
                 center = as.integer(start + length/2))]

    # cluster CNVs
    tmp <- data.frame(X = scale(cnvs$center), Y = scale(cnvs$length))
    # better than nothing
    if (sum(is.na(tmp)) > 0)
      tmp <- data.frame(X = cnvs$center, Y = cnvs$length)

    if (nrow(tmp) > 3) {
      maxk <- nrow(tmp)-1
      rr <- c(1, seq(2, 54, by = 4))
      r1 <- rr[1:8]

      # it's hard to imagine more than 24 CNVRs per window
      pam_cl <- fpc::pamk(tmp, krange = r1[r1 < maxk])
      nc <- pam_cl$nc
      # just in case
      if (nc == 24) pam_cl <- fpc::pamk(tmp, krange = rr[9:15])
      # fine tune the actual number of clusters
      pam_cl <- fpc::pamk(tmp, krange = nc-3:nc+3)
      cl <- as.integer(pam_cl$pamobject$clustering)
    }
    else {
      cl <- rep(1, nrow(tmp))
    }

    # assign each CNV to his cluster
    cnvs[, cnvr := cl]
    # create CNVRs table
    cnvrs <- data.table()
    # cluster 0 are the outliers, skip them for now
    for (i in 1:max(cl)) {
      # smooth the boundaries, median is less sensible to extreme values
      st <- as.integer(median(cnvs[cnvr == i, start]))
      en <- as.integer(median(cnvs[cnvr == i, end]))
      cnvrs <- rbind(cnvrs, data.table(cnvr = i, chr = cnvs[1, chr] ,start = st, end = en, freq = nrow(cnvs[cnvr == i, ])))
    }

    cnvs[, cnvr := paste0(l, "_", cnvr)]
    cnvrs[, cnvr := paste0(l, "_", cnvr)]

    cnvs_all <- rbind(cnvs_all, cnvs)
    cnvrs_all <- rbind(cnvrs_all, cnvrs)
  }

  return(list(cnvs = cnvs_all, cnvrs = cnvrs_all))

}



