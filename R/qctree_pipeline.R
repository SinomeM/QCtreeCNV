#' QCtree pipeline in a convenient and automated form
#'
#' All steps from the PennCNV raw data to the objects required by
#' `qctree()`
#'
#' Please note that at the moment having the input in the required format
#' is still something the user needs to take care of. Only defaults values
#' for each step are supported at the moment (will change soon)

# ideally some step will be parallelized
# also the idea is that cnvrs can be passed already computed, however
# in that case the CNV calls must also be already the equivalent of put_cnvs


qctree_pipeline <- function(loci, calls, pennqc, samples_list, rds_path, cnvrs = NA,
                            hg_version = c("gh18", "hg19", "hg38"), snppos = NA) {

  check_inputs(loci, calls, pennqc, samples_list)

  # check and correct chr & GT/CN
  loci <- chr_uniform(loci)
  calls <- chr_uniform(calls)
  calls <- uniform_GT_CN(calls)

  # add length
  loci[, length := end - start + 1]
  calls[, length := end - start + 1]

  # select and stitch calls in the required loci
  put_cnvs <- select_stitch_calls(calls, loci)

  # compute CNVRs (long step)
  if (is.na(cnvrs)) {
    #     if (hg_version == "hg18") arms <- [...]
    #     if (hg_version == "hg19") arms <- [...]
    #     if (hg_version == "hg38") arms <- [...]
    cnvrs <- cnvrs_create(put_cnvs, arms)
    put_cnvs <- cnvrs[[2]]
  }

###  # extract RDS (very long step because of I/O)
###  extractRDS(loci, samples_list, rds_path)

  # create final QC table
  qc <- extractMetrics(loci, put_cnvs, pennqc, samples_list, snppos)

  # if there is more than one call per locus in a sample, keep the largest one
  if (rm_dup) {
    duprm <- data.table()
    for (l in unique(put_cnvs[, locus])) {
      a <- put_cnvs[locus == l, ]
      b <- a[sample_ID %in% a$sample_ID[duplicated(a$sample_ID)], ]
      for (s in unique(b$sample_ID)) {
        lk <- max(b[sample_ID == s, length])
        duprm <- rbind(duprm, b[sample_ID == s & length != lk, ])
      }
    }
    put_cnvs <- fsetdiff(put_cnvs, duprm)
  }

  return(list(put_cnvs, cnvrs[[1]], qc, loci))
}

check_inputs <- function(loci, calls, pennqc, samples_list) {
  # check column names for all main objects
  if (!all(c("locus", "chr", "start", "end") %in% colnames(loci)))
    stop("Some required columns are missing from loci object")

  if (!all(c("sample_ID", "chr", "start", "end", "CN", "numsnp") %in% colnames(calls)))
    stop("Some required columns are missing from CNV calls object")

  if (!all(c("sample_ID", "BAFdrift", "LRRSD", "GCWF") %in% colnames(pennqc)))
    stop("Some required columns are missing from PennCNV QC object")

  if (!all(c("sample_ID", "file_path", "file_path_tabix") %in% colnames(samples_list)))
    stop("Some required columns are missing from samples list object")

  if(!is.na(cnvrs)) {
    if (!all(c() %in% colnames(cnvrs)))
      stop("Some required columns are missing from cnvrs object")
    if (!all(c("CNVR_ID", "locus") %in% colnames(calls))) # other????
      stop("cnvrs are passed but no CNVR_ID | locus column is found in the cnv object")
  }

  # error if samples are not unique in the QC file
  if (length(pennqc$sample_ID) != length(unique(pennqc$sample_ID)) )
    stop("At least one sample has more than one line in the QC file")

  # warning if calls form samples not in the QC are present
  if (!all(calls$sample_ID %in% pennqc$sample_ID))
    warning("Some calls are from samples not present in the QC object")

  # error if calls not in the samplee list are present
  if (!all(calls$sample_ID %in% samples_list$sample_ID))
    stop("Some calls are from samples not present in the samples list object")
}
