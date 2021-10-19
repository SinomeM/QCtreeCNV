#' QCtree pipeline in a convenient and automated form
#'
#' All steps from the PennCNV raw data to the objects required by
#' `qctree()`
#'
#' Please note that at the moment having the input in the required format
#' is still something the user needs to take care of. Only defaults values
#' for each step are supported at the moment (will change soon)


qctree_pipeline <- function(loci, calls, pennqc, samples_list, rds_path,
                            hg_version = c("gh18", "hg19", "hg38")) {
  # check column names for all main objects
  if (!all(c("locus", "chr", "start", "end") %in% colnames(loci)))
    stop("Some columns are missing from loci object")

  if (!all(c() %in% colnames(calls)))
    stop("Some columns are missing from CNV calls object")

  if (!all(c() %in% colnames(pennqc)))
    stop("Some columns are missing from PennCNV QC object")

  if (!all(c() %in% colnames(samples_list)))
    stop("Some columns are missing from samples list object")

  # error if samples are not unique in the QC file
  if (length(pennqc$sample_ID) != length(unique(pennqc$sample_ID)) )
    stop("At least one sample has more than one line in the QC file")

  # warning if calls form samples not in the QC are present
  if (!all(calls$sample_ID %in% pennqc$sample_ID))
    warning("Some calls are from samples not present in the QC object")

  # warning if calls not in the samplee list are present
  if (!all(calls$sample_ID %in% samples_list$sample_ID))
    warning("Some calls are from samples not present in the samples list object")

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
  if (hg_version == "hg18") arms <- [...]
  if (hg_version == "hg19") arms <- [...]
  if (hg_version == "hg38") arms <- [...]
  cnvrs <- cnvrs_create(put_cnvs, arms)
  put_cnvs <- cnvrs[[2]]

  # extract RDS (very long step because of I/O)
  extractRDS(loci, samples_list, rds_path)

  # create final QC table
  ...

  # if there is more than one call per locus in a sample, keep the largest one
  duprm <- data.table()
  for (l %in% unique(put_cnvs[, locus])) {
    a <- put_cnvs[locus == l, ]
    b <- a[sample_ID %in% a$sample_ID[duplicated(a$sample_ID)], ]
    for (s in unique(b$sample_ID)) {
      lk <- max(b[sample_ID == s, length])
      duprm <- rbind(duprm, b[sample_ID == s & length != lk, ])
    }
  }
  put_cnvs <- fsetdiff(put_cnvs, duprm)



}
