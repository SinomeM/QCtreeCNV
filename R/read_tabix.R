#' Load raw data from tabix indexed intensity file
#'
#' @import data.table

get_region_tabix <- function(chr, start, end, f_path, snppos = NA) {
  a <- fread(cmd = paste("tabix", f_path, chr, ":", start, "-", end, sep = " "), header = F)
  colnames(a) <- c("chr", "postion", "end", "LRR", "BAF", "LRRadj")

  # Turn NaN into standard NAs
  a[LRR == "NaN", LRR := NA][BAF == "NaN", BAF := NA]

  # filter SNPs if required
  if (!is.na(snppos)) a <- a[snp_name %in% snp_pos$snp_name, ]

  setorder(a, position)
  return(a)
}
