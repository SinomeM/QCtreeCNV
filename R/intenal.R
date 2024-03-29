
## TODO!
# 3. Add micro documentation in each function

get_region <- function(my_line, prop = 1) {
  # given a line consisting of a single CNV, returns a vector constaing chr,
  # start, end , length * prop as integers
  chr <- as.integer(my_line$chr)
  st <- as.integer(my_line$start)
  en <- as.integer(my_line$end)
  len <- (en - st +1) * prop
  reg <- c(chr, st, en, len)

  return(reg)
}

# -----------------------------------------------------------------------------

get_region_with_CNVR_ID <- function(my_line, prop = 1) {
  # same as get_region but returns also r_ID if present, in that case reg is a
  # list of vectors, in this way reg[[1]] remain a numeric vector
  chr <- as.integer(my_line$chr)
  st <- as.integer(my_line$start)
  en <- as.integer(my_line$end)
  len <- (en - st +1) * prop
  reg <- c(chr, st, en, len)

  # attach also CNVR_ID if present
  if ("CNVR_ID" %in% colnames(my_line))
    reg <- list(reg, my_line$CNVR_ID)

  return(reg)
}

# -----------------------------------------------------------------------------

get_regions_list <- function(my_lines, prop = 1) {
  # same as get_region but with multiple lines, returns a list of vectors
  chr <- as.integer(my_lines$chr)
  st <- as.integer(my_lines$start)
  en <- as.integer(my_lines$end)
  len <- (en - st +1) * prop

  # list of list when CNVR_ID or cnvrs information are present (character/numeric)
  if ("CNVR_ID" %in% colnames(my_lines))
    reg <- list(chr, st, en, len, my_lines$CNVR_ID)
  else if ("cnvr" %in% colnames(my_lines))
      reg <- list(chr, st, en, len, my_lines$cnvr)
  else
    reg <- list(chr, st, en, len)

  return(reg)
}

# -----------------------------------------------------------------------------

check_overlap <- function(cnvs, my_reg, prop) {
  # search reciprocal overlap between "my_reg" and any entry in "cnvs", if found
  # returns 1, 0 otherwise.
  res <- 0
  for (n in 1:nrow(cnvs)) {
    tmp_reg <- get_region(cnvs[n], prop)
    overl <- min(my_reg[3], tmp_reg[3]) - max(my_reg[2], tmp_reg[2]) + 1
    if (overl >= my_reg[4] & overl >= tmp_reg[4]) {
      res <- 1
      break # unnecessary?
    }
  }
  return(res)
}

# -----------------------------------------------------------------------------

chr_uniform <- function(DT_in) {

  # standardise col name, within certain limits
  if ("Chr" %in% colnames(DT_in)) setnames(DT_in, "Chr", "chr")
  if ("Chromosome" %in% colnames(DT_in))
    setnames(DT_in, "Chromosome", "chr")
  if ("chromosome" %in% colnames(DT_in))
    setnames(DT_in, "chromosome", "chr")

  if (!"chr" %in% colnames(DT_in))
    stop("No 'chr' columns found!\n")
  if (!is.data.table(DT_in))
    stop("Input must be a data.table1\n")

  DT_in[, chr := tolower(gsub(" ", "", chr))]
  if (substr(DT_in$chr[1], 1, 3) == "chr")
    DT_in[, chr := substring(chr, 4)]

  DT_in[chr == "x", chr := "23"][chr == "y", chr := "24"]

  # drop calls not in chrs 1:22, X, Y
  DT_in <- DT_in[chr %in% as.character(1:24), ]
}

# -----------------------------------------------------------------------------

uniform_GT_CN <- function(DT_in) {
  if (!any(c("CN", "type", "Type") %in% colnames(DT_in)))
    stop("No Copy Number column found!\n")
  if (!is.data.table(DT_in))
    stop("'DT_in' must be a data.table1\n")

  # standardise col name, within certain limits
  if ("type" %in% colnames(DT_in)) setnames(DT_in, "type", "CN")
  if ("Type" %in% colnames(DT_in)) setnames(DT_in, "Type", "CN")

  # standardise CN
  DT_in[, CN := tolower(gsub(" ", "", CN))]
  # +/-
  if ("+" %in% unique(DT_in$CN) | "-" %in% unique(DT_in$CN)) {
    # Autosomes
    DT_in[chr %in% as.character(1:22) & CN != "+" & CN != "-", CN := "2"][
      chr %in% as.character(1:22) & CN == "+", CN := "3"][
        chr %in% as.character(1:22) & CN == "-", CN := "1"]
  }
  # deletion/duplication, it also accept them with capital letters
  if ("deletion" %in% unique(DT_in$CN) | "duplication" %in% unique(DT_in$CN)) {
    # Autosomes
    DT_in[chr %in% as.character(1:22) & CN != "duplication" & CN != "deletion",
          CN := "2"][chr %in% as.character(1:22) & CN == "duplication", CN := "3"][
                     chr %in% as.character(1:22) & CN == "deletion", CN := "1"]
  }

  # compute GT
  # only chrs 1:22
  DT_in[chr %in% as.character(1:22) & CN == 2, GT := 0][
          chr %in% as.character(1:22) & CN < 2, GT := 1][
          chr %in% as.character(1:22) & CN > 2, GT := 2]
}


# ---------------------

# return sample_ID, locus, GT and CN as a vector of <chr>
# data.table is smart enough to convert the values back to
# <int> when needed (GT, CN, st/en etc.)

# BE CAREFUL, not all functions are smart enough to do this!!!
# change it tot a list ????

getline_cnv <- function(cnv) {
  if (!"CNVR_ID" %in% colnames(cnv)) 

  st <- cnv$start
  en <- cnv$end
  if ("length" %in% colnames(cnv)) l <- cnv$length
  else l <- en - st + 1
  cen <- st + l/2
  if (!"CNVR_ID" %in% colnames(cnv))
    #           1              2          3       4       5            6   7   8  9
    myline <- c(cnv$sample_ID, cnv$locus, cnv$GT, cnv$CN, cnv$CNVR_ID, st, en, l, cen)
  else
    myline <- c(cnv$sample_ID, cnv$locus, cnv$GT, cnv$CN, st, en, l, cen)
  return(myline)
}


# return start end and freq

getline_cnvr <- function(cnvr) {
  st <- cnvr$start
  en <- cnvr$end
  if ("length" %in% colnames(cnvr)) l <- cnvr$length
  else l <- en - st + 1
  cen <- st + l/2
  #           1   2   3          4         5             6  7
  myline <- c(st, en, cnvr$freq, cnvr$chr, cnvr$CNVR_ID, l, cen)
  return(myline)
}

# same but for loci
getline_locus <- function(locus) {
  st <- locus$start
  en <- locus$end
  if ("length" %in% colnames(locus)) l <- locus$length
  else l <- en - st + 1
  cen <- st + l/2
  #           1            2          3   4   5  6
  myline <- c(locus$locus, locus$chr, st, en, l, cen)
  return(myline)
}

getline_locus2 <- function(locus) {
  st <- locus$start
  en <- locus$end
  if ("length" %in% colnames(locus)) l <- locus$length
  else l <- en - st + 1
  cen <- st + l/2
  #                      1          2   3   4  5
  myline <- as.integer(c(locus$chr, st, en, l, cen))
  return(myline)
}
