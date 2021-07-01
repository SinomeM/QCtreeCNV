
## TODO!
# 1. UPDATE DATA FORMAT (e.g. end -> stop)!!!
# 2. UPDATE chruniform !!!!!
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

get_region_with_rID <- function(my_line, prop = 1) {
  # same as get_region but returns also r_ID if present, in that case reg is a
  # list of vectors, in this way reg[[1]] remain a numeric vector
  chr <- as.integer(my_line$chr)
  st <- as.integer(my_line$start)
  en <- as.integer(my_line$end)
  len <- (en - st +1) * prop
  reg <- c(chr, st, en, len)

  # attach also r_ID if present
  if ("r_ID" %in% colnames(my_line))
    reg <- list(reg, my_line$r_ID)

  return(reg)
}

# -----------------------------------------------------------------------------

get_regions_list <- function(my_lines, prop = 1) {
  # same as get_region but with multiple lines, returns a list of vectors
  chr <- as.integer(my_lines$chr)
  st <- as.integer(my_lines$start)
  en <- as.integer(my_lines$end)
  len <- (en - st +1) * prop

  # list of list when R_ID or cnvrs information are present (character/numeric)
  if ("r_ID" %in% colnames(my_lines))
    reg <- list(chr, st, en, len, my_lines$r_ID)
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
  if (!"chr" %in% colnames(DT_in))
    stop("No 'chr' columns found!\n")
  if (!is.data.table(DT_in))
    stop("'DT_in' must be a dat.table1\n")

  DT_in[, chr := tolower(gsub(" ", "", chr))]
  # this won't work if the the notation in the column is not coherent, like the
  # results of biomaRt::getBM()
  if (substr(DT_in$chr[1], 1, 3) == "chr")
    DT_in[, chr := substring(chr, 4)]
  # or sub(".+(\\d+|x|y)", "\\1", chr)
  DT_in[chr == "x", chr := "23"][chr == "y", chr := "24"]
  # drop calls not in chrs 1:22, X, Y
  DT_in <- DT_in[chr %in% as.character(1:24), ]

  return(DT_in)
}

# -----------------------------------------------------------------------------

uniform_GT_CN <- function(DT_in) {
  if (!any(c("CN", "type", "Type") %in% colnames(DT_in)))
    stop("No Copy Number column found!\n")
  if (!is.data.table(DT_in))
    stop("'DT_in' must be a data.table1\n")

  # standardise col name
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

  return(DT_in)
}
