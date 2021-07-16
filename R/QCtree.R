#' CNVs QC Filtering Tree
#'
#' The main functions of the \code{QCtreeCNV} package
#'
#' @param cnvs lorem ipsum
#' @param cnvrs lorem ipsum
#' @param qsdt lorem ipsum
#' @param loci lorem ipsum
#' @param maxLRRSD lorem ipsum
#' @param maxBAFdrift lorem ipsum
#' @param maxGCWF lorem ipsum
#' @param minGCWF lorem ipsum
#' @param commonCNVRsMinFreq lorem ipsum
#' @param st4minlogr1 lorem ipsum
#' @param st4maxlogr1 lorem ipsum
#' @param st4maxBAFcDEL lorem ipsum
#' @param st4maxBAFcDUP lorem ipsum
#' @param st4maxBAFbDEL lorem ipsum
#' @param st5maxLRRSDlocus lorem ipsum
#' @param st5maxBAFcDEL lorem ipsum
#' @param st5maxBAFcDUP lorem ipsum
#' @param st5maxBAFbDEL lorem ipsum
#' @param st5maxBAFbDUP lorem ipsum
#' @param clean_out lorem ipsum
#'
#' @export
#'
#' @import data.table

## TODO!

# 1. Some documentation, keep in mind most of it will live in more discoursive
#    manner in the vignette.
# 2. Test and check w/ Andres default values.
# 3. Each step as a function for easy testing
# 4. Initial checks


# Each step create a columns 1/0 that mean the line has passed the check (YES/NO)
# the successive step will be performed only on the relevant lines using the previous
# columns for filtering.
# the initialized columns will have the value -1.
# NOTE that pass = YES (1) is not necessarily the "good" outcome

# USE THE SCHEME IN THE VIGNETTE TO FOLLOW THE FLOW

# remember BAFb should not be used for duplications

# function parameters are organized by row:
# - inputs objects
# - CNVRs
# - step 1 (QC)
# - step 4
qctree <- function(cnvs, cnvrs, qsdt, loci,
                   maxLRRSD=.3, maxBAFdrift=.01,maxGCWF=.015, minGCWF=-.015,
                   commonCNVRsMinFreq = NA,
                   st4minlogr1=-.35,st4maxlogr1=.4, st4maxBAFcDEL=.03,
                   st4maxBAFcDUP=.125, st4maxBAFbDEL=0.075,
                   st5maxLRRSDlocus=0.35, st5maxBAFcDEL=.05, st5maxBAFcDUP=.15,
                   st5maxBAFbDEL=.1, st5maxBAFbDUP=NA,
                   clean_out = T) {

  ### Pre process ###

  ### initial checks ###
  # # TODO!

  message("# -------------------------- #\n",
          "Step 0, pre-process")

  # add qs measures, i.e. merge the two tables
  setkey(cnvs, c("sample_ID", "locus"))
  setkey(qsdt, c("sample_ID", "locus"))
  # join tables using data.table keys
  cnvsOUT <- cnvs[qsdt]
  # add excl column to keep track of the process
  cnvsOUT[, excl := -1]


  ### STEP 1, QC outliers ###

  # a mild (in the default settings) filter on the three
  # classic measures, LRRSD, BAFdrift and GCWF
  message("# -------------------------- #\n",
          "Step 1, QC outliers removal")
  cnvsOUT <- step1(cnvsOUT, maxLRRSD, maxBAFdrift, minGCWF, maxGCWF)


  ### STEP 2 & 3 CNVRs ###
  message("# -------------------------- #\n",
          "Steps 2 and 3, CNVRs")
  cnvrs_groups <- sortCNVRs(cnvs, cnvrs, loci, commonCNVRsMinFreq)

  # STEP 2
  cnvsOUT <- step2(cnvsOUT, cnvrs_groups)

  # STEP 3
  cnvsOUT <- step3(cnvsOUT, cnvrs_groups)



  ### STEP 4, logr1, BAFc & BAFb checks 1 ###
  message("# -------------------------- #\n",
          "Step 4, logr1, BAFc & BAFb for small CNVRs calls")
  cnvsOUT<- step4(cnvsOUT, st4minlogr1, st4maxlogr1,
                  st4maxBAFcDEL, st4maxBAFcDUP, st4maxBAFbDEL)


  ### STEP 5, LRRSDlocus, mLRRlocus ###
  message("# -------------------------- #\n",
          "Step 5, LRRSD and mean LRR checks inside locus of interest\n",
          "plus logr1, BAFc & BAFb check for all calls")
  cnvsOUT <- step5(cnvsOUT, 0, -0.3, st5maxLRRSDlocus,
                   st5maxBAFcDEL, st5maxBAFcDUP, st5maxBAFbDEL,
                   st4maxBAFcDEL, st4maxBAFbDEL)

  ### RETURN ###
  message("# -------------------------- #\n",
          "Pipeline complete!")
  # check all CNVs were evaluated
  if (nrow(cnvsOUT[excl == -1,] == 0)) message("All CNVs were evaluated")

  # remove temporary columns
  if (clean_out) cnvsOUT[, c("st1", "st2", "st3", "st4", "st5") := NULL]

  # return good and excluded as a list
  return(list(cnvsOUT[excl == 0, ], cnvsOUT[excl == 1, ]))
}


step1 <- function(cnvs, mlrrsd, mbafd, mingc, maxgc) {
  cnvs[, st1 := -1]
  # these pass the check
  cnvs[LRRSD <= mlrrsd & BAFdrift <= mbafd &
       between(GCWF, mingc, maxgc,incbounds=T), st1 := 0]
  # these are excluded
  cnvs[st1 == -1, st1 := 1][st1 == 1, excl := 1]
  #   message("# ", nrow(cnvs[excl == 1, ]),
  #           "calls are from QC outlier samples.\n")
  return(cnvs)
}

step2 <- function(cnvs, cnvrs) {
  cnvs[, st2 := -1]
  # CNVs from step 1 == 0 that are in a cnvrA will pass step 2 (to good CNVs)
  cnvs[st1 == 0 & CNVR_ID %in% cnvrs[[1]], `:=` (st2 = 1, excl = 0)]
  # CNVs from step 1 == 0 that are in a cnvrA will fail step 2 (to step3)
  cnvs[st1 == 0 & !CNVR_ID %in% cnvrs[[1]], st2 := 0]

  # check all CNVs from step 1 are assigned
  if (!all(cnvs[st1 == 0, st2] %in% c(1,0)))
    stop("There is a problem in step 2")
  return(cnvs)
}


step3 <- function(cnvs, cnvrs) {
  cnvs[, st3 := -1]
  # CNVs from step 2 == 0 that are in a cnvrB will go to step 4
  # CNVs from step 2 == 0 that are in a cnvrC will go to step 5
  cnvs[st2 == 0 & CNVR_ID %in% cnvrs[[2]], st3 := 1]
  cnvs[st2 == 0 & CNVR_ID %in% cnvrs[[3]], st3 := 0]

  # check all CNVs from step 2 are assigned
  if (!all(cnvs[st2 == 0, st3] %in% c(1,0)))
    stop("There is a problem in step 3")
  return(cnvs)
}

step4 <- function(cnvs, minlogr1, maxlogr1, maxbafcdel, maxbafcdup, maxbafbdel) {
  cnvs[, st4 := -1]
  # Deletions and duplications are treated differently here
  # Takes calls from st3 == 1
  # If log1 or at least one between BAFc and BAFb are whithin limits, it goes
  # to step 5, i.e. st4 = 0
  # If logr1 AND BAFc | BAFb are out of limits (i.e. consistent with small CNVR)
  # exclude them
  cnvs[st3 == 1 & GT == 1 & (between(logr1, minlogr1, maxlogr1) |
          (BAFc <= maxbafcdel | BAFb <= maxbafbdel)), st4 := 0]
  cnvs[st3 == 1 & GT == 1 & st4 == -1, `:=` (st4 = 1, excl = 1)]

  cnvs[st3 == 1 & GT == 2 & (between(logr1, minlogr1, maxlogr1) |
          BAFc <= maxbafcdup), st4 := 0]
  cnvs[st3 == 1 & GT == 2 & st4 == -1, `:=` (st4 = 1, excl = 1)]

  # check all CNVs from step 3 are assigned
  if (nrow(cnvs[st3 == 0,] != nrow(cnvs[st4 %in% c(1,0), ])))
    stop("There is a problem in step 4")
  return(cnvs)
}

step5 <- function(cnvs, maxmLRRdel, minmLRRdup, maxlrrsd,
                  maxbafcdel, maxbafcdup, maxbafbdel,
                  maxbafcdel2, maxbafbdel2) {
  cnvs[, st5 := -1]
  # Takes calls from st3 = 0 and st4 = 0

  # mLRRlocus particularly out
  cnvs[GT == 1 & st3 == 0 & st4 == 0 & mLRRlocus > maxmLRRdel, `:=` (st5 = 1, excl = 1)]
  cnvs[GT == 2 & st3 == 0 & st4 == 0 & mLRRlocus < minmLRRdup, `:=` (st5 = 1, excl = 1)]

  # In calls with high LRRSDlocus, check BAFc and BAFb, if at least one of them
  # well out of range or both with lower thresholds then it's excluded
  # Deletions
  cnvs[GT == 1 & st3 == 0 & st4 == 0 & LRRSDlocus > maxlrrsd &
          ((BAFc <= maxbafcdel | BAFb <= maxbafbdel) |
           (BAFc <= maxbafcdel2 & BAFb <= maxbafbdel2)), `:=` (st5 = 1, excl = 1)]
  # Duplications
  cnvs[GT == 2 & st3 == 0 & st4 == 0 & LRRSDlocus > maxlrrsd &
          BAFc <= maxbafcdup, `:=` (st5 = 1, excl = 1)]

  # LRRSDlocus extremely high
  cnvs[st3 == 0 & st4 == 0 & LRRSDlocus > 0.55, `:=` (st5 = 1, excl = 1)]

  # all the other
  cnvs[st3 == 0 & st4 == 0 & st5 == -1, `:=` (st5 = 0, excl = 0)]

  # check all CNVs from step 3 are assigned
  if (!all(cnvs[st3 == 0 & st4 == 0, st5] %in% c(1,0)))
    stop("There is a problem in step 5")
  return(cnvs)
}


sortCNVRs <- function(cnvs, loci, cnvrs, minFreq) {

  message("Checking CNVRs frequency and overlaps with the loci")
  # Sort CNVRs in the different categories, done per locus
  # Three categories:
  # - A: rare or large, alls CNVs are good
  # - B: small and common, CNVs to to step 4
  # - C: something in between, CNVs to step 5

  # If the CNVR is rare (freq <= 1%) or large (overlap >= 0.75*llen)
  # all CNVs in it are good candidates, category A
  # If the CNVRs is frequent (freq >= 5% & freq < min_feq ,e.g. 25) and
  # small (overlap <= 55%*llen), then it is category B
  # If it has no category it's C
  cnvrsA <- c()
  cnvrsB <- c()
  cnvrsC <- c()

  for (loc in unique(loci$locus)) {

    nCNVS <- nrow(cnvs[locus == loc, ])

    # the two CNVRs frequency thresholds, 1% (rare) and 5% (common)
    th1 <- nCNVS*0.01
    th5 <- nCNVS*0.05
    # IF setted commonCNVRsMinFreq override th5, if th5 is smaller for the
    # specific locus. This is useful to avoid excluding calls in very rare loci
    if (!is.na(minFreq)) {
      if (minFreq > th5) th5 <- minFreq
    }

    loc_line <- getline_locus(loci[locus == loc, ])
    # all overlapping CNVRs, it is not a perfect methods but the
    # CNVRs that can give problems are the very large ones, that
    # would be saved anyway
    lcnvrs <- cnvrs[chr == loc_line[2] &
                    start <= loc_line[4] & end >= loc_line[3], ]
    # compute overlap prop
    lcnvrs[, op := (pmax(start,loc_line[3])-
                    pmin(end,loc_line[4])+1) / loc_line[5]]

    cnvrsA <- c(cnvrsA, lcnvrs[freq <= th1 & op >= 0.75, CVNR_ID])
    # the user should be able to change these values
    cnvrsB <- c(cnvrsB, lcnvrs[freq >= th5 & op <= 0.55, CVNR_ID])
    cnvrsC <- c(cnvrC, lcnvrs[, CVNR_ID][!lcnvrs[, CNVR_ID] %in% c(cnvrsA,cnvrsB)])
  }
  return(list(cnvrsA, cnvrsB, cnvrsC))
}
