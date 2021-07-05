
## TODO!
# 1. Move here most of the code from exlcudeCNVs.R, all that can be reused
# 2. Integrate additional steps in the pipeline
# 3. Update the suggested values for each step, reflecting the new run in GDK
#    (after initial QC outliers removal)
# 4. Some documentation, keep in mind most of it will live in more discoursive
#    manner in the vignette.


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
                   st4minlogr1=-.35,st4maxlogr1=.4,
                   st4maxBAFcDEL=.03, st4maxBAFcDUP=.125,
                   st4maxBAFbDEL=0.075, st4maxBAFbDUP=NA,
                   st5maxLRRSDlocus,
                   st5maxBAFcDEL=.05, st5maxBAFcDUP=.15,
                   st5maxBAFbDEL=.1, st5maxBAFbDUP=NA) {
  ### initial checks ###
  # # TODO!

  ### Pre process ###
  message("# -------------------------- #\n",
          "Step 0, preprocess")

  # add index to cnvs
  cnvs[, ix := 1:nrow(cnvs)]
  # add qs measures, i.e. merge the two tables
  setkey(cnvs, c("sample_ID", "locus")
  setkey(qsdt, c("sample_ID", "locus")
  # join tables using data.table keys
  cnvsOUT <- cnvs[qsdt]
  # add excl column to keep track of the process
  cnvsOUT[, excl := -1]


  ### STEP 1, QC outliers ###

  # a mild (in the default settings) filter on the three
  # classic measures, LRRSD, BAFdrift and GCWF
  message("# -------------------------- #\n",
          "Step 1, QC outliers removal")
  cnvsOUT[, st1 := -1]
  # these pass the check
  cnvsOUT[LRRSD <= maxLRRSD & BAFdrift <= maxBAFdrift &
          between(GCWF, minGCWF, maxGCWF,incbounds=T), st1 := 1]
  # these are excluded
  cnvsOUT[st1 == -1, st1 := 0][st1 == 0, excl := 1]
  message("# ", nrow(cnvsOUT[excl == 1, ]),
          "calls are from QC outlier samples.\n")


  ### STEP 2 & 3 CNVRs ###
  message("# -------------------------- #\n",
          "Steps 2 and 3, CNVRs")

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

    nCNVS <- nrow(cnvsOUT[locus == loc, ])

    # the two CNVRs frequency thresholds, 1% (rare) and 5% (common)
    th1 <- nCNVS*0.01
    th5 <- nCNVS*0.05
    # IF setted commonCNVRsMinFreq override th5, if th5 is smaller for the
    # specific locus. This is useful to avoid excluding calls in very rare loci
    if (!is.na(commonCNVRsMinFreq)) {
      if (commonCNVRsMinFreq > th5) th5 <- commonCNVRsMinFreq
    }

    loc_line <- getline_locus(loci, loc)
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

  # STEP 2
  cnvsOUT[, st2 := -1]
  # CNVs from step 1 == 0 that are in a cnvrA will pass step 2 (to good CNVs)
  cnvsOUT[st1 == 0 & CNVR_ID %in% cnvrsA, `:=` (st2 = 1, excl = 0)]
  # CNVs from step 1 == 0 that are in a cnvrA will fail step 2 (to step3)
  cnvsOUT[st1 == 0 & !CNVR_ID %in% cnvrsA, st2 := 0]
  # check all CNVs from step 1 are assigned
  if (nrow(cnvsOUT[st1 == 0,] != nrow(cnvsOUT[st2 %in% c(1,0), ]))
    stop("There is a problem in step 2")

  # STEP 3
  cnvsOUT[, st3 := -1]
  # CNVs from step 2 == 0 that are in a cnvrB will go to step 4
  # CNVs from step 2 == 0 that are in a cnvrC will go to step 5
  cnvsOUT[st2 == 0 & CNVR_ID %in% cnvrB, st3 := 1]
  cnvsOUT[st2 == 0 & CNVR_ID %in% cnvrC, st3 := 0]

  # check all CNVs from step 2 are assigned
  if (nrow(cnvsOUT[st2 == 0,] != nrow(cnvsOUT[st3 %in% c(1,0), ]))
    stop("There is a problem in step 3")



  ### STEP 4, logr1, BAFc & BAFb checks 1 ###
  message("# -------------------------- #\n",
          "Step 4, logr1, BAFc & BAFb for small CNVRs calls")
  cnvsOUT[, st4 := -1]
  # Deletions and duplications are treated differently here
  # Takes calls from st3 == 1
  # If log1 or at least one between BAFc and BAFb are whithin limits, it goes
  # to step 5, i.e. st4 = 0
  # If logr1 AND BAFc | BAFb are out of limits (i.e. consistent with small CNVR)
  # exclude them
  cnvsOUT[st3 == 1 & GT == 1 & (between(logr1, st4minlogr1, st4maxlogr1) |
          (BAFc <= st4maxBAFcDEL | BAFb <= st4maxBAFbDEL)), st4 := 0]
  cnvsOUT[st3 == 1 & GT == 1 & st4 == -1, `:=` (st4 = 1, excl = 1)]

  cnvsOUT[st3 == 1 & GT == 2 & (between(logr1, st4minlogr1, st4maxlogr1) |
          BAFc <= st4maxBAFcDUP), st4 := 0]
  cnvsOUT[st3 == 1 & GT == 2 & st4 == -1, `:=` (st4 = 1, excl = 1)]

  # check all CNVs from step 3 are assigned
  if (nrow(cnvsOUT[st3 == 1,] != nrow(cnvsOUT[st4 %in% c(1,0), ]))
    stop("There is a problem in step 4")


  ### STEP 5, LRRSDlocus, mLRRlocus ###
  message("# -------------------------- #\n",
          "Step 5, LRRSD and mean LRR checks inside locus of interest\n",
          "plus logr1, BAFc & BAFb check for all calls")
  cnvsOUT[, st5 := -1]
  # Takes calls from st3 = 0 and st4 = 0

  # mLRRlocus particularly out
  cnvsOUT[GT == 1 & st3 == 0 & st4 == 0 & mLRRlocus > 0, `:=` (st5 = 1, excl = 1)]
  cnvsOUT[GT == 2 & st3 == 0 & st4 == 0 & mLRRlocus < -0.3, `:=` (st5 = 1, excl = 1)]

  # In calls with high LRRSDlocus, check BAFc and BAFb, if at least one of them
  # well out of range or both with lower thresholds then it's excluded
  # Deletions
  cnvsOUT[GT == 1 & st3 == 0 & st4 == 0 & LRRSDlocus > st5maxLRRSDlocus &
          ((BAFc <= st5maxBAFcDEL | BAFb <= st5maxBAFbDEL) |
           (BAFc <= st4maxBAFcDEL & BAFb <= st4maxBAFbDEL)), `:=` (st5 = 1, excl = 1)]
  # Duplications
  cnvsOUT[GT == 2 & st3 == 0 & st4 == 0 & LRRSDlocus > st5maxLRRSDlocus &
          BAFc <= st5maxBAFcDUP, `:=` (st5 = 1, excl = 1)]

  cnvsOUT[GT == 1 & st3 == 0 & st4 == 0 & st5 == -1, `:=` (st5 = 0, excl = 0)]

  ### RETURN ###
  message("# -------------------------- #\n",
          "Pipeline complete!")
  # check all CNVs were evaluated
  if (nrow(cnvsOUT[excl == -1,] == 0) message("All CNVs were evaluated")

  # remove temporary columns
  cnvsOUT[, c("st1", "st2", "st3", "st4", "st5") := NULL]
  # return good and excluded
  return(list(cnvsOUT[excl == 0, ], cnvsOUT[excl == 1, ]))
}


# ---------------------

# return sample_ID, locus, GT and CN as a vector of <chr>
# data.table is smart enough to convert the values back to
# <int> when needed (of GT and CN)

getline_cnv <- function(cnvs, i) {
  ll <- cnvs[i]
  myline <- c(ll$sample_ID, ll$locus, ll$GT, ll$CN, ll$CNVR_ID)
  return(myline)
}


# return start end and freq

getline_cnvr <- (cnvrs, cnvr_id) {
  ll <- cnvrs[CNVR_ID == cnvr_id, ]
  myline <- c(ll$start, ll$end, ll$freq)
  return(myline)
}

# same but for loci
getline_locus <- function(loci, locus) {
  myline <- c(locus, loci[locus == locus, chr], loci[locus == locus, start],
              loci[locus == locus, end], loci[locus == locus, length])
  return(myline)
}
