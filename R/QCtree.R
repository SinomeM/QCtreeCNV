
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

qctree <- function(cnvs, cnvrs, qsdt,
                   maxLRRSD=.3, maxBAFdrift=.01,maxGCWF=.015, minGCWF=-.015,
                   ) {
  ### initial checks ###
  # # TODO!

  ### Pre process ###

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
  message("Step 1, QC outliers removal")
  cnvsOUT[, st1 := -1]
  # these pass the check
  cnvsOUT[LRRSD <= maxLRRSD & BAFdrift <= maxBAFdrift &
          between(GCWF, minGCWF, maxGCWF,incbounds=T), st1 := 1]
  # these are excluded
  cnvsOUT[st1 == -1, st1 := 0][st1 == 0, excl := 1]


  ### STEP 2 & 3 CNVRs ###
  message("Steps 2 and 3, CNVRs")

  # Sort CNVRs in the different categories, done per locus


  ### STEP 4, logr1, BAFc & BAFb checks 1 ###
  message("Step 4, logr1, BAFc & BAFb for small CNVRs calls")


  ### STEP 5, LRRSDlocus, mLRRlocus ###
  message("Step 5, LRRSD and mean LRR checks inside locus of interest")


  ### STEP 6, logr1, BAFc & BAFb checks 2 ###
  message("Step 4, logr1, BAFc & BAFb for bad LRRSD/mLRRlocus calls")


  ### RETURN ###

  # remove temporary columns
  cnvsOUT[, c("st1", "st2", "st3", "st4", "st5", "st6") := NULL]
  # return good and excluded
  return(list(cnvsOUT[excl == F, ], cnvsOUT[excl == T, ]))
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
