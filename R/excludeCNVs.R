# loci
# chr start stop length

# cnvs
# chr start stop cnvr locus [ logr1 BAFc ]

# cnvrs
# chr start stop r_ID freq length

# chr is like chr1 chrX
# length is stop-start+1

cnv_excl <- function(loci,cnvrs,cnvs,showT=F,printdt=T,printthr=T,
                     cnvr_len=0.5,cnv_len=0.35,cnvr_freq1=0.025,cnvr_freq2=25,
                     minlog1=-0.175,maxlog1=0.175,maxlocusBAFc=0.05) {

  if (printthr)
    message("Exclusion results based on the following thresholds\n",
            "(logr1 and locusBAFc are not mandatory):\n",
            "- max CNVR overlap = ", cnvr_len*100,"% of locus\n",
            "- max CNV overlap  = ", cnv_len*100, "% of locus\n",
            "- min CNVR freq    = ", cnvr_freq1*100, "% of tot CNVs in locus\n",
            "- min CNVR freq    = ", cnvr_freq2, " CNVs (absolute number)\n",
            "- logr1 boundaries = ", minlog1," - ", maxlog1, "\n",
            "- max locusBAFc    = ", maxlocusBAFc)

  cnvrsExcl <- data.table()
  cnvsExcl <- data.table()
  excltrue <- data.table()

  for (i in 1:nrow(loci)) {

    cc <- copy(cnvs)
    rr <- copy(cnvrs)

    lll <- loci$locus[i]
    lcc <- loci$chr[i]
    lst <- loci$start[i]
    lsp <- loci$stop[i]
    lle <- loci$length[i]

    # ALL overlapping CNVRs
    lcnvRs <- rr[chr == lcc & start < lsp & stop > lst,]
    lff <- nrow(cc[locus == lll,])

    # knob 2,3 here, FREQ
    lcnvRsExcl <- lcnvRs[freq >= lff*cnvr_freq1 & freq >= cnvr_freq2,]
    # knob 1 here
    lcnvRsExcl[, a := pmax(lst, lcnvRsExcl$start)][
                         , b:= pmin(lsp,lcnvRsExcl$stop)][
                         , overlap := b-a+1][, c("a","b") := NULL]
    lcnvRsExcl <- lcnvRsExcl[overlap <= lle*cnvr_len,]
    # if the overlap is less than 10% the region must belong to another locus,
    # do not even consider it
    lcnvRsExcl <- lcnvRsExcl[overlap >= 0.1*lle,]
    lcnvRsExcl[, locus := lll]

    # all CVNs in locus & in selected regions
    lcnvsExcl <- cc[locus == lll & cnvr %in% lcnvRsExcl$r_ID,]

    # knob 4 here
    lcnvsExcl[, a := pmax(lst, lcnvsExcl$start)][
                         , b:= pmin(lsp,lcnvsExcl$stop)][
                         , overlap := b-a+1][, c("a","b") := NULL]
    lcnvsExcl <- lcnvsExcl[overlap <= lle*cnv_len,]

    # logr1, calls in those boundaries clusterize with the TRUE and should not be
    # excluded
    if (!is.na(minlog1) & !is.na(maxlog1))
      lcnvsExcl <- lcnvsExcl[logr1<=minlog1 | logr1>=maxlog1,]
    # locusBAFc, calls with a very low value are likely TRUE and should not be excluded
    if (!is.na(maxlocusBAFc))
      lcnvsExcl <- lcnvsExcl[BAFc >= maxlocusBAFc,]

    if (i == 1) {
      # initialize tables
      dt <- data.table(ix=i,locus=lll,excl=nrow(lcnvsExcl),
                       exclPerc=nrow(lcnvsExcl)/lff*100)
      # add Exclude TRUE if information is present
      if ("eval" %in% colnames(cnvs))
        dt[ix == i, exclTRUE := nrow(lcnvsExcl[eval == 1,])]
    } else {
      # add Exclude TRUE if information is present
      if ("eval" %in% colnames(cnvs))
        dt <- rbind(dt, data.table(ix=i,locus=lll,excl=nrow(lcnvsExcl),
                                   exclPerc=nrow(lcnvsExcl)/lff*100,
                                   exclTRUE = nrow(lcnvsExcl[eval == 1,])))
      else
        dt <- rbind(dt, data.table(ix=i,locus=lll,excl=nrow(lcnvsExcl),
                                   exclPerc=nrow(lcnvsExcl)/lff*100))
    } # end if else

    cnvrsExcl <- rbind(cnvrsExcl, lcnvRsExcl)
    cnvsExcl <- rbind(cnvsExcl, lcnvsExcl)
    if ("eval" %in% colnames(cnvs))
      excltrue <- rbind(excltrue, lcnvsExcl[eval == 1,])

  } # end for

  if (printdt) print(dt)
  if (showT) print(excltrue)

  if ("eval" %in% colnames(cnvs))
    return(list(dt,cnvsExcl,cnvrsExcl,excltrue))
  else
    return(list(dt,cnvsExcl,cnvrsExcl))
}
