
#' Extratc QC metrics
#'
#' \code{extractMetrics()} extract the QC metrics from the raw data in RDS
#' format.
#'
#' @export
#'
#' @import data.table


## TODO!
# 1. move here the code from the RMD file from GDK and standardize it
# 2. minimal documentation, in particular on the input files!

extractMetrics <- function(loci, cnvs, pennQC, int_rds_path, tmp_rds_path) {

  # initial checks
  # # TODO

  ids <- unique(cnvs$sample_ID)

  for (l in 1:nrow(loci)) {
    loc <- getline_locus(loci[i])
    dt <- data.table(sample_ID = ids)
    message("Locus #", l, ": ", ll)

    for (s in ids) {
      tmp <-
        readRDS(file.path(int_rds_path,s,".rds"))[Chr == loc[2] &
                                                  between(Position, loc[3], loc[4]),]

      bafc <- nrow(tmp[between(`B Allele Freq`, 0.4, 0.6, incbounds=T), ]) /
                nrow(tmp)
      bafb <- nrow(tmp[between(`B Allele Freq`, 0.2, 0.4, incbounds=F) |
                       between(`B Allele Freq`, 0.6, 0.8, incbounds=F) |
                       `B Allele Freq` %in% c(0.2, 0.8), ]) / nrow(tmp)

      # locus measures (compute them regardless of the presence of a call)
      dt[sample_ID == s, `:=` (mLRRlocus = mean(tmp[, `Log R Ratio`], na.rm=T),
                         LRRSDlocus = sd(tmp[, `Log R Ratio`], na.rm=T),
                         BAFc = bafc, BAFb = bafb,
                         centLocus = lst+(lsp-lst+1)/2 ,lenLocus=lsp-lst+1)]

      # check if this sample has a call in the locus
      put <- cnvs[Locus == ll & sample_ID == s,]
      if (nrow(put) == 0) {
        dt[sample_ID == s, `:=` (mLRRcall = NA_real_, LRRSDcall = NA_real_,
                                 centCall = NA_real_, lenCall = NA_real_)]
      } else {
        # A sample can have only one call per locus
        if (nrow(put) > 1) stop(paste0("Sample ", s, " has more than one call in",
                                       " locus ", loc[1]))
        putline <- getline_cnv(put)
        tmp1 <- tmp[between(Position, putline[6], putline[7]),]
        dt[sample_ID == s, `:=` (mLRRcall = mean(tmp1[, `Log R Ratio`], na.rm=T),
                                 LRRSDcall = sd(tmp1[, `Log R Ratio`], na.rm=T),
                                 centCall = putline[9], lenCall = putline[8])]
      }
    }
    # rind all loci together
    dtOUT <- rbind(dtOUT, dt)
  }

  return(dtOUT)
}

# TODO:
# Integrate the code below in the main function,
# the less step required the better.
# Do not compute stuff not needed.
# Remember that the eval results are only in our case

{
dt <- data.table()

for (samp in c("2012", "2015i")) {
  message(samp)
  files <- list.files(paste0("dts/",samp,"/"))

  for (f in files) {
    message(f)
    tmp <- readRDS(paste0("dts/",samp,"/",f))
    ll <- gsub(".rds", "", f)
    tmp$locus <- ll
    tmp$sample <- samp
    if (samp == "2012") {
      put <- put2012[Locus == ll,]
    } else
        put <- put2015i[Locus == ll,]
    # vectors of pids
    pput <- put[, pid]
    pputT <- put[cons_eval == 1, pid]
    pputF <- put[cons_eval == 2, pid]
    pputU <- put[cons_eval == 3, pid]
    pputNA <- put[is.na(cons_eval), pid]

    # update tmp
    tmp[pid %in% pput, putCarrier := T][!pid %in% pput, putCarrier := F][
          pid %in% pputT, CNV := 1][pid %in% pputF, CNV := 2][
          pid %in% pputU, CNV := 3][pid %in% pputNA, CNV := NA]

    lst <- loci30[locus == ll, start]
    lsp <- loci30[locus == ll, stop]
    overlaps <- pmin(put[match(pput, pid),stop],lsp) - max(put[match(pput, pid),start],lst) + 1
    tmp[match(pput, pid), overlap:=overlaps]

    dt <- rbind(dt, tmp)
  }
}

dt[, r1:=abs(mLRRcall / mLRRlocus)][, logr1 := log(r1 + 0.000001)][
     , centDistProp := abs(centLocus-centCall)/lenLocus]
dt[, CNVedit := CNV][is.na(CNV) & putCarrier ==T, CNVedit := 2][CNV == 3, CNVedit := 2]
dt[, r2 := abs(LRRSDcall/LRRSDlocus)]
dt[, overlapProp:= overlap/lenLocus]
dt[, r3 := abs(mBAFcall/mBAFlocus)][, logr3 := log(r3 + 0.000001)]

dt[,id := paste0(sample,locus,pid)]
put2012[, id := paste0("2012",Locus,pid)]
put2015i[, id := paste0("2015i",Locus,pid)]
setorder(dt, id) ; setorder(put2012, id) ; setorder(put2015i, id)
dt[id %in% put2012$id,  `:=` (type = put2012$Type, conf = put2012$conf)][
     id %in% put2015i$id,  `:=` (type = put2015i$Type, conf = put2015i$conf)]

dt <- dt[, .(pid,locus,sample,putCarrier,CNV,CNVedit,type,conf,
             mLRRlocus,mLRRcall,r1,logr1,LRRSDlocus,LRRSDcall,r2,
             BAFc,BAFb,mBAFlocus,mBAFcall,r3,logr3,
             centCall,centLocus,centDistProp,lenLocus,overlapProp)]
}
