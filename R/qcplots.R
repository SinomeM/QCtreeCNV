
#' Create Quality Control plots for the filtering pipeline
#'
#'
#'
#' @export
#'
#' @import data.table ggplot2


# ALL columns
# chr, start, end, GT, numsnp, , numsnp, densnp, overlap, sample_ID, locus, putCarrier, LRRSD, BAFdrift, GCWF, logr1, LRRSDlocus, BAFc, BAFb, centDistProp, mLRRlocus

# FILTERS
# - select_stitch_calls() filters: numsnp, overlap
# - sample-wise filters: LRRSD, BAFdrift, GCWF
# - CNVRS filter is difficult to test, but calls belonging to class B CNVRs should not have any true CNV
# - logr1 and BAFb and BAFc are more complex to evaluate

# PLOTS
# - dot plot true CNV prevalence across 3 LRRSD chunks (DEL/DUP separated, but in the same plot).
#   Same for BAFdrift and GCWF?
# - prevalence vs numsnp and same for overlap. a smoothed line of some sort (DEL/DUP separated,
#   but in the same plot)

# I think the more general approach possibler is to create a folder a put all plots in it.
# For plotting only one (or a selection) of loci one should filter the table before
# passing it to the function, and give a meaningful name to the plot folder.


qc_plots_cnvs <- function(cnvs, folder_name, qc,
                          maxLRRSD=.35, maxBAFdrift=.01, maxGCWF=.02, minGCWF=-.02) {

  # apply the same filter to compute the prevalences
  qc <- qc[LRRSD <= maxLRRSD & BAFdrift <= maxBAFdrift &
             between(GCWF, minGCWF, maxGCWF,incbounds=T), kept := T]

  # these groups are created ad-hoc but should be quit general
  qc[between(LRRSD, 15, 25, incbounds = F), lrrsd_group := "medium"][
      LRRSD <= 25, lrrsd_group := "low"][LRRSD >= 15, lrrsd_group := "high"]
  qc[between(BAFdrift, 0.001, 0.002, incbounds = F), bafd_group := "medium"][
      BAFdrift <= 0.001, bafd_group := "low"][BAFdrift >= 0.002, bafd_group := "high"]
  qc[between(abs(GCWF), 0.005, 0.01, incbounds = F), gcwf_group := "medium"][
      abs(GCWF) <= 0.005, gcwf_group := "low"][abs(GCWF) >= 0.01, gcwf_group := "high"]

  # plot1, LRRSD chunks.
  dtlrr <- data.table(LRRSD_group = rep(c("low", "medium", "high"), 2), GT = rep(c(1,2),3),
                   prevalence = NA, CImin = NA, CImax = NA)

  a <- qc[, .N, by = .(lrrsd_group, GT)]; setnames(a, "N", "n")
  b <- qc[sample_ID %in% cnvs[visual_output == 1, sample_ID], .N, by = .(lrrsd_group, GT)]; setnames(a, "N", "x")

  dtlrr <- merge(dtlrr, merge(a,b))
  dtlrr[, prevalence := (x/n)*100]
  ## BETA distribution to compute the Confidence Intervals? ##
  ##
  pl1 <- ggplot(aes(y = prevalence, x = LRRSD_group, colour = GT), data = dtlrr) +
           geom_point(position = "dodge") +
           geom_errorbar(aes(ymin = CImin, ymax = CImax, colour = GT)) + theme_bw()

  # plot 2&3, numsnp and overlap density
  pl2 <- ggplot(aes(x = numsnp, colour = GT, linetype = visual_output), data = dt) +
           geom_density() + theme_bw()

  pl3 <- ggplot(aes(x = overlap, colour = GT, linetype = visual_output), data = dt) +
           geom_density() + theme_bw()

  # similar to pl1 but for BAFdrift and GCWF
  pl4 <- ggplot()
  pl5 <- ggplot()

}
