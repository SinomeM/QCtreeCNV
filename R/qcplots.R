
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
                          min_numsnp = 20, min_overlap = 0.2,
                          maxLRRSD=.35, maxBAFdrift=.01, maxGCWF=.02, minGCWF=-.02) {

  # keep only validated CNVs
  cnvs <- cnvs[visual_output == 1, ]

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
  dtlrr <- data.table(lrrsd_group = rep(c("low", "medium", "high"), 2), GT = rep(c(1,2),3),
                   prevalence = NA, CImin = NA, CImax = NA)

  a <- qc[, .N, by = .(lrrsd_group, GT)]; setnames(a, "N", "n")
  b <- qc[sample_ID %in% cnvs[visual_output == 1, sample_ID], .N, by = .(lrrsd_group, GT)]
  setnames(a, "N", "x")

  dtlrr <- merge(dtlrr, merge(a,b))
  dtlrr[, prevalence := round((x/n)*100, digits=3)]
  # BETA distribution to compute the Confidence Intervals
  # computed applying  the Bayesian credible interval using the Jeffreys prior
  # as in ‘Brown, Lawrence D., T. Tony Cai, and Anirban DasGupta. “Interval Estimation
  # for a Binomial Proportion.” Statistical Science 16, no. 2 (2001): 101–17.
  # http://www.jstor.org/stable/2676784.’
  dtlrr[, CImin := round(100*qbeta(c(0.05/2,1-0.05/2), x+0.05, n-x+0.05)[1], digits = 3)][,
          CImax := round(100*qbeta(c(0.05/2,1-0.05/2), x+0.05, n-x+0.05)[2], digits = 3)]

  pl1 <- ggplot(aes(y = prevalence, x = lrrsd_group, colour = as.character(GT)), data = dtlrr) +
           geom_point(position = position_dodge(0.3), size = 3) +
           geom_errorbar(aes(ymin = CImin, ymax = CImax, colour = as.character(GT)),
                         position = position_dodge(0.3), size = 0.5, width = 0.3) +
           scale_colour_discrete(name = "Type", labels = c("Del", "Dup")) + theme_bw()

  # plot 2&3, numsnp and overlap density
  pl2 <- ggplot(aes(x = numsnp, colour = as.character(GT)), data = cnvs) +
             geom_density() +
             scale_colour_discrete(name = "Type", labels = c("Del", "Dup")) +
             geom_vline(xintercept = min_numsnp, linetype="dashed") +
             theme_bw()

  pl3 <- ggplot(aes(x = overlap, colour = as.character(GT)), data = cnvs) +
             geom_density() +
             scale_colour_discrete(name = "Type", labels = c("Del", "Dup")) +
             geom_vline(xintercept = min_overlap, linetype="dashed") +
             theme_bw()


  # similar to pl1 but for BAFdrift and GCWF
  dtbaf <- data.table(bafd_group = rep(c("low", "medium", "high"), 2), GT = rep(c(1,2),3),
                   prevalence = NA, CImin = NA, CImax = NA)
  a <- qc[, .N, by = .(bafd_group, GT)]; setnames(a, "N", "n")
  b <- qc[sample_ID %in% cnvs[visual_output == 1, sample_ID], .N, by = .(bafd_group, GT)]
  setnames(a, "N", "x")
  dtbaf <- merge(dtlrr, merge(a,b))
  dtbaf[, prevalence := round((x/n)*100, digits=3)]
  dtbaf[, CImin := round(100*qbeta(c(0.05/2,1-0.05/2), x+0.05, n-x+0.05)[1], digits = 3)][,
          CImax := round(100*qbeta(c(0.05/2,1-0.05/2), x+0.05, n-x+0.05)[2], digits = 3)]

  pl4 <- ggplot(aes(y = prevalence, x = bafd_group, colour = as.character(GT)), data = dtbaf) +
           geom_point(position = position_dodge(0.3), size = 3) +
           geom_errorbar(aes(ymin = CImin, ymax = CImax, colour = as.character(GT)),
                         position = position_dodge(0.3), size = 0.5, width = 0.3) +
           scale_colour_discrete(name = "Type", labels = c("Del", "Dup")) + theme_bw()

  dtgc <- data.table(gcwf_group = rep(c("low", "medium", "high"), 2), GT = rep(c(1,2),3),
                   prevalence = NA, CImin = NA, CImax = NA)
  a <- qc[, .N, by = .(gcwf_group, GT)]; setnames(a, "N", "n")
  b <- qc[sample_ID %in% cnvs[visual_output == 1, sample_ID], .N, by = .(gcwf_group, GT)]
  setnames(a, "N", "x")
  dtgc <- merge(dtlrr, merge(a,b))
  dtgc[, prevalence := round((x/n)*100, digits=3)]
  dtgc[, CImin := round(100*qbeta(c(0.05/2,1-0.05/2), x+0.05, n-x+0.05)[1], digits = 3)][,
         CImax := round(100*qbeta(c(0.05/2,1-0.05/2), x+0.05, n-x+0.05)[2], digits = 3)]

  pl5 <- ggplot(aes(y = prevalence, x = gcwf_group, colour = as.character(GT)), data = dtgc) +
           geom_point(position = position_dodge(0.3), size = 3) +
           geom_errorbar(aes(ymin = CImin, ymax = CImax, colour = as.character(GT)),
                         position = position_dodge(0.3), size = 0.5, width = 0.3) +
           scale_colour_discrete(name = "Type", labels = c("Del", "Dup")) + theme_bw()

  # save plots in PDF
  dir.create(folder_name)
  pdf(system.file(folder_name, "plot1.pdf"))
  print(pl1)
  pdf(system.file(folder_name, "plot2.pdf"))
  print(pl2)
  pdf(system.file(folder_name, "plot3.pdf"))
  print(pl3)
  pdf(system.file(folder_name, "plot4.pdf"))
  print(pl4)
  pdf(system.file(folder_name, "plot5.pdf"))
  print(pl5)
  dev.off()
}
