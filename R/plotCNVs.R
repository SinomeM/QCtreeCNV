#' Save CNVs plots for Visual Inspection
#'
#' Plot LRR BAF plots and save them to be used in visual inspection
#'
#' @param cnvs lorem ipsum
#' @param loci lorem ipsum
#' @param plots_path lorem ipsum
#' @param samples_list lorem ipsum
#' @param reg_len lorem ipsum
#'
#' @export
#'
#' @import data.table
#' @import ggplot2

saveCNVplots <- function(cnvs, plots_path, loci, samples_file, reg_len = 2000000, overwrite = F) {

  for (i in 1:nrow(cnvs)) {

    if (overwrite)
      if (file.exists(pl_path) next

    s <- cnvs[i, sample_ID]
    loc <- cnvs[i, locus]
    ll <- getline_locus(loci[locus == loc,])
    r_st <- as.integer(ll[3]) - reg_len
    r_en <- as.integer(ll[3]) + reg_len

    fp <- samples_file[sample_ID == s, file_path][1]
    if (length(fp) == 0) stop("Can't find sample ", s,
                              " in samples list file")
    raw <- fread(fp, skip = "Sample ID")[
             Chr == ll[2] & between(Position, r_st, r_en),
             .(Chr, Position, `Log R Ratio`, `B Allele Freq`)]
    colnames(raw) <- c("Chr", "Position", "LRR", "BAF")
    raw[between(as.integer(Position), as.integer(ll[3]), as.integer(ll[4])),
        core := T][is.na(core), core := F]

    m <- 1000000
    cc <- ifelse(cnvs[i, GT] == 1, "red", "green")
    vv <- ifelse(cnvs[i, GT] == 1, -1.5, 1.5)
    pl_lrr <-
      ggplot() +
        geom_point(data = raw, aes(Position/m, LRR, colour=core)) +
        scale_color_manual(values = c("black", "blue")) +
        geom_segment(aes(x=r_st/m, xend=r_en/m, y=0, yend=0),
                     linetype = "dashed", size=0.15, alpha=0.75) +
        xlab("Position (Mbp)") +
        guides(colour = "none") +
        geom_segment(data=cnvs[i], aes(y=vv, yend=vv, x=start/m, xend=end/m),
                     colour=cc, arrow=arrow(angle=90, ends="both",
                                            length=unit(0.05, "inches"))) +
        ylim(-2, 2) +
        xlim(r_st/m, r_en/m) +
        theme_classic()
    pl_baf <-
      ggplot() +
        geom_point(data = raw, aes(Position/m, BAF, colour=core)) +
        scale_color_manual(values = c("black", "blue")) +
        geom_segment(aes(x=r_st/m, xend=r_en/m, y=0.5, yend=0.5),
                     linetype = "dashed", size=0.15, alpha=0.75) +
        xlab("Position (Mbp)") +
        guides(colour = "none") +
        ylim(0, 1) +
        xlim(r_st/m, r_en/m) +
        theme_classic()

      pl <- cowplot::plot_grid(pl_baf, pl_lrr, nrow = 2)

      pl_path <- paste0(plots_path, cnvs[i, paste0("/",locus, "_",sample_ID, ".png")])

      cnvs[i, Plotname := pl_path]

      png(pl_path, 520, 480)
      print(pl)
      dev.off()
  }
  return(cnvs)
}
