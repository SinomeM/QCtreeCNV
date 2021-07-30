
select_stitch_calls <- function(cnvs, loci, minsnp = 20,
                               maxgap = 0.5, minoverlap = 0.2) {
  l <- 1

  for (i in 1:nrow(loci)) {
    lloc <- loci$locus[i]
    lchr <- loci$chr[i]
    lst <- loci$start[i]
    lsp <- loci$end[i]
    message("Locus ", lloc)

    lcnvs <- cnvs[chr == lchr & start <= lsp &
                  stop >= lst & numsnp >= minsnp, ]
    lcnvs[, densnp := round(length / numsnp , digits=0)]
    setorder(lcnvs, sample, type, start)
    lcnvs[, `:=` (stitch = 0, remove = F, gap = NA_real_)]

    tmpround <- 1
    tmpadd <- 1

    while (tmpadd > 0) {
      inrows <- nrow(lcnvs)

      for (j in 1:(nrow(lcnvs)-1)) {
      k <- j+1

      if (lcnvs$sample[j] == lcnvs$sample[k] & lcnvs$type[j] == lcnvs$type[k]) {
        ggap <- (lcnvs[k, start] - lcnvs[j, stop]) / (lcnvs[k, stop] - lcnvs[j, start])

        if (ggap < maxgap) {
          lcnvs[k, `:=` (stitch = 1, gap = ggap, start = lcnvs[j, start],
                        length = stop - start + 1, numsnp = numsnp + lcnvs[j, numsnp],
                        densnp = round(length / numsnp, digits = 0))]
          lcnvs[j, remove := T]
          }
        }
      }
      lcnvs <- lcnvs[remove == F,]
      outrows <- nrow(lcnvs)
      tmpadd <- inrows - outrows
      if (tmpadd > 0) message("round", tmpround)
      else message("done")
      tmpround <- tmpround+1
    }
    lcnvs[, remove := NULL]

    # compute overlap & add locus
    lcnvs[, `:=` (overlap = (pmin(stop, lsp) - pmax(start, lst) + 1) / (lsp-lst+1),
                  locus = lloc)]
    lcnvs[overlap >= minoverlap]

    if (!"overlap" %in% colnames(lcnvs)) lcnvs$overlap <- NA_real_
    
    if (nrow(lcnvs) > 0) {
      if (l == 1) {
        cnvsOUT <- lcnvs
        l <- 2
        }
      else cnvsOUT <- rbind(cnvsOUT, lcnvs)
    }
  }
  return(cnvsOUT)
}