
devtools::build_vignettes()

devtools::document()

devtools::test()

devtools::load_all()

devtools::check()

devtools::build()

devtools::install()

browseVignettes("QCtreeCNV")

library(data.table)

f <-  QCtreeCNV:::step3
debugonce(f)

rhub::check_on_linux("../")
rhub::check_on_windows("../")
rhub::check_for_cran("../")


# good (GT 1, deletions) and bad (GT 2, duplications) examples plots 1
dt1 <- rbind(data.table(prevalence=1, CImin=0.85, CImax=1.15, GT=1, LRRSD_group="low"),
             data.table(prevalence=1.4, CImin=1.25, CImax=1.55, GT=2, LRRSD_group="low"),
             data.table(prevalence=1.05, CImin=0.85, CImax=1.3, GT=1, LRRSD_group="medium"),
             data.table(prevalence=1.4, CImin=1.15, CImax=1.6, GT=2, LRRSD_group="medium"),
             data.table(prevalence=0.95, CImin=0.55, CImax=1.45, GT=1, LRRSD_group="high"),
             data.table(prevalence=0.9, CImin=0.75, CImax=1.05, GT=2, LRRSD_group="high"))
# LRRSD as a factor
dt1[, LRRSD_group := factor(LRRSD_group, levels = c("low", "medium", "high"))]

pla<- ggplot(aes(y = prevalence, x = LRRSD_group, colour = as.character(GT)), data = dt1) +
          geom_point(position = position_dodge(0.3), size = 3) +
          geom_errorbar(aes(ymin = CImin, ymax = CImax, colour = as.character(GT)),
                        position = position_dodge(0.3), size = 0.5, width = 0.3) +
          scale_colour_discrete(name = "Type", labels = c("Del", "Dup")) + theme_bw()


dt2 <- rbind(data.table(numsnp = rnorm(100, 41, 5), GT = 1),
             data.table(numsnp = rnorm(100, 62, 10), GT = 1),
             data.table(numsnp = rnorm(100, 33, 11), GT = 2),
             data.table(numsnp = rnorm(100, 58, 7), GT = 2))

plb <- ggplot(aes(x = numsnp, colour = as.character(GT)), data = dt2[numsnp >= 20, ]) +
           geom_density() +
           scale_colour_discrete(name = "Type", labels = c("Del", "Dup")) +
           geom_vline(xintercept = 20, linetype="dashed") +
           theme_bw()

cowplot::plot_grid(pla,plb, nrow = 2, labels="AUTO")


# simulate a CNV in a larger region marked by a smaller call
pos <- seq(from = 1, by = 1000, length.out = 1000)
loc_pos <- 301:700 ; cnv_pos <- c(pos[350],pos[550])
baf0 <- rep(0, times= 1000) + abs(rnorm(1000, 0, 0.06))
baf05 <- rep(0.5, times= 1000) + ( sample(c(1, -1), 1000, T) * abs(rnorm(1000, 0, 0.075)) )
baf1 <- rep(1, times= 1000) - abs(rnorm(1000, 0, 0.06))
lrr0 <- rep(0, times= 1000) + ( sample(c(1, -1), 1000, T) * abs(rnorm(1000, 0, 0.1)) )
lrrDEL <- rep(-0.5, times= 1000) + ( sample(c(1, -1), 1000, T) * abs(rnorm(1000, 0, 0.1)) )

baf <- c( sample(c(baf0, baf05, baf1), 300), sample(c(baf0, baf1), 400),
          sample(c(baf0, baf05, baf1), 300) )
lrr <- c( sample(c(lrr0), 300), sample(c(lrrDEL), 400), sample(c(lrr0), 300) )

dt3 <- data.table(position = pos, BAF = baf, LRR = lrr, locus = F)
dt3[loc_pos, locus := T]

pl1 <- ggplot(aes(x = position, y = LRR, colour = locus), data = dt3) +
         geom_point() + theme_bw() + theme(legend.position = "none") + ylim(-1, 1) +
         geom_segment(aes(x = cnv_pos[1], y = 0, xend = cnv_pos[2], yend = 0),
                      arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")), colour = "black")

pl2 <- ggplot(aes(x = position, y = BAF, colour = locus), data = dt3) +
         geom_point() + theme_bw() + theme(legend.position = "none") + ylim(0, 1) +
         geom_segment(aes(x = cnv_pos[1], y = 0.5, xend = cnv_pos[2], yend = 0.5),
                      arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")), colour = "black")
cowplot::plot_grid(pl1,pl2, nrow = 2, labels="AUTO")
