#' Compute Copy Number Variable Regions (CNVRs)
#'
#' Faster CNVRs computation algorithm
#'
#' It is based on pamk clustering algorithm. It is a much simpler process that
#' takes advantage of the fact that this package is meant to deal only with groups
#' of calls overlapping a specific locus.
#'
#' @param cnvs lorem ipsum
#' @param chr_arms lorem ipsum
#' @param final_merge lorem ipsum
#'
#' @return a \code{list} of two elements. The first element is a \code{data.table}
#'   that contains the actual CNVR information, genomic location and frequency in
#'   the cohort. The second element is the \code{CNVresults}
#'
#' @export
#'
#' @import data.table


cnvr_fast <- function(put_cnvs) {

  cnvs_all <- data.table()
  cnvrs_all <- data.table()

  for (l in unique(put_cnvs$locus)) {
    cnvs <- put_cnvs[locus == l, ]
    cnvs[, `:=` (start = as.integer(start), end = as.integer(end),
                 length = as.integer(end - start + 1),
                 center = as.integer(start + length/2))]

    # cluster CNVs
    tmp <- data.frame(X = scale(cnvs$center), Y = scale(cnvs$length))
    # better than nothing
    if (sum(is.na(tmp)) > 0)
      tmp <- data.frame(X = cnvs$center, Y = cnvs$length)

    if (nrow(tmp) > 3) {
      maxk <- nrow(tmp)-1
      rr <- c(1, seq(2, 54, by = 4))
      r1 <- rr[1:8]

      # it's hard to imagine more than 24 CNVRs per window
      pam_cl <- fpc::pamk(tmp, krange = r1[r1 < maxk])
      nc <- pam_cl$nc
      # just in case
      if (nc == 24) pam_cl <- fpc::pamk(tmp, krange = rr[9:15])
      # fine tune the actual number of clusters
      pam_cl <- fpc::pamk(tmp, krange = nc-3:nc+3)
      cl <- as.integer(pam_cl$pamobject$clustering)
    }
    else {
      cl <- rep(1, nrow(tmp))
    }

    # assign each CNV to his cluster
    cnvs[, cnvr := cl]
    # create CNVRs table
    cnvrs <- data.table()
    # cluster 0 are the outliers, skip them for now
    for (i in 1:max(cl)) {
      # smooth the boundaries, median is less sensible to extreme values
      st <- as.integer(median(cnvs[cnvr == i, start]))
      en <- as.integer(median(cnvs[cnvr == i, end]))
      cnvrs <- rbind(cnvrs, data.table(cnvr = i, start = st, end = en))
    }

    cnvs[, cnvr := paste0(l, "_", cnvr)]
    cnvrs[, cnvr := paste0(l, "_", cnvr)]

    cnvs_all <- rbind(cnvs_all, cnvs)
    cnvrs_all <- rbind(cnvrs_all, cnvrs)
  }

  return(list(cnvs = cnvs_all, cnvrs = cnvrs_all))

}


#' Compute Copy Number Variable Regions (CNVRs)
#'
#' \code{cnvrs_create} compute CNVRs from a list of CNVs.
#'
#' Copy Number Variable Regions (CNVRs) are defined as groups of reciprocal overlapping
#' CNVs. This function first try to assign every call to a CNVR (or create a new
#' one if it is not possible), then check if adjacent CNVRs cam be merged, and
#' finally recheck all the CNVs in each CNVRs and un-assign them if the reciprocal
#' overlap is no longer fulfilled with all the member of the CNVR. If any event
#' is touched by this last step, a new cycle begins, this until no more CNVs can
#' be removed from the assigned CNVR.
#'
#' @param cnvs lorem ipsum
#' @param chr_arms lorem ipsum
#' @param prop reciprocal overlap proportion, default 0.3 (30%).
#' @param final_merge lorem ipsum
#'
#' @return a \code{list} of two elements. The first element is a \code{data.table}
#'   that contains the actual CNVR information, genomic location and frequency in
#'   the cohort. The second element is the \code{CNVresults}
#'
#' @export
#'
#' @import data.table


## TODO!
# 1. CNVRs should be computed per LOCUS here!!!!
# 1. add last step, merge 75% overlap CNVRs + update CNVRs + update CNVs
# 2. update documentation and initial checks

cnvrs_create <- function(cnvs, chr_arms, prop = 0.3, final_merge = T) {
  # check input
  if (!is.data.table(cnvs) | !is.data.table(chr_arms))
    stop("Inputs must be data.table!\n")

  # data.table "set" and ":=" functions act by reference, I create a copy to
  # avoid modifying the original object (perhaps there is a better option?)
  cnvs_cp <- copy(cnvs)
  rm(cnvs)

  # check input formats, in particular for chr_arms
  chr_arms[, `:=` (start = as.integer(start), end = as.integer(end))]
  chr_uniform(chr_arms)

  # sort per chr, start & end
  setorder(cnvs_cp, chr, start, end)
  setorder(chr_arms, chr, start)
  cnvs_cp[, CNVR_ID := NA_character_]

  # cnvrs: arm_ID, chr, start, end
  # this create a line with all "NA", I remove it later
  cnvrs <- data.table("CNVR_ID" = NA_character_, "chr" = NA_character_,
                      "start" = NA_integer_, "end" = NA_integer_)
  res <- data.table()

  # proceed per chromosomal arms
  for (arm in chr_arms$arm_ID) {
    message("\n-- Computing CNVRs in chromosome arm:", arm, "--\n")

    # cnvrs for this arm
    cnvrs_tmp <- data.table("CNVR_ID" = NA_character_, "chr" = NA_character_,
                            "start" = NA_integer_, "end" = NA_integer_)
    # arm related variables
    reg_arm <- get_region(chr_arms[arm_ID == arm,])

    # subset compatible cnvs
    DT <- cnvs_cp[chr == reg_arm[1] &
                    between(start, reg_arm[2], reg_arm[3]) &
                    between(end, reg_arm[2], reg_arm[3]), ]
    DT[, ix := 1:nrow(DT)]
    # if no CNVs are present in the arm, just skip
    if (nrow(DT) == 0) next

    # keep track of CNVR and loop number
    n_loop <- 1
    n <- 1

    # for each arm keep running until no more cnv is excluded
    while (any(is.na(DT$CNVR_ID)) & n_loop <= 50) {

      # create CNVRs/fill CNVRs
      message("\nCreating/filling CNVRs, loop #", n_loop, "...\n")
      indexes <- DT[is.na(CNVR_ID), ix]
      message(length(indexes), " CNVs to be assigned")

      tmp <- create_fill_CNVR(cnvrs_tmp, DT, n, prop, indexes, arm, reg_arm)
      DT <- tmp[[2]]
      cnvrs_tmp <- tmp[[1]]
      cnvrs_tmp <- cnvrs_tmp[!is.na(CNVR_ID) & CNVR_ID %in% DT$CNVR_ID, ]
      n <- tmp[[3]]


      # Check if some CNVRs can be merged
      message("\nRe-checking CNVRs 1...")
      if (nrow(cnvrs_tmp > 1)) {
        # first run
        tmp <- check_cnvrs(cnvrs_tmp, DT, n, prop, arm, reg_arm)

        # if CNVRs are updated tmp[[4]] is TRUE
        while(tmp[[4]]) {
          tmp <- check_cnvrs(tmp[[1]], tmp[[2]], tmp[[3]], prop, arm, reg_arm)
        }
        cnvrs_tmp <- tmp[[1]]
        DT <- tmp[[2]]
        n <- tmp[[3]]
      }

      # Check if some CNVRs can be merged
      # the function can get stuck here so the final 5 loop are run with an even
      # more stringent threshold
      message("\nRe-checking CNVRs 2...")
      if (nrow(cnvrs_tmp > 1)) {
        if (n_loop <= 25) ppp <- 0.9
        else ppp <- 0.95
        # first run
        tmp <- merge_cnvrs(cnvrs_tmp, DT, n, prop=ppp, arm, reg_arm)

        # if CNVRs are updated tmp[[4]] is TRUE
        while(tmp[[4]]) {
          tmp <- merge_cnvrs(tmp[[1]], tmp[[2]], tmp[[3]], prop=ppp, arm, reg_arm)
          #           Sys.sleep(0.5)
        }
        cnvrs_tmp <- tmp[[1]]
        DT <- tmp[[2]]
        n <- tmp[[3]]
      }

      # Recheck CNVs
      message("\nRe-checking CNVs ...")

      # # Move CVNs if a better overlap is found
      # for (i in DT[, ix]) {
      #
      # }

      # Remove CNVs if necessary
      if (n_loop < 50) DT <- remove_cnvs(DT, prop)

      # remove CNVRs if necessary
      cnvrs_tmp <- cnvrs_tmp[CNVR_ID %in% unique(DT$CNVR_ID), ]

      n_loop <- n_loop + 1
    }

    # Recreate the output
    res <- rbind(res, DT[, ix := NULL])
    cnvrs <- rbind(cnvrs, cnvrs_tmp)
  }

  setorder(cnvrs, end, start)
  # Check if there are duplicated CNVRs (same start & stop) for some reason
  message("\nCNVRs final check...")
  if (nrow(cnvrs > 1)) {
    # first run
    tmp <- dupl_cnvrs(cnvrs, res)

    # if CNVRs are updated tmp[[4]] is TRUE
    while(tmp[[3]]) {
      tmp <- dupl_cnvrs(tmp[[1]], tmp[[2]])
    }
    cnvrs <- tmp[[1]]
    res <- tmp[[2]]
  }

  # CNVRs can no longer have CNVs, clean those ones
  cnvrs <- cnvrs[CNVR_ID %in% unique(res$CNVR_ID), ]
  # count the final CNVs frequency per CNVRs
  freqs <- res[, .N, by = CNVR_ID]
  cnvrs[, freq := freqs[match(cnvrs$CNVR_ID, freqs$CNVR_ID), N]]

  # Merge overlapping CNVRs regardless of the CNVs contained
  if (final_merge) {
    message("Merging 75% overlapping CNVRs regardless of the CNVs")
    tmp <- merge75(cnvrs, res)
    cnvrs <- tmp[[1]]
    res <- tmp[[2]]
  }

  return(list(cnvrs[!is.na(chr), ], res))
}


## SUBFUNCTIONS

create_fill_CNVR <- function(cnvrs, DT, n, prop, ixs, arm, reg_arm) {

  for (i in ixs) {
    bb <- FALSE
    my_reg <- get_region(DT[ix == i, ], prop)
    # possible cnvrs
    if (nrow(cnvrs[!is.na(CNVR_ID), ]) != 0) {
      cnvr_m <- cnvrs[!is.na(CNVR_ID),][
        between(start, my_reg[2], my_reg[3], incbounds = TRUE) |
          between(end, my_reg[2], my_reg[3], incbounds = TRUE), ]
    } # this if else is needed for the initial CNVR (of each chr)
    else cnvr_m <- cnvrs[!is.na(CNVR_ID), ]

    # no match, initialize new cnvr
    if (nrow(cnvr_m) == 0) {
      cnvrs <- rbind(cnvrs, data.table("CNVR_ID" = paste0(arm, "-", n),
                                       "chr" = reg_arm[1], "start" = my_reg[2],
                                       "end" = my_reg[3]))
      DT$CNVR_ID[i] <- paste0(arm, "-", n)
      n <- n + 1
    }
    # at least one cnvr overlap with cnv under examination
    else {
      for (r in cnvr_m$CNVR_ID) {
        reg_list <- get_regions_list(DT[CNVR_ID == r, ], prop)
        overlaps <- pmin(my_reg[3], reg_list[[3]]) - pmax(my_reg[2], reg_list[[2]]) + 1
        # reciprocal overlaps
        if (all(overlaps >= my_reg[4]) &
            ### sapply -> vapply
            all(sapply(overlaps, function(x) x >= reg_list[[4]]))) {
          DT$CNVR_ID[i] <- r
          # update cnvrs boundaries
          cnvrs[CNVR_ID == r, `:=` (start = min(start, my_reg[2]), end = max(end, my_reg[3]))]
          bb <- TRUE
        }
      if (bb) break
      }

      # cnvr not set mean no match with candidates cnvrs
      if (!bb) {
        cnvrs <- rbind(cnvrs, data.table("CNVR_ID" = paste0(arm, "-", n),
                                         "chr" = reg_arm[1], "start" = my_reg[2],
                                         "end" = my_reg[3]))
        DT$CNVR_ID[i] <- paste0(arm, "-", n)
        n <- n + 1
      }
    }
  }

  return(list(cnvrs, DT, n))
}

merge_cnvrs <- function(cnvrs, DT, n, prop, arm, reg_arm) {

  for (i in 1:nrow(cnvrs)) {
    b <- FALSE
    my_reg <- get_region_with_CNVR_ID(cnvrs[i], prop)
    # search compatible cnvrs
    #     cnvrs_m <-
    #       cnvrs[start <= my_reg[[1]][3] & end >= my_reg[[1]][2], ]
    cnvrs_m <-
      cnvrs[between(start, my_reg[[1]][2], my_reg[[1]][3], incbounds = FALSE) |
            between(end, my_reg[[1]][2], my_reg[[1]][3], incbounds = FALSE), ]

    # if there are no hits skip
    if (nrow(cnvrs_m) == 0 ) next

    # if there is a match, check if two CNVRs can be combined, what comes first get merged
    # in any case there will be at least another loop
    for (ii in 1:nrow(cnvrs_m)) {

      my_reg_m <- get_region_with_CNVR_ID(cnvrs_m[ii], prop)
      overl <-
        min(my_reg[[1]][3], my_reg_m[[1]][3]) - max(my_reg[[1]][2], my_reg_m[[1]][2])

      if (overl >= my_reg[[1]][4] & overl >= my_reg_m[[1]][4]) {

        # if two CNVRs have the desired overlap, update CNVR, remove the old ones and update DT
        cnvrs <- rbind(cnvrs,
                       data.table("CNVR_ID" = paste0(arm, "-", n), "chr" = reg_arm[1],
                                  "start" = min(my_reg[[1]][2], my_reg_m[[1]][2]),
                                  "end" = max(my_reg[[1]][3], my_reg_m[[1]][3])))
        cnvrs <- cnvrs[!CNVR_ID %in% c(my_reg[[2]], my_reg_m[[2]])]
        DT[CNVR_ID %in% c(my_reg[[2]], my_reg_m[[2]]), CNVR_ID := paste0(arm, "-", n)]
        n <- n + 1
        message("CNVR updated")
        # at this point all the for loops are interrupted
        b <- TRUE
      } # fi
      if (b) break
    }
    if (b) break
  }
  return(list(cnvrs, DT, n, b))
}

dupl_cnvrs <- function(cnvrs, DT) {

  for (i in 1:nrow(cnvrs)) {
    b <- FALSE
    my_reg <- get_region_with_CNVR_ID(cnvrs[i])
    # search compatible cnvrs
    cnvrs_m <- cnvrs[chr == my_reg[[1]][1] & CNVR_ID != my_reg[[2]] &
                     start == my_reg[[1]][2] & end == my_reg[[1]][3], ]

    # if there are no hits skip
    if (nrow(cnvrs_m) == 0 ) next

    # if there is a match, the CNVRs can be combined
    cnvrs <- cnvrs[!CNVR_ID %in% cnvrs_m$CNVR_ID,]
    DT[CNVR_ID %in% cnvrs_m$CNVR_ID, CNVR_ID := my_reg[[2]]]
    message("CNVRs updated")
    b <- TRUE
    break
  }
  return(list(cnvrs, DT, b))
}

check_cnvrs <- function(cnvrs, DT, n, prop, arm, reg_arm) {

  for (i in 1:nrow(cnvrs)) {
    b <- FALSE
    my_reg <- get_region_with_CNVR_ID(cnvrs[i], prop)
    # search compatible cnvrs
    cnvrs_m <-
      cnvrs[between(start, my_reg[[1]][2], my_reg[[1]][3], incbounds = FALSE) |
            between(end, my_reg[[1]][2], my_reg[[1]][3], incbounds = FALSE), ]
          #     cnvrs_m <-
          #       cnvrs[start <= my_reg[[1]][3] & end >= my_reg[[1]][2], ]

    # if there are no hits skip
    if (nrow(cnvrs_m) == 0 ) next

    # if there is a match, check the relative cnvs
    for (ii in 1:nrow(cnvrs_m)) {

      my_reg_m <- get_region_with_CNVR_ID(cnvrs_m[ii], prop)
      overl <-
        min(my_reg[[1]][3], my_reg_m[[1]][3]) - max(my_reg[[1]][2], my_reg_m[[1]][2])
      pass_check <- TRUE

      if (overl >= my_reg[[1]][4] & overl >= my_reg_m[[1]][4]) {

        reg_list <-
          get_regions_list(DT[CNVR_ID %in% c(my_reg[[2]], my_reg_m[[2]]), ], prop)
        cnvs_tmp <- DT[CNVR_ID %in% c(my_reg[[2]], my_reg_m[[2]]), ]

        # test the CNVs one against all the other
        for (iii in 1:nrow(cnvs_tmp)) {
          overlaps <- pmin(cnvs_tmp$end[iii], cnvs_tmp$end) -
            pmax(cnvs_tmp$start[iii], cnvs_tmp$start) + 1

          # reciprocal overlaps
          if (!(all(overlaps >= (cnvs_tmp$end[iii]-cnvs_tmp$start[iii]+1) * prop) &
                all(sapply(overlaps, function(x) x >= reg_list[[4]])))) {
            pass_check <- FALSE
            break
          }
        }
        if (pass_check) {
          # if all cnvs are compatible, update CNVR, remove the old ones and update DT
          nnew <- paste0(arm, "-", n)
          cnvrs <- rbind(cnvrs,
                         data.table("CNVR_ID" = nnew, "chr" = reg_arm[1],
                                    "start" = min(my_reg[[1]][2], my_reg_m[[1]][2]),
                                    "end" = max(my_reg[[1]][3], my_reg_m[[1]][3])))
          cnvrs <- cnvrs[!CNVR_ID %in% c(my_reg[[2]], my_reg_m[[2]]),]
          DT[CNVR_ID %in% c(my_reg[[2]], my_reg_m[[2]]), CNVR_ID := nnew]
          n <- n + 1
          message("CNVR updated")
          # at this point all the for loops are interrupted
          b <- TRUE
        }
        if (b) break
      } # fi
      if (b) break
    }
    if (b) break
  }

  return(list(cnvrs, DT, n, b))
}

move_cnvs <- function() {

}

remove_cnvs <- function(DT, prop) {

  for (i in DT[, ix]) {
    my_reg <- get_region(DT[ix == i, ], prop)
    # all cnvs in the cnvr
    reg_list <- get_regions_list(DT[CNVR_ID == DT$CNVR_ID[i]], prop)

    overlaps <- pmin(my_reg[3], reg_list[[3]]) - pmax(my_reg[2], reg_list[[2]]) + 1

    if (!(all(overlaps >= my_reg[4]) &
          ### sapply -> vapply
          all(sapply(overlaps, function(x) x >= reg_list[[4]]))))
      DT$CNVR_ID[i] <- NA_character_
  }

  message(length(is.na(DT$CNVR_ID)[is.na(DT$CNVR_ID) == TRUE]),
      " CNVs removed from the assigned CNVR")

  return(DT)
}

merge75 <- function(cnvrs, cnvs) {

  cnvrs_edit <- copy(cnvrs)
  for (i in 1:nrow(cnvrs)) {
    # c(chr, st, en, len)
    rid <- cnvrs$CNVR_ID[i]
    if (nrow(cnvrs_edit[CNVR_ID == rid,]) == 0 ) next

    reg <- get_region(cnvrs[i])
    cc <- cnvrs_edit[chr == reg[1] & start < reg[3] & end > reg[2] &
                    CNVR_ID !=rid,]

    if (nrow(cc) < 1) next
    comp <- c(cnvrs$CNVR_ID[i])

    for (n in 1:nrow(cc)) {
      overlap <- min(reg[3], cc$end[n]) - max(reg[2], cc$start[n]) +1
      if (overlap >= 0.75*reg[4] & overlap >= 0.75*cc[n, end-start+1])
        comp <- c(comp, cc[n,CNVR_ID])
    }

    # update CNVR_ID
    if (length(comp) >1) {
      if (nrow(cnvrs_edit[CNVR_ID %in% comp,]) < 2) next
      nst <- min(cnvrs_edit[CNVR_ID %in% comp, start])
      nsp <- max(cnvrs_edit[CNVR_ID %in% comp, end])
      nfr <- sum(cnvrs_edit[CNVR_ID %in% comp, freq])
      newline <- cnvrs_edit[CNVR_ID == comp[1],]
      newline[,`:=` (start = nst, end = nsp, freq = nfr)]
      cnvrs_edit <- rbind(cnvrs_edit[!CNVR_ID %in% comp,], newline)
      # update also CNVs
      cnvs[CNVR_ID %in% comp, CNVR_ID := comp[1]]
    }
    if (cnvrs[, sum(freq)] != cnvrs_edit[, sum(freq)]) {
      message("Ouch!")
      break
    }
  }

  return(list(cnvrs_edit, cnvs))
}
