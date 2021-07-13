# PRE-PROCESSING #

# ...

# STEP 1 #

# default values
a <- 0.3
b <- 0.01
c <- -0.015
d <- -c
# for the function step1() only few columns are necessary
dt <- data.table(locus = c("a", "b", "c", "a", "b", "c", "a", "b", "c"),
                 LRRSD = c(a+0.001, a, rep.int(0, 7)),
                 BAFdrift = c(0, 0, b+0.001, b, rep.int(0, 5)),
                 GCWF = c(0, 0, 0, 0, c-0.001, c, d, d+0.001, 0),
                 excl = -1)
dt <- QCtreeCNV:::step1(dt, a, b, c, d)

test_that("Step 1 behaves correctly.", {
# LRRSD higher than the threshold result in st1 col value = 1
            expect_equal(dt$st1[1], 1)
# ("BAFdrift higher than the threshold result in st1 col value = 1"
            expect_equal(dt$st1[3], 1)
# ("GCWF outside boundaries result in st1 col value = 1"
            expect_equal(dt$st1[c(5,8)], c(1,1))
# ("All three inside the respective boundaries result in st1 col value = 0"
            expect_equal(dt$st1[9], 0)
# ("Borders values results in st1 = 0"
            expect_equal(dt$st1[c(2,4,6,7)], c(0,0,0,0))
# ("When st1 = 1 also excl = 1"
            expect_equal(unique(dt[st1 == 1, excl]), 1)
# ("When st1 = 0 excl = -1"
            expect_equal(unique(dt[st1 == 0, excl]), -1)
})


# STEP 2 #

# only locus columns is needed in the CNV object for this function
dtc <- data.table(locus = c(rep.int("a", 10), rep.int("b", 15),
                            rep.int("c", 20), rep.int("d", 30)))
# simplified case where each CNVR is in a separate chromosome
dtl <- data.table(locus = c("a", "b", "c", "d"), chr = c(1, 2, 3, 4),
                  start = c(1, 10, 100, 1000), end = c(1001, 2010, 3100, 5000))
dtr <- data.table(CNVR_ID = c("a1", "b1", "b2", "c1", "c2", "d1", "d2", "d3", "d4"),
                  chr= c("..."))

test_that("CNVRs sorting behave as expected", {
            # ...
})


# - st1 must be either 1 or 0, only those with st1 = 0 will be considered
# - cnvs belonging to CVNR in the first group, rr[[1]] will have st2 = 1 and exlc = 0,
#   the others will have st2 = 0 and excl = -1
# - cnvs with st = 1 will have st2 = -1
rr <- list(c("a", "b"), "c", "d")
dt <- data.table(st1 = c(0, 0, 0, rep.int(1,7)),
                 CNVR_ID = c(rr[[1]], rep_len(c(rr[[2]],rr[[3]]), 7),
                             rr[[1]][1]),
                 excl = -1)
dt <- QCtreeCNV:::step2(dt, rr)

test_that("Step 2 behave as expected", {
            # st1 = 0 , CNV_ID in rr[[1]]
            expect_equal((unique(dt$st2[1:2])), 1)
            expect_equal((unique(dt$excl[1:2])), 0)
            # st1 = 1
            expect_equal((unique(dt$st2[4:10])), 0)
            expect_equal((unique(dt$excl[4:10])), -1)
            # st1 = 0 , CNV_ID !in rr[[1]]
            expect_equal(dt$st2[3], 0)
            expect_equal(dt$excl[3], 0)
})



# STEP 3 #

# - st2 must be either 1 or 0, only those with st2 = 0 will be considered
# - cnvs belonging to CNVR in the second group will have st3 = 1
# - cnvs belonging to CNVR in the third group will have st3 = 0
# - cnvs with st2 = 1 will have st3 = -1
# - all cnvs with st2 = 0 will have excl = -1

rr <- list("a", c("b", "c"), "d")
dt <- data.table(st2 = c(0, 0, 0, 0, rep.int(1,6)),
                 CNVR_ID = c(rr[[2]], rr[[3]],
                             rep_len(rr[[1]], 6), rr[[3]]))
dt <- QCtreeCNV:::step3(dt, rr)

test_that("Step 3 behave as expected", {
            # st2 = 0, CNVR_ID in rr[[2]]
            expect_equal((unique(dt$st2[1:2])), 1)
            # st2 = 0, CNVR_ID in rr[[3]]
            expect_equal(dt$st2[3], 1)
            # st2 = 0, CNVR_ID in rr[[1]]
            expect_equal(dt$st2[4], 0)
            # st2 = 1
            expect_equal((unique(dt$st2[5:10])), -1)
            # st2 = 0
            expect_equal((unique(dt$excl[1:5])), -1)
})


# STEP 4



# STEP 5


