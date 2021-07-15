
# getline_locus() #

# no length column
dtl <- data.table(locus = c("a", "b", "c", "d"), chr = c(10, 2, 3, 4),
                  start = c(1, 10, 100, 1000), end = c(1001, 2010, 3100, 5000))
ln1 <- QCtreeCNV:::getline_locus(dtl[1])
ln2 <- QCtreeCNV:::getline_locus(dtl[2])
# wrong length
dtl <- data.table(locus = c("a", "b", "c", "d"), chr = c(10, 2, 3, 4),
                  start = c(1, 10, 100, 1000), end = c(1001, 2010, 3100, 5000),
                  length = c(100, 100, 100, 100))
ln3 <- QCtreeCNV:::getline_locus(dtl[3])
ln4 <- QCtreeCNV:::getline_locus(dtl[4])

test_that("getline_locus internal function" , {
            # ln[1] is locus, ln[2] is chr, ln[3] is start, ln[4] is end
            expect_equal(ln1[1:4], c(dtl$locus[1], dtl$chr[1], dtl$start[1], dtl$end[1]))
            #ln[5] is length.    Length = end - start + 1
            expect_equal(as.numeric(ln1[5]), 1001)
            # if length columns is present is up to the user compute it right
            expect_equal(as.numeric(ln3[5]), 100)
            # ln[6] is center. Center = start + length/2
            expect_equal(as.numeric(ln2[6]), 10 + as.numeric(ln2[5])/2)
            # if length columns is present and not correct also center will
            # be affected
            expect_equal(as.numeric(ln4[6]), 1050)
})


# getline_cnv() #


# getline_cnvr() #


# chr_uniform() #
# this should take care of column name and content (whitespace included)
test_that("chr_uniform() works fine", {

            dt <- QCtreeCNV:::chr_uniform(data.table(Chr = "chr1"))
            expect_equal(colnames(dt), "chr")
            expect_equal(dt$chr, "1")

            dt <- QCtreeCNV:::chr_uniform(data.table(Chromosome = "Chr 1"))
            expect_equal(colnames(dt), "chr")
            expect_equal(dt$chr, "1")

            dt <- QCtreeCNV:::chr_uniform(data.table(chromosome = " chr 11"))
            expect_equal(colnames(dt), "chr")
            expect_equal(dt$chr, "11")

            dt <- QCtreeCNV:::chr_uniform(data.table(chr = "chrX"))
            expect_equal(colnames(dt), "chr")
            expect_equal(dt$chr, "23")

            dt <- QCtreeCNV:::chr_uniform(data.table(chr = "chry"))
            expect_equal(colnames(dt), "chr")
            expect_equal(dt$chr, "24")
})



# uniform_GT_CN() #
# this should take care of CN column name and content, plus create columns GT
test_that("uniform_GT_CN() works fine", {
            #
})

