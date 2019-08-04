gi <- GenomicInteractions(anchor1 = GRanges(seqnames = "chr1", 
                                            ranges = IRanges(start = c(1, 1, 11, 11, 21, 21, 31, 31),
                                                             end = c(10, 10, 20, 20, 30, 20, 40, 40))),
                          anchor2 = GRanges(seqnames = rep(c("chr1", "chr2"), 4),
                                            ranges = IRanges(start = c(11, 11, 1, 11, 1, 21, 21, 1),
                                                             end = c(20, 20, 10, 20, 10, 30, 30, 10)))
                                            )

vp1 <- GenomicInteractions(anchor1 = GRanges(seqnames = "chr1", 
                                             ranges = IRanges(start = rep(1, 7),
                                                              end = rep(20, 7))),
                           anchor2 = GRanges(seqnames = c(rep("chr1", 5), rep("chr2", 2)),
                                             ranges = IRanges(start = c(1, 1, 11, 11, 21, 11, 11),
                                                              end = c(10, 10, 20, 20, 30, 20, 20))))

vp2 <- GenomicInteractions(anchor1 = GRanges(seqnames = "chr1", 
                                             ranges = IRanges(start = rep(1, 5),
                                                              end = rep(20, 5))),
                           anchor2 = GRanges(seqnames = c(rep("chr1", 5)),
                                             ranges = IRanges(start = c(1, 1, 11, 11, 21),
                                                              end = c(10, 10, 20, 20, 30))),
                           regions = GRanges(seqnames = c(rep("chr1", 4), "chr2"),
                                            ranges = IRanges(start = c(1, 1, 11,21, 11),
                                                             end = c(10, 20, 20, 30, 20))))

pos <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1, end = 20))
region <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1, end = 40))

test_that("Viewpoints are as expected", {
  expect_equal(vp1, viewPoint(gi, bait = pos))
  expect_equal(vp2, viewPoint(gi, bait = pos, region))
})

