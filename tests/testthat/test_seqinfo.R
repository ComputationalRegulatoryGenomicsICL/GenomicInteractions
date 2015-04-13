make_gi <- function() {
  new("GenomicInteractions",
      metadata = list(experiment_name="test", description = "this is a test"),
      anchor_one = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                           ranges = IRanges(1:10, width = 10:1),
                           strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                           seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""), 
                                             seqlengths = rep(1000, 3))),
      anchor_two = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                           ranges = IRanges(101:110, width = 10:1),
                           strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                           seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""),
                           seqlengths = rep(1000, 3))),
      counts = as.integer(10:1),
      elementMetadata = new("DataFrame", nrows = as.integer(10))
  )
}

test_that("seqinfo and associated methods work", {
  gi <- make_gi()
  
  expect_seqinfo <- Seqinfo(seqnames = paste("chr", 1:3, sep=""), 
                            seqlengths = rep(1000, 3))
  expect_seqlevels <- paste("chr", 1:3, sep="")
  expect_seqlengths <- rep(1000, 3)
  names(expect_seqlengths) <- paste("chr", 1:3, sep="")

  expect_equal(seqinfo(gi), expect_seqinfo)
  expect_equal(seqlevels(gi), expect_seqlevels)
  expect_equal(seqlengths(gi), expect_seqlengths)
})