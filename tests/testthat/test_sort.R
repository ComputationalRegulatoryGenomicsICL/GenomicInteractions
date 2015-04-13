
test_that("Sort returns correctly reversed anchors", {
  gi <- new("GenomicInteractions",
            metadata = list(experiment_name="test", description = "this is a test"),
            anchor_one = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1")), c(1, 2, 2)),
                                 ranges = IRanges(1:5, width = 10:6),
                                 strand = S4Vectors::Rle(strand(c("+", "+", "-", "-", "-"))),
                                 seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep=""))),
            anchor_two = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1")), c(2, 2, 1)),
                                 ranges = IRanges(7:3, width = 10:6),
                                 strand = S4Vectors::Rle(strand(c("-", "+", "-", "-", "-"))),
                                 seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep=""))),
            counts = as.integer(1:5),
            elementMetadata = new("DataFrame", nrows = as.integer(5))
  )
  
  expect <- new("GenomicInteractions",
                metadata = list(experiment_name="test", description = "this is a test"),
                anchor_one = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1")), c(2, 1, 2)),
                                     ranges = IRanges(start = c(1,6,3,4,3), end = c(10,14,10,10,8)),
                                     strand = S4Vectors::Rle(strand(c("+", "+", "-", "-", "-"))),
                                     seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep=""))),
                anchor_two = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1")), c(1, 3, 1)),
                                     ranges = IRanges(start = c(7,2,5,4,5), end = c(16,10,12,10,10)),
                                     strand = S4Vectors::Rle(strand(c("-", "+", "-", "-", "-"))),
                                     seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep=""))),
                counts = as.integer(1:5),
                elementMetadata = new("DataFrame", nrows = as.integer(5))
  )
  
  expect_equal(sort(gi, order.interactions = FALSE), expect)
  
})


test_that("Sort returns correctly reversed and ordered anchors", {
  gi <- new("GenomicInteractions",
            metadata = list(experiment_name="test", description = "this is a test"),
            anchor_one = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1")), c(1, 2, 2)),
                                 ranges = IRanges(1:5, width = 10:6),
                                 strand = S4Vectors::Rle(strand(c("+", "+", "-", "-", "-"))),
                                 seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep=""))),
            anchor_two = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1")), c(2, 2, 1)),
                                 ranges = IRanges(7:3, width = 10:6),
                                 strand = S4Vectors::Rle(strand(c("-", "+", "-", "-", "-"))),
                                 seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep=""))),
            counts = as.integer(1:5),
            elementMetadata = new("DataFrame", nrows = as.integer(5))
  )
  
  expect <- new("GenomicInteractions",
                metadata = list(experiment_name="test", description = "this is a test"),
                anchor_one = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2")), c(4, 1)),
                                     ranges = IRanges(start = c(1,6,3,4,3), end = c(10,14,8,10,10)),
                                     strand = S4Vectors::Rle(strand(c("+", "+", "-", "-", "-"))),
                                     seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep=""))),
                anchor_two = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr2")), c(1,1,1,2)),
                                     ranges = IRanges(start = c(7,2,5,4,5), end = c(16,10,10,10,12)),
                                     strand = S4Vectors::Rle(strand(c("-", "+", "-", "-", "-"))),
                                     seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep=""))),
                counts = as.integer(c(1,2,5,4,3)),
                elementMetadata = new("DataFrame", nrows = as.integer(5))
  )
  
  expect_equal(sort(gi, order.interactions = TRUE), expect)
  
})