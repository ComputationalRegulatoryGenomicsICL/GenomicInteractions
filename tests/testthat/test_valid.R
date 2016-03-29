## test validObject method

make_gi <- function() {
  new("GenomicInteractions",
      metadata = list(experiment_name="test", description = "this is a test"),
      anchor1 = as.integer(c(1,  7,  8,  9,  2,  3, 13, 14, 15, 16)),
      anchor2 = as.integer(c(4, 10, 11, 12,  5,  6, 17, 18, 19, 20)),
      regions = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr3")), c(6,6,8)),
                        ranges = IRanges(start = c(1L, 5L, 6L, 101L, 105L, 106L, 2L, 3L, 4L, 102L, 103L, 104L, 
                                                   7L, 8L, 9L, 10L, 107L, 108L, 109L, 110L),
                                         width = c(10L, 6L, 5L, 10L, 6L, 5L, 9L, 8L, 7L, 9L, 8L, 7L, 4L, 3L, 2L, 
                                                   1L, 4L, 3L, 2L, 1L)),
                        strand = S4Vectors::Rle(strand(c("-", "+", "-")), c(6, 6, 8)),
                        seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
      elementMetadata = DataFrame(counts = 10:1))
}

make_unequal_anchors <- function() {
        new("GenomicInteractions",
            metadata = list(experiment_name="test", description = "this is a test"),
            anchor1 = as.integer(c(1,  7,  8,  9,  2,  3)),
            anchor2 = as.integer(c(4, 10, 11, 12,  5,  6, 17, 18, 19, 20)),
            regions = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr3")), c(6,6,8)),
                              ranges = IRanges(start = c(1L, 5L, 6L, 101L, 105L, 106L, 2L, 3L, 4L, 102L, 103L, 104L, 
                                                         7L, 8L, 9L, 10L, 107L, 108L, 109L, 110L),
                                               width = c(10L, 6L, 5L, 10L, 6L, 5L, 9L, 8L, 7L, 9L, 8L, 7L, 4L, 3L, 2L, 
                                                         1L, 4L, 3L, 2L, 1L)),
                              strand = S4Vectors::Rle(strand(c("-", "+", "-")), c(6, 6, 8)),
                              seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
            elementMetadata = DataFrame(counts = 10:1))
  }

make_bad_counts <- function() {
  new("GenomicInteractions",
      metadata = list(experiment_name="test", description = "this is a test"),
      anchor1 = as.integer(c(1,  7,  8,  9,  2,  3, 13, 14, 15, 16)),
      anchor2 = as.integer(c(4, 10, 11, 12,  5,  6, 17, 18, 19, 20)),
      regions = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr3")), c(6,6,8)),
                        ranges = IRanges(start = c(1L, 5L, 6L, 101L, 105L, 106L, 2L, 3L, 4L, 102L, 103L, 104L, 
                                                   7L, 8L, 9L, 10L, 107L, 108L, 109L, 110L),
                                         width = c(10L, 6L, 5L, 10L, 6L, 5L, 9L, 8L, 7L, 9L, 8L, 7L, 4L, 3L, 2L, 
                                                   1L, 4L, 3L, 2L, 1L)),
                        strand = S4Vectors::Rle(strand(c("-", "+", "-")), c(6, 6, 8)),
                        seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
      elementMetadata = DataFrame(counts = 20:1))
}


make_neg_counts <- function() {
  new("GenomicInteractions",
      metadata = list(experiment_name="test", description = "this is a test"),
      anchor1 = as.integer(c(1,  7,  8,  9,  2,  3, 13, 14, 15, 16)),
      anchor2 = as.integer(c(4, 10, 11, 12,  5,  6, 17, 18, 19, 20)),
      regions = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr3")), c(6,6,8)),
                        ranges = IRanges(start = c(1L, 5L, 6L, 101L, 105L, 106L, 2L, 3L, 4L, 102L, 103L, 104L, 
                                                   7L, 8L, 9L, 10L, 107L, 108L, 109L, 110L),
                                         width = c(10L, 6L, 5L, 10L, 6L, 5L, 9L, 8L, 7L, 9L, 8L, 7L, 4L, 3L, 2L, 
                                                   1L, 4L, 3L, 2L, 1L)),
                        strand = S4Vectors::Rle(strand(c("-", "+", "-")), c(6, 6, 8)),
                        seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
      elementMetadata = DataFrame(counts = 8:-1))
}

test_that("Invalid objects are invalid", {
  expect_error(make_unequal_anchors(), "'x@anchor2' is not parallel to 'x'")
  expect_error(make_bad_counts())
  expect_error(make_neg_counts(), "'counts' must be finite positive integers")
  
})
