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

test_that("Anchors are returned correctly", {
  gi <- make_gi()
  expect_anchor_one <- GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                               ranges = IRanges(1:10, width = 10:1),
                               strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                               seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")))
  expect_anchor_two <- GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                               ranges = IRanges(101:110, width = 10:1),
                               strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                               seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")))
  
  expect_equal(anchorOne(gi), expect_anchor_one)
  expect_equal(anchorTwo(gi), expect_anchor_two)
  
})

test_that("Name/description are returned correctly", {
  gi <- make_gi()
  expect_name <- "test"
  expect_description <- "this is a test"
  
  expect_equal(name(gi), expect_name)
  expect_equal(description(gi), expect_description)
})

test_that("counts are returned correctly", {
  gi <- make_gi()
  expect_counts <- as.integer(10:1)
  
  expect_equal(interactionCounts(gi), expect_counts)
})

test_that("Annotation features are returned correctly", {
  gi <- make_gi()
  
  expect_equal(annotationFeatures(gi), NA_character_)
  
  gi <- new("GenomicInteractions",
      metadata = list(experiment_name="test", description = "this is a test"),
      anchor1 = as.integer(c(1,  7,  8,  9,  2,  3, 13, 14, 15, 16)),
      anchor2 = as.integer(c(4, 10, 11, 12,  5,  6, 17, 18, 19, 20)),
      regions = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr3")), c(6,6,8)),
                        ranges = IRanges(start = c(1L, 5L, 6L, 101L, 105L, 106L, 2L, 3L, 4L, 102L, 103L, 104L, 
                                                   7L, 8L, 9L, 10L, 107L, 108L, 109L, 110L),
                                         width = c(10L, 6L, 5L, 10L, 6L, 5L, 9L, 8L, 7L, 9L, 8L, 7L, 4L, 3L, 2L, 
                                                   1L, 4L, 3L, 2L, 1L)),
                        strand = S4Vectors::Rle(strand(c("-", "+", "-")), c(6, 6, 8)),
                        seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
                        node.class = rep(c("distal", "distal", "promoter", "distal", "promoter"), 4)),
      elementMetadata = DataFrame(counts = 10:1))
  expect <- "promoter|distal"
  expect_match(annotationFeatures(gi), expect)
})