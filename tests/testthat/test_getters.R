make_gi <- function() {
  new("GenomicInteractions",
      metadata = list(experiment_name="test", description = "this is a test"),
      anchor_one = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                           ranges = IRanges(1:10, width = 10:1),
                           strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                           seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
      anchor_two = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                           ranges = IRanges(101:110, width = 10:1),
                           strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                           seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
      counts = as.integer(10:1),
      elementMetadata = new("DataFrame", nrows = as.integer(10))
  )
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
            anchor_one = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                                 ranges = IRanges(1:10, width = 10:1),
                                 strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                                 seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
                                 node.class = rep(c("distal", "distal", "promoter", "promoter", "promoter"),2)),
            anchor_two = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                                 ranges = IRanges(101:110, width = 10:1),
                                 strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                                 seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
                                 node.class = rep(c("distal", "distal", "promoter", "distal", "promoter"), 2)),
            counts = as.integer(10:1),
            elementMetadata = new("DataFrame", nrows = as.integer(10))
  )
  expect <- "promoter|distal"
  expect_match(annotationFeatures(gi), expect)
})