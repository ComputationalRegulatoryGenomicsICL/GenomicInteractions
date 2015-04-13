## test validObject method

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

#make_no_anchors <- function() {
  #new("GenomicInteractions",
      #metadata = list(experiment_name="test", description = "this is a test"),
      #anchor_one = GRanges(seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
      #anchor_two = GRanges(seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
      #counts = as.integer(10:1),
      #elementMetadata = new("DataFrame")
  #)
#}

make_unequal_anchors <- function() {
    new("GenomicInteractions",
        metadata = list(experiment_name="test", description = "this is a test"),
        anchor_one = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 5)),
                             ranges = IRanges(1:11, width = 11:1),
                             strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 3)),
                             seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
        anchor_two = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                             ranges = IRanges(101:110, width = 10:1),
                             strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                             seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
        counts = as.integer(10:1),
        elementMetadata = new("DataFrame", nrows = as.integer(11))
    )
  }

make_diff_seqinfo <- function() {
  new("GenomicInteractions",
      metadata = list(experiment_name="test", description = "this is a test"),
      anchor_one = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr2")), c(1, 3, 2, 4)),
                           ranges = IRanges(1:10, width = 10:1),
                           strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                           seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep=""))),
      anchor_two = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                           ranges = IRanges(101:110, width = 10:1),
                           strand = S4Vectors::Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                           seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep=""))),
      counts = as.integer(10:1),
      elementMetadata = new("DataFrame", nrows = as.integer(10))
  )
}

make_bad_counts <- function() {
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
      counts = as.integer(20:1),
      elementMetadata = new("DataFrame", nrows = as.integer(10))
  )
}


make_neg_counts <- function() {
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
      counts = as.integer(8:-1),
      elementMetadata = new("DataFrame", nrows = as.integer(10))
  )
}




test_that("Invalid objects are invalid", {
  #expect_error(make_no_anchors(), "anchor one cannot be of length 0")
  expect_error(make_unequal_anchors(), "length of anchor one and anchor two do not match")
  expect_error(make_diff_seqinfo(), "seqinfo must be indentical for both GRanges")
  expect_error(make_bad_counts(), "length of counts must be the same as the anchors")
  expect_error(make_neg_counts(), "Counts should contain only non-negative integers")
  
})
