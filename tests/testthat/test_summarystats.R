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

promoters <- GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2")), c(2, 1)),
                     ranges = IRanges(c(9, 14, 11), width = 4:6),
                     strand = S4Vectors::Rle(strand(c("+", "-", "-"))),
                     seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep="")),
                     id = paste0("P", 1:3))

## Annotation by overlaps
annotateInteractions(gi, list(promoter = promoters))

test_that("categoriseInteraction returns correct result", {
  res <- data.frame(category = c("promoter-promoter", "promoter-distal", "distal-distal"),
                    count = c(1, 2, 2), 
                    stringsAsFactors = FALSE)
  expect_equal(categoriseInteractions(gi), res)
  expect_equal(categoriseInteractions(gi, viewpoints = "promoter"), res[1:2,])
  expect_equal(categoriseInteractions(gi, node.classes = "promoter"), res[1,])
})

test_that("plotSummaryStats works as previously", {
  expect_equal_to_reference(plotSummaryStats(gi), 
                            file = "plotsummarystats.rds")
  resetAnnotations(gi)
  expect_equal_to_reference(plotSummaryStats(gi), 
                            file = "plotsummarystats_reset.rds")
})

