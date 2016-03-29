gi <- new("GenomicInteractions",
          metadata = list(experiment_name="test", description = "this is a test"),
          anchor1 = as.integer(c(1, 7, 8, 4, 5)),
          anchor2 = as.integer(c(6, 2, 10, 9, 3)),
          regions = GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2")), c(6, 4)),
                            ranges = IRanges(start = c(1,6,3,4,5,7,2,3,4,5),
                                             width = c(10,9,6,7,6,10,9,8,7,8)),
                            strand = S4Vectors::Rle(c("+", "-", "+", "-"), c(2,4,1,3)),
                            seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep=""))),
          elementMetadata = DataFrame(counts = 1:5))
promoters <- GRanges(seqnames = S4Vectors::Rle(factor(c("chr1", "chr2")), c(2, 1)),
                     ranges = IRanges(c(9, 14, 11), width = 4:6),
                     strand = S4Vectors::Rle(strand(c("+", "-", "-"))),
                     seqinfo = Seqinfo(seqnames = paste("chr", 1:2, sep="")),
                     id = paste0("P", 1:3))

## Annotation by overlaps
suppressMessages(annotateInteractions(gi, list(promoter = promoters)))

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

