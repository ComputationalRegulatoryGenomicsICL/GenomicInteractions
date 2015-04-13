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

test_that("Annotate interactions returns expected results", {
  expect_node_class_one <- c("promoter", rep("distal", 4))
  expect_node_class_two <- c(rep("promoter", 3), rep("distal", 2))
  
  expect_ids_one <- c("P1", as.list(rep(NA, 4)))
  expect_ids_two <- c(c("P2", "P1","P3", as.list(rep(NA, 2))))
  
  expect_equal(anchorOne(gi)$node.class, expect_node_class_one)
  expect_equal(anchorTwo(gi)$node.class, expect_node_class_two)
  
  expect_equal(anchorOne(gi)$promoter.id, expect_ids_one)
  expect_equal(anchorTwo(gi)$promoter.id, expect_ids_two)
  })

## Annotation by setting mcols

annotateAnchors(gi, 1, "signal", c(100, 102, 2, 320, 40))
annotateAnchors(gi, 2, "signal", c(150, 62, 300, 0, 56))

test_that("Annotate anchors returns expected results", {
  expect_ann_one <- c(100, 102, 2, 320, 40)
  expect_ann_two <-  c(150, 62, 300, 0, 56)

  expect_true("signal" %in% names(mcols(anchorOne(gi))))
  expect_true("signal" %in% names(mcols(anchorTwo(gi))))
  expect_equal(anchorOne(gi)$signal, expect_ann_one)
  expect_equal(anchorTwo(gi)$signal, expect_ann_two)

})
## Resetting annotation
resetAnnotations(gi)

test_that("resetting annotations removes mcols", {
  expect_false("signal" %in% names(mcols(anchorOne(gi))))
  expect_false("signal" %in% names(mcols(anchorTwo(gi))))
  
  expect_false("node.class" %in% names(mcols(anchorOne(gi))))
  expect_false("node.class" %in% names(mcols(anchorTwo(gi))))
  
  expect_true(length(mcols(anchorOne(gi)))==0)
  expect_true(length(mcols(anchorTwo(gi)))==0)
})

## Distance calculations?

