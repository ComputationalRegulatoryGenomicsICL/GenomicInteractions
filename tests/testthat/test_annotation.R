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
  
  expect_error(annotateAnchors(gi, 3, "signal", c(100, 102, 2, 320, 40)),
               "anchor is neither 1 or 2")

})

## test summarisation

test_that("summariseByFeatures results are the same as before", {
  expect_equal_to_reference(summariseByFeatures(gi, promoters, feature.name = "promoter"),
                            file = "summariseByFeatures.rds")
})

test_that("summariseByFeaturePairs results are the same as before", {
  expect_equal_to_reference( summariseByFeaturePairs(gi, features.one = promoters, features.two = promoters, 
                                                     feature.name.one = "promoter", feature.name.two = "promoter"),
                            file = "summariseByFeaturePairs.rds")
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

## Distance calculations

test_that("distance calculations work as expected", {
  expect_equal(calculateDistances(gi),
               c(5, NA, 1, NA, 1))
  expect_equal(calculateDistances(gi, method = "midpoint"),
               c(5, NA, 1, NA, 1))
  expect_equal(calculateDistances(gi, method = "inner"),
               c(0, NA, 0, NA, 0))
  expect_warning(calculateDistances(gi, method = "inner"),
                 "setting negative distances to 0, this is due to the presence of overlapping anchors in your dataset")
  expect_equal(calculateDistances(gi, method = "outer"),
               c(14, NA, 8, NA, 6))
  expect_error(calculateDistances(gi, method = "my_method"))
})

one.df <- as.data.frame(anchorOne(gi))
two.df <- as.data.frame(anchorTwo(gi))

test_that("distance calculations on dataframes work as expected", {
  expect_equal(GenomicInteractions:::.calculateDistances.df(one.df, two.df),
               calculateDistances(gi))
  expect_equal(GenomicInteractions:::.calculateDistances.df(one.df, two.df, method = "midpoint"),
               calculateDistances(gi, method = "midpoint"))
  expect_equal(GenomicInteractions:::.calculateDistances.df(one.df, two.df, method = "inner"),
               calculateDistances(gi, method = "inner"))
  expect_equal(GenomicInteractions:::.calculateDistances.df(one.df, two.df, method = "outer"),
               calculateDistances(gi, method = "outer"))
  })

