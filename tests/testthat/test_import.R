test_that("Constructor function throws errors with invalid inputs", {
  fn <- file.path(system.file(package="GenomicInteractions"), "extdata", "k562.rep1.cluster.pet3+.txt")
  
  expect_error(GenomicInteractions(fn, type="chiapet.tool", experiment_name="k562", description="k562 pol2 8wg16"), 
               "Anchors must be GRanges objects")
  
  anchor_one <- GRanges()
  anchor_two <- GRanges()
  
  expect_error(GenomicInteractions(anchor_one, anchor_two, counts = 1.5), 
               "counts must contain integer values")
})