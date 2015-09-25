test_that("Constructor function throws errors with invalid inputs", {
  fn <- file.path(system.file(package="GenomicInteractions"), "extdata", "k562.rep1.cluster.pet3+.txt")
  
  expect_error(GenomicInteractions(fn, type="chiapet.tool", experiment_name="k562", description="k562 pol2 8wg16"), 
               "Anchors must be GRanges objects")
  
  anchor_one <- GRanges()
  anchor_two <- GRanges()
  
  expect_error(GenomicInteractions(anchor_one, anchor_two, counts = 1.5), 
               "counts must contain integer values")
})

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

drop_strand <- function(gi){
  strand(gi@anchor_one) <- "*"
  strand(gi@anchor_two) <- "*"
  return(gi)
  }

test_that("bed12 export/import is consistent", {
  tmp <- tempfile()
  export.bed12(gi, fn = tmp)
  #attributes dropped on export
  expect_equivalent(sort(drop_strand(gi)), #bed12 cannot store strand so is dropped
               sort(makeGenomicInteractionsFromFile(tmp, type = "bed12")))
  unlink(tmp)
})

test_that("bedpe export/import is consistent", {
  tmp <- tempfile()
  export.bedpe(gi, fn = tmp)
  #attributes dropped on export
  expect_equivalent(sort(gi), sort(makeGenomicInteractionsFromFile(tmp, type = "bedpe")))
  unlink(tmp)
  #fails, anchor strands don't match after import!
})


test_that("homer import is consistent with previous behaviour", {
  fn <- system.file("extdata", "Seitan2013_WT_100kb_interactions.txt", 
                    package="GenomicInteractions")
  
  expect_equal_to_reference(makeGenomicInteractionsFromFile(fn, type = "homer"), 
                            file = "importhomer.rds")
})

test_that("chiapet tool import is consistent with previous behaviour", {
  fn <- system.file("extdata/k562.rep1.cluster.pet3+.txt", 
                    package="GenomicInteractions")
  
  expect_equal_to_reference(makeGenomicInteractionsFromFile(fn, type = "chiapet.tool"), 
                            file = "importchiapettool.rds")
})