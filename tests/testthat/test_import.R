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

drop_strand <- function(gi){
  strand(gi@regions) <- "*"
  return(gi)
  }

test_that("asBED works with just cis and just trans interactions", {
  expect_silent(asBED(gi))
  b <- asBED(gi)
  expect_equal(asBED(gi[is.cis(gi)]), b[c(2,4,6)])
  expect_equal(asBED(gi[is.trans(gi)]), b[c(1,3,5,7)])
})

test_that("bed12 export/import is consistent", {
  tmp <- tempfile()
  export.bed12(gi, fn = tmp)
  #attributes dropped on export
  expect_equivalent(sort(drop_strand(gi)), #bed12 cannot store strand so is dropped
               sort(makeGenomicInteractionsFromFile(tmp, type = "bed12")))
  unlink(tmp)
})

test_that("bed12 export/import via rtracklayer is consistent", {
  tmp <- tempfile()
  export(asBED(gi), tmp, format="bed")
  #attributes dropped on export
  expect_equivalent(sort(asBED(gi)), #bed12 cannot store strand so is dropped
               sort(rtracklayer::import(tmp, format="bed")))
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
