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

make_gr <- function() {
  GRangesList(nomatch = GRanges(seqnames = "chr1",
                                ranges = IRanges(start=5, end=10),
                                strand = "+"),
              onematch_one = GRanges(seqnames = "chr3",
                                 ranges = IRanges(start=2, end=7),
                                 strand = "-"),
              twomatch_one = GRanges(seqnames = "chr1",
                                 ranges = IRanges(start=1, end=5),
                                 strand = "-"),
              onematch_two = GRanges(seqnames = "chr1",
                                     ranges = IRanges(start=100, end=105),
                                     strand = "-"),
              matchboth = GRanges(seqnames = "chr1",
                                  ranges = IRanges(start=10, end=105),
                                  strand = "-"))
}

test_that("No overlaps returns list of empty matches (gi, gr)", {
  gi <- make_gi()
  gr <- make_gr()
  expect <- list(one = new("Hits", queryHits = integer(0), subjectHits = integer(0),
                          queryLength = 10L, subjectLength = 1L),
                two = new("Hits", queryHits = integer(0), subjectHits = integer(0),
                          queryLength = 10L, subjectLength = 1L))
  ## select = "all"
  for (type in c("any", "start", "end")){
    expect_equal(findOverlaps(gi, gr$nomatch, type = type, select = "all"), expect)
  }
  
  expect <- list(one = rep(NA_integer_, 10),
                 two = rep(NA_integer_, 10))
  ## select = "first"
  for (type in c("any", "start", "end")){
    expect_equal(findOverlaps(gi, gr$nomatch, type = type, select = "first"), expect)
  }
  
})

test_that("No overlaps returns list of empty matches (gr, gi)", {
  gi <- make_gi()
  gr <- make_gr()
  expect <- list(one = new("Hits", queryHits = integer(0), subjectHits = integer(0),
                           queryLength = 1L, subjectLength = 10L),
                 two = new("Hits", queryHits = integer(0), subjectHits = integer(0),
                           queryLength = 1L, subjectLength = 10L))
  ## select = "all"
  for (type in c("any", "start", "end")){
    expect_equal(findOverlaps(gr$nomatch, gi, type = type, select = "all"), expect) 
  }
  
  ## select = "first"
  expect <- list(one = rep(NA_integer_, 1),
                 two = rep(NA_integer_, 1))
  for (type in c("any", "start", "end")){
    expect_equal(findOverlaps(gr$nomatch, gi, type = type, select = "first"), expect)  
  }
  
  
})

test_that("No / one anchor overlaps returns list of empty matches (gi, gi)", {
  gi1 <- make_gi()
  gi2 <- gi1  
  gi2@anchor_one <- shift(anchorOne(gi2), 1000L)
  gi2@anchor_two <- shift(anchorTwo(gi2), 1000L)
  gi3 <- gi1
  gi3@anchor_two <- shift(anchorTwo(gi3), 1000L)
  
  expect <- new("Hits", queryHits = integer(0), subjectHits = integer(0),
                           queryLength = 10L, subjectLength = 10L)
  
  expect_equal(findOverlaps(gi1, gi2), expect)
  
  #first anchor overlaps, second does not
  expect_equal(findOverlaps(gi1, gi3), expect)

})

test_that("Empty GRanges as query/subject returns list of empty matches",{
  gi <- make_gi()
  ##empty gr
  gr <- GRanges()
  expect <- list(one = new("Hits", queryHits = integer(0), subjectHits = integer(0),
                           queryLength = 0L, subjectLength = 10L),
                 two = new("Hits", queryHits = integer(0), subjectHits = integer(0),
                           queryLength = 0L, subjectLength = 10L))

  for (type in c("any", "start", "end")){
    expect_equal(findOverlaps(gr, gi, type = type, select = "all"), expect)
  }
  
  expect <- list(one = new("Hits", queryHits = integer(0), subjectHits = integer(0),
                           queryLength = 10L, subjectLength = 0L),
                 two = new("Hits", queryHits = integer(0), subjectHits = integer(0),
                           queryLength = 10L, subjectLength = 0L))
  for (type in c("any", "start", "end")){
    expect_equal(findOverlaps(gi, gr, type = type, select = "all"), expect)
  }
})
