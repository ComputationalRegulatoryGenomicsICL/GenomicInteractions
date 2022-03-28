# test InteractionTrack methods
anchor.one <- GRanges(c('chr1', 'chr1', 'chr1', 'chr1'),  IRanges(c(10, 20, 30, 20), width=5))
anchor.two <- GRanges(c('chr1', 'chr1', 'chr1', 'chr2'), IRanges(c(100, 200, 300, 50), width=5))
test <- GenomicInteractions(anchor.one, anchor.two, counts=1:4)

interaction_track1 <- InteractionTrack(name='Test', test, chromosome='chr1')
interaction_track2 <- InteractionTrack(name='Test', test, chromosome='chr1', start = 5, end = 250)
interaction_track3 <- InteractionTrack(name='Test', test[c(1,2,4)], chromosome='chr1')

test_that("InteractionTrack start is correct", {
  expect_equal(start(interaction_track1), 10)
  expect_equal(start(interaction_track2), 5)
})

test_that("InteractionTrack end is correct", {
  expect_equal(end(interaction_track1), 304)
  expect_equal(end(interaction_track2), 250)
})

test_that("InteractionTrack chromosome is correct", {
  expect_equal(chromosome(interaction_track1), "chr1")
})

test_that("InteractionTrack subset works as expected", {
  expect_equal(subset(interaction_track1, from = 5, to = 25, chromosome = "chr1"),
               interaction_track3)
})


