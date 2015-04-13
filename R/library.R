.loadGenome = function(genome){
  require(genome, character.only=TRUE)
  genome.name = unlist(strsplit(genome, split='\\.'))
  return(get(genome.name[2]))
}

#' Capitalize first letter of string
#'
#' This function will capitalize the first letter of each string in
#' a character vector, and lowercase following letters.
#'
#' @param x A character vector
#'
#' @return a string with the first letter capitalised
capitalize = function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
}


.pasteAnchor = function(x) {
    paste(paste(seqnames(x), start(x), sep=":"), end(x), sep="..")
}

.generateInteractionName = function(anchorString1, anchorString2, counts){
    paste(paste(anchorString1, anchorString2, sep="-"), counts, sep=",")
}

.isSortedStart = function(x) {
    (start(anchorOne(x)) < start(anchorTwo(x))) | is.trans(x)
}

.isSortedChrom = function(x) {
    as.numeric(seqnames(anchorOne(x))) <= as.numeric(seqnames(anchorTwo(x)))
}

.isSorted = function(x) {
    .isSortedStart(x) & .isSortedChrom(x)
}

.isEqualSeqInfo = function(one, two) {
    seqinfo.one = seqinfo(one)
    seqinfo.two = seqinfo(two)
    if (length(seqinfo.one) != length(seqinfo.two)) {
        value = FALSE
    } else {
        chr = all(seqlevels(seqinfo.one) == seqlevels(seqinfo.two), na.rm=TRUE)
        len = all(seqlengths(seqinfo.one) == seqlengths(seqinfo.two), na.rm=TRUE)
        cir = all(isCircular(seqinfo.one) == isCircular(seqinfo.two), na.rm=TRUE)
        gen = all(genome(seqinfo.one) == genome(seqinfo.two), na.rm=TRUE)
        value = chr & len & cir & gen
    }
    return(value)
}

