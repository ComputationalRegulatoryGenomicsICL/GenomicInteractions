# Capitalize first letter of string
#
# This function will capitalize the first letter of each string in
# a character vector, and lowercase following letters.
#
# @param x A character vector
#
# @return a string with the first letter capitalised
capitalize = function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
}

#' @importFrom GenomeInfoDb seqinfo seqlevels 'seqlevels<-' seqlengths genome isCircular 

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

