## Definition of GenomicInteractions

#' A S4 class to represent interactions between genomic regions.
#'
#'  @slot metadata List, defaults to "experiment_name" and "description", inherited from S4Vectors::Vector
#'  @slot anchor_one,anchor_two GRanges. Set of anchors of interactions.
#'  @slot counts integer vector, contains raw counts
#'  @slot elementMetadata DataFrame
#'
#' This class is used to store information on which genomic regions are
#' interacting with each other. Objects of this class contain information of
#' the genomic coordinates of the interacting regions and the strength of these
#' interactions, and associated metadata such as the name of the dataset and a
#' brief description of the dataset.  Interacting regions are stored as a pair
#' of GenomicRanges: each set of anchor regions is stored as a separate
#' GenomicRanges object, accessed by \code{getAnchorOne} and
#' \code{getAnchorTwo}.
#'
#' @examples
#'
#' showClass("GenomicInteractions")
#' library(GenomicRanges)
#'
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5)))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5)))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, experiment_name="test", 
#'                            description="this is a test", counts=interaction_counts)
#'
#' @import GenomicRanges
#' @import S4Vectors
#'
#' @export GenomicInteractions
setClass("GenomicInteractions",
    representation(anchor_one = "GRanges",
                   anchor_two = "GRanges",
                   counts = "integer",
                   elementMetadata="DataFrame"),
    prototype(anchor_one = GRanges(),
              anchor_two = GRanges(),
              counts = integer(0),
              elementMetadata = DataFrame() ),
    contains="Vector",
    validity = function(object){
        if(length(object@anchor_one) != length(object@anchor_two)) {
            return("length of anchor one and anchor two do not match")
        } else if(!.isEqualSeqInfo(object@anchor_one, object@anchor_two)) {
            return("seqinfo must be indentical for both GRanges") # this is order-dependent which is not desireable
        } else if(length(object@counts) != length(object@anchor_one)) {
            return("length of counts must be the same as the anchors")
        } else if(any(object@counts < 0)) {
            return("Counts should contain only non-negative integers")
        } else{ return(TRUE)}}
)

#' Function to create a GenomicInteraction object
#'
#' Create GenomicInteraction objects from two GRanges ojects.
#'
#' @param anchor_one, anchor_two GRanges objects.
#' @param counts An integer vector, defaults to 1.
#' @param experiment_name Experiment name.
#' @param description Description of experiment.
#' @param ... Additional data to be added to mcols
#' @return a GenomicInteractions object
#'
#' @examples
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5), seqlengths=seqlengths(Mmusculus))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5), seqlengths=seqlengths(Mmusculus))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, experiment_name="test", 
#'                            description="this is a test", counts=interaction_counts)
#'
#' @import GenomeInfoDb
#' @export
GenomicInteractions = function(anchor_one=GRanges(), anchor_two=GRanges(), 
                               counts=integer(), experiment_name=NULL, description=NULL, ...) {
    if (class(anchor_one)!="GRanges" || class(anchor_two)!="GRanges"){
      stop("Anchors must be GRanges objects")
    }
    if (!all(counts == floor(counts)))
        stop("counts must contain integer values")
    if (length(counts) == 1)
        counts = rep(counts, length(anchor_one))
    mcols = DataFrame(...)
    if (ncol(mcols) == 0L)
        mcols = new("DataFrame", nrows = length(anchor_one))
    if (nrow(mcols) == 1L)
        mcols = mcols[rep(1, length(anchor_one)), ]
	if (!.isEqualSeqInfo(anchor_one, anchor_two)) {
        seqinfo_both = merge(seqinfo(anchor_one), seqinfo(anchor_two))
        seqlevels(anchor_one) = seqlevels(seqinfo_both)
        seqinfo(anchor_one) = seqinfo_both
        seqlevels(anchor_two) = seqlevels(seqinfo_both)
        seqinfo(anchor_two) = seqinfo_both
    }
    new("GenomicInteractions",
        metadata=list(experiment_name=experiment_name, description=description),
        anchor_one=anchor_one,
        anchor_two=anchor_two,
        counts=as.integer(counts),
        elementMetadata=mcols)
}

#' Get the length of a GenomicInteractions GIObject
#'
#' @param x GenomicInteractions GIObject
#' @return A numeric vector containing the length of the GIObject
#' @docType methods
#' @export
setMethod(length, "GenomicInteractions", function(x) length(x@anchor_one))

# Quick access to mcols()

setMethod("$", "GenomicInteractions",
    function(x, name) mcols(x)[[name]]
)

setReplaceMethod("$", "GenomicInteractions",
    function(x, name, value) {mcols(x)[[name]] <- value; x}
)

