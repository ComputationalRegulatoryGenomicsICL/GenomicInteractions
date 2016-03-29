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
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, counts=interaction_counts)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors DataFrame setValidity2
#' @import InteractionSet
#' @import methods
#'
#' @exportClass GenomicInteractions

setClass("GenomicInteractions", 
         contains="GInteractions",
         representation(
           anchor1="integer",
           anchor2="integer",
           regions="GRanges",
           NAMES="characterORNULL"
         ),
         prototype(
           anchor1=integer(0),
           anchor2=integer(0),
           regions=GRanges(),
           elementMetadata = DataFrame(counts = integer(0)),
           NAMES=NULL
         )
)

setValidity2("GenomicInteractions", function(object) {
  if(!("counts" %in% names(object@elementMetadata))){
    stop("Valid GenomicInteractions object must contain counts")
  }
  if (!is.integer(object@elementMetadata$counts) || 
      !(all(is.finite(object@elementMetadata$counts))) ||
      !all(object@elementMetadata$counts >=0)){
    stop("'counts' must be finite positive integers")
  }
  return(TRUE)
})


#' Function to create a GenomicInteractions object
#'
#' Create GenomicInteractions objects from two GRanges ojects.
#'
#' @param anchor1,anchor2 GRanges objects.
#' @param counts An integer vector, defaults to 1.
#' @param ... Additional data to be added to mcols
#' @return a GenomicInteractions object
#' @examples
#' 
#' library(GenomicRanges)
#'
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, counts=interaction_counts)
#'
#' @export
setGeneric("GenomicInteractions", function(anchor1, anchor2, counts, ...){ standardGeneric("GenomicInteractions")})

#' @rdname GenomicInteractions
setMethod("GenomicInteractions", c("GRanges", "GRanges", "numeric"), 
          function(anchor1, anchor2, counts, ...){
            out <- GInteractions(anchor1, anchor2, counts = counts, ...)
            class(out) <- "GenomicInteractions"
            out
          })

#' @rdname GenomicInteractions
setMethod("GenomicInteractions", c("GInteractions"), 
          function(anchor1){
            out <- anchor1
            
            if(!("counts" %in% names(out@elementMetadata))){
              out@elementMetadata$counts <- rep(as.integer(1), length(out))
            } else {
              out@elementMetadata$counts <- as.integer(out@elementMetadata$counts)
            }
            class(out) <- "GenomicInteractions"
            out
          })

#' @rdname GenomicInteractions
setMethod("GenomicInteractions", c("GInteractions", "numeric"), 
          function(anchor1, anchor2){
            out <- anchor1
            out$counts <- as.integer(anchor2)
            class(out) <- "GenomicInteractions"
            out
          })

#' @rdname GenomicInteractions
setMethod("GenomicInteractions", c("numeric", "numeric", "GRanges"),
          function(anchor1, anchor2, counts, ...){
            GenomicInteractions(GInteractions(anchor1, anchor2, counts, ...))
          })

#' @rdname GenomicInteractions
setMethod("GenomicInteractions", c("GRanges", "GRanges", "GenomicRangesORmissing"),
          function(anchor1, anchor2, counts, ...){
            if(missing(counts)){
              GenomicInteractions(GInteractions(anchor1, anchor2, ...))
            } else {
              GenomicInteractions(GInteractions(anchor1, anchor2, counts, ...))
            }
          })

#' @rdname GenomicInteractions
setMethod("GenomicInteractions", c("missing", "missing", "GenomicRangesORmissing"),
          function(anchor1, anchor2, counts, ...){
            if(missing(counts)){
              GenomicInteractions(GInteractions())
            } else {
              GenomicInteractions(GInteractions(integer(0), integer(0), counts))
            }
          })

#' @rdname GenomicInteractions
setMethod("GenomicInteractions", c("ANY", "ANY", "ANY"), 
          function(anchor1, anchor2, counts, ...){
            GenomicInteractions(GInteractions(anchor1, anchor2, counts, ...))
          })

###############################################################################
#' updateObject method for GenomicInteractions 1.3.7 and earlier
#' 
#' @inheritParams BiocGenerics::updateObject
#' @return A GenomicInteractions object
#' @importFrom Biobase updateObject
#' @export

setMethod("updateObject", signature(object="GenomicInteractions"),
          function(object, ..., verbose = FALSE){
            if (verbose){message("updating GenomicInteractions object")}
            
            anchor1 <- object@anchor_one
            anchor2 <- object@anchor_two
            all_anchors <- unique(c(object@anchor_one, object@anchor_two))
            mcols(anchor1) <- NULL
            mcols(anchor2) <- NULL
            
            em <- DataFrame(counts = object@counts, object@elementMetadata)
            
            newobj <- GInteractions(anchor1, anchor2, 
                                    all_anchors,
                                    metadata = object@metadata,
                                    elementMetadata = em)
            
            class(newobj) <- "GenomicInteractions"
            names(mcols(newobj)) <- gsub("elementMetadata.", "", names(mcols(newobj)))
            newobj
          })
