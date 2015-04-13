#' Subset a GenomicInteractions object by features
#'
#' Subsets interactions for which at least one of the anchors overlaps with a given GRanges object.
#' Alternatively, subsets interactions based on annotated feature IDs for a particular feature.
#'
#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @docType methods
#' @param GIObject A GenomicInteractions object
#' @param features A GRanges or GRangesList object, or a character vector containing
#'                      IDs of annotated features, e.g. promoter IDs.
#' @param feature.class If `features' is a character vector, the corresponding feature name, e.g. "promoter".
#' @return a subsetted GenomicInteractions object
#' @export
setGeneric("subsetByFeatures",function(GIObject, features, feature.class=NULL){standardGeneric ("subsetByFeatures")})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @import GenomicRanges
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "GRanges", "missing"), function(GIObject, features, feature.class=NULL){
    i = unique(c(subjectHits(findOverlaps(features, GIObject@anchor_one)), subjectHits(findOverlaps(features, GIObject@anchor_two))))
    GIObject[i]
})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @import GenomicRanges
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "GRangesList", "missing"), function(GIObject, features, feature.class=NULL){
    i = unique(c(subjectHits(findOverlaps(features, GIObject@anchor_one)), subjectHits(findOverlaps(features, GIObject@anchor_two))))
    GIObject[i]
})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @import GenomicRanges
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "character", "character"), function(GIObject, features, feature.class){
    if(!"node.class" %in% names(elementMetadata(GIObject@anchor_one)) & feature.class %in% unique(c(GIObject@anchor_one$node.class, GIObject@anchor_two$node.class)))
        stop(paste(feature.class," has not been annotated on this GenomicInteractions object"))
    i = sapply(elementMetadata(GIObject@anchor_one)[[paste(feature.class, "id", sep=".")]],
               function(x){ features %in% x }) | sapply(elementMetadata(GIObject@anchor_two)[[paste(feature.class, "id", sep=".")]], function(x){ features %in% x })
    GIObject[i]
})

#' Standard subsetting methods for GenomicInteractions objects
#'
#' @name [
#' @param x A genomicInteractions object
#' @param i A numeric, logical or Rle vector
#'
#' @return A GenomicInteractions object containing only the features specified by `i`.
#' @rdname GenomicInteractions-subset-methods
NULL

#' @name [
#' @aliases [,GenomicInteractions-method
#' @docType methods
#' @rdname GenomicInteractions-subset-methods
#' @import BiocGenerics
#' @export
setMethod(f="[", "GenomicInteractions", function(x, i, j, drop) {
          if (!missing(i)) {
            ans_anchor_one = x@anchor_one[i]
            ans_anchor_two = x@anchor_two[i]
            ans_counts = x@counts[i]
            ans_mcols = mcols(x)[i, ,drop=FALSE]
            x = BiocGenerics:::updateS4(x, anchor_one=ans_anchor_one,
                                        anchor_two=ans_anchor_two,
                                        counts=ans_counts,
                                        elementMetadata=ans_mcols)
        }
        if (!missing(j))
            mcols(x) = mcols(x)[ , j, drop=FALSE]
        return(x)
} )

#' Combine GenomicInteractions Methods
#'
#' This method will fail if the seqlengths of the objects to be combined do not match.
#' If some chromosomes appear in one set of seqinfo but not the other, the seqinfo will
#' be merged.
#' 
#' @name c
#' @param x, ... GenomicInteractions objects to be concatenated
#' @param ignore.mcols Logical, default FALSE, remove mcols in combined object.
#' @param recursive Not supported
#' @return A GenomicInteractions object.
#'
#' @aliases c,GenomicInteractions
#' @docType methods
#' @rdname GenomicInteractions-combine-methods
#' @export
setMethod(f="c", signature="GenomicInteractions", function(x, ..., ignore.mcols=FALSE, recursive=FALSE) {
          if (!identical(recursive, FALSE))
              stop("'recursive' argument not supported")
          if (missing(x))
              args = unname(list(...))
          else
              args = unname(list(x, ...))
          ans_anchor_one=do.call(c, lapply(args, anchorOne)) # does this implicitly check for seqinfo?
          ans_anchor_two=do.call(c, lapply(args, anchorTwo))
          ans_counts=do.call(c, lapply(args, interactionCounts))
          if (ignore.mcols)
              ans_mcols = new("DataFrame", nrows=length(ans_anchor_one))
          else
              ans_mcols = DataFrame(do.call(rbind, lapply(args, mcols, FALSE)))
          new("GenomicInteractions",
              metadata = list(experiment_name="", description=""), # users can set this later
              anchor_one=ans_anchor_one, # does this implicitly check for seqinfo?
              anchor_two=ans_anchor_two,
              counts=ans_counts,
              elementMetadata=ans_mcols)
} )

