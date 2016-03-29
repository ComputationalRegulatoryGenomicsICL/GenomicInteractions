#' Functions to set data held in a GInteractions object.
#'
#' Use these functions to set data stored in each of the slots of a
#' GInteractions object.
#'
#' @name setters
#' @param GIObject A GenomicInteractions object
#' @param value A vector to replace a slot in the object
#' @return GenomicInteractions object
#' @rdname setters
#' @examples
#'
#' library(GenomicRanges)
#'
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, experiment_name="test", 
#'                            description="this is a test", counts=interaction_counts)
#'
#' name(test) <- "Mouse test"
#' name(test)
#'
#' description(test) <- "This is a test using the mouse genome"
#' description(test)
#'
#' interactionCounts(test) <- c(2,3,8,5)
#' interactionCounts(test)
#'
## GENERICS
#' @rdname setters
#' @aliases name<-
#' @export
setGeneric("name<-",function(GIObject, value){standardGeneric ("name<-")})

#' @rdname setters
#' @export
setGeneric("interactionCounts<-",function(GIObject, value){standardGeneric ("interactionCounts<-")})

## METHODS

#' @rdname setters
#' @export
setReplaceMethod("name", "GInteractions", function(GIObject, value){
    GIObject@metadata$experiment_name = value
    GIObject
    })

#' @rdname setters
#' @inheritParams Biobase::'description<-'
#' @importMethodsFrom Biobase 'description<-'
#' @export
setReplaceMethod("description", "GInteractions", function(object, value){
    object@metadata$description = value
    object
})

#' @rdname setters
#' @import BiocGenerics
#' @export
setReplaceMethod("interactionCounts", "GInteractions", function(GIObject, value){
    if (!all(value == floor(value)))
        stop("value must contain integer values")
    value = as.integer(value)
    if (length(value) == 1)
        value = rep(value, length(GIObject))
    GIObject@elementMetadata$counts <- value
    GIObject
})
