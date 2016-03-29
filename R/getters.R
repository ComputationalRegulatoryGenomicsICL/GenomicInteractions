#' Functions to access data held in a GenomicInteractions object.
#'
#' Use these functions to access data stored in each of the slots of a
#' GenomicInteractions object.
#'
#' @name getters
#' @param GIObject A Gnteractions object
#' @rdname getters
#' 
#' @return For 'anchorOne' and 'anchorTwo', a GRanges. For 'interactionCounts', 
#' a numeric vector with counts for each interaction in the object. For
#'   'description' and 'name',  a character vector with
#'   length 1. For 'annotationFeatures', a character vector of features with
#'   which the object was previously annotated, or 'NA' if the object is unannotated.
#'
#' @examples
#'
#' library(GenomicRanges)
#'
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, counts=interaction_counts)
#'
#' name(test)
#' description(test)
#' anchorOne(test)
#' anchorTwo(test)
#' interactionCounts(test)
#'
## GENERICS

#' @rdname getters
#' @export
setGeneric("name",function(GIObject){standardGeneric ("name")})

#' @rdname getters
#' @export
setGeneric("anchorOne",function(GIObject){standardGeneric ("anchorOne")})

#' @rdname getters
#' @export
setGeneric("anchorTwo",function(GIObject){standardGeneric ("anchorTwo")})

#' @rdname getters
#' @export
setGeneric("interactionCounts",function(GIObject){standardGeneric ("interactionCounts")})

#' @rdname getters
#' @export
setGeneric("annotationFeatures",function(GIObject){standardGeneric ("annotationFeatures")})

## METHODS

#' @rdname getters
#' @export
#' @aliases name
setMethod("name", "GInteractions", function(GIObject){ 
  return(GIObject@metadata$experiment_name) } )

#' @rdname getters
#' @inheritParams Biobase::description
#' @importMethodsFrom Biobase description
#' @export
setMethod("description", c("GInteractions"), function(object){ 
  return(object@metadata$description) } )

#' @rdname getters
#' @export
setMethod("anchorOne", "GInteractions", function(GIObject){ 
  return(anchors(GIObject, type = "first")) } )

#' @rdname getters
#' @export
setMethod("anchorTwo", "GInteractions", function(GIObject){ 
  return(anchors(GIObject, type = "second")) } )

## N.B. this may not apply to all GInteractions objects
#' @rdname getters
#' @export
setMethod("interactionCounts", "GInteractions", function(GIObject){
  if (!"counts" %in% names(elementMetadata(GIObject))){
    warning("'counts' not in mcols of object; will return NULL")
  }
  return(GIObject@elementMetadata$counts) })

#' @rdname getters
#' @export
setMethod("annotationFeatures", "GInteractions", function(GIObject){
  if( "node.class" %in% names(elementMetadata(GIObject@regions))) {
    annotation = unique(GIObject@regions$node.class)
  } else { annotation = NA_character_ }
  return(annotation)
} )
