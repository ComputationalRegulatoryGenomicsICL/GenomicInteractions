#' Functions to access data held in a GenomicInteractions object.
#'
#' Use these functions to access data stored in each of the slots of a
#' GenomicInteractions object.
#'
#' @name getters
#' @param GIObject A GenomicInteractions object
#'
#' @return For 'anchorOne' and 'anchorTwo', a GRanges. For 'interactionCounts', 
#' a numeric vector with counts for each interaction in the object. For
#'   'description' and 'name',  a character vector with
#'   length 1. For 'annotationFeatures', a character vector of features with
#'   which the object was previously annotated, or 'NA' if the object is unannotated.
#'
#'  @examples
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5), seqlengths=seqlengths(Mmusculus))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5), seqlengths=seqlengths(Mmusculus))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, experiment_name="test", 
#'                            description="this is a test", counts=interaction_counts)
#'
#' name(test)
#' description(test)
#' anchorOne(test)
#' anchorTwo(test)
#' interactionCount(test)
#'
## GENERICS

#' @rdname getters
#' @export
setGeneric("name",function(GIObject){standardGeneric ("name")})

#' @rdname getters
#' @export
setGeneric("description",function(GIObject){standardGeneric ("description")})

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
setMethod("name", "GenomicInteractions", function(GIObject){ return(GIObject@metadata$experiment_name) } )

#' @rdname getters
#' @export
setMethod("description", "GenomicInteractions", function(GIObject){ return(GIObject@metadata$description) } )

#' @rdname getters
#' @export
setMethod("anchorOne", "GenomicInteractions", function(GIObject){ return(GIObject@anchor_one) } )

#' @rdname getters
#' @export
setMethod("anchorTwo", "GenomicInteractions", function(GIObject){ return(GIObject@anchor_two) } )

#' @rdname getters
#' @export
setMethod("interactionCounts", "GenomicInteractions", function(GIObject){ return(GIObject@counts) })

#' @rdname getters
#' @export
setMethod("annotationFeatures", "GenomicInteractions", function(GIObject){
  if( "node.class" %in% names(elementMetadata(GIObject@anchor_one))) {
    annotation = unique(c(GIObject@anchor_one$node.class, GIObject@anchor_two$node.class))
  } else { annotation = NA_character_ }
  return(annotation)
} )
