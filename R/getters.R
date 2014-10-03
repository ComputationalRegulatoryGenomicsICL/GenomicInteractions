#' Functions to access data held in a GenomicInteractions object.
#' 
#' Use these functions to access data stored in each of the slots of a 
#' GenomicInteractions object.
#' 
#' @name getters
#' @param GIObject A GenomicInteractions object
#' 
#' @return For 'anchorOne' and 'anchorTwo', a GRanges. For 'counts', 
#'  'normalisedCount', pValue', 'FDR', a numeric vector with counts, 
#'  normalised counts, p-values or FDRs for each interaction in the object. For
#'   'description','name', and 'genomeName', a character vector with 
#'   length 1. For 'annotationFeatures', a character vector of features with
#'   which the object was previously annotated, or 'NA' if the object is unannotated.
#'  
#'  @examples
#' library(BSgenome.Mmusculus.UCSC.mm9)
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5), seqlengths=seqlengths(Mmusculus))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5), seqlengths=seqlengths(Mmusculus))
#' test <- new("GenomicInteractions", experiment_name="test", description="this is a test", 
#'                  genome_name="BSgenome.Mmusculus.UCSC.mm9", anchor_one = anchor.one, 
#'                  anchor_two = anchor.two, counts=as.integer(c(2,1,2,3)), pvalue=c(0.1, 0.3, 0.1, 0.08))
#'  
#' name(test)
#' description(test)
#' anchorOne(test)
#' anchorTwo(test)
#' count(test)
#' pValue(test)
#' genomeName(test)
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
setGeneric("count",function(GIObject){standardGeneric ("count")})

#' @rdname getters
#' @export
setGeneric("pValue",function(GIObject){standardGeneric ("pValue")})

#' @rdname getters
#' @export
setGeneric("FDR",function(GIObject){standardGeneric ("FDR")})

#' @rdname getters
#' @export
setGeneric("normalisedCount",function(GIObject){standardGeneric ("normalisedCount")})

#' @rdname getters
#' @export
setGeneric("description",function(GIObject){standardGeneric ("description")})

#' @rdname getters
#' @export
setGeneric("genomeName",function(GIObject){standardGeneric ("genomeName")})

#' @rdname getters
#' @export
setGeneric("annotationFeatures",function(GIObject){standardGeneric ("annotationFeatures")})

#' Get the length of a GenomicInteractions GIObject
#'
#' @param x GenomicInteractions GIObject
#' @return A numeric vector containing the length of the GIObject
#' @docType methods
#' @export
setMethod("length", "GenomicInteractions", function(x){ return(length(x@anchor_one)) } )

#' Return the total number of interactions in a GenomicInteractions GIObject
#'
#' @param x GenomicInteractions GIObject
#' @return The sum of the counts in GIObject 
#' @docType methods
#' @export
setMethod("sum", "GenomicInteractions", function(x){ return( sum(x@counts)) })

#' @rdname getters
#' @export
#' @aliases name
setMethod("name", "GenomicInteractions", function(GIObject){ return(GIObject@experiment_name) } )

#' @rdname getters
#' @export
setMethod("anchorOne", "GenomicInteractions", function(GIObject){ return(GIObject@anchor_one) } )

#' @rdname getters 
#' @export
setMethod("anchorTwo", "GenomicInteractions", function(GIObject){ return(GIObject@anchor_two) } )

#' @rdname getters 
#' @export
setMethod("count", "GenomicInteractions", function(GIObject){ return(GIObject@counts) } )

#' @rdname getters
#' @export
setMethod("normalisedCount", "GenomicInteractions", function(GIObject){ return(GIObject@normalised_counts) } )

#' @rdname getters 
#' @export
setMethod("pValue", "GenomicInteractions", function(GIObject){ return(GIObject@pvalue) } )

#' @rdname getters 
#' @export
setMethod("FDR", "GenomicInteractions", function(GIObject){ return(GIObject@fdr) } )

#' @rdname getters 
#' @export
setMethod("description", "GenomicInteractions", function(GIObject){ return(GIObject@description) } )

#' @rdname getters 
#' @export
setMethod("genomeName", "GenomicInteractions", function(GIObject){ return(GIObject@genome_name) } )


#' @rdname getters 
#' @export
setMethod("annotationFeatures", "GenomicInteractions", function(GIObject){
  if( "node.class" %in% names(elementMetadata(GIObject@anchor_one))) {
    annotation = unique(c(GIObject@anchor_one$node.class, GIObject@anchor_two$node.class))
  } else { annotation = "NA" }
  return(annotation)
} )
