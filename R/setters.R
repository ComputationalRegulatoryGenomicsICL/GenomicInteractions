#' Functions to set data held in a GenomicInteractions object.
#' 
#' Use these functions to set data stored in each of the slots of a 
#' GenomicInteractions object.
#' 
#' @name setters
#' @param GIObject A GenomicInteractions object
#' @param value A vector to replace a slot in the object
#' @return GenomicInteractions object
#' 
#'  @examples
#' library(BSgenome.Mmusculus.UCSC.mm9)
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5), seqlengths=seqlengths(Mmusculus))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5), seqlengths=seqlengths(Mmusculus))
#' test <- new("GenomicInteractions", experiment_name="test", description="this is a test", 
#'                  genome_name="BSgenome.Mmusculus.UCSC.mm9", anchor_one = anchor.one, 
#'                  anchor_two = anchor.two, counts=as.integer(c(2,1,2,3)) )
#' 
#' name(test) <- "Mouse test"
#' name(test)
#' 
#' description(test) <- "This is a test using the mouse genome"
#' description(test)
#' 
#' pValue(test) = c(0.1, 0.3, 0.1, 0.08)
#' pValue(test)
#' 
#' FDR(test) = p.adjust(pValue(test), "bonferroni")
#' FDR(test)
#' 
#' normalisedCount(test) = count(test) / sum(test)
#' normalisedCount(test)
#' 
## GENERICS
#' @rdname setters
#' @aliases name<-
#' @export
setGeneric("name<-",function(GIObject, value){standardGeneric ("name<-")})

#' @rdname setters
#' @export
setGeneric("description<-",function(GIObject, value){standardGeneric ("description<-")})

#' @rdname setters
#' @export
setGeneric("pValue<-",function(GIObject, value){standardGeneric ("pValue<-")})

#' @rdname setters
#' @export
setGeneric("FDR<-",function(GIObject, value){standardGeneric ("FDR<-")})

#' @rdname setters
#' @export
setGeneric("normalisedCount<-",function(GIObject, value){standardGeneric ("normalisedCount<-")})

## METHODS

#' @rdname setters
#' @export
setReplaceMethod("name", "GenomicInteractions", function(GIObject, value){ 
    GIObject@experiment_name = value
    validObject(GIObject)
    GIObject
    })

#' @rdname setters
#' @export
setReplaceMethod("normalisedCount", "GenomicInteractions", function(GIObject, value){ 
    GIObject@normalised_counts = value
    validObject(GIObject)
    GIObject
    })

#' @rdname setters 
#' @export
setReplaceMethod("pValue", "GenomicInteractions", function(GIObject, value){ 
    GIObject@pvalue = value
    validObject(GIObject)
    GIObject
    })

#' @rdname setters 
#' @export
setReplaceMethod("FDR", "GenomicInteractions", function(GIObject, value){ 
    GIObject@fdr = value
    validObject(GIObject)
    GIObject
    })

#' @rdname setters 
#' @export
setReplaceMethod("description", "GenomicInteractions", function(GIObject, value){ 
    GIObject@description = value
    validObject(GIObject)
    GIObject
})
