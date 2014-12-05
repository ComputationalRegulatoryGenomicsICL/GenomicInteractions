#generic definitions of functions defined in this file
#' @rdname InteractionHelpers
#' @export
setGeneric("is.pp",function(GIObject){standardGeneric ("is.pp")})
#' @rdname InteractionHelpers
#' @export
setGeneric("is.pd",function(GIObject){standardGeneric ("is.pd")})
#' @rdname InteractionHelpers
#' @export
setGeneric("is.pt",function(GIObject){standardGeneric ("is.pt")})
#' @rdname InteractionHelpers
#' @export
setGeneric("is.dd",function(GIObject){standardGeneric ("is.dd")})
#' @rdname InteractionHelpers
#' @export
setGeneric("is.dt",function(GIObject){standardGeneric ("is.dt")})
#' @rdname InteractionHelpers
#' @export
setGeneric("is.tt",function(GIObject){standardGeneric ("is.tt")})
#' @rdname InteractionHelpers
#' @export
setGeneric("isInteractionType",function(GIObject, x, y){standardGeneric ("isInteractionType")})
#' @rdname InteractionHelpers
#' @param x,y Names of annotated node classes
#' @export
setGeneric("is.trans",function(GIObject){standardGeneric ("is.trans")})
#' @rdname InteractionHelpers
#' @export
setGeneric("is.cis",function(GIObject){standardGeneric ("is.cis")})

#' Interaction Type Helpers
#' 
#' Functions to classify interactions within GenomicInteractions objects.
#' \itemize{
#'     \item "isInteractionType" takes two character arguments which are
#'           annotated node classes and returns interactions between them.
#'     \item "is.pp", "is.pd" etc. are bindings for common annotations:
#'     \describe{ \item{p}{promoter}
#'                \item{d}{distal}
#'                \item{t}{terminator} }
#'     \item "is.trans" & "is.cis" select trans-chromosomal and 
#'           intra-chromosomal interactions, respectively }
#' @param GIObject A GenomicInteractions object
#' @return A logical vector
#' @name InteractionHelpers
#' @rdname InteractionHelpers
NULL

#' @rdname InteractionHelpers
#' @export
setMethod("is.pp", "GenomicInteractions", 
            function(GIObject){ 
                return( GIObject@anchor_one$node.class == "promoter" & GIObject@anchor_two$node.class == "promoter") 
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.pd", "GenomicInteractions", 
            function(GIObject){ 
                return( (GIObject@anchor_one$node.class == "distal" & GIObject@anchor_two$node.class == "promoter" ) | 
                        (GIObject@anchor_one$node.class == "promoter" & GIObject@anchor_two$node.class == "distal" )) 
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.pt", "GenomicInteractions", 
            function(GIObject){ 
                return( (GIObject@anchor_one$node.class == "distal" & GIObject@anchor_two$node.class == "terminator" ) | 
                        (GIObject@anchor_one$node.class == "terminator" & GIObject@anchor_two$node.class == "distal" )) 
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.dd", "GenomicInteractions", 
            function(GIObject){ 
                return( GIObject@anchor_one$node.class == "distal" & GIObject@anchor_two$node.class == "distal") 
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.dt", "GenomicInteractions", 
            function(GIObject){ 
                return( (GIObject@anchor_one$node.class == "distal" & GIObject@anchor_two$node.class == "terminator" ) | 
                        (GIObject@anchor_one$node.class == "distal" & GIObject@anchor_two$node.class == "terminator" )) 
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.tt", "GenomicInteractions", 
            function(GIObject){ 
                return( GIObject@anchor_one$node.class == "terminator" & GIObject@anchor_two$node.class == "terminator") 
            })

#' @rdname InteractionHelpers
#' @export
setMethod("isInteractionType", "GenomicInteractions", 
            function(GIObject, x, y){ 
                return( (GIObject@anchor_one$node.class %in% x & GIObject@anchor_two$node.class  %in% y) | 
                        (GIObject@anchor_one$node.class %in% y & GIObject@anchor_two$node.class %in% x ))  
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.trans", "GenomicInteractions", 
            function(GIObject){ 
                return( as.character(seqnames(GIObject@anchor_one)) != as.character(seqnames(GIObject@anchor_two)) ) 
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.cis", "GenomicInteractions", 
            function(GIObject){ 
                return( as.character(seqnames(GIObject@anchor_one)) == as.character(seqnames(GIObject@anchor_two)) )
            })

#' Find overlaps between a GRanges and a GenomicInteractions object
#'
#' This function calls findOverlaps separately on each anchor and
#' returns a list. See 'findOverlaps' in the GenomicRanges package for detailed
#' documentation for this function.
#'
#' @param query GenomicInteractions or GRanges
#' @param subject GRanges or GenomicInteractions
#' @param maxgap,minoverlap,type,select See 'findOverlaps' in the IRanges package.
#' 
#' @return A list containing Hits objects for anchors one and two.
#' 
#' @rdname GenomicInteractions-overlaps-methods
#' @docType methods
#' @name findOverlaps
#' @import GenomicRanges 
#' @export
NULL

#' @rdname GenomicInteractions-overlaps-methods
#' @export
setMethod("findOverlaps", c("GenomicInteractions", "GRanges"), function(query, subject,  maxgap = 0L, minoverlap = 1L,
                                                                        type = c("any", "start", "end", "within", "equal"),
                                                                        select = c("all", "first", "last", "arbitrary"),
                                                                        algorithm = c("nclist", "intervaltree")){
    algorithm <- match.arg(algorithm)
    return(list(one=findOverlaps(anchorOne(query), subject, maxgap=maxgap, minoverlap=minoverlap, type=type, select=select, algorithm=algorithm), 
                two=findOverlaps(anchorTwo(query), subject, maxgap=maxgap, minoverlap=minoverlap, type=type, select=select, algorithm=algorithm)))  
})

#' @rdname GenomicInteractions-overlaps-methods
#' @export
setMethod("findOverlaps", c("GRanges", "GenomicInteractions"), function(query, subject, maxgap = 0L, minoverlap = 1L,
                                                                        type = c("any", "start", "end", "within", "equal"),
                                                                        select = c("all", "first", "last", "arbitrary"),
                                                                        algorithm = c("nclist", "intervaltree")){
  algorithm <- match.arg(algorithm)
  return(list(one=findOverlaps(query, anchorOne(subject), maxgap=maxgap, minoverlap=minoverlap, type=type, select=select, algorithm=algorithm), 
              two=findOverlaps(query, anchorTwo(subject), maxgap=maxgap, minoverlap=minoverlap, type=type, select=select, algorithm=algorithm)
              ))  
})

#' Print function for GenomicInteractions
#'
#' @param x GenomicInteractionsObject
#' @return invisible(1)
#' @docType methods
#' @export
setMethod("print", "GenomicInteractions", function(x){
    cat("GenomicInteractions\n")
    cat("\tName: ", x@experiment_name, "\n", sep="")
    cat("\tDescription: ", x@description, "\n", sep="")
    cat("\tGenome: ", x@genome_name,"\n", sep="")
    cat("\tNumber of individual interactions: ",  length(x), "\n", sep="")
    cat("\tNumber of interactions: ",  sum(x), "\n", sep="")
    cat("\t", "Annotated: ", ifelse( "node.class" %in% names(elementMetadata(x@anchor_one)), "yes", "no"), "\n")
        cat("\t\t", "Annotated with: ",  ifelse( "node.class" %in% names(elementMetadata(x@anchor_one)), 
                                                paste(unique(c(x@anchor_one$node.class, x@anchor_two$node.class)), collapse=", "),
                                                "N/A"), 
            "\n", sep="")
    cat("\tInteractions:\n")
    cat("\t\t", paste(.pasteAnchor(as.character(seqnames(x@anchor_one))[1:min(10, length(x))], 
                                   start(x@anchor_one)[1:min(10, length(x))], 
                                   end(x@anchor_one)[1:min(10, length(x))]),
                        .pasteAnchor(as.character(seqnames(x@anchor_two))[1:min(10, length(x))], 
                                     start(x@anchor_two)[1:min(10, length(x))], 
                                     end(x@anchor_two)[1:min(10, length(x))]),
                        sep="\t-----\t", collapse="\n\t\t"), sep="")
    cat(ifelse(length(x)>10, "\n\t\t....\n", ""))
    cat("\n")
    return(invisible(1))
})

#' Representation function for GenomicInteractions
#'
#' @param object A GenomicInteractionsObject
#' @return invisible(1)
#' @docType methods
#' @export
setMethod("show", "GenomicInteractions", function(object){ 
    cat("GenomicInteractions\n")
    cat("\tName: ", object@experiment_name, "\n", sep="")
    cat("\tDescription: ", object@description, "\n", sep="")
    cat("\tGenome: ", object@genome_name,"\n", sep="")
    cat("\tNumber of individual interactions: ",  length(object), "\n", sep="")
    cat("\tNumber of interactions: ",  sum(object), "\n", sep="")
    cat("\t", "Annotated: ", ifelse( "node.class" %in% names(elementMetadata(object@anchor_one)), "yes", "no"), "\n", sep="")
    cat("\t\t", "Annotated with: ",  ifelse( "node.class" %in% names(elementMetadata(object@anchor_one)), 
                                                paste(unique(c(object@anchor_one$node.class, object@anchor_two$node.class)), collapse=", "),
                                                "N/A"), 
               "\n", sep="")
    cat("\tInteractions:\n")
    cat("\t\t", paste(.pasteAnchor(as.character(seqnames(object@anchor_one))[1:min(10, length(object))], 
                                   start(object@anchor_one)[1:min(10, length(object))], 
                                   end(object@anchor_one)[1:min(10, length(object))]),
                             .pasteAnchor(as.character(seqnames(object@anchor_two))[1:min(10, length(object))], 
                                          start(object@anchor_two)[1:min(10, length(object))], 
                                          end(object@anchor_two)[1:min(10, length(object))]),
                             sep="\t-----\t", collapse="\n\t\t"), sep="")
    cat(ifelse(length(object)>10, "\n\t\t....\n", ""))
    cat("\n")
    return(invisible(1))
})
       
