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
#' Functions to classify interactions within GInteractions objects.
#' \itemize{
#'     \item "isInteractionType" takes two character arguments which are
#'           annotated node classes and returns interactions between them.
#'     \item "is.pp", "is.pd" etc. are bindings for common annotations:
#'     \describe{ \item{p}{promoter}
#'                \item{d}{distal}
#'                \item{t}{terminator} }
#'     \item "is.trans" & "is.cis" select trans-chromosomal and
#'           intra-chromosomal interactions, respectively }
#' @param GIObject A GInteractions object
#' @return A logical vector
#' @examples 
#' data(hic_example_data)
#' table(is.cis(hic_example_data))
#' sum(interactionCounts(hic_example_data))
#' 
#' @name InteractionHelpers
#' @rdname InteractionHelpers
NULL

#' @rdname InteractionHelpers
#' @export
setMethod("is.pp", "GInteractions",
            function(GIObject){
                return(GIObject@regions$node.class[GIObject@anchor1] == "promoter" 
                      & GIObject@regions$node.class[GIObject@anchor2] == "promoter")
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.pd", "GInteractions",
            function(GIObject){
                return( (GIObject@regions$node.class[GIObject@anchor1] == "distal" & 
                           GIObject@regions$node.class[GIObject@anchor2] == "promoter" ) |
                        (GIObject@regions$node.class[GIObject@anchor1] == "promoter" & 
                           GIObject@regions$node.class[GIObject@anchor2] == "distal" ))
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.pt", "GInteractions",
            function(GIObject){
                return( (GIObject@regions$node.class[GIObject@anchor1] == "terminator" & 
                           GIObject@regions$node.class[GIObject@anchor2] == "promoter" ) |
                        (GIObject@regions$node.class[GIObject@anchor1] == "promoter" & 
                           GIObject@regions$node.class[GIObject@anchor2] == "terminator" ))
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.dd", "GInteractions",
            function(GIObject){
                return( GIObject@regions$node.class[GIObject@anchor1] == "distal" & 
                          GIObject@regions$node.class[GIObject@anchor2] == "distal")
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.dt", "GInteractions",
            function(GIObject){
                return( (GIObject@regions$node.class[GIObject@anchor1] == "distal" & 
                           GIObject@regions$node.class[GIObject@anchor2] == "terminator" ) |
                        (GIObject@regions$node.class[GIObject@anchor1] == "terminator" & 
                           GIObject@regions$node.class[GIObject@anchor2] == "distal" ))
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.tt", "GInteractions",
            function(GIObject){
                return( GIObject@regions$node.class[GIObject@anchor1] == "terminator" 
                        & GIObject@regions$node.class[GIObject@anchor2] == "terminator")
            })

#' @rdname InteractionHelpers
#' @export
setMethod("isInteractionType", "GInteractions",
            function(GIObject, x, y){
                return( (GIObject@regions$node.class[GIObject@anchor1] %in% x 
                         & GIObject@regions$node.class[GIObject@anchor2]  %in% y) |
                        (GIObject@regions$node.class[GIObject@anchor1] %in% y 
                         & GIObject@regions$node.class[GIObject@anchor2] %in% x ))
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.trans", "GInteractions",
            function(GIObject){
                return( as.character(seqnames(GIObject@regions[GIObject@anchor1])) != 
                          as.character(seqnames(GIObject@regions[GIObject@anchor2])) )
            })

#' @rdname InteractionHelpers
#' @export
setMethod("is.cis", "GInteractions",
            function(GIObject){
                return( as.character(seqnames(GIObject@regions[GIObject@anchor1])) == 
                          as.character(seqnames(GIObject@regions[GIObject@anchor2])) )
            })

#' Return the total number of interactions in a GInteractions GIObject
#'
#' @param x GInteractions GIObject
#' @return The sum of the counts in GIObject
#' @docType methods
#' @export
setMethod("sum", "GInteractions", function(x){ return( sum(interactionCounts(x))) })
