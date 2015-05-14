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
                return( (GIObject@anchor_one$node.class == "terminator" & GIObject@anchor_two$node.class == "promoter" ) |
                        (GIObject@anchor_one$node.class == "promoter" & GIObject@anchor_two$node.class == "terminator" ))
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
                        (GIObject@anchor_one$node.class == "terminator" & GIObject@anchor_two$node.class == "distal" ))
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

#' Return the total number of interactions in a GenomicInteractions GIObject
#'
#' @param x GenomicInteractions GIObject
#' @return The sum of the counts in GIObject
#' @docType methods
#' @export
setMethod("sum", "GenomicInteractions", function(x){ return( sum(interactionCounts(x))) })

#' Find overlaps between GRanges and GenomicInteractions objects
#'
#' When called with a GRanges and a GenomicInteractions object, this function
#' calls findOverlaps separately on each anchor and returns a list. countOverlaps
#' and overlapsAny return a list of integer vectors and logical vectors respectively.
#'
#' When
#'
#' See
#' 'findOverlaps' in the GenomicRanges package for detailed
#' documentation for this function.
#'
#' @param query GenomicInteractions or GRanges
#' @param subject GRanges or GenomicInteractions
#' @param maxgap,minoverlap,type,select,ignore.strand See 'findOverlaps' in the IRanges package.
#'
#' @return A Hits object or a list containing Hits objects for both anchors.
#'
#' @rdname GenomicInteractions-overlaps-methods
#' @docType methods
#' @name findOverlaps
#' @import GenomicRanges
#' @export
NULL

#' @rdname GenomicInteractions-overlaps-methods
#' @export
setMethod("findOverlaps", c("GenomicInteractions", "GRanges"), 
          function(query, subject,  maxgap = 0L, minoverlap = 1L,
                  type = c("any", "start", "end", "within", "equal"),
                  select = c("all", "first", "last", "arbitrary"),
                  ignore.strand = FALSE){
    return(list(one=findOverlaps(anchorOne(query), subject, 
                                 maxgap=maxgap, minoverlap=minoverlap, type=type, 
                                 select=select, ignore.strand = ignore.strand),
                two=findOverlaps(anchorTwo(query), subject, 
                                 maxgap=maxgap, minoverlap=minoverlap, type=type, 
                                 select=select, ignore.strand = ignore.strand)
                ))
})

#' @rdname GenomicInteractions-overlaps-methods
#' @export
setMethod("findOverlaps", c("GRanges", "GenomicInteractions"), 
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                  type = c("any", "start", "end", "within", "equal"),
                  select = c("all", "first", "last", "arbitrary"),
                  ignore.strand=FALSE){
  return(list(one=findOverlaps(query, anchorOne(subject), 
                               maxgap=maxgap, minoverlap=minoverlap, type=type, 
                               select=select, ignore.strand = ignore.strand),
              two=findOverlaps(query, anchorTwo(subject), 
                               maxgap=maxgap, minoverlap=minoverlap, type=type, 
                               select=select, ignore.strand = ignore.strand)
              ))
})

#' @rdname GenomicInteractions-overlaps-methods
#' @export
setMethod("countOverlaps", c("GenomicInteractions", "GRanges"), 
          function(query, subject,  maxgap = 0L, minoverlap = 1L,
                   select = c("all", "first", "last", "arbitrary"),
                   type = c("any", "start", "end", "within", "equal"),
                   ignore.strand=FALSE){
    return(list(one=countOverlaps(anchorOne(query), subject, 
                                  maxgap=maxgap, minoverlap=minoverlap, type=type, 
                                  select=select, ignore.strand = ignore.strand),
                two=countOverlaps(anchorTwo(query), subject, 
                                  maxgap=maxgap, minoverlap=minoverlap, type=type, 
                                  select=select, ignore.strand = ignore.strand)
                ))
})

#' @rdname GenomicInteractions-overlaps-methods
#' @export
setMethod("countOverlaps", c("GRanges", "GenomicInteractions"), 
            function(query, subject, maxgap = 0L, minoverlap = 1L,
                    type = c("any", "start", "end", "within", "equal"),
                    select = c("all", "first", "last", "arbitrary"),
                    ignore.strand=FALSE){
  return(list(one=countOverlaps(query, anchorOne(subject), 
                                maxgap=maxgap, minoverlap=minoverlap, type=type, 
                                select=select, ignore.strand = ignore.strand),
              two=countOverlaps(query, anchorTwo(subject), 
                                maxgap=maxgap, minoverlap=minoverlap, type=type, 
                                select=select, ignore.strand = ignore.strand)
              ))
})

#' @rdname GenomicInteractions-overlaps-methods
#' @export
setMethod("overlapsAny", c("GenomicInteractions", "GRanges"), 
          function(query, subject,  maxgap = 0L, minoverlap = 1L,
                  type = c("any", "start", "end", "within", "equal"),
                  select = c("all", "first", "last", "arbitrary"),
                  ignore.strand=FALSE){
    return(list(one=overlapsAny(anchorOne(query), subject,
                                maxgap=maxgap, minoverlap=minoverlap, type=type, 
                                select=select, ignore.strand = ignore.strand),
                two=overlapsAny(anchorTwo(query), subject, 
                                maxgap=maxgap, minoverlap=minoverlap, type=type, 
                                select=select, ignore.strand = ignore.strand)
                ))
})

#' @rdname GenomicInteractions-overlaps-methods
#' @export
setMethod("overlapsAny", c("GRanges", "GenomicInteractions"), 
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                  type = c("any", "start", "end", "within", "equal"),
                  select = c("all", "first", "last", "arbitrary"),
                  ignore.strand=FALSE){
  return(list(one=overlapsAny(query, anchorOne(subject), 
                              maxgap=maxgap, minoverlap=minoverlap, type=type, 
                              select=select, ignore.strand = ignore.strand),
              two=overlapsAny(query, anchorTwo(subject), 
                              maxgap=maxgap, minoverlap=minoverlap, type=type, 
                              select=select, ignore.strand = ignore.strand)
              ))
})

#' @rdname GenomicInteractions-overlaps-methods
#' @export
setMethod("findOverlaps", c("GenomicInteractions", "GenomicInteractions"), function(query, subject) {
    if (!all(.isSorted(query)) & !all(.isSorted(subject)))
        stop("GenomicInteractions object must be sorted")
    overlaps.one = as.matrix(findOverlaps(anchorOne(query), anchorOne(subject)))
    overlaps.two = as.matrix(findOverlaps(anchorTwo(query), anchorTwo(subject)))
    # this is afaik the best way to find duplicates in R
    overlaps.all = rbind(overlaps.one, overlaps.two)
    overlaps.common = overlaps.all[duplicated(overlaps.all), ]
    hits_object = new("Hits", queryHits = overlaps.common[,"queryHits"], subjectHits=overlaps.common[,"subjectHits"], queryLength=length(query), subjectLength=length(subject))
    return(hits_object)
})

#' Acessing/modifying sequence information for a GenomicInteractions object
#'
#' Allows access/modification of seqinfo for GenomicInteractions objcets. When
#' used with "force=True", interactions with either (or both) anchors on invalid
#' chromosomes will be removed.
#'
#' For more information see ?seqinfo in the GenomeInfoDb
#' package.
#'
#' @param x A GenomicInteractions object
#' @return A seqinfo object,
#' @docType methods
#' @rdname seqinfo-GenomicInteractions-method
#' @import GenomicRanges
#' @export
setMethod("seqinfo", "GenomicInteractions", function(x) {
    if (!.isEqualSeqInfo(anchorOne(x), anchorTwo(x))) {
        objName = deparse(substitute(x))
        stop(paste("Seqinfo differs between anchors in", objName))
    }
    return(seqinfo(anchorOne(x)))
})

#' @param new2old Mapping between new and old seqnames. See ?seqinfo in GenomeInfoDb for details.
#' @param force A logical indicating whether or not to drop invalid levels.
#' @param value A replacement seqinfo object
#' @rdname seqinfo-GenomicInteractions-method
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import BiocGenerics
#' @export
setReplaceMethod("seqinfo", "GenomicInteractions", function(x, new2old=NULL, force=FALSE, value) {
    if (!is(value, "Seqinfo"))
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    if (!.isEqualSeqInfo(anchorOne(x), anchorTwo(x))) {
        objName = deparse(substitute(x))
        stop(paste("Seqinfo differs between anchors in", objName))
    }
    dangling_seqlevels_one <- GenomeInfoDb:::getDanglingSeqlevels(anchorOne(x), new2old=new2old, force=force, seqlevels(value))
    dangling_seqlevels_two <- GenomeInfoDb:::getDanglingSeqlevels(anchorTwo(x), new2old=new2old, force=force, seqlevels(value))
    dangling_seqlevels = unique(c(dangling_seqlevels_one, dangling_seqlevels_two))
    if (length(dangling_seqlevels) != 0L) {
        x.one = !(seqnames(anchorOne(x)) %in% dangling_seqlevels)
        x.two = !(seqnames(anchorTwo(x)) %in% dangling_seqlevels)
        x = x[x.one & x.two]
    }
    old_seqinfo = seqinfo(x)
    new_seqnames_one <- GenomeInfoDb:::makeNewSeqnames(anchorOne(x), new2old=new2old, seqlevels(value))
    new_seqnames_two <- GenomeInfoDb:::makeNewSeqnames(anchorTwo(x), new2old=new2old, seqlevels(value))
    x@anchor_one = BiocGenerics:::updateS4(anchorOne(x), seqnames=new_seqnames_one, seqinfo=value, check=FALSE)
    x@anchor_two = BiocGenerics:::updateS4(anchorTwo(x), seqnames=new_seqnames_two, seqinfo=value, check=FALSE)
    # check for object validity? need valid.GenomicInteractions.seqinfo method ?
    return(x)
})


#' Sort GenomicInteractions Object
#'
#' This method will sort a GenomicInteractions object by first
#' arranging all interactions start on the lower-ordered anchor;
#' for trans-chromosomal interactions this is the anchor on the
#' lower ordered chromomsome (defined by the seqlevels factor);
#' and then by the position of the first anchor. See 
#' ?"GenomicRanges-comparison" for ordering rules.
#'
#' @param x GenomicInteractions Object
#' @param decreasing A logical indicating sort order
#' @param order.interactions A logical indicating if interactions should be reordered,
#'                           or only rearranged so that start(anchorOne) < start(anchorTwo)
#' @import GenomicRanges
#' @return A sorted GenomicInteractions object
#' @docType methods
#' @export
setMethod("sort", "GenomicInteractions", function(x, decreasing=FALSE, order.interactions=TRUE) {
          if (any(names(mcols(anchorOne(x))) != names(mcols(anchorTwo(x)))))
              stop("Metadata differs between anchors and will be lost in sort.")
          anchor.one = anchorOne(x)
          anchor.two = anchorTwo(x)
          reversed = !.isSorted(x)
          if (decreasing==TRUE) {reversed = !reversed}
          one.rev = anchor.one[reversed]
          anchor.one[reversed] = anchor.two[reversed]
          anchor.two[reversed] = one.rev
          x <- BiocGenerics:::updateS4(x, anchor_one=anchor.one, anchor_two=anchor.two)
          if (order.interactions==TRUE) {
            i = order(anchor.one, decreasing=decreasing)
            x = x[i]
          }
          return(x)
})

#' Trim a GenomicInteractions object
#'
#' This will remove any interactions with an anchor falling outside
#' of the seqlengths in a GenomicInteractions object, and trim ranges
#' which cross the ends of chromosomes.
#'
#' @param x A GenomicInteractions object
#' @param minAnchorSize The minimum size anchor to allow when trimming ranges.
#' @param ... any additional arguments to trim
#' @return A trimmed GenomicInteractions object
#' @importFrom IRanges trim
#' @docType methods
#' @export
setMethod("trim", "GenomicInteractions", function(x, minAnchorSize=1, ...) {
          one.valid = start(x@anchor_one) < seqlengths(x@anchor_one)[as.character(seqnames(x@anchor_one))] - minAnchorSize
          two.valid = start(x@anchor_two) < seqlengths(x@anchor_two)[as.character(seqnames(x@anchor_two))] - minAnchorSize
          suppressWarnings({ x = x[one.valid & two.valid] }) # otherwise warns if x needs trimming
          x@anchor_one = trim(anchorOne(x), ...)
          x@anchor_two = trim(anchorTwo(x), ...)
          return(x)
})

.makeNakedMatFromGenomicInteractions = function(x) {
    lx <- length(x)
    nc <- ncol(mcols(x))
    ans <- cbind("Anchor One"=.pasteAnchor(x@anchor_one),
                 "   "=rep.int("---", lx),
                 "Anchor Two"=.pasteAnchor(x@anchor_two),
                 "Counts"=x@counts)
    if (nc > 0L) {
        tmp <- do.call(data.frame, c(lapply(mcols(x), showAsCell), list(check.names=FALSE)))
        ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
    }
    ans
}


showGenomicInteractions = function(x, margin="", print.seqinfo=FALSE) {
    lx <- length(x)
    nc <- ncol(mcols(x))
    cat(class(x), " object with ",
        lx, " ", ifelse(lx == 1L, "interaction", "interactions"),
        " and ",
        nc, " metadata ", ifelse(nc == 1L, "column", "columns"),
        ":\n", sep="")
    if (!is.null(name(x))) cat("\tName: ", name(x), "\n", sep="")
    if (!is.null(description(x))) cat("\tDescription: ", description(x), "\n", sep="")
    cat("\tSum of interactions: ",  sum(x), "\n", sep="")
    is_atd = ifelse( "node.class" %in% names(elementMetadata(x@anchor_one)), "yes", "no")
    cat("\tAnnotated: ", is_atd, "\n")
    if (is_atd == "yes") {
        annotations = paste(unique(c(x@anchor_one$node.class, x@anchor_two$node.class)), collapse=", ")
        cat("\t\t", "Annotated with: ",  annotations, "\n")
    }
    cat("\tInteractions:\n")
    out = S4Vectors:::makePrettyMatrixForCompactPrinting(x, .makeNakedMatFromGenomicInteractions)
    if (nrow(out) != 0L)
        rownames(out) = paste0(margin, rownames(out))
    print(out, quote=FALSE, right=TRUE, max=length(out))
    #cat(ifelse(lx>10, "\n\t\t....\n", ""))
    cat("\n")
    if (print.seqinfo) {
        cat(margin, "-------\n")
        cat(margin, "seqinfo: ", summary(seqinfo(x)), "\n", sep="")
    }
}

#' Print function for GenomicInteractions
#'
#' @param x GenomicInteractionsObject
#' @return invisible(1)
#' @docType methods
#' @export
setMethod("print", "GenomicInteractions", function(x){
    showGenomicInteractions(x, margin="  ", print.seqinfo=TRUE)
})

#' Representation function for GenomicInteractions
#'
#' @param object A GenomicInteractionsObject
#' @return invisible(1)
#' @docType methods
#' @export
setMethod("show", "GenomicInteractions", function(object){
    showGenomicInteractions(object, margin="  ", print.seqinfo=TRUE)
})


.duplicated.GenomicInteractions <- function(x, fromLast=FALSE, dropMetadata = FALSE)
{
  if (dropMetadata == TRUE){
    dat <- cbind(as.data.frame(unname(anchorOne(x))), #duplicated names ok for GRanges but not for df 
                as.data.frame(unname(anchorTwo(x))), 
                interactionCounts(x))
  } else if (dropMetadata == FALSE){
    dat <- cbind(as.data.frame(unname(anchorOne(x))), #duplicated names ok for GRanges but not for df 
                 as.data.frame(unname(anchorTwo(x))), 
                 interactionCounts(x),
                 mcols(x))
  } else {
    stop("dropMetadata must be TRUE or FALSE")
  }
  return(duplicated(dat, incomparables = FALSE, 
                    fromLast = fromLast))
}

#' duplicated, GenomicInteractions-method
#' 
#' Finds duplicated interactions in a GenomicInteractions object. 
#' 
#' Uniqueness is based on anchor positions and metadata, interaction counts, and
#' interaction metadata.
#'
#' @param x A GenomicInteractions object
#' @param fromLast Whether to identify duplicates starting from last item in the 
#'        Genomicinteractions object or not. Default: FALSE.
#' @param dropMetadata Logical, default FALSE. Whether to drop interaction mcols 
#'        when considering unique interactions.
#' @return A vector containing indices of duplicated interactions
#' @docType methods
#' @export
setMethod("duplicated", "GenomicInteractions", .duplicated.GenomicInteractions)

.unique.GenomicInteractions <- function(x, dropMetadata = FALSE)
{
  idx <- !duplicated(x, dropMetadata = dropMetadata)
  return(x[idx])
}

#' unique, GenomicInteractions-method
#' 
#' Finds unique interactions in a GenomicInteractions object. 
#' 
#' Uniqueness is based on anchor positions and metadata, interaction counts, and
#' interaction metadata (unless dropMetadata is TRUE)
#'
#' @param x GenomicInteractionsObject
#' @param dropMetadata Logical, default FALSE. Whether to drop interaction mcols 
#'        when considering unique interactions.
#' @return A GenomicInteractions object
#' @docType methods
#' @export
#' @examples
#' 
#' library(GenomicInteractions)
#' 
#' data(hic_example_data)
#' unique(hic_example_data[c(1:4, 1:5)])

setMethod("unique", "GenomicInteractions", .unique.GenomicInteractions)

