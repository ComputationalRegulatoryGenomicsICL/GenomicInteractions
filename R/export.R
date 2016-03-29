#' Export interactions in BED12 format.
#'
#' @param GIObject  A GInteractions object.
#' @param fn        A filename to write the object to
#' @param score     Which metadata column to export as score
#' @param drop.trans Logical indicating whether to drop trans interactions.
#'
#' Exports a GInteractions object to BED12 format, and writes to a specified file. If filename is not specified,
#' then a data.frame containing the information is returned.
#'
#' Bed12 files provide a method for visualising interactions, it is not a good format for storing all of the data associated
#' with an interaction dataset, particularly for trans-chromosomal interactions, which can only be stored in the bed12 names
#' field.
#'
#' @return invisible(1) if outputting to file or a data.frame containing all of the corresponding information
#' @export
#' @examples 
#' data(hic_example_data)
#' export.bed12(hic_example_data, fn = tempfile(), score = "counts", drop.trans=TRUE)
#' @docType methods
#' @rdname export.bed12

setGeneric("export.bed12",function(GIObject, fn=NULL, score="counts", drop.trans=c(FALSE, TRUE)){standardGeneric ("export.bed12")})
#' @rdname export.bed12
#' @export
#' @importFrom utils write.table
setMethod("export.bed12", c("GInteractions"),
        function(GIObject, fn=NULL, score="counts", drop.trans=c(FALSE, TRUE)){
           
          
           # drop.trans not used???
            GIObject = sort(GIObject) # to do: check this sorting is okay

            is_trans = is.trans(GIObject)
            cis = GIObject[!is_trans]
            trans = GIObject[is_trans]

            len = length(cis)
            s1 = as.vector(strand(anchorOne(cis)))
            s2 = as.vector(strand(anchorTwo(cis)))

            score_vector = .getScore(cis, score)
            if (is.null(score_vector)) stop("Supplied score field not in element metadata.")
            
            output_cis = data.frame(
                chr=as.character(seqnames(anchorOne(cis))),
                start=start(anchorOne(cis))-1,
                end=end(anchorTwo(cis)),
                name=.exportName(cis),
                score=score_vector,
                strand=ifelse(s1 == s2 & s1 %in% c("+", "-"), s1, "."), # avoid case where strand == "*"
                thickStart=start(anchorOne(cis))-1,
                thickEnd=end(anchorTwo(cis)),
                itemRgb=rep("255,0,0", len),
                blockCount=2,
                blockSizes=paste(as.character(width(anchorOne(cis))),
                                 as.character(width(anchorTwo(cis))), sep=","),
                                 blockStarts=paste(0, start(anchorTwo(cis)) - start(anchorOne(cis)), sep=","),
                stringsAsFactors = FALSE)

            score_vector = .getScore(trans, score)
            if (is.null(score_vector)) stop("Supplied score field not in element metadata.")
            output_trans = data.frame(
                chr=c(as.character(seqnames(anchorOne(trans))),
                      as.character(seqnames(anchorTwo(trans)))),
                start=c(as.character(start(anchorOne(trans))),
                        as.character(start(anchorTwo(trans)))),
                end=c(as.character(end(anchorOne(trans))),
                      as.character(end(anchorTwo(trans)))),
                name=rep(.exportName(trans), 2),
                score=rep(score_vector, 2),
                strand=c(as.character(strand(anchorOne(trans))),
                         as.character(strand(anchorTwo(trans)))),
                thickStart=c(as.character(start(anchorOne(trans))),
                             as.character(start(anchorTwo(trans)))),
                thickEnd=c(as.character(end(anchorOne(trans))),
                           as.character(end(anchorTwo(trans)))),
                itemRgb=rep("255,0,0", length(trans)),
                blockCount=1,
                blockSizes=c(as.character(width(anchorOne(trans))), as.character(width(anchorTwo(trans)))),
                blockStarts=0,
                stringsAsFactors = FALSE)

            if (!is.null(fn)) {
                write.table(output_cis, fn, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE )
                write.table(output_trans, fn, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE, append=TRUE)
                return(invisible(1))
			} else {
                return(rbind(output_cis, output_trans))
			}
})

.exportName = function(gi) {
    paste0(
        seqnames(anchorOne(gi)), ":",
        start(anchorOne(gi)) - 1 , "..",
        end(anchorOne(gi)), "-",
        seqnames(anchorTwo(gi)), ":",
        start(anchorTwo(gi)) - 1, "..",
        end(anchorTwo(gi)), ",",
        interactionCounts(gi))
}

#' Export interactions in BED Paired-End format.
#'
#' #' Exports a GInteractions object to BED-PE format, and writes to a specified file. If filename is not specified,
#' then a data.frame containing the information is returned. The value of the score parameter defines which field is used
#' to populate the score field.
#'
#'
#' @param GIObject A GInteractions object.
#' @param fn	   A filename to write the interactions data to
#' @param score    Which metadata column to use as score
#' @return invisible(1) if outputting to file or a data.frame containing all of the corresponding information
#'
#' @export
#' @docType methods
#' @rdname export.bedpe
#' @export
#' @examples 
#' data(hic_example_data)
#' export.bedpe(hic_example_data, fn = tempfile(), score = "counts")
setGeneric("export.bedpe", function(GIObject, fn=NULL, score="counts"){ standardGeneric("export.bedpe")} )
#' @rdname export.bedpe
#' @export
setMethod("export.bedpe", c("GInteractions"), function(GIObject, fn=NULL, score="counts"){
    score_vector = .getScore(GIObject, score)
    if (is.null(score_vector)) stop("Supplied score field not in element metadata.")
    output = cbind(as.character(seqnames(anchorOne(GIObject))),
                   start(anchorOne(GIObject))-1,
                   end(anchorOne(GIObject)),
                   as.character(seqnames(anchorTwo(GIObject))),
                   start(anchorTwo(GIObject))-1,
                   end(anchorTwo(GIObject)),
                   paste("interaction:", 1:length(GIObject), sep=""),
                   score_vector,
                   as.character(strand(anchorOne(GIObject))),
                   as.character(strand(anchorTwo(GIObject))))

    if(!is.null(fn)){
        write.table(output, fn, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE )
    }else{
        return(output)
    }

    return(invisible(1))
})


#' Export interactions in a BEDPE-like format for use with ChiaSig
#'
#' Exports a GInteractions object to BEDPE like format, (anchor specifications and a column for reads connecting them)
#' and writes to a specified file. If filename is not specified,
#' then a data.frame containing the information is returned. The value of the score parameter defines which field is used
#' to populate the score field.
#'
#'
#' @param GIObject A GInteractions object.
#' @param fn     A filename to write the interactions data to
#' @param score    Which metadata column to use as the score: counts or normalised
#' @return invisible(1) if outputting to file or a data.frame containing all of the corresponding information
#'
#' @export
#' @docType methods
#' @rdname export.chiasig
#' @export
#' @examples 
#' data(hic_example_data)
#' export.chiasig(hic_example_data, fn = tempfile(), score = "counts")
setGeneric("export.chiasig", function(GIObject, fn=NULL, score="counts"){ standardGeneric("export.chiasig")} )
#' @rdname export.chiasig
#' @export
setMethod("export.chiasig", c("GInteractions"), function(GIObject, fn=NULL, score="counts"){
    score_vec = .getScore(GIObject, score)
    if (is.null(score_vec)) stop("Supplied score field not in element metadata.")
    output = cbind(as.character(seqnames(anchorOne(GIObject))),
                    start(anchorOne(GIObject))-1,
                    end(anchorOne(GIObject)),
                    as.character(seqnames(anchorTwo(GIObject))),
                    start(anchorTwo(GIObject))-1,
                    score_vec)

    if(!is.null(fn)){
        write.table(output, fn, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE )
        return(invisible(1))
    }else{
        return(output)
    }
})

.getScore = function(x, score) {
    if (score=="counts")
        ans = interactionCounts(x)
    else
        ans = mcols(x)[[score]]
    ans
}

.getNames = function(x) {
    if ("name" %in% colnames(mcols(x)))
        names = mcols(x)[["name"]]
    else
        names = paste0("interaction_", 1:length(x))
    names
}

