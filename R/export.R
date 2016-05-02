#' Export interactions in BED12 format.
#'
#' @param GIObject  A GInteractions object.
#' @param fn        A filename to write the object to
#' @param score     Which metadata column to export as score
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
#' export.bed12(hic_example_data, fn = tempfile(), score = "counts")
#' @docType methods
#' @rdname export.bed12

setGeneric("export.bed12",function(GIObject, fn=NULL, score="counts"){standardGeneric ("export.bed12")})
#' @rdname export.bed12
#' @export
#' @importFrom utils write.table
#' @importFrom rtracklayer export
#' @importFrom grDevices col2rgb
setMethod("export.bed12", c("GInteractions"),
        function(GIObject, fn=NULL, score="counts"){
            bed = asBED(GIObject, score)
            if (!is.null(fn)) {
                export(bed, fn, format="bed")
                return(invisible(1))
			} else {
			    blocks = bed$blocks
			    bed$blocks = NULL
			    strand_info = as.character(strand(bed))
			    strand_info[strand_info == "*"] = NA
			    df = data.frame(
			        seqnames=seqnames(bed),
			        start=start(bed)-1,
			        end=end(bed),
			        names=bed$name,
			        score=bed$score,
			        strand=strand_info,
			        thickStart=bed$thickStart,
			        thickEnd=bed$thickEnd,
			        itemRgb=apply(col2rgb(bed$itemRgb), 2, paste, collapse=","),
			        blockCount=elementNROWS(blocks),
			        blockStarts=unlist(lapply(start(blocks), paste, collapse = ","), use.names=FALSE),
			        blockSizes=unlist(lapply(width(blocks), paste, collapse = ","), use.names=FALSE))
                return(bed)
			}
})

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
    if (is.null(ans))
        ans = rep(0, length(x))
    ans
}

.getNames = function(x) {
    if ("name" %in% colnames(mcols(x)))
        names = mcols(x)[["name"]]
    else
        names = paste0("interaction_", 1:length(x))
    names
}

#' Coerce to BED structure
#'
#' Coerce the structure of an object to one following BED-like
#' conventions, i.e., with columns for blocks and thick regions.
#'
#' @param x Generally, a tabular object to structure as BED
#' @param keep.mcols logical whether to keep non-BED12 columns in final
#'                   output (may cause problems with some parsers).
#' @param score character, which field to export as "score" in BED12.
#'              Defaults to "auto" which will choose score, then counts,
#'              if present, or fill column with zeros.
#'
#' @param ... Arguments to pass to methods
#'
#'      The exact behavior depends on the class of `object`.
#'
#'      `GRangesList` This treats `object` as if it were a list of
#'           transcripts, i.e., each element contains the exons of a
#'           transcript. The `blockStarts` and `blockSizes` columns are
#'           derived from the ranges in each element. Also, add `name`
#'           column from `names(object)`.
#'
#' @return A `GRanges`, with the metadata columns `name`, `blockStarts` and
#'         `blockSizes` added.
#'
#' @importFrom rtracklayer asBED
#' @importFrom IRanges PartitioningByWidth
#' @importFrom S4Vectors elementNROWS
#' @export
#' @docType methods
#' @examples
#' data(hic_example_data)
#' asBED(hic_example_data)
setMethod("asBED", c("GInteractions"),
    function(x, keep.mcols=FALSE, score="score") {
        if (!is.null(names(x)) || !is.null(x$name))
            warning("Names will be dropped during BED12 export")

        x = swapAnchors(x, mode="order")

        is_trans = as.vector(seqnames(regions(x))[x@anchor1] != seqnames(regions(x))[x@anchor2])
        a1_cis = x@anchor1[!is_trans]
        a2_cis = x@anchor2[!is_trans]
        a1_trans = x@anchor1[is_trans]
        a2_trans = x@anchor2[is_trans]

        scores = .getScore(x, score)

        if (is.null(x$color)) {
            x$color = "#000000"
        }

        names = .exportName(x, scores)

        cis_blocks = relist(IRanges(
                start=c(rbind(rep(1L, length(a1_cis)), start(x@regions)[a2_cis] - start(x@regions)[a1_cis] + 1L)),
                width=c(rbind(width(x@regions)[a1_cis], width(x@regions)[a2_cis]))),
            PartitioningByWidth(rep(2, length(a1_cis))))

        output_cis = GRanges(
            seqnames=as.character(seqnames(x@regions)[a1_cis]),
            IRanges(start=start(x@regions)[a1_cis],
                    end=end(x@regions)[a2_cis]),
            name=names[!is_trans],
            score=scores[!is_trans],
            strand=ifelse(
                strand(x@regions)[a1_cis] == strand(x@regions)[a2_cis] &
                       as.vector(strand(x@regions)[a1_cis]) %in% c("+", "-"),
                as.vector(strand(x@regions)[a1_cis]),
                "*"),
            thickStart=start(x@regions)[a1_cis],
            thickEnd=end(x@regions)[a2_cis],
            itemRgb=x$color[!is_trans],
            blocks=cis_blocks
        )

        trans_blocks = relist(IRanges(
                start=rep(1, 2*length(a1_trans)),
                width=c(width(x@regions)[a1_trans], width(x@regions)[a2_trans])),
            PartitioningByWidth(rep(1, 2*length(a1_trans))))

        output_trans = GRanges(
            seqnames=c(as.character(seqnames(x@regions)[a1_trans]),
                    as.character(seqnames(x@regions)[a2_trans])),
            IRanges(start=c(start(x@regions)[a1_trans],
                            start(x@regions)[a2_trans]),
                    end=c(end(x@regions)[a1_trans],
                            end(x@regions)[a2_trans])),
            name=rep(names[is_trans], 2),
            score=rep(scores[is_trans], 2),
            strand=c(as.character(strand(x@regions)[a1_trans]),
                        as.character(strand(x@regions)[a2_trans])),
            thickStart=c(start(x@regions)[a1_trans],
                         start(x@regions)[a2_trans]),
            thickEnd=c(end(x@regions)[a1_trans],
                       end(x@regions)[a2_trans]),
            itemRgb=rep(x$color[is_trans], 2),
            blocks=trans_blocks
        )

        extra_cols = setdiff(colnames(mcols(x)), c("score", "name"))

        if(length(extra_cols) && keep.mcols==TRUE) {
            mcols(output_cis) = cbind(mcols(output_cis), mcols(x)[!is_trans, extra_cols,drop=FALSE])
            mcols(output_trans) =
                cbind(mcols(output_trans),
                      rep(mcols(x)[is_trans, extra_cols,drop=FALSE], 2))
        }

        return(sort(c(output_cis, output_trans)))
})

.exportName = function(gi, score=0) {
    paste0(
        seqnames(gi@regions)[gi@anchor1], ":",
        start(gi@regions)[gi@anchor1] - 1 , "..",
        end(gi@regions)[gi@anchor1], "-",
        seqnames(gi@regions)[gi@anchor2], ":",
        start(gi@regions)[gi@anchor2] - 1, "..",
        end(gi@regions)[gi@anchor2], ",",
        score)
}
