#' Export interactions in BED12 format.
#'
#' @param GIObject  A GenomicInteractions object.
#' @param fn        A filename to write the object to
#' @param score     Which metadata column to export as score
#' @param drop.trans Logical indicating whether to allow dropping trans interactions (otherwise an error is produced)
#'
#' Exports a GenomicInteractions object to BED12 format, and writes to a specified file. If filename is not specified,
#' then a data.frame containing the information is returned. Please note some large datasets may take a long time to export.
#'
#' @return invisible(1) if outputting to file or a data.frame containing all of the corresponding information
#' @export
#' @docType methods
#' @rdname export.bed12
#' @export
setGeneric("export.bed12",function(GIObject, fn=NULL, score="counts", drop.trans=c(FALSE, TRUE)){standardGeneric ("export.bed12")})
#' @rdname export.bed12
#' @export
setMethod("export.bed12", c("GenomicInteractions"),
        function(GIObject, fn=NULL, score="counts", drop.trans=c(FALSE, TRUE)){

            is_trans = is.trans(GIObject)

            if (any(is_trans)) {
                if (drop.trans==FALSE) stop("trans interactions present in object and drop.trans == FALSE")
            }

            GIObject = sort(GIObject[!is_trans], order.interactions=FALSE)

            lx = length(GIObject)
            names = .getNames(GIObject)
            score = .getScore(GIObject, score)
            s1 = strand(anchorOne(GIObject))
            s2 = strand(anchorTwo(GIObject))

            output = data.frame(chr=as.character(seqnames(anchorOne(GIObject))),
                                start=start(anchorOne(GIObject))-1,
                                end=end(anchorTwo(GIObject)),
                                name=names,
                                score=score,
                                strand=ifelse(s1 == s2 & s1 %in% c("+", "-"), s1, "."), # avoid case where strand == "*"
                                thickStart=start(anchorOne(GIObject))-1,
                                thickEnd=end(anchorTwo(GIObject)),
                                itemRgb=rep("255,0,0", lx),
                                blockCount=2,
                                blockSizes=paste(as.character(width(anchorOne(GIObject))), as.character(width(anchorTwo(GIObject))), sep=","),
                                blockStarts=paste(0, start(anchorTwo(GIObject)) - start(anchorOne(GIObject)), sep=","),
                                stringsAsFactors=FALSE)

            if(!is.null(fn)){
			    write.table(output, fn, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE )
			}else{
			    return(output)
			}

			return(invisible(1))
})

#' Export interactions in BED Paired-End format.
#'
#' #' Exports a GenomicInteractions object to BED-PE format, and writes to a specified file. If filename is not specified,
#' then a data.frame containing the information is returned. The value of the score parameter defines which field is used
#' to populate the score field.
#'
#'
#' @param GIObject A GenomicInteractions object.
#' @param fn	   A filename to write the interactions data to
#' @param score    Which metadata column to use as score
#' @return invisible(1) if outputting to file or a data.frame containing all of the corresponding information
#'
#' @export
#' @docType methods
#' @rdname export.bedpe
#' @export
setGeneric("export.bedpe", function(GIObject, fn=NULL, score="counts"){ standardGeneric("export.bedpe")} )
#' @rdname export.bedpe
#' @export
setMethod("export.bedpe", c("GenomicInteractions"), function(GIObject, fn=NULL, score="counts"){
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
#' Exports a GenomicInteractions object to BEDPE like format, (anchor specifications and a column for reads connecting them)
#' and writes to a specified file. If filename is not specified,
#' then a data.frame containing the information is returned. The value of the score parameter defines which field is used
#' to populate the score field.
#'
#'
#' @param GIObject A GenomicInteractions object.
#' @param fn     A filename to write the interactions data to
#' @param score    Which metadata column to use as the score: counts or normalised
#' @return invisible(1) if outputting to file or a data.frame containing all of the corresponding information
#'
#' @export
#' @docType methods
#' @rdname export.chiasig
#' @export
setGeneric("export.chiasig", function(GIObject, fn=NULL, score="counts"){ standardGeneric("export.chiasig")} )
#' @rdname export.chiasig
#' @export
setMethod("export.chiasig", c("GenomicInteractions"), function(GIObject, fn=NULL, score="counts"){
    score_vec = .getScore(GIObject, score)
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

#' Export interactions to an igraph object.
#'
#' Exports a GenomicInteractions object to graph.data.frame for use by igraph package. This uses unique anchors
#' as nodes and generates edges between them. For the resulting graph to be easily interpretable, anchors
#' should be non-overlapping. This should already be the case for HiC data (either binned or restriction
#' fragments), however ChIA-PET data can contain overlapping anchors, which may need to be reduced to
#' non-overlapping regions before graph export.
#' @param GIObject A GenomicInteractions object.
#'
#' @return a graph.data.frame representation of the GenomicInteractions object
#' @importFrom igraph graph.data.frame
#'
#' @export
#' @docType methods
#' @rdname export.igraph
#' @export
setGeneric("export.igraph",function(GIObject){standardGeneric ("export.igraph")})
#' @rdname export.igraph
#' @export
setMethod("export.igraph", "GenomicInteractions", function(GIObject){
    nodes = unique( c(anchorOne(GIObject), c(anchorTwo(GIObject))))
    gi.nodes.ol = findOverlaps(GIObject, nodes, type="equal")

    mapping = cbind(subjectHits(gi.nodes.ol[["one"]]), subjectHits(gi.nodes.ol[["two"]]), queryHits(gi.nodes.ol[["one"]]))
    names =  paste(paste(paste(as.character(seqnames(nodes))), as.character(start(nodes)), sep=":"), as.character(end(nodes)), sep="..")
    edge_mcols = mcols(GIObject)
    edge_mcols$count = interactionCounts(GIObject)
    edges = cbind(data.frame(from=names[mapping[,1]], to=names[mapping[,2]]), edge_mcols)

    verts = data.frame(name=names)
    if("node.class" %in% unique(c(names(elementMetadata(anchorOne(GIObject))), names(elementMetadata(anchorTwo(GIObject)))))){
        node.class = sapply(1:length(nodes), function(x){ paste(unique(c(anchorOne(GIObject)$node.class[unique(queryHits(gi.nodes.ol[["one"]][ subjectHits(gi.nodes.ol[["one"]]) == x ]))],
                                                                anchorTwo(GIObject)$node.class[unique(queryHits(gi.nodes.ol[["two"]][ subjectHits(gi.nodes.ol[["two"]]) == x ]))])),
                                                                collapse=",")})
        verts = data.frame(names=names, nodeclass=node.class)
        potential.node.classes = unique(c(GIObject@anchor_one$node.class, GIObject@anchor_two$node.class))
        potential.node.classes = potential.node.classes[ potential.node.classes != "distal" ]
        for(i in potential.node.classes){
            verts[, paste(i, "id", sep=".")] = sapply(1:length(nodes),
                                                    function(x){ paste(unique(c(
                                                            unlist(elementMetadata(anchorOne(GIObject))[[paste(i, "id", sep=".")]][
                                                                    unique(queryHits(gi.nodes.ol[["one"]][ subjectHits(gi.nodes.ol[["one"]]) == x ]))
                                                                    ]),
                                                            unlist(elementMetadata(anchorTwo(GIObject))[[paste(i, "id", sep=".")]][
                                                                    unique(queryHits(gi.nodes.ol[["two"]][ subjectHits(gi.nodes.ol[["two"]]) == x ]))
                                                                    ])
                                                            )), collapse=",")
                                                            })
        }
    }
    return(graph.data.frame(edges, directed=FALSE, vertices = verts))
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

