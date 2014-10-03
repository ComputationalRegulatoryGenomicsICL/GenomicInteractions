#' Export interactions in BED12 format.
#' 
#' @param GIObject  A GenomicInteractions object.
#' @param fn        A filename to write the object to
#' 
#' Exports a GenomicInteractions object to BED12 format, and writes to a specified file. If filename is not specified, 
#' then a data.frame containing the information is returned. Please note some large datasets may take a long time to export.
#' 
#' @return invisible(1) if outputting to file or a data.frame containing all of the corresponding information
#' @export
#' @docType methods
#' @rdname export.bed12
#' @export
setGeneric("export.bed12",function(GIObject, fn=NULL){standardGeneric ("export.bed12")})
#' @rdname export.bed12
#' @export
setMethod("export.bed12", c("GenomicInteractions"), 
        function(GIObject, fn=NULL){
            
            cis.ind = is.cis(GIObject)
            trans.ind = is.trans(GIObject)
            cis.length = sum(cis.ind)
            trans.length = sum(trans.ind)
            obj.length = cis.length + ( 2 * trans.length) 

            output = data.frame(chr=character(obj.length), start=integer(obj.length), end=integer(obj.length), name=character(obj.length),
                                score=integer(obj.length), strand=rep(".", obj.length), thickStart=integer(obj.length), thickEnd=integer(obj.length),
                                itemRgb = rep("255,0,0", obj.length), blockCount=integer(obj.length), blockSizes = character(obj.length), blockStarts=character(obj.length),
                                #signalValue=character(obj.length), pValue=character(obj.length), qValue=character(obj.length), 
                                stringsAsFactors=FALSE)
            output[1:cis.length, "chr"] = as.character(seqnames(anchorOne(GIObject)))[cis.ind] 
            output[1:cis.length, "start"] = apply(cbind(as.numeric(start(anchorOne(GIObject)))[cis.ind], as.numeric(start(anchorTwo(GIObject)))[cis.ind]), 1, min) -1
            output[1:cis.length, "end"] = apply(cbind(as.numeric(end(anchorOne(GIObject)))[cis.ind], as.numeric(end(anchorTwo(GIObject)))[cis.ind] ), 1, max)
            output[1:cis.length, "name"] = .generateInteractionName(
                                                .pasteAnchor(as.character(seqnames(anchorOne(GIObject)))[cis.ind],
                                                             as.character(as.numeric(start(anchorOne(GIObject)))[cis.ind]-1),
                                                             as.character(as.numeric(end(anchorOne(GIObject)))[cis.ind])),
                                                .pasteAnchor(as.character(seqnames(anchorTwo(GIObject)))[cis.ind],
                                                             as.character(as.numeric(start(anchorOne(GIObject)))[cis.ind]-1),
                                                             as.character(as.numeric(end(anchorOne(GIObject)))[cis.ind])),
                                                as.numeric(count(GIObject)[cis.ind])
                                                )
            output[1:cis.length, "score"] = ifelse(count(GIObject)[cis.ind] > 10, 1000, 100 * count(GIObject)[cis.ind])
            output[1:cis.length, "thickStart"] = apply(cbind(as.numeric(start(anchorOne(GIObject)))[cis.ind], as.numeric(start(anchorTwo(GIObject)))[cis.ind]), 1, min) -1
            output[1:cis.length, "thickEnd"] = apply(cbind(as.numeric(end(anchorOne(GIObject)))[cis.ind], as.numeric(end(anchorTwo(GIObject)))[cis.ind] ), 1, max)
            output[1:cis.length, "blockCount"] = 2
            output[1:cis.length, "blockSizes"] = ifelse(as.numeric(start(anchorOne(GIObject)))[cis.ind] < as.numeric(start(anchorTwo(GIObject)))[cis.ind], 
                                                    paste(as.character(width(anchorOne(GIObject))[cis.ind]), as.character(width(anchorTwo(GIObject))[cis.ind]), sep=","), 
                                                    paste(as.character(width(anchorTwo(GIObject))[cis.ind]), as.character(width(anchorOne(GIObject))[cis.ind]), sep=","))
            output[1:cis.length, "blockStarts"] = ifelse(as.numeric(start(anchorOne(GIObject)))[cis.ind] < as.numeric(start(anchorTwo(GIObject)))[cis.ind], 
                                                    paste(0, as.character(as.numeric(start(anchorTwo(GIObject)))[cis.ind] - as.numeric(start(anchorOne(GIObject)))[cis.ind]), sep=","),
                                                    paste(0, as.character(as.numeric(start(anchorOne(GIObject)))[cis.ind] - as.numeric(start(anchorTwo(GIObject)))[cis.ind]), sep=","))
            # now output trans-interactions
            output[(cis.length+1):nrow(output), "chr"] = c(as.character(seqnames(anchorOne(GIObject)))[trans.ind] , as.character(seqnames(anchorTwo(GIObject)))[trans.ind] ) 
            output[(cis.length+1):nrow(output), "start"] = c(as.numeric(start(anchorOne(GIObject)))[trans.ind]-1, as.numeric(start(anchorTwo(GIObject)))[trans.ind] -1)
            output[(cis.length+1):nrow(output), "end"] = c(as.numeric(end(anchorOne(GIObject)))[trans.ind], as.numeric(end(anchorTwo(GIObject)))[trans.ind])
            output[(cis.length+1):nrow(output), "name"] = rep(.generateInteractionName(
                                                                .pasteAnchor(as.character(seqnames(anchorOne(GIObject)))[trans.ind],
                                                                            as.character(as.numeric(start(anchorOne(GIObject)))[trans.ind]-1),
                                                                            as.character(as.numeric(end(anchorOne(GIObject)))[trans.ind])),
                                                                .pasteAnchor(as.character(seqnames(anchorTwo(GIObject)))[trans.ind],
                                                                            as.character(as.numeric(start(anchorOne(GIObject)))[trans.ind]-1),
                                                                            as.character(as.numeric(end(anchorOne(GIObject)))[trans.ind])),
                                                                as.numeric(count(GIObject)[trans.ind])), 2)
            output[(cis.length+1):nrow(output), "score"] = rep(ifelse(count(GIObject)[trans.ind] > 10, 1000, 100 * count(GIObject)[trans.ind]), 2)
            output[(cis.length+1):nrow(output), "thickStart"] = c(as.numeric(start(anchorOne(GIObject)))[trans.ind]-1, as.numeric(start(anchorTwo(GIObject)))[trans.ind] -1)
            output[(cis.length+1):nrow(output), "thickEnd"] = c(as.numeric(end(anchorOne(GIObject)))[trans.ind], as.numeric(end(anchorTwo(GIObject)))[trans.ind])
            output[(cis.length+1):nrow(output), "blockCount"] = 1
            output[(cis.length+1):nrow(output), "blockSizes"] = c(as.character(width(anchorOne(GIObject))[trans.ind]), as.character(width(anchorTwo(GIObject))[trans.ind]))
            output[(cis.length+1):nrow(output), "blockStarts"] = 0
            #output[(cis.length+1):nrow(output), "signalValue"] = rep(count(GIObject)[trans.ind], 2)
                    
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
#' @param score    Which metadata column to use as the score: counts, pvalue, fdr, normalised
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
    output = cbind(as.character(seqnames(anchorOne(GIObject))), start(anchorOne(GIObject))-1, end(anchorOne(GIObject)),
                    as.character(seqnames(anchorTwo(GIObject))), start(anchorTwo(GIObject))-1, end(anchorTwo(GIObject)),paste("interaction:", 1:length(GIObject), sep=""))
    if(score == "counts"){
        output = cbind(output, count(GIObject))
    }else if(score=="fdr"){
        output = cbind(output, FDR(GIObject))
    }else if(score=="pvalue"){ 
        output = cbind(output, pValue(GIObject))
    }else if(score=="normalised"){
        output = cbind(output, normalisedCount(GIObject))
    }else{
        output = cbind(output, ".")
    }
    
    output = cbind(output, as.character(strand(anchorOne(GIObject))), as.character(strand(anchorTwo(GIObject))))
    if(!is.null(fn)){
        write.table(output, fn, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE )
    }else{
        return(output)
    }
    return(invisible(1))
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
    edges = data.frame(from=names[mapping[,1]], to=names[mapping[,2]], 
                       weight=GIObject@counts[mapping[,3]], 
                       count=GIObject@counts[mapping[,3]], 
                       pval=GIObject@pvalue[mapping[,3]], 
                       fdr=GIObject@fdr[mapping[,3]]) 

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

