.calcCoords = function(coordinates, plot.width, block.start, block.width){
  return(plot.width * (( coordinates - block.start) / block.width))
}

.calcYCoords = function(score, max.score, bottom, top){
  return( (score / max.score) * (top - bottom))
}

.calcMdpt = function(start, end){
  foo = ifelse( end > start, start + ((end - start) / 2),  end + (( start - end) / 2 ))
  return(foo)	 
}

.bezier.curve = function(p0x, p1x, p2x, n){
  t = seq(0,1,1/n)
  return( ((1 - t)^2 * p0x) + (2 * ( 1 - t ) * t * p1x ) + ( t^2 * p2x ) )
}

.specialrect = function(x.start, y.start, x.end, y.end, strand, col){
    for(i in 1:length(x.start)){
        if(as.character(strand[i])=="*"){
            rect(x.start[i], y.start[i], x.end[i], y.end[i], col=col)
        }else if(as.character(strand[i])=="+"){
            feature.width = x.end[i] - x.start[i]
            feature.height = y.end[i] - y.start[i]
            turn = x.start[i] + 0.5 * sqrt(3) * (0.9 * feature.width)
            xs = c(x.start[i], turn, x.end[i], turn, x.start[i])
            ys = c(y.start[i], y.start[i], y.start[i] + feature.height/2 , y.end[i], y.end[i])
            polygon(xs, ys, col=col)
        }else if(as.character(strand[i])=="-"){
            feature.width = x.end[i] - x.start[i]
            feature.height = y.end[i] - y.start[i]
            turn = x.end[i] - 0.5 * sqrt(3) * (0.9 * feature.width)
            xs = c(x.end[i], turn, x.start[i], turn, x.end[i])
            ys = c(y.start[i], y.start[i], y.start[i] + feature.height/2 , y.end[i], y.end[i]) 
            polygon(xs, ys, col=col)
        }   
    }
}



#' Plot interactions within a specified region. 
#' 
#' This is function allows the plotting of an interactions between annotated features in a specified area. 
#' The resulting plot shows unique interactions as curves between interaction anchor points with the number 
#' of counts supporting that interaction proportional to the thickness of that line. It is also possible to 
#' add cis-interactions which are not within the window/region and to also plot regions that are involved
#' in trans-interactions. Plotting the data this way makes it possible to examine a region and easily examine
#' which regions are highly interacting with each other. It is not recommended to use this style of plot to examine
#' regions larger than 5Mb. 
#' 
#' 
#' @param GIObject GenomicInteractions object
#' @param region A GRanges specifying the genomic region to plot 
#' @param annotation.features a list of GRanges specifying the features within the region to plot
#' @param annotation.cols a named vector specifying which colour to plot the individual tracks
#' @param reduce.anchors a logical specifying whether to reduce the anchor GRanges
#' @param plot.trans a logical specifying whether to show trans-interactions
#' @param plot.cis a logical specifying whether to plot cis-interactions that are outside of the specified region
#' @param order.cis logical specifying whether to order cis-interactions by their distances to the specified region
#' @param plot.cis.names a logical specifying whether to plot textual information on the other anchor region of cis-interactions
#' @param plot.header a logical specifying whether to plot a header describing the genomic region plotted
#' @param plot.lines a logical specifying whether to plot dashed lines indicating anchor regions and associated features
#' @param plot.ids a logical specifying whether to plot ids for features in annotation.features. Looks for the presence of an id or name
#' @param anchor.col colour for anchor regions
#' @return invisible(1)
#' @importFrom plotrix ablineclip
#' 
#' @docType methods
#' @rdname .plotRegion
setGeneric("plotRegion",function(GIObject, region, annotation.features, annotation.cols=NULL, reduce.anchors=TRUE,
                                 plot.trans=TRUE, plot.cis=TRUE, order.cis=TRUE, plot.cis.names=TRUE, plot.header=TRUE,
                                 plot.lines=TRUE, anchor.col="darkred", plot.ids=FALSE){standardGeneric ("plotRegion")})
#' @rdname .plotRegion
setMethod("plotRegion", 
    signature=c("GenomicInteractions", "GRanges", "list"), function(GIObject, region, annotation.features, annotation.cols=NULL, 
                       reduce.anchors=TRUE, plot.trans=TRUE, plot.cis=TRUE, order.cis=TRUE, plot.cis.names=TRUE, plot.header=TRUE,
                       plot.lines=TRUE, anchor.col = "darkred", plot.ids=FALSE
                       ){
    plot.width=1000
    spacer = 0.05
    block.start = start(region) 
    block.end = end(region)
    block.width=width(region)
    plot.height = length(annotation.features) + 4
    
    feature.locations =c( 0:(length(annotation.features)+1) )
    names(feature.locations) = c(rev(names(annotation.features)), "interactions")
    feature.sizes = c( header=0.5, interactions=3, anchor=0.15 )
    
    if(is.null(annotation.cols)){
        feature.cols =  rainbow(length(annotation.features))
    }else{
        feature.cols = annotation.cols
    }
    names(feature.cols) = names(annotation.features)
    
    par(oma=c(0.1, 0, 0.1, 0), mar=c(0, 0, 0, 0 ))
    plot.new()
    plot.window(c(-50,plot.width+50), c(0, plot.height))
      
    for(featureSet in names(annotation.features)){ 
        tmp.features.ol = NULL
        if(class(annotation.features[[featureSet]])=="GRanges"){ # deal with strand
            tmp.features.ol = restrict(annotation.features[[featureSet]][ unique(subjectHits(findOverlaps(region, annotation.features[[featureSet]])))], start(region), end(region))
            tmp.features.ol = unique(tmp.features.ol)
        }else if(class(annotation.features[[featureSet]])=="GRangesList"){
            tmp.features.ol = annotation.features[[featureSet]][ unique(subjectHits(findOverlaps(region, annotation.features[[featureSet]])))]
            tmp.granges = unlist(tmp.features.ol)
            elementMetadata(tmp.granges) = NULL
            ids = unlist(sapply(names(tmp.features.ol), function(x){rep(x, length(tmp.features.ol[[x]]))}))
            tmp.granges$id = ids 
            tmp.features.ol = restrict(tmp.granges, start(region), end(region))
        }else{
            stop("annotation.features must contain GRanges or GRangesLists")
        }
        if(length(tmp.features.ol) > 0){
            xcoord.start = .calcCoords(start(tmp.features.ol), plot.width, block.start, block.width)
            xcoord.end = .calcCoords(end(tmp.features.ol), plot.width, block.start, block.width)
            overlapping = countOverlaps(tmp.features.ol, tmp.features.ol, ignore.strand=TRUE)
            if(max(overlapping) > 1){
                x = reduce(tmp.features.ol, with.revmap=TRUE, ignore.strand=TRUE)
                ycoord.start = rep(0, length(tmp.features.ol))
                    
                multiples = sapply(1:length(tmp.features.ol), 
                                        function(y){unlist(lapply(x$revmap, 
                                            function(z){ if(y %in% unlist(z)){
                                                        return(which( y == unlist(z)))
                                                        }}
                                        ))}
                                    )
                ycoord.start = feature.locations[featureSet]  + ((multiples - 1) * (1 / max(overlapping))) + spacer
                ycoord.end = ycoord.start + 1 / max(overlapping)
            }else{
                ycoord.start = feature.locations[featureSet] + spacer
                ycoord.end = ycoord.start + 1 / max(overlapping)
            }
            .specialrect(xcoord.start, ycoord.start, xcoord.end, ycoord.end, as.character(strand(tmp.features.ol)), col=feature.cols[featureSet])
            # add identifiers
            if(plot.ids & !is.null(names(tmp.features.ol))){
                text( xcoord.start + (xcoord.end - xcoord.start)/2, ycoord.start + ((ycoord.end - ycoord.start)/2), names(tmp.features.ol), cex=0.6, offset=2 )
            }else if(plot.ids & "id" %in% names(elementMetadata(tmp.features.ol))){
                text( xcoord.start + (xcoord.end - xcoord.start)/2, ycoord.start + ((ycoord.end - ycoord.start)/2), tmp.features.ol$id, cex=0.6, offset=2 )
            }
        }else{
            text( plot.width/2, feature.locations[featureSet]+ 0.5, "N/A", cex=0.6)
        }
        text(-5, feature.locations[featureSet]+ 1/2, featureSet, pos=2, cex=0.6)
    }

    interaction.indexes = unique(c(subjectHits(findOverlaps(region, anchorOne(GIObject))), subjectHits(findOverlaps(region, anchorTwo(GIObject)))))
    tmp.interactions = GIObject[interaction.indexes]

    # determine whether interactions are trans, within the region or cis-interacting on the left or right side of the region
    one.status = ifelse(start(anchorOne(tmp.interactions)) < block.start & as.character(seqnames(anchorOne(tmp.interactions))) == as.character(seqnames(region)), "cis-left", 
                    ifelse(end(anchorOne(tmp.interactions)) > block.end & as.character(seqnames(anchorOne(tmp.interactions))) == as.character(seqnames(region)), "cis-right",
                        ifelse(as.character(seqnames(anchorOne(tmp.interactions))) != as.character(seqnames(region)), "trans", "inside")))
    two.status = ifelse(start(anchorTwo(tmp.interactions)) < block.start & as.character(seqnames(anchorTwo(tmp.interactions))) == as.character(seqnames(region)), "cis-left", 
                    ifelse(end(anchorTwo(tmp.interactions)) > block.end & as.character(seqnames(anchorTwo(tmp.interactions))) == as.character(seqnames(region)), "cis-right",
                        ifelse(as.character(seqnames(anchorTwo(tmp.interactions))) != as.character(seqnames(region)), "trans", "inside")))
  
    cis.left.reduced = reduce(c(anchorOne(tmp.interactions)[one.status=="cis-left"],	anchorTwo(tmp.interactions)[two.status=="cis-left"]))
    if(order.cis){ cis.left.reduced = cis.left.reduced[ order(end(cis.left.reduced), decreasing=TRUE) ] }
    cis.right.reduced = reduce(c(anchorOne(tmp.interactions)[one.status=="cis-right"],	anchorTwo(tmp.interactions)[two.status=="cis-right"]))
    if(order.cis){ cis.right.reduced = cis.right.reduced[ order(end(cis.right.reduced), decreasing=TRUE) ] }
    trans.reduced = reduce(c(anchorOne(tmp.interactions)[one.status=="trans"], anchorTwo(tmp.interactions)[two.status=="trans"]))
  
    all.anchors = restrict(c(anchorOne(GIObject)[ interaction.indexes ], anchorTwo(GIObject)[interaction.indexes]), start(region), end(region))
    if(reduce.anchors){
        all.anchors = restrict(reduce(all.anchors[  unique(queryHits(findOverlaps(all.anchors, region)))]), start(region), end(region))
    }else{ 
        all.anchors = restrict(all.anchors[  unique(queryHits(findOverlaps(all.anchors, region)))], start(region), end(region))
    }
    
    left.y.ends = seq( feature.locations["interactions"] + feature.sizes["interactions"], plot.height-feature.sizes["header"]-1, length.out=length(cis.left.reduced)) 
    right.y.ends = seq( feature.locations["interactions"] + feature.sizes["interactions"], plot.height-feature.sizes["header"]-1, length.out=length(cis.right.reduced)) 
  
    # plot trans-interactions 
    if(plot.trans){
        trans.interaction.indexes = which( one.status == "trans" | two.status == "trans")
        for(i in trans.interaction.indexes){
            a.ol = findOverlaps(all.anchors[i], tmp.interactions) 
            a.1.ol = a.ol[["one"]]
            a.2.ol = a.ol[["two"]]
            if( sum(two.status[ unique(subjectHits(a.1.ol)) ] == "trans")>0 ){
                a.1.counts = sum( tmp.interactions[ two.status[ unique(subjectHits(a.1.ol)) ] == "trans" ])
            }else{ 
                a.1.counts = 0 
            }
            if( sum(one.status[ unique(subjectHits(a.2.ol)) ] == "trans")>0 ){
                a.2.counts = sum( tmp.interactions[ one.status[ unique(subjectHits(a.2.ol)) ] == "trans" ])
            }else{
                a.2.counts = 0
            }
            sum.counts = a.1.counts + a.2.counts
            if(sum.counts > 0){
                p0x = .calcCoords(.calcMdpt(start(all.anchors[ i ]), end(all.anchors[ i ])), plot.width, block.start, block.width)
                if(p0x < (plot.width/2)){
                    p2x = 0
                    p1x = p0x * (2/3)   
                }else{
                    p2x = 1000
                    p1x = (( plot.width - p0x ) / 3) + p0x
                }
                x.lines = .bezier.curve(p0x, p1x, p2x, 1000) 
                y.lines = .bezier.curve(feature.locations["interactions"] + (feature.sizes["anchor"]/2) + spacer, 
                                        plot.height - ifelse(plot.header, feature.sizes["header"], 0), 
                                        plot.height - ifelse(plot.header, feature.sizes["header"], 0) - 1, 1000)
                lines(x.lines, y.lines, type="l", col="lightgray", lwd=log(sum.counts))
            }
        }
    }    
    
    if(plot.cis){
        cis.interaction.indexes = which(one.status == "cis-left" | two.status == "cis-left" | one.status == "cis-right" | two.status == "cis-right")
        for(i in cis.interaction.indexes){
            if(((one.status[i] == "inside" & two.status[i] == "cis-left") | (one.status[i] == "cis-left" & two.status[i] == "inside"))){
                p2x = 0
                if(one.status[i] == "inside"){
                    p0x = .calcCoords(.calcMdpt(start(anchorOne(tmp.interactions)[ i ]), end(anchorOne(tmp.interactions)[ i ])), plot.width, block.start, block.width)
                    p2y = left.y.ends[  subjectHits(findOverlaps( anchorTwo(tmp.interactions)[i], cis.left.reduced )) ] 
                }else if(two.status[i] == "inside"){
                    p0x = .calcCoords(.calcMdpt(start(anchorTwo(tmp.interactions)[ i ]), end(anchorTwo(tmp.interactions)[ i ])), plot.width, block.start, block.width)											   
                    p2y = left.y.ends[  subjectHits(findOverlaps( anchorOne(tmp.interactions)[i], cis.left.reduced )) ]
                }
                p1x = p0x * (2/3)
                x.lines = .bezier.curve(p0x, p1x, p2x, 1000) 
                y.lines = .bezier.curve(feature.locations["interactions"] + (feature.sizes["anchor"]/2) + spacer, 
                                        plot.height - ifelse(plot.header, feature.sizes["header"], 0) - 1.5, 
                                        p2y, 1000)	 													
                lines(x.lines, y.lines, type="l", col="darkslategray", lwd=log(count(tmp.interactions)[i]))
            }else if(((one.status[i] == "inside" & two.status[i] == "cis-right") | (one.status[i] == "cis-right" & two.status[i] == "inside"))){
                p2x = 1000
                if(one.status[i] == "inside"){
                    p0x = .calcCoords(.calcMdpt(start(anchorOne(tmp.interactions)[ i ]), end(anchorOne(tmp.interactions)[ i ])), plot.width, block.start, block.width)
                    p2y = right.y.ends[  subjectHits(findOverlaps( anchorTwo(tmp.interactions)[i], cis.right.reduced )) ] 
                }else if(two.status[i] == "inside"){
                    p0x = .calcCoords(.calcMdpt(start(anchorTwo(tmp.interactions)[ i ]), end(anchorTwo(tmp.interactions)[ i ])), plot.width, block.start, block.width)
                    p2y = right.y.ends[  subjectHits(findOverlaps( anchorOne(tmp.interactions)[i], cis.right.reduced )) ] 
                }
                p1x = (( plot.width - p0x ) / 3) + p0x
                x.lines = .bezier.curve(p0x, p1x, p2x, 1000) 
                y.lines = .bezier.curve(feature.locations["interactions"] + (feature.sizes["anchor"]/2) + spacer, 
                                        plot.height - ifelse(plot.header, feature.sizes["header"], 0) - 1.5, 
                                        p2y, 1000) 
                lines(x.lines, y.lines, type="l", col="darkslategray", lwd=log(count(tmp.interactions)[i]))
            }
        }
    }
    
    internal.interaction.indexes = which( one.status == "inside" & two.status == "inside")
    for(i in internal.interaction.indexes){
        if(one.status[i] == "inside" & two.status[i] == "inside"){
            p0x = .calcCoords(.calcMdpt(start(anchorOne(tmp.interactions)[ i ]), end(anchorOne(tmp.interactions)[ i ])), plot.width, block.start, block.width)
            p2x = .calcCoords(.calcMdpt(start(anchorTwo(tmp.interactions)[ i ]), end(anchorTwo(tmp.interactions)[ i ])), plot.width, block.start, block.width)
            p1x = .calcMdpt(p0x, p2x)
            x.lines = .bezier.curve(p0x, p1x, p2x, 1000) 
            y.lines = .bezier.curve(feature.locations["interactions"] + (feature.sizes["anchor"]/2) + spacer, 
                                    feature.locations["interactions"] + feature.sizes["interactions"], 
                                    feature.locations["interactions"] + (feature.sizes["anchor"]/2) + spacer, 1000)
            lines(x.lines, y.lines, type="l", col="black", lwd=log(count(tmp.interactions)[i]))
        }
    }
    #if(reduce.anchors==TRUE){
    xcoord.start = .calcCoords(start(all.anchors), plot.width, block.start, block.width)
    xcoord.end = .calcCoords(end(all.anchors), plot.width, block.start, block.width)
    rect(xcoord.start, feature.locations["interactions"] + spacer, xcoord.end, feature.locations["interactions"] + feature.sizes["anchor"] + spacer, col=anchor.col)
    if(plot.lines){ablineclip(v=unique(xcoord.start), y1=0, y2=plot.height-ifelse(plot.header, feature.sizes["header"], 0), lty=2,lwd=0.5)}
    if(plot.lines){ablineclip(v=unique(xcoord.end), y1=0, y2=plot.height-ifelse(plot.header, feature.sizes["header"], 0), lty=2, lwd=0.5)}
    #}else{
    #    xcoord.start = c(.calcCoords(start(anchorOne(tmp.interactions)[one.status == "inside"]), plot.width, block.start, block.width), 
    #                     .calcCoords(start(anchorTwo(tmp.interactions)[two.status == "inside"]), plot.width, block.start, block.width))
    #    xcoord.end = c(.calcCoords(end(anchorOne(tmp.interactions)[one.status == "inside"]), plot.width, block.start, block.width), 
    #                   .calcCoords(end(anchorTwo(tmp.interactions)[two.status == "inside"]), plot.width, block.start, block.width))
    #    rect(xcoord.start, feature.locations["interactions"] + spacer, xcoord.end, feature.locations["interactions"]+ feature.sizes["anchor"] + spacer, col=anchor.col)
    #    if(plot.lines){ablineclip(v=unique(xcoord.start), y1=0, y2=plot.height-ifelse(plot.header, feature.sizes["header"], 0), lty=2,lwd=0.5)}
    #    if(plot.lines){ablineclip(v=unique(xcoord.end), y1=0, y2=plot.height-ifelse(plot.header, feature.sizes["header"], 0), lty=2, lwd=0.5)}
    #}
  
    if(plot.cis & plot.cis.names){
        if(length(cis.left.reduced) > 0){
            text( -4, left.y.ends, labels=paste(seqnames(cis.left.reduced), ":", start(cis.left.reduced),  "-", end(cis.left.reduced), sep=""), pos=2, cex=0.45) 
        }	   					   					   	      	     
        if(length(cis.right.reduced) > 0){	
            text( 1004, right.y.ends, labels=paste(seqnames(cis.right.reduced), ":", start(cis.right.reduced),  "-", end(cis.right.reduced), sep=""), pos=4, cex=0.45) 
        } 
    }
    if(plot.header){
        rect(0, plot.height-feature.sizes["header"], plot.width, plot.height, col="gray")
        text(plot.width/2, plot.height-( feature.sizes["header"] /2), paste(seqnames(region), ":", start(region), "-", end(region), sep=""))
    } 
    return(invisible(1))
})
