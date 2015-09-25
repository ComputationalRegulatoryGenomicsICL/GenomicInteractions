#' Virtual 4C viewpoint
#'
#' This function creates a GenomicInteractions object representing interactions 
#' originating at a given viewpoint ("bait"), or set of viewpoints. This is 
#' similar to the idea of a virtual 4C experiment where you are interested in 
#' interactions with a specific region. 
#' 
#' The object returned has the "bait" as anchor one, and the interacting regions 
#' as anchor two. By default this is genome wide. If you only want to consider 
#' interactions within a certain distance around the bait, you can specify a 
#' region to consider.
#' 
#' Multiple baits can be given, e.g. to find all interactions around promoters.
#' 
#' You may want to visualise the resulting interactions in a genome browser - 
#' you can do this by creating coverage over anchor two of the object and exporting 
#' as a wig or bedgraph file.
#' 
#' 
#' @param x A GenomicInteractions object.
#' @param bait A GRanges object describing bait regions.
#' @param region If present, a GenomicInteractions object specifying the
#'               region to look for bait interactions in.
#' @param ... additional arguments to findoverlaps
#' 
#' @return A GenomicInteractions object.
#' 
#' @import GenomicRanges
#' @export
#' @examples
#' \dontrun{
#' data(hic_data)
#' library(GenomicRanges)
#' pos <- GRanges(seqnames="chr5", ranges=IRanges(start=115938063, end=115941352))
#' region <- GRanges(seqnames="chr5", ranges=IRanges(start=115838063, end=116041352))
#' viewPoint(hic_data, pos, region)
#' }
viewPoint = function(x, bait, region=NULL, ...) {
    if (any(names(mcols(anchorOne(x))) != names(mcols(anchorTwo(x))))) {
        warning("Metadata differs between anchors and will be lost.")
        mcols(x@anchor_one) = NULL
        mcols(x@anchor_two) = NULL
    }
    hits = findOverlaps(x, bait, ...)
    if (length(hits) == 0) return(GenomicInteractions())
    vp = GenomicInteractions(anchor_one=bait[c(subjectHits(hits$one), 
                                          subjectHits(hits$two))],
                             anchor_two=c(x@anchor_two[queryHits(hits$one)], 
                                          x@anchor_one[queryHits(hits$two)]),
                             counts=c(x@counts[queryHits(hits$one)], 
                                      x@counts[queryHits(hits$two)]))
    if (!is.null(region)) { vp = x[overlapsAny(x@anchor_two, region, ...)] }
    ord = order(vp@anchor_one, vp@anchor_two)
    return(vp[ord])
}

#' Plot coverage around a virtual 4C viewpoint
#' 
#' Plots coverage of interactions around a given viewpoint. This function requires 
#' the output of `viewPoint()` as input. You should additionally specify the total 
#' region you wish to plot. 
#'
#' @param x a GenomicInteractions object which is output from viewPoint
#' @param region The genomic region to plot
#' @param ylab Y axis label.
#' @param xlab X axis label. By default this is the chromosome of the region 
#' that is being plotted.
#' @param ... additional arguments to plot
#' 
#' @return Coverage that is plotted (invisibly)
#'
#' @import GenomicRanges
#' @export
plotViewpoint = function(x, region, ylab="Signal", xlab=NULL, ...) {
    if (length(region) > 1) stop("region must be a single range")
    x = x[overlapsAny(x@anchor_two, region, type="within")]
    cov = as(coverage(x@anchor_two, weight = interactionCounts(x))[region], "GRanges")
    points_x = c(start(cov), end(cov)) + start(region)
    points_y = rep.int(cov$score, 2)
    ord = order(points_x)
    
    if(is.null(xlab)){ 
    xlab <- as.character(seqnames(region))
    }
    
    plot(points_x[ord], points_y[ord], t="l", 
             ylab = ylab, xlab = xlab, ...)
    
    return(invisible(coverage(cov, weight = "score")))
}

#' Plot coverage around a set of virtual 4C viewpoints
#' 
#' Plots summarised coverage of interactions around a set of viewpoints, 
#' e.g. promoters. This function requires the output of `viewPoint()` as input. 
#'
#' @param x A GenomicInteractions object which is output from viewPoint
#' @param left_dist Distance 'left' of interactions to consider, in bp.
#' @param right_dist Distance 'right' of interactions to consider, in bp.
#' @param fix One of "center", "start", "end". Passed to `resize`. Interaction 
#' distances are calculated relative to this part of the bait. 
#' @param ylab Y axis label.
#' @param xlab X axis label. 
#' @param ... additional arguments to plot
#' 
#' @return  Coverage that is plotted (invisibly)
#'
#' @import GenomicRanges
#' @export
plotAvgViewpoint = function(x, left_dist = 100000, right_dist = 100000, ylab="Average signal", 
                            xlab="Relative position", fix = "center",...) {
  
  cov <- .makeRelativeVP(x, fix = fix, left_dist = left_dist, right_dist = right_dist )
  
  points_x = c(start(cov), end(cov)) - left_dist
  points_y = rep.int(runValue(cov), 2)
  ord = order(points_x)
  
  p = plot(points_x[ord], points_y[ord], t="l", 
           ylab = ylab, xlab = xlab, ...)
  
  return(invisible(cov))
}

.makeRelativeVP <- function(GIObject, fix = "center", left_dist = 100000, right_dist = 100000){
  #make ranges relative to bait
  bait <- resize(anchorOne(GIObject), width = 1, fix = fix)
  ints <- ranges(anchorTwo(GIObject))
  ints <- GenomicRanges::shift(ints, shift = -(start(bait)))
  ints <- GenomicRanges::shift(ints, shift = left_dist)
  
  #make coverage and adjust to mean coverage per bait
  cov <- coverage(ints, weight = interactionCounts(GIObject), 
                  width = right_dist+left_dist)
  cov <- cov / length(unique(bait))
  
  return(cov)
  #get points to plot
  
}