#' Virtual 4C viewpoint
#'
#' Some explanation needed
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
#' @param x a GenomicInteractions object which is output from viewPoint
#' @param region The genomic region to plot
#' @param ... additional arguments to plot
#' 
#' @return output of plot()
#'
#' @import GenomicRanges
#' @export
plotViewpoint = function(x, region, ...) {
    if (length(region) > 1) stop("region must be a single range")
    x = x[overlapsAny(x@anchor_two, region, type="within")]
    cov = as(coverage(x@anchor_two)[region], "GRanges")
    points_x = c(start(cov), end(cov)) + start(region)
    points_y = rep.int(x$score, 2)
    ord = order(points_x)
    p = plot(points_x[ord], points_y[ord], t="l", ...)
    return(p)
}

