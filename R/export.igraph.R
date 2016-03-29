
#' Export interactions to an igraph object.
#'
#' Exports a GInteractions object to graph.data.frame for use by igraph package. This uses unique anchors
#' as nodes and generates edges between them. For the resulting graph to be easily interpretable, anchors
#' should be non-overlapping. This should already be the case for HiC data (either binned or restriction
#' fragments), however ChIA-PET data can contain overlapping anchors, which may need to be reduced to
#' non-overlapping regions before graph export.
#' @param GIObject A GInteractions object.
#'
#' @examples
#' data(hic_example_data)
#' ig <- export.igraph(hic_example_data)
#' 
#' @return a graph.data.frame representation of the GInteractions object
#' @importFrom igraph graph_from_data_frame
#'
#' @export
#' @docType methods
#' @rdname export.igraph
#' @export
#' 
setGeneric("export.igraph",function(GIObject){standardGeneric ("export.igraph")})
#' @rdname export.igraph
#' @export
setMethod("export.igraph", "GInteractions", function(GIObject){
  
  if(!is.null(names(regions(GIObject)))){
    nodes <- names(regions(GIObject))
  } else {
    nodes <- as.character(regions(GIObject))
  }
  nodes <- data.frame(name = nodes, mcols(regions(GIObject)))
  
  edges <- data.frame(from = nodes$name[GIObject@anchor1],
                      to = nodes$name[GIObject@anchor2])
  edges <- cbind(edges, mcols(GIObject))
  
  return(graph_from_data_frame(edges, directed=FALSE, vertices = nodes))
})
