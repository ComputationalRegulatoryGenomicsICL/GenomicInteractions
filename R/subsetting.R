#' Subset a GInteractions object by features
#'
#' Subsets interactions for which at least one of the anchors overlaps with a given GRanges object.
#' Alternatively, subsets interactions based on annotated feature IDs for a particular feature.
#'
#' @rdname GInteractions-subsetByFeatures-methods
#' @docType methods
#' @param GIObject A GInteractions object
#' @param features A GRanges or GRangesList object, or a character vector containing
#'                      IDs of annotated features, e.g. promoter IDs.
#' @param feature.class If `features' is a character vector, the corresponding feature name, e.g. "promoter".
#' @return a subsetted GInteractions object
#' @export
#' 
#' @examples 
#' data("hic_example_data")
#' data("mm9_refseq_promoters")
#' annotateInteractions(hic_example_data, list(promoter = mm9_refseq_promoters))
#' ids <- names(mm9_refseq_promoters[1:10])
#' subsetByFeatures(hic_example_data, ids, "promoter")
setGeneric("subsetByFeatures",function(GIObject, features, feature.class=NULL){standardGeneric ("subsetByFeatures")})

#' @rdname GInteractions-subsetByFeatures-methods
#' @export
setMethod("subsetByFeatures", c("GInteractions", "GRanges", "missing"), function(GIObject, features, feature.class=NULL){
    i <- overlapsAny(GIObject, features)
    GIObject[i]
})

#' @rdname GInteractions-subsetByFeatures-methods
#' @export
setMethod("subsetByFeatures", c("GInteractions", "GRangesList", "missing"), function(GIObject, features, feature.class=NULL){
  i <- overlapsAny(GIObject, features)
  GIObject[i]
})

#' @rdname GInteractions-subsetByFeatures-methods
#' @export
setMethod("subsetByFeatures", c("GInteractions", "character", "character"), 
          function(GIObject, features, feature.class){
    if(!"node.class" %in% names(elementMetadata(regions(GIObject))) & 
       feature.class %in% unique(mcols(regions(GIObject))$node.class)){
      stop(paste(feature.class," has not been annotated on this object"))
    }
    #get regions which are annotated with given feature IDs
    region_idx <- which(sapply(mcols(regions(GIObject))[[paste(feature.class, "id", sep=".")]], 
                         function(x){any(features %in% x)}))
    #get object index for region idx
    gi_idx <- (anchors(GIObject, type = "first", id = TRUE) %in% region_idx) | 
              (anchors(GIObject, type = "second", id = TRUE) %in% region_idx)
    
    GIObject[gi_idx]
})
