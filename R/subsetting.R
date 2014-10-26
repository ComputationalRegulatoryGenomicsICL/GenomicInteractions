#' Subset a GenomicInteractions object by features
#'
#' Subsets interactions for which at least one of the anchors overlaps with a given GRanges object.
#' Alternatively, subsets interactions based on annotated feature IDs for a particular feature.
#'
#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @docType methods
#' @param GIObject A GenomicInteractions object
#' @param features A GRanges or GRangesList object, or a character vector containing
#'                      IDs of annotated features, e.g. promoter IDs.
#' @param feature.class If `features' is a character vector, the corresponding feature name, e.g. "promoter".
#' @return a subsetted GenomicInteractions object
#' @export
setGeneric("subsetByFeatures",function(GIObject, features, feature.class=NULL){standardGeneric ("subsetByFeatures")})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @import GenomicRanges
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "GRanges", "missing"), function(GIObject, features, feature.class=NULL){
    indic = unique(c(subjectHits(findOverlaps(features, GIObject@anchor_one)), subjectHits(findOverlaps(features, GIObject@anchor_two))))
    return( new("GenomicInteractions", 
        experiment_name = GIObject@experiment_name, 
        description = GIObject@description, 
        genome_name = GIObject@genome_name, 
        anchor_one=GIObject@anchor_one[indic], 
        anchor_two=GIObject@anchor_two[indic], 
        counts=GIObject@counts[indic], 
        pvalue=GIObject@pvalue[indic], 
        fdr=GIObject@fdr[indic]) )
})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @import GenomicRanges
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "GRangesList", "missing"), function(GIObject, features, feature.class=NULL){
    indic = unique(c(subjectHits(findOverlaps(features, GIObject@anchor_one)), subjectHits(findOverlaps(features, GIObject@anchor_two))))
    return( new("GenomicInteractions", 
        experiment_name = GIObject@experiment_name, 
        description = GIObject@description, 
        genome_name = GIObject@genome_name, 
        anchor_one=GIObject@anchor_one[indic], 
        anchor_two=GIObject@anchor_two[indic], 
        counts=GIObject@counts[indic], 
        pvalue=GIObject@pvalue[indic], 
        fdr=GIObject@fdr[indic]) )
})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @import GenomicRanges
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "character", "character"), function(GIObject, features, feature.class){
    if("node.class" %in% names(elementMetadata(GIObject@anchor_one)) & feature.class %in% unique(c(GIObject@anchor_one$node.class, GIObject@anchor_two$node.class))){
        indic = sapply(elementMetadata(GIObject@anchor_one)[[paste(feature.class, "id", sep=".")]], 
                    function(x){ features %in% x }) | sapply(elementMetadata(GIObject@anchor_two)[[paste(feature.class, "id", sep=".")]], function(x){ features %in% x })
        return( new("GenomicInteractions", 
            experiment_name = GIObject@experiment_name, 
            description = GIObject@description, 
            genome_name = GIObject@genome_name, 
            anchor_one=GIObject@anchor_one[indic], 
            anchor_two=GIObject@anchor_two[indic], 
            counts=GIObject@counts[indic], 
            pvalue=GIObject@pvalue[indic], 
            fdr=GIObject@fdr[indic]) )
    }else{
        stop(paste(feature.class," has not been annotated on this GenomicInteractions object"))
    }                                                                         
  
})

#' Standard subsetting methods for GenomicInteractions objects
#'
#' @name [
#' @param x A genomicInteractions object
#' @param i A numeric, logical or Rle vector
#' 
#' @return A GenomicInteractions object containing only the features specified by `i`.
#' @rdname GenomicInteractions-subset-methods
NULL

#' @name [
#' @aliases [,GenomicInteractions,logical,missing-method
#' @docType methods
#' @rdname GenomicInteractions-subset-methods
#' @export
setMethod(f="[", signature=c("GenomicInteractions", "logical", "missing"), definition=function(x,i){
  return( new("GenomicInteractions", 
              experiment_name = x@experiment_name, 
              description = x@description, 
              genome_name = x@genome_name, 
              anchor_one=x@anchor_one[i], 
              anchor_two=x@anchor_two[i], 
              counts=x@counts[i], 
              pvalue=x@pvalue[i], 
              fdr=x@fdr[i]) )
} )

#' @name [
#' @aliases [,GenomicInteractions,numeric,missing-method
#' @docType methods
#' @rdname GenomicInteractions-subset-methods
#' @export
setMethod(f="[", signature=c("GenomicInteractions", "numeric", "missing"), definition=function(x,i){
  return( new("GenomicInteractions", 
              experiment_name = x@experiment_name, 
              description = x@description, 
              genome_name = x@genome_name, 
              anchor_one=x@anchor_one[i], 
              anchor_two=x@anchor_two[i], 
              counts=x@counts[i], 
              pvalue=x@pvalue[i], 
              fdr=x@fdr[i]) )
} )

#' @name [
#' @aliases [,GenomicInteractions,Rle,missing-method
#' @docType methods
#' @rdname GenomicInteractions-subset-methods
#' @export
setMethod(f="[", signature=c("GenomicInteractions", "Rle", "missing"), definition=function(x,i){
  i = as.vector(i) # bit of a hack but unsure how to access Rle directly
  return( new("GenomicInteractions", 
              experiment_name = x@experiment_name, 
              description = x@description, 
              genome_name = x@genome_name, 
              anchor_one=x@anchor_one[i], 
              anchor_two=x@anchor_two[i], 
              counts=x@counts[i], 
              pvalue=x@pvalue[i], 
              fdr=x@fdr[i]) )
} )

#' @name [
#' @aliases [,GenomicInteractions,rle,missing-method
#' @docType methods
#' @rdname GenomicInteractions-subset-methods
#' @export
setMethod(f="[", signature=c("GenomicInteractions", "rle", "missing"), definition=function(x,i){
  i = as.vector(i) # bit of a hack but unsure how to access Rle directly
  return( new("GenomicInteractions", 
              experiment_name = x@experiment_name, 
              description = x@description, 
              genome_name = x@genome_name, 
              anchor_one=x@anchor_one[i], 
              anchor_two=x@anchor_two[i], 
              counts=x@counts[i], 
              pvalue=x@pvalue[i], 
              fdr=x@fdr[i]) )
} )


