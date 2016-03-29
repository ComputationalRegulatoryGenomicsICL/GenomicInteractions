
#' Calculate interaction distances
#'
#' This function takes a GInteractions object and calculates the distances
#' between the anchors according to the value of \code{method}. The distances returned
#' follow the same convention as distance(x, y) in GenomicRanges where the
#' distance between adjacent regions is 0. Note that if anchors are overlapping
#' this method will print a warning and return the distance as 0.
#' 
#' @param GIObject A GInteractions object
#' @param method Character vector indicating how to calculate distances, must
#'        be one of `midpoint', `outer', `inner'.
#' @param floor A logical specifying whether to round down distances to nearest 
#' base pair or not. Default TRUE.
#'
#' @return An vector containing the distances between anchors/GRanges,
#'         NA if on different chromosomes, rounded down to the nearest bp.
#' @examples
#'
#' library(GenomicRanges)
#' 
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, experiment_name="test", 
#'                            description="this is a test", counts=interaction_counts)
#' calculateDistances(test, method="midpoint")
#'         
#' @docType methods
#' @rdname calculateDistances
#' @export
setGeneric("calculateDistances", function(GIObject, method="midpoint", floor=TRUE){
  standardGeneric ("calculateDistances")})

#' @rdname calculateDistances
#' @export
setMethod("calculateDistances", c("GInteractions"), 
          function(GIObject, method="midpoint", floor=TRUE){ 
            if(method=="midpoint"){
              distances <- pairdist(GIObject, type = "mid")
            }else if(method=="outer"){
              distances <- pairdist(GIObject, type = "span")
            }else if(method=="inner"){
              distances <- pairdist(GIObject, type = "gap")
            }else{
              distances <- pairdist(GIObject, type = method)
            }
            
            if(floor){
              return(floor(distances))
            }else{
              return(distances)
            }
          })


setGeneric(".calculateDistances.df", function(object1, object2, method="midpoint", floor=TRUE){standardGeneric (".calculateDistances.df")})

setMethod(".calculateDistances.df", c("data.frame", "data.frame"), 
          function(object1, object2, method="midpoint", floor=TRUE){ 
            if(method=="midpoint"){
              midpt.one = (object1$start + object1$end) / 2
              midpt.two = (object2$start + object2$end) / 2
              distances = ifelse(object1$seqnames==object2$seqnames,
                                 ifelse(midpt.two > midpt.one,
                                        midpt.two - midpt.one, 
                                        midpt.one - midpt.two), NA)
            }else if(method=="outer"){
              distances = ifelse(object1$seqnames==object2$seqnames,
                                 ifelse(object2$start > object1$start, 
                                        object2$end - object1$start + 1,
                                        object1$end - object2$start + 1), NA)
            }else if(method=="inner"){ 
              distances = ifelse(object1$seqnames==object2$seqnames,
                                 ifelse(object2$start > object1$start, 
                                        object2$start - object1$end - 1 ,
                                        object1$start - object2$end - 1 ), NA)
            }else{
              stop( "method must be one of c(\"midpoint\", \"outer\", \"inner\")" )
            }
            if(floor){
              return(floor(distances))
            }else{
              return(distances)
            }
          })

