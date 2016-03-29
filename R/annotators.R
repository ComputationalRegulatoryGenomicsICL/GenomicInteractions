#' Reset annotations made to a GInteractions object
#'
#' This function removes all annotations from a GInteractions object by
#' deleting  all of the metadata columns associated with both anchors.
#'
#' @param GIObject An annotated GInteractions object
#' @return invisible(1)
#' @docType methods
#' @rdname resetAnnotations
#' @export
#' @examples
#' data(hic_example_data)
#' resetAnnotations(hic_example_data)
setGeneric("resetAnnotations", function(GIObject){standardGeneric ("resetAnnotations")})
#' @rdname resetAnnotations
#' @export
setMethod("resetAnnotations", c("GInteractions"), function(GIObject){ 
      objName = deparse(substitute(GIObject))
      elementMetadata(GIObject@regions)=NULL
      assign(objName, GIObject, envir = parent.frame())
      return(invisible(1))
})

#' Annotate anchors - DEPRECATED
#'
#' This function is deprecated and will be removed in the next release of Bioconductor.
#' Use `annotateRegions` instead.
#'
#' @param GIObject A GenomicInteractions object
#' @param oneOrTwo An integer indicating which anchor to annotate
#' @param name Character. Will be used as a column name for the elementMetadata 
#' of the annotated anchor. 
#' @param dat Vector of the same length as the GenomicInteractions object,
#' containing data with which to annotate the object. 
#' @return invisible(1)
#' 
#' @rdname annotateAnchors
#' @docType methods
#' @export
#' 
setGeneric("annotateAnchors",function(GIObject, oneOrTwo, name, dat){
  standardGeneric ("annotateAnchors")})
#' @rdname annotateAnchors
#' @export
setMethod("annotateAnchors", c("GenomicInteractions", "numeric", "character", "vector"), 
            function(GIObject, oneOrTwo, name, dat){
              .Deprecated("annotateRegions")
              
                # need to check the validity of the arguments
                objName = deparse(substitute(GIObject))
                elementMetadata(GIObject@regions)[[name]] <- NA
                if(oneOrTwo == 1 ){
                  if(length(dat) > length(GIObject@anchor1)){
                    warning("Annotation data longer than anchors! Try `annotateRegions` instead.")
                  }          
                  elementMetadata(GIObject@regions)[[name]][GIObject@anchor1] <- dat
                }else if(oneOrTwo == 2){
                  if(length(dat) > length(GIObject@anchor2)){
                    warning("Annotation data longer than anchors! Try `annotateRegions` instead.")
                  } 
                  elementMetadata(GIObject@regions)[[name]][GIObject@anchor2] <- dat
                  }else{
                    stop("anchor is neither 1 or 2")
                  }
                assign(objName, GIObject, envir = parent.frame())
                return(invisible(1)) 
          })

#' Annotate regions
#'
#' Use this function to add metadata parallel to the `regions` slot of a 
#' GenomicInteractions or GInteractions object.
#' 
#' @param GIObject A GenomicInteractions or GInteractions object
#' @param name Character. Will be used as a column name.
#' @param dat Vector of the same length as the GInteractions object,
#' containing data with which to annotate the object. 
#' @return invisible(1)
#' 
#' @rdname annotateRegions
#' @docType methods
#' @export
#' 
#' @examples 
#' data(hic_example_data)
#' chip <- runif(n = length(regions(hic_example_data)), max = 1000)
#' annotateRegions(hic_example_data, "chip", chip)
setGeneric("annotateRegions", function(GIObject, name, dat){standardGeneric ("annotateRegions")})
#' @rdname annotateRegions
#' @export
setMethod("annotateRegions", c("GInteractions", "character", "vector"),
  function(GIObject, name, dat){
    objName = deparse(substitute(GIObject))
    GIObject@regions@elementMetadata[[name]] <- dat
    assign(objName, GIObject, envir = parent.frame())
    return(invisible(1))
})

#' Annotate the interactions in a GInteractions object
#'
#' This function will annotate both anchors with a list of named GRanges
#' objects. Each metadata column is labeled "name.id" and contains the id of
#' the genomic interval(s) it overlaps. Anonymous lists will be given names
#' "FEATURE#.id" where # is the position in the list.
#'
#' For each anchor a "node.class" metadata column will also be added, containing
#' the name of the list element which was \emph{first} annotated to each range.
#' Ranges with no overlaps will be classified as "distal". The identifiers for each 
#' individual feature/annotation are taken from either the name of the list item in the 
#' case of a GRangesList or from either the names of a the provided GRanges or an id column 
#' in its associated metadata.
#'
#' @param GIObject A GInteractions object to be annotated
#' @param annotations A list containing GRanges (or GRangesList) objects with which to annotate
#'             the GInteractions object.
#' @return invisible(1)
#' @rdname annotateInteractions
#' @docType methods
#' @export
#' 
#' @examples
#' 
#' library("GenomicRanges")
#' data(hic_example_data)
#' data(mm9_refseq_promoters)
#' mm9_refseq_grl = split(mm9_refseq_promoters, mm9_refseq_promoters$id)
#' annotateInteractions(hic_example_data, list(promoter=mm9_refseq_grl))
setGeneric("annotateInteractions",function(GIObject, annotations){standardGeneric ("annotateInteractions")})
#' @rdname annotateInteractions
#' @export
#' @importFrom S4Vectors queryHits subjectHits
setMethod("annotateInteractions", c("GInteractions", "list"), 
            function(GIObject, annotations){
                objName = deparse(substitute(GIObject))
                mcols.reg <- mcols(GIObject@regions)
                mcols.reg$node.class <- NA
                if (is.null(names(annotations))){
                  names(annotations) <- paste0("FEATURE", 1:length(annotations))
                }

                feature_names_list = lapply(annotations, .get_gr_names)
                if (any(vapply(feature_names_list, function(x) any(duplicated(x)), logical(1)))) {
                    warning("Some features contain duplicate IDs which will result in duplicate annotations")
                }

                for(name in names(annotations)){
                    message(paste("Annotating with", name, "..."))
                    field_name = paste(name, "id", sep=".")
                    feature_names = feature_names_list[[name]]
                    mcols.reg[[field_name]] = NA
                    reg.ol = findOverlaps(GIObject@regions, annotations[[name]])
                    mcols.reg[[field_name]][ unique(queryHits(reg.ol)) ] <- split(feature_names[subjectHits(reg.ol)], 
                                                                                  queryHits(reg.ol) )
                    mcols.reg$node.class = ifelse(is.na(mcols.reg$node.class) & !is.na(mcols.reg[[field_name]]), 
                                                  name, mcols.reg$node.class)
                }

                mcols.reg$node.class = ifelse(is.na(mcols.reg$node.class), 
                                              "distal", mcols.reg$node.class)
                mcols(GIObject@regions) = mcols.reg
                assign(objName, GIObject, envir = parent.frame())
                return(invisible(1))
})

.get_gr_names = function(x) {
    if (class(x)=="GRanges") { 
        if (!is.null(names(x))) {
            value = names(x)
        } else if("id" %in% names(mcols(x))) {
            value = x$id
        } else {
            stop("annotations requires an id column in elementMetadata or names to be non-null")
        }
    } else if(class(x)=="GRangesList") {
        value = names(x)
    } else {
        stop("annotations must be GRanges or GRangesList objects")
    }
    as.character(value)
}
