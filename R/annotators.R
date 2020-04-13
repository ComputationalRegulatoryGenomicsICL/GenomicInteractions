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
setGeneric("resetAnnotations", function(GIObject) {
    standardGeneric("resetAnnotations")
})
#' @rdname resetAnnotations
#' @export
setMethod("resetAnnotations", c("GInteractions"), function(GIObject) {
    objName <- deparse(substitute(GIObject))
    elementMetadata(GIObject@regions) <- NULL
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
#' annotateRegions(hic_example_data, 'chip', chip)
setGeneric("annotateRegions", function(GIObject, name, dat) {
    standardGeneric("annotateRegions")
})
#' @rdname annotateRegions
#' @export
setMethod("annotateRegions", c("GInteractions", "character", "vector"), function(GIObject, name, dat) {
    objName <- deparse(substitute(GIObject))
    GIObject@regions@elementMetadata[[name]] <- dat
    assign(objName, GIObject, envir = parent.frame())
    return(invisible(1))
})

#' Annotate the interactions in a GInteractions object
#'
#' This function annotates the regions of a GInteractions object according to
#' their overlaps with `annotations`, a list of named GRanges (or GRangesList)
#' objects.
#'
#' Metadata columns will be added to the regions of the GInteractions object,
#' named according to the names of `annotations` and containing the id(s) of the
#' corresponding overlapping genomic interval(s). If `annotations` is not named,
#' the metadata columns will be named 'FEATURE#.id' where # is the position in
#' the list.
#'
#' IDs for features in each element of annotations will be extracted according
#' to the following rules: if the annotation features are a GRanges object and a
#' metadata column is specified using `id.col`, this column will be used as
#' feature IDs. Otherwise, if the annotation features are named, their names
#' will be used. If neither of these are present, they will be given numeric
#' IDs.
#'
#' For each anchor a 'node.class' metadata column will also be added, containing
#' the name of the list element which was \emph{first} annotated to each range.
#' Ranges with no overlaps will be classified as 'distal'.
#'
#' @param GIObject A GInteractions object to be annotated
#' @param annotations A list containing GRanges (or GRangesList) objects with
#'   which to annotate the GInteractions object.
#' @param id.col Metadata column of GRanges objects to use as feature ID.
#'   Default: NULL.
#' @return invisible(1)
#' @rdname annotateInteractions
#' @docType methods
#' @export
#'
#' @examples
#'
#' library('GenomicRanges')
#' data(hic_example_data)
#' data(mm9_refseq_promoters)
#' mm9_refseq_grl = split(mm9_refseq_promoters, mm9_refseq_promoters$id)
#' # This adds a `promoter.id` metadata column to regions(hic_example_data)
#' # containing IDs of overlapping promoters for each region, taken from
#' # names(mm9_refseq_grl)
#' annotateInteractions(hic_example_data, list(promoter=mm9_refseq_grl))
setGeneric("annotateInteractions", function(GIObject, annotations, id.col = NULL) {
    standardGeneric("annotateInteractions")
})
#' @rdname annotateInteractions
#' @export
#' @importFrom S4Vectors queryHits subjectHits
setMethod("annotateInteractions", c("GInteractions", "list"), function(GIObject, annotations, id.col = NULL) {
    objName <- deparse(substitute(GIObject))
    mcols.reg <- mcols(GIObject@regions)
    mcols.reg$node.class <- NA
    if (is.null(names(annotations))) {
        names(annotations) <- paste0("FEATURE", seq_along(annotations))
    }
    
    feature_names_list <- lapply(annotations, .get_gr_names, id.col = id.col)
    if (any(vapply(feature_names_list, function(x) any(duplicated(x)), logical(1)))) {
        warning("Some features contain duplicate IDs which will result in duplicate annotations")
    }
    
    for (name in names(annotations)) {
        message(paste("Annotating with", name, "..."))
        field_name <- paste(name, "id", sep = ".")
        feature_names <- feature_names_list[[name]]
        mcols.reg[[field_name]] <- NA
        reg.ol <- findOverlaps(GIObject@regions, annotations[[name]])
        
        mcols.reg[[field_name]][unique(queryHits(reg.ol))] <- split(feature_names[subjectHits(reg.ol)], queryHits(reg.ol))
        
        mcols.reg$node.class <- ifelse(is.na(mcols.reg$node.class) & !is.na(mcols.reg[[field_name]]), name, mcols.reg$node.class)
    }
    
    mcols.reg$node.class <- ifelse(is.na(mcols.reg$node.class), "distal", mcols.reg$node.class)
    mcols(GIObject@regions) <- mcols.reg
    assign(objName, GIObject, envir = parent.frame())
    return(invisible(1))
})

.get_gr_names <- function(x, id.col = NULL) {
    if (is(x, "GRanges")) {
        if (!is.null(id.col)) {
            if (!id.col %in% colnames(mcols(x))) {
                stop("id.col `", id.col, "` not found in metadata columns")
            }
            value <- mcols(x)[[id.col]]
        } else if (!is.null(names(x))) {
            value <- names(x)
        } else {
            warning("Ranges are not named and no ID column specified; using numeric IDs")
            value <- seq_along(x)
        }
    } else if (is(x, "GRangesList")) {
        if (!is.null(names(x))) {
            value <- names(x)
        } else {
            warning("GRangesList is not named; using numeric IDs")
            value <- seq_along(x)
        }
    } else {
        stop("annotations must be GRanges or GRangesList objects")
    }
    return(as.character(value))
}
