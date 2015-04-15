#' Reset annotations made to a GenomicInteractions object
#'
#' This function removes all annotations from a GenomicInteractions object by
#' deleting  all of the metadata columns associated with both anchors.
#'
#' @param GIObject An annotated GenomicInteractions object
#' @return invisible(1)
#' @docType methods
#' @rdname resetAnnotations
#' @export
setGeneric("resetAnnotations", function(GIObject){standardGeneric ("resetAnnotations")})
#' @rdname resetAnnotations
#' @export
setMethod("resetAnnotations", c("GenomicInteractions"), function(GIObject){ 
      objName = deparse(substitute(GIObject))
      elementMetadata(GIObject@anchor_one)=NULL
      elementMetadata(GIObject@anchor_two)=NULL
      assign(objName, GIObject, envir = parent.frame())
      return(invisible(1))
})

#' Annotate anchors
#'
#' This function directly annotates a single set of anchors using the GRanges
#' elementMetadata.
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
#' @import GenomicRanges
#' @export
setGeneric("annotateAnchors",function(GIObject, oneOrTwo, name, dat){standardGeneric ("annotateAnchors")})
#' @rdname annotateAnchors
#' @export
setMethod("annotateAnchors", c("GenomicInteractions", "numeric", "character", "vector"), 
            function(GIObject, oneOrTwo, name, dat){
                # need to check the validity of the arguments
                objName = deparse(substitute(GIObject))
                if(oneOrTwo == 1 ){
                    tmp.anchor = anchorOne(GIObject)
                    elementMetadata(tmp.anchor)[name] = dat
                    GIObject@anchor_one = tmp.anchor
                }else if(oneOrTwo == 2){
                    tmp.anchor = anchorTwo(GIObject)
                    elementMetadata(tmp.anchor)[name] = dat    
                    GIObject@anchor_two = tmp.anchor
                  }else{
                    stop("anchor is neither 1 or 2")
                  }
                assign(objName, GIObject, envir = parent.frame())
                return(invisible(1)) 
          })


#' Calculate interaction distances
#'
#' This function takes a GenomicInteractions object and calculates the distances
#' between the anchors according to the value of \code{method}. The distances returned
#' follow the same convention as distance(x, y) in GenomicRanges where the
#' distance between adjacent regions is 0. Note that if anchors are overlapping
#' this method will print a warning and return the distance as 0.
#' 
#' @param GIObject A GenomicInteractions object
#' @param method Character vector indicating how to calculate distances, must
#'        be one of `midpoint', `outer', `inner'.
#' @param floor A logical specifying whether to round down distances to nearest 
#' base pair or not. Default TRUE.
#'
#' @return An vector containing the distances between anchors/GRanges,
#'         NA if on different chromosomes, rounded down to the nearest bp.
#' @import GenomicRanges 
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
setGeneric("calculateDistances", function(GIObject, method="midpoint", floor=TRUE){standardGeneric ("calculateDistances")})

#' @rdname calculateDistances
#' @export
setMethod("calculateDistances", c("GenomicInteractions"), 
        function(GIObject, method="midpoint", floor=TRUE){ 
            if(method=="midpoint"){
                midpt.one = (start(GIObject@anchor_one) + end(GIObject@anchor_one)) / 2
                midpt.two = (start(GIObject@anchor_two) + end(GIObject@anchor_two)) / 2
                distances = ifelse(as.character(seqnames(GIObject@anchor_one))==as.character(seqnames(GIObject@anchor_two)),
                                   ifelse(midpt.two > midpt.one,
                                          midpt.two - midpt.one - 1, 
                                          midpt.one - midpt.two - 1), NA)
            }else if(method=="outer"){
                distances = ifelse(as.character(seqnames(GIObject@anchor_one))==as.character(seqnames(GIObject@anchor_two)),
                                ifelse(start(GIObject@anchor_two) > start(GIObject@anchor_one), 
                                    end(GIObject@anchor_two) - start(GIObject@anchor_one) -1 ,
                                    end(GIObject@anchor_one) - start(GIObject@anchor_two) -1 ), NA)
            }else if(method=="inner"){
                distances = ifelse(as.character(seqnames(GIObject@anchor_one))==as.character(seqnames(GIObject@anchor_two)),
                                ifelse(start(GIObject@anchor_two) > start(GIObject@anchor_one), 
                                    start(GIObject@anchor_two) - end(GIObject@anchor_one) - 1 ,
                                    start(GIObject@anchor_one) - end(GIObject@anchor_two) - 1 ), NA)
            }else{
                stop( "method must be one of c(\"midpoint\", \"outer\", \"inner\")" )
            }
            
            if(sum(distances[!is.na(distances)] < 0) > 0){
                warning("setting negative distances to 0, this is due to the presence of overlapping anchors in your dataset")
                distances = ifelse(distances < 0, 0, distances)
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
                                    midpt.two - midpt.one - 1, 
                                    midpt.one - midpt.two - 1), NA)
            }else if(method=="outer"){
                distances = ifelse(object1$seqnames==object2$seqnames,
                                ifelse(object2$start > object1$start, 
                                    object2$end - object1$start -1 ,
                                    object1$end) - object2$start -1, NA)
            }else if(method=="inner"){ 
                distances = ifelse(object1$seqnames==object2$seqnames,
                                ifelse(object2$start > object1$start, 
                                        object2$start - object1$end - 1 ,
                                        object1$start - object2$end - 1 ), NA)
            }else{
                stop( "method must be one of c(\"midpoint\", \"outer\", \"inner\")" )
            }
            if(sum(distances[!is.na(distances)] < 0) > 0){
                #warning("setting negative distances to 0, this is due to the presence of overlapping anchors in your dataset")
                distances = ifelse(distances < 0, 0, distances)
            }
            if(floor){
                return(floor(distances))
            }else{
                return(distances)
            }
})


#' Annotate the interactions in a GenomicInteractions object
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
#' @param GIObject A GenomicInteractions object to be annotated
#' @param annotations A list containing GRanges (or GRangesList) objects with which to annotate
#'             the GenomicInteractions object.
#' @return invisible(1)
#' @rdname annotateInteractions
#' @docType methods
#' @import GenomicRanges
#' @export
#' 
#' @examples
#' 
#' data(hic_example_data)
#' data(mm9_refseq_promoters)
#' \dontrun{
#' annotateInteractions(hic_example_data, list(promoter=mm9_refseq_promoters))
#' }
setGeneric("annotateInteractions",function(GIObject, annotations){standardGeneric ("annotateInteractions")})
#' @rdname annotateInteractions
#' @export
setMethod("annotateInteractions", c("GenomicInteractions", "list"), 
            function(GIObject, annotations){
                objName = deparse(substitute(GIObject))
                GIObject@anchor_one$node.class = NA
                GIObject@anchor_two$node.class = NA
                
                if (is.null(names(annotations))){
                  names(annotations) <- paste0("FEATURE", 1:length(annotations))
                }

                for(name in names(annotations)){
                    message(paste("Annotating with", name, "..."))
                    elementMetadata(GIObject@anchor_one)[[paste(name, "id", sep=".")]] = NA
                    elementMetadata(GIObject@anchor_two)[[paste(name, "id", sep=".")]] = NA
                    one.ol = findOverlaps( GIObject@anchor_one, annotations[[name]])
                    two.ol = findOverlaps( GIObject@anchor_two, annotations[[name]])
                    
                    if(class(annotations[[name]])=="GRanges"){ 
                    # GRanges for annotation must have an id column or have names set
                        if(!is.null(names(annotations[[name]]))){
                            elementMetadata( GIObject@anchor_one)[[ paste(name, "id", sep=".")]][ unique(queryHits(one.ol)) ] = 
                                lapply(unique(queryHits(one.ol)), function(x){ unique(names(annotations[[name]])[unique(subjectHits(one.ol[ queryHits(one.ol) == x]))])})
                            elementMetadata( GIObject@anchor_two)[[ paste(name, "id", sep=".")]][ unique(queryHits(two.ol)) ] = 
                                lapply(unique(queryHits(two.ol)), function(x){ unique(names(annotations[[name]])[unique(subjectHits(two.ol[ queryHits(two.ol) == x]))])})
                        }else if("id" %in% names(mcols(annotations[[name]]))){
                            elementMetadata( GIObject@anchor_one)[[ paste(name, "id", sep=".")]][ unique(queryHits(one.ol)) ] = 
                                lapply(unique(queryHits(one.ol)), function(x){ unique(annotations[[name]]$id[unique(subjectHits(one.ol[ queryHits(one.ol) == x]))])})
                            elementMetadata( GIObject@anchor_two)[[ paste(name, "id", sep=".")]][ unique(queryHits(two.ol)) ] = 
                                lapply(unique(queryHits(two.ol)), function(x){ unique(annotations[[name]]$id[unique(subjectHits(two.ol[ queryHits(two.ol) == x]))])})
                        }else{
                            stop("annotations requires an id column in elementMetadata or names to be non-null")
                        }
                    }else if(class(annotations[[name]])=="GRangesList"){
                      # GRangesList names must be ids for the regions
                      elementMetadata( GIObject@anchor_one)[[ paste(name, "id", sep=".")]][ unique(queryHits(one.ol)) ] = 
                          lapply(unique(queryHits(one.ol)), function(x){ names(annotations[[name]])[unique(subjectHits(one.ol[ queryHits(one.ol) == x]))]})
                      elementMetadata( GIObject@anchor_two)[[ paste(name, "id", sep=".")]][ unique(queryHits(two.ol)) ] = 
                          lapply(unique(queryHits(two.ol)), function(x){ names(annotations[[name]])[unique(subjectHits(two.ol[ queryHits(two.ol) == x]))]})
                    }  
                    GIObject@anchor_one$node.class = ifelse(is.na(GIObject@anchor_one$node.class) & !is.na(elementMetadata( GIObject@anchor_one)[[ paste(name, "id", sep=".")]]), 
                                                            name, 
                                                            GIObject@anchor_one$node.class)
                    GIObject@anchor_two$node.class = ifelse(is.na(GIObject@anchor_two$node.class) & !is.na(elementMetadata( GIObject@anchor_two)[[ paste(name, "id", sep=".")]]), 
                                                            name, 
                                                            GIObject@anchor_two$node.class)  
                }
                GIObject@anchor_one$node.class = ifelse(is.na(GIObject@anchor_one$node.class), "distal", GIObject@anchor_one$node.class)
                GIObject@anchor_two$node.class = ifelse(is.na(GIObject@anchor_two$node.class), "distal", GIObject@anchor_two$node.class)
                assign(objName, GIObject, envir = parent.frame())
                return(invisible(1))
})

#' Summary statistics of interactions for a given feature set
#'
#' This function will calculate summary statistics for each element in the
#' given feature set, including the number of interactions (the sum of all
#' interaction counts), number of unique interactions and number of trans-
#' (interchromosomal) interations.  It also returns some statistics for the 
#' distances of interactions for all interactions of the feature, and for the 
#' different interaction types e.g. promoter-distal.
#'
#' @param GIObject An annotated GenomicInteractions object
#' @param features A GRanges object containing the feature set
#' @param feature.name The name of the feature set
#' @param distance.method Method for calculating distances between anchors, see
#'                        ?calculateDistances
#' @param annotate.self Logical. Indicates whether to annotate self interactions,
#' i.e. where a feature in `features` overlaps both anchors of an interaction.
#'  Default: FALSE.
#'
#' @return A data frame with one line for each range in `features'
#' @rdname summariseByFeatures
#' @docType methods
#' @import GenomicRanges
#' @export
setGeneric("summariseByFeatures",function(GIObject, features, feature.name, distance.method="midpoint", annotate.self=FALSE){standardGeneric ("summariseByFeatures")})
#' @rdname summariseByFeatures
#' @export
setMethod("summariseByFeatures", "GenomicInteractions", 
          function(GIObject, features, feature.name, distance.method="midpoint", annotate.self=FALSE){  
                if(!("node.class" %in% names(elementMetadata(GIObject@anchor_one))) | !("node.class" %in% names(elementMetadata(GIObject@anchor_two)))){stop("GIObject has not been annotated")}
                    potential.node.classes = unique(c(GIObject@anchor_one$node.class, GIObject@anchor_two$node.class))
                    summary.df = data.frame(matrix(0, ncol=(5 + (length(potential.node.classes) * 2) + ifelse(annotate.self, (length(potential.node.classes)-1) * 2, 0) + 5), nrow=length(features)))
                    summary.names = c(paste(capitalize(feature.name), "id", sep="."), 
                    paste("numberOf", capitalize(feature.name), "Interactions", sep=""),
                    paste("numberOf", capitalize(feature.name), "UniqueInteractions", sep=""),
                    paste("numberOf", capitalize(feature.name), "InterChromosomalInteractions", sep=""),
                    paste("numberOf", capitalize(feature.name), "UniqueInterChromosomalInteractions", sep=""),
                    paste("numberOf", capitalize(feature.name), capitalize(potential.node.classes), "Interactions", sep=""),
                    paste("numberOfUnique", capitalize(feature.name), capitalize(potential.node.classes), "Interactions", sep=""),
                    paste(capitalize(feature.name), "DistanceMedian", sep=""),
                    paste(capitalize(feature.name), "DistanceMean", sep=""),
                    paste(capitalize(feature.name), "DistanceMinimum", sep=""),
                    paste(capitalize(feature.name), "DistanceMaximum", sep=""),
                    paste(capitalize(feature.name), "DistanceWeightedMedian", sep=""))
                if(annotate.self){    
                    pc= potential.node.classes[ -which(potential.node.classes == "distal")]
                    summary.names = append(summary.names, c( paste("numberOfSelf", capitalize(feature.name), capitalize(pc), "Interactions", sep=""),
                                                             paste("numberOfSelfUnique", capitalize(feature.name), capitalize(pc), "Interactions", sep="")
                                                            ))
                }
                names(summary.df) = summary.names

                if(class(features) == "GRangesList"){
                    row.names(summary.df) = names(features)
                    summary.df[,paste(capitalize(feature.name), "id", sep=".")] = names(features)
                }else if(class(features) == "GRanges"){
                    if(!is.null(names(features))){
                        row.names(summary.df) = names(features)
                        summary.df[,paste(capitalize(feature.name), "id", sep=".")] = names(features)
                    }else if(is.null(names(features)) & "id" %in% names(elementMetadata(features))){
                        row.names(summary.df) = features$id
                        summary.df[,paste(capitalize(feature.name), "id", sep=".")] = features$id
                    }else{
                        stop("features missing names or an id column")    
                    }
                }else{
                    stop("features must be GRanges object")
                }
                                                                              
                one.ol = findOverlaps(features, GIObject@anchor_one)
                two.ol = findOverlaps(features, GIObject@anchor_two)
                one.indexes = queryHits(one.ol)[GIObject@anchor_one[ subjectHits(one.ol) ]$node.class == feature.name]  
                two.indexes = queryHits(two.ol)[GIObject@anchor_two[ subjectHits(two.ol) ]$node.class == feature.name] 
                features.with.interactions.indexes = unique(c(one.indexes, two.indexes))
                features.with.interactions.indexes = features.with.interactions.indexes[order(features.with.interactions.indexes)]
                anchor_one.df = data.frame(seqnames=as.character(seqnames(anchorOne(GIObject))),
                                                            start=start(anchorOne(GIObject)),
                                                             end=end(anchorOne(GIObject)),
                                                             width=width(anchorOne(GIObject)),
                                                             strand=as.character(strand(anchorOne(GIObject))),stringsAsFactors=FALSE)
                anchor_two.df = data.frame(seqnames=as.character(seqnames(anchorTwo(GIObject))),
                                                         start=start(anchorTwo(GIObject)),
                                                         end=end(anchorTwo(GIObject)),
                                                         width=width(anchorTwo(GIObject)),
                                                         strand=as.character(strand(anchorTwo(GIObject))),stringsAsFactors=FALSE)
                for(i in features.with.interactions.indexes){
                    interactions = unique(c(subjectHits(one.ol[ queryHits(one.ol)==i]), subjectHits(two.ol[ queryHits(two.ol)==i])))
                    interactions.one = unique(subjectHits(one.ol[ queryHits(one.ol)==i]))
                    interactions.two = unique(subjectHits(two.ol[ queryHits(two.ol)==i]))
                
                    numberOfInteractions = sum(GIObject@counts[ interactions ])
                    numberOfUniqueInteractions = length(interactions)
                    summary.df[i,paste("numberOf", capitalize(feature.name), "Interactions", sep="")] = numberOfInteractions
                    summary.df[i,paste("numberOf", capitalize(feature.name), "UniqueInteractions", sep="")] = numberOfUniqueInteractions
                    
                    intercis = interactions[as.character(seqnames(GIObject@anchor_one[interactions])) != as.character(seqnames( GIObject@anchor_two[interactions]))]
                    if(length(intercis)>0){
                      numberOfInterChromosomalInteractions = sum(GIObject@counts[intercis])
                      summary.df[i,paste("numberOf", capitalize(feature.name), "InterChromosomalInteractions", sep="")] = numberOfInterChromosomalInteractions
                      numberOfUniqueInterChromosomalInteractions = length(intercis)
                      summary.df[i,paste("numberOf", capitalize(feature.name), "UniqueInterChromosomalInteractions", sep="")] = numberOfUniqueInterChromosomalInteractions
                    }
                    
                    for(nc in potential.node.classes){
                      nc1Counts = 0
                      nc1Indexes = c()
                      nc2Counts = 0
                      nc2Indexes = c()  
                      if(length(interactions.one)>0){
                        nc1Indexes = interactions.one[which( GIObject@anchor_one$node.class[interactions.one] == feature.name &  GIObject@anchor_two$node.class[interactions.one]  == nc)]
                        nc1Counts = sum(GIObject@counts[nc1Indexes])
                      }
                      if(length(interactions.two)>0){
                        nc2Indexes = interactions.two[which(GIObject@anchor_one$node.class[interactions.two] == nc & GIObject@anchor_two$node.class[interactions.two] == feature.name)]
                        if(nc == feature.name ){ 
                            if(length(nc2Indexes[ (nc2Indexes %in% nc1Indexes) ])>0){
                                summary.df[i,paste("numberOfSelf", capitalize(feature.name), capitalize(nc), "Interactions", sep="")] = sum(GIObject@counts[nc2Indexes[(nc2Indexes %in% nc1Indexes)]])
                                summary.df[i,paste("numberOfSelfUnique", capitalize(feature.name), capitalize(nc), "Interactions", sep="")] = length(nc2Indexes[(nc2Indexes %in% nc1Indexes)])
                            }
                            nc2Indexes = nc2Indexes[!(nc2Indexes %in% nc1Indexes)] # stop double counting of interactions where both anchors are in the same feature
                        } 
                        nc2Counts = sum(GIObject@counts[nc2Indexes])
                      }
                      
                      if((nc1Counts > 0) | (nc2Counts > 0)){
                        numberOfNCInteractions = nc1Counts + nc2Counts
                        summary.df[i,paste("numberOf", capitalize(feature.name), capitalize(nc), "Interactions", sep="")] = numberOfNCInteractions
                        numberOfUniqueNCInteractions = length(c(nc1Indexes, nc2Indexes))
                        summary.df[i,paste("numberOfUnique", capitalize(feature.name), capitalize(nc), "Interactions", sep="")] = numberOfUniqueNCInteractions
                      }
                    }
                    # annotate.self .. ie id is the same id at both ends for different classes
                    if(annotate.self){
                        for(nc in potential.node.classes){
                            if(nc != feature.name & nc != "distal"){
                                ncs1Counts = 0 
                                ncs1Indexes = c()
                                ncs2Counts = 0
                                ncs2Indexes = c()
                                if(length(interactions.one)>0){
                                    ncs1Indexes = interactions.one[which(GIObject@anchor_one$node.class[interactions.one] == feature.name &  GIObject@anchor_two$node.class[interactions.one]  == nc)] 
                                    ncs1Indexes = names(features)[i] %in% elementMetadata(GIObject@anchor_two)[[ paste(nc, "id", sep=".") ]][ncs1Indexes]
                                    ncs1Counts = sum(GIObject@counts[ncs1Indexes])
                                }
                                if(length(interactions.two)>0){
                                    nc2sIndexes = interactions.two[which(GIObject@anchor_one$node.class[interactions.two] == nc & GIObject@anchor_two$node.class[interactions.two] == feature.name)]
                                    ncs2Indexes = names(features)[i]  %in%  elementMetadata(GIObject@anchor_one)[[ paste(nc, "id", sep=".") ]][ncs2Indexes]
                                    ncs2Counts = sum(GIObject@counts[ncs2Indexes])
                                }
                                if((ncs1Counts > 0) | (ncs2Counts > 0)){
                                    numberOfNCSInteractions = ncs1Counts + ncs2Counts
                                    summary.df[i,paste("numberOfSelf", capitalize(feature.name), capitalize(nc), "Interactions", sep="")] = numberOfNCInteractions
                                    numberOfUniqueNCSInteractions = length(c(ncs1Indexes, ncs2Indexes))
                                    summary.df[i,paste("numberOfSelfUnique", capitalize(feature.name), capitalize(nc), "Interactions", sep="")] = numberOfUniqueNCInteractions
                                }
                            }
                        }
                    }
                    distances = .calculateDistances.df(anchor_one.df[interactions,], anchor_two.df[interactions,], distance.method)
                    dis = distances[!is.na(distances)]
                    wdis = rep(distances, GIObject@counts[interactions])
                    wdis = wdis[!is.na(wdis)]
                    
                    median.distance = median(dis)
                    median.distance = ifelse(is.infinite(median.distance) | is.nan(median.distance), NA, median.distance)
                    mean.distance = median(dis)
                    mean.distance = ifelse(is.infinite(mean.distance) | is.nan(mean.distance), NA, mean.distance)
                    min.distance = min(dis)
                    min.distance = ifelse(is.infinite(min.distance) | is.nan(min.distance), NA, min.distance)
                    max.distance = max(dis)
                    max.distance = ifelse(is.infinite(max.distance) | is.nan(max.distance), NA, max.distance)
                    wmedian.distance = median(wdis)
                    wmedian.distance = ifelse(is.infinite(wmedian.distance) | is.nan(wmedian.distance), NA, wmedian.distance)
                    
                    summary.df[i,paste(capitalize(feature.name), "DistanceMedian", sep="")] = median.distance
                    summary.df[i,paste(capitalize(feature.name), "DistanceMean", sep="")] = mean.distance
                    summary.df[i,paste(capitalize(feature.name), "DistanceMinimum", sep="")] = min.distance
                    summary.df[i,paste(capitalize(feature.name), "DistanceMaximum", sep="")] = max.distance
                    summary.df[i,paste(capitalize(feature.name), "DistanceWeightedMedian", sep="")] = wmedian.distance
                }
              return(summary.df)                                                                                                                                       
})
