#' Summary statistics of interactions for a given feature set
#'
#' This function will calculate summary statistics for each element in the
#' given feature set, including the number of interactions (the sum of all
#' interaction counts), number of unique interactions and number of trans-
#' (interchromosomal) interations.  It also returns some statistics for the 
#' distances of interactions for all interactions of the feature, and for the 
#' different interaction types e.g. promoter-distal.
#'
#' @param GIObject An annotated GInteractions object
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
#' @importFrom stats median
#' @export
#' @examples 
#' data("hic_example_data")
#' data("mm9_refseq_promoters")
#' annotateInteractions(hic_example_data, list(promoter = mm9_refseq_promoters))
#' summariseByFeatures(hic_example_data, mm9_refseq_promoters[1:10], "promoter")
setGeneric("summariseByFeatures",function(GIObject, features, feature.name, distance.method="midpoint", annotate.self=FALSE){standardGeneric ("summariseByFeatures")})
#' @rdname summariseByFeatures
#' @export
setMethod("summariseByFeatures", "GInteractions", 
          function(GIObject, features, feature.name, distance.method="midpoint", annotate.self=FALSE){  
            if(!("node.class" %in% names(elementMetadata(regions(GIObject))))){
              stop("GIObject has not been annotated")}
              
              potential.node.classes = unique(GIObject@regions$node.class)
              
              feature.names.full <- .get_gr_names(features)
              feature.names <- unique(feature.names.full)
              
              summary.df = data.frame(matrix(0, 
                                             ncol=(5 + (length(potential.node.classes) * 2) + ifelse(annotate.self, (length(potential.node.classes)-1) * 2, 0) + 5), 
                                             nrow=length(feature.names)))
              
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
              
              summary.df[,paste(capitalize(feature.name), "id", sep=".")] = feature.names
              
              #hack
              anchor_one <- anchorOne(GIObject)
              anchor_two <- anchorTwo(GIObject)
              
              one.ol = findOverlaps(features, anchor_one)
              two.ol = findOverlaps(features, anchor_two)
              one.indexes = queryHits(one.ol)[anchor_one[ subjectHits(one.ol) ]$node.class == feature.name]  
              two.indexes = queryHits(two.ol)[anchor_two[ subjectHits(two.ol) ]$node.class == feature.name] 
              features.with.interactions.indexes = unique(c(one.indexes, two.indexes))
              features.with.interactions.indexes = features.with.interactions.indexes[order(features.with.interactions.indexes)]
              
              feature.names.with.interactions = feature.names.full[ features.with.interactions.indexes]
              
              anchor_one.df = data.frame(seqnames=as.character(seqnames(anchor_one)),
                                         start=start(anchor_one),
                                         end=end(anchor_one),
                                         width=width(anchor_one),
                                         strand=as.character(strand(anchor_one)),stringsAsFactors=FALSE)
              anchor_two.df = data.frame(seqnames=as.character(seqnames(anchor_two)),
                                         start=start(anchor_two),
                                         end=end(anchor_two),
                                         width=width(anchor_two),
                                         strand=as.character(strand(anchor_two)),stringsAsFactors=FALSE)
              #for(i in features.with.interactions.indexes){
              for(fn in unique(feature.names.with.interactions)){
                #print(fn)
                i = which(summary.df[,paste(capitalize(feature.name), "id", sep=".")]==fn)
                iss = which(feature.names.full==fn)
                #print(i)
                #if(length(iss)>1){
                #  print(fn)
                #  print(i)
                #  print(iss)
                #}
                interactions = unique(c(subjectHits(one.ol[ queryHits(one.ol) %in% iss]), subjectHits(two.ol[ queryHits(two.ol) %in% iss])))
                interactions.one = unique(subjectHits(one.ol[ queryHits(one.ol) %in% iss]))
                interactions.two = unique(subjectHits(two.ol[ queryHits(two.ol) %in% iss]))
                
                numberOfInteractions = sum(interactionCounts(GIObject)[ interactions ])
                numberOfUniqueInteractions = length(interactions)
                summary.df[i,paste("numberOf", capitalize(feature.name), "Interactions", sep="")] = numberOfInteractions
                summary.df[i,paste("numberOf", capitalize(feature.name), "UniqueInteractions", sep="")] = numberOfUniqueInteractions
                
                intercis = interactions[as.character(seqnames(anchor_one[interactions])) != as.character(seqnames( anchor_two[interactions]))]
                if(length(intercis)>0){
                  numberOfInterChromosomalInteractions = sum(interactionCounts(GIObject)[intercis])
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
                    nc1Indexes = interactions.one[which( anchor_one$node.class[interactions.one] == feature.name &  anchor_two$node.class[interactions.one]  == nc)]
                    nc1Counts = sum(interactionCounts(GIObject)[nc1Indexes])
                  }
                  if(length(interactions.two)>0){
                    nc2Indexes = interactions.two[which(anchor_one$node.class[interactions.two] == nc & anchor_two$node.class[interactions.two] == feature.name)]
                    if(nc == feature.name ){ 
                      if(length(nc2Indexes[ (nc2Indexes %in% nc1Indexes) ])>0){
                        summary.df[i,paste("numberOfSelf", capitalize(feature.name), capitalize(nc), "Interactions", sep="")] = sum(interactionCounts(GIObject)[nc2Indexes[(nc2Indexes %in% nc1Indexes)]])
                        summary.df[i,paste("numberOfSelfUnique", capitalize(feature.name), capitalize(nc), "Interactions", sep="")] = length(nc2Indexes[(nc2Indexes %in% nc1Indexes)])
                      }
                      nc2Indexes = nc2Indexes[!(nc2Indexes %in% nc1Indexes)] # stop double counting of interactions where both anchors are in the same feature
                    } 
                    nc2Counts = sum(interactionCounts(GIObject)[nc2Indexes])
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
                        ncs1Indexes = interactions.one[which(anchor_one$node.class[interactions.one] == feature.name &  anchor_two$node.class[interactions.one]  == nc)] 
                        ncs1Indexes = names(features)[i] %in% elementMetadata(anchor_two)[[ paste(nc, "id", sep=".") ]][ncs1Indexes]
                        ncs1Counts = sum(interactionCounts(GIObject)[ncs1Indexes])
                      }
                      if(length(interactions.two)>0){
                        nc2sIndexes = interactions.two[which(anchor_one$node.class[interactions.two] == nc & anchor_two$node.class[interactions.two] == feature.name)]
                        ncs2Indexes = names(features)[i]  %in%  elementMetadata(anchor_one)[[ paste(nc, "id", sep=".") ]][ncs2Indexes]
                        ncs2Counts = sum(interactionCounts(GIObject)[ncs2Indexes])
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
                wdis = rep(distances, interactionCounts(GIObject)[interactions])
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


#' Summarise the number of interactions between two sets of features.
#' 
#' This function will calculate the number of observed interactions between
#' two sets of features provided by the end-user. This allows the summarisation
#' of the number of features of a specific type a particular region is involved in 
#' and how many interactions exist between them. 
#'
#' @param GIObject An annotated GInteractions object
#' @param features.one A GRanges object containing the feature set of interest
#' @param feature.name.one The name of the first feature set of interest
#' @param features.two A GRanges object containing the second feature set of interest
#' @param feature.name.two The name of the second feature set of interest
#'
#' @return A data frame with one line for each range in `features'
#' @rdname summariseByFeaturePairs
#' @docType methods
#' @import GenomicRanges
#' @export
#' @examples
#' data("hic_example_data")
#' data("mm9_refseq_promoters")
#' data("thymus_enhancers")
#' annotateInteractions(hic_example_data, list(promoter = mm9_refseq_promoters, enhancer = thymus_enh))
#' # can be slow so subset of features used for examples
#' p <- unique(unlist(head(regions(hic_example_data)$promoter.id)))
#' e <- unique(unlist(head(regions(hic_example_data)$enhancer.id)))
#' p <- p[!is.na(p)]
#' p <- mm9_refseq_promoters[p]
#' e <- e[!is.na(e)]
#' e <- thymus_enh[e]
#' ep_summary <- summariseByFeaturePairs(hic_example_data, p, "promoter", e, "enhancer")
setGeneric("summariseByFeaturePairs",function(GIObject, features.one, feature.name.one, features.two, feature.name.two){standardGeneric ("summariseByFeaturePairs")})
#' @rdname summariseByFeaturePairs
#' @export
setMethod("summariseByFeaturePairs", "GInteractions", 
          function(GIObject, features.one, feature.name.one, features.two, feature.name.two){  
            if(!("node.class" %in% names(elementMetadata(regions(GIObject))))){
              stop("GIObject has not been annotated")}
              
              potential.node.classes = unique(GIObject@regions$node.class)
              
              feature.one.names.full <- .get_gr_names(features.one)
              feature.one.names <- unique(feature.one.names.full)
              
              feature.two.names.full <- .get_gr_names(features.two)
              feature.two.names <- unique(feature.two.names.full)
              
              x.gi  = subsetByFeatures(GIObject, features.one)
              x.gi  = subsetByFeatures(x.gi, features.two)
              
              anchor_one <- anchorOne(x.gi)
              anchor_two <- anchorTwo(x.gi)
              
              one.one.ol = findOverlaps(features.one, anchor_one)
              two.one.ol = findOverlaps(features.one, anchor_two)
              
              f1.one.f2.two = subjectHits(one.one.ol)[
                which(anchor_one[ subjectHits(one.one.ol) ]$node.class == feature.name.one & 
                      anchor_two[ subjectHits(one.one.ol) ]$node.class == feature.name.two)]
              
              f1.two.f2.one = subjectHits(two.one.ol)[
                which(anchor_two[ subjectHits(two.one.ol) ]$node.class == feature.name.one & 
                      anchor_one[ subjectHits(two.one.ol) ]$node.class == feature.name.two)]
              
              results = NULL
              # this now results in a GI object only contain feature.name.one:feature.name.two interactions
              x.gi = x.gi[unique(c(f1.one.f2.two, f1.two.f2.one)),]
              
              one.one.ol = findOverlaps(features.one, anchor_one)
              two.one.ol = findOverlaps(features.one, anchor_two)
              
              one.two.ol = findOverlaps(features.two, anchor_one)
              two.two.ol = findOverlaps(features.two, anchor_two)
              
              one.indexes = queryHits(one.one.ol)[anchor_one[ subjectHits(one.one.ol) ]$node.class == feature.name.one]  
              two.indexes = queryHits(two.one.ol)[anchor_two[ subjectHits(two.one.ol) ]$node.class == feature.name.one] 
              features.one.with.interactions.indexes = unique(c(one.indexes, two.indexes))
              
              all.features.one.with.interactions = features.one[features.one.with.interactions.indexes]
              feature.ones.id = feature.one.names.full[features.one.with.interactions.indexes]
              for( fn in unique(feature.ones.id)){
                #print(fn)
                iss = which(feature.one.names.full ==fn)
                
                interactions = unique(c(subjectHits(one.one.ol[ queryHits(one.one.ol) %in% iss]), subjectHits(two.one.ol[ queryHits(two.one.ol) %in% iss])))
                interactions.one = unique(subjectHits(one.one.ol[ queryHits(one.one.ol) %in% iss]))
                interactions.two = unique(subjectHits(two.one.ol[ queryHits(two.one.ol) %in% iss]))
                
                features.two.involved.one = unique(unlist(elementMetadata(anchorOne(x.gi))[[paste(feature.name.two, "id", sep=".")]][interactions.two]))
                features.two.involved.two = unique(unlist(elementMetadata(anchorTwo(x.gi))[[paste(feature.name.two, "id", sep=".")]][interactions.one]))
                #print(unique(c(features.two.involved.one, features.two.involved.two)))
                for( fn.two in unique(c(features.two.involved.one, features.two.involved.two))){
                  #print(feature.two)
                  iss.two = which(feature.two.names.full ==fn.two)
                  #print(iss.two)
                  indexes = unique(intersect(interactions, unique(c(subjectHits(one.two.ol[ queryHits(one.two.ol) %in% iss.two]), 
                                                                    subjectHits(two.two.ol[ queryHits(two.two.ol) %in% iss.two])))))
                  counts = sum(interactionCounts(x.gi)[indexes])
                  #print(counts)
                  results = rbind(results, c(fn, fn.two, counts))
                }
              }
              results = data.frame(results[,1], results[,2], as.numeric(results[,3]))
              colnames(results) = c(paste(capitalize(feature.name.one), "id", sep="."), 
                                    paste(capitalize(feature.name.two), "id", sep="."), 
                                    "counts")
              return(results)
          })
            
            
            
            