#' Plot coverage of interactions originating at a given viewpoint.
#' 
#' @param pos A single region in GRanges format. 
#' @param GIObject A GenomicInteractions object.
#' @param leftflank An integer; flank size in bp upstream of pos centre.
#' @param rightflank An integer; flank size in bp downstream of pos centre.
#' @param plot Logical. Whether to plot vector (default), or just return vector of coverage.
#' 
#' @return Plot of coverage or Rle-vector of coverage, depending on plot parameter.
#' 
#' @importFrom IRanges as.vector subjectHits queryHits
#' @export
#' @examples
#' \dontrun{
#' data(hic_data)
#' pos <- GRanges(seqnames="chr5", ranges=IRanges(start=115938063, end=115941352))
#' viewPoint(pos,hic_data, 100000, 100000, plot=FALSE)
#' }

viewPoint <- function(pos, GIObject, leftflank, rightflank, plot=TRUE){
	#pos should be a single region in GRanges format
	#GIObject should be a GenomicInteractions object
	#leftflank and rightflank should be integers giving flank size left and right of pos centre
  #plot: whether to plot (default), or just return vector of coverage 

	if (length(pos) > 1){
		stop("pos should be a single position, you have given ", length(pos), " positions")
	}

	#get all interactions that overlap pos
	fO <- findOverlaps(pos, GIObject)
	#find the other anchor for interactions where one anchor overlaps with pos
	#make new GRanges and create coverage
	#weight by counts
  signal_gr <- c(GIObject@anchor_two[subjectHits(fO$one)], GIObject@anchor_one[subjectHits(fO$two)])
  
	signal_cvg <- coverage(signal_gr, weight= GIObject@counts[c(subjectHits(fO$one), subjectHits(fO$two))])
  
  region_start_actual <- start(resize(pos, width=1, fix="center"))-leftflank
  region_end_actual <- start(resize(pos, width=1, fix="center"))+rightflank
  chr_end <- seqlengths(GIObject@anchor_one)[[as.character(seqnames(pos))]]
  
  #in case region is more than chromosome length
  region_start <- max(1,region_start_actual)
  region_end <- min(chr_end, region_end_actual)
	
  signal_vec <- as.numeric(signal_cvg[[as.character(seqnames(pos))]][region_start : region_end])
  
  #pad with zeroes if necessary
  if (length(signal_vec) < (leftflank+rightflank+1)){
    length_diff = (leftflank+rightflank+1) - length(signal_vec)
    if (region_start_actual < 1){
      signal_vec <- c(rep(0, length_diff), signal_vec)
    }
    if (region_end_actual > chr_end){
      signal_vec <- c(signal_vec, rep(0, length_diff))
    }
  }
  #make spline for quicker plotting
  if (length(signal_gr)==0){
    warning(paste("no interactions here:", seqnames(pos),":", start(pos), "-", end(pos)))
    spline_n <- max(100,ceiling(length(signal_vec/1000)))
  }else{
    spline_n <- ceiling(length(signal_vec) / min(width(signal_gr)))
  }
	
  if (plot){
    #plotting
    plot(spline(signal_vec, n = spline_n), type="l", xaxt="n", ylab="Signal", xlab=as.character(seqnames(pos)))
    
    labs <- seq(from=region_start, to=region_end, length.out=5)
    ats <- seq(from=1, to=length(signal_vec), length.out=length(labs))
    
    axis(side=1, at=ats, labels=labs)
    
    abline(v=leftflank-(width(pos)/2), lty=2)
    abline(v=leftflank+(width(pos)/2), lty=2)
  }else{
    return(signal_vec)
  }
}

#' Plot coverage of interactions originating at a given set of viewpoints.
#' 
#' @param positions A set of regions in GRanges format. 
#' @param GIObject A GenomicInteractions object.
#' @param leftflank An integer; flank size in bp upstream of position centre.
#' @param rightflank An integer; flank size in bp downstream of position centre.
#' @param plot Logical. Whether to plot vector (default), or just return vector of coverage.
#' @param flip Logical. If TRUE, flip window of coverage around a position on the minus strand. Use e.g. to plot around promoters. Default: FALSE.
#' 
#' @return Plot of mean coverage around all positions or Rle-vector of mean coverage, depending on plot parameter.
#' 
#' @importFrom IRanges as.vector subjectHits queryHits
#' @export
viewPointAverage <- function(positions, GIObject, leftflank, rightflank, plot=TRUE, flip=FALSE){
  #positions should be regions in GRanges format
  #GIObject should be a GenomicInteractions object
  #leftflank and rightflank should be integers giving flank size left and right of positions centre
  #plot: whether to plot (default), or just return matrix of coverage
  #flip: whether to flip coverage vectors around a position on the - strand: use for e.g. signal around promoters.
  
  if (!flip){
    #do not flip according to strand
    #pre-allocate matrix (takes a long time if matrix is very large)
    #do not recommend flanks larger than ~200kb
    results_vec <- rep(0, times=leftflank+rightflank+1)
    
    #fill matrix
    for (i in 1:length(positions)){
      results_vec <- results_vec + as.numeric(viewPoint(positions[i], GIObject, leftflank, rightflank, plot=FALSE))
    }
  }else{
    pos_plus <- positions[as.character(strand(positions)) %in% c("+", "*")]
    pos_minus <- positions[as.character(strand(positions)) =="-"]
    
    results_plus <- rep(0, times=leftflank+rightflank+1)
    #fill matrix
    for (i in 1:length(pos_plus)){
      results_plus <- results_plus + as.vector(viewPoint(pos_plus[i], GIObject, leftflank, rightflank, plot=FALSE))
    }
    
    results_minus <- rep(0, times=leftflank+rightflank+1)
    #fill matrix
    for (i in 1:length(pos_minus)){
      #flank lengths swapped here
      results_minus <- results_minus + rev(as.vector(viewPoint(pos_minus[i], GIObject, rightflank, leftflank, plot=FALSE)))
    }
    #rownames(results_plus) <- which(strand(positions) %in% c("+", "*"))
    #rownames(results_minus) <- which(strand(positions)=="-")
    
    #combine
    results_vec <- results_minus + results_plus
    #reorder to match original order of positions
    #results_matrix <- results_matrix[order(as.numeric(rownames(results_matrix))),]
  }

  #plot
 if (plot){
   plot(spline(results_vec/length(positions), n=max(c(100,(length(results_vec)/1000)))), 
        type="l", ylab="Signal", xlab="Distance from midpoint", xaxt="n")
   
   labs <- seq(from=-leftflank, to=rightflank, length.out=5)
   ats <- seq(from=1, to=length(results_vec), length.out=length(labs))
   
   axis(side=1, at=ats, labels=labs)
   
 }else{
   return(results_vec/length(positions))
 }
  
}  
