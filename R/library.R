.loadGenome = function(genome){
  require(genome, character.only=TRUE)
  genome.name = unlist(strsplit(genome, split='\\.'))
  return(get(genome.name[2]))
}

#' Capitalize first letter of string
#'
#' This function will capitalize the first letter of each string in
#' a character vector, and lowercase following letters.
#'
#' @param x A character vector
#' 
#' @return a string with the first letter capitalised
capitalize = function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
}


.pasteAnchor = function(seqname, start, end){
    paste(paste(seqname, start, sep=":"), end, sep="..")
}

.generateInteractionName = function(anchorString1, anchorString2, counts){
    paste(paste(anchorString1, anchorString2, sep="-"), counts, sep=",")
}