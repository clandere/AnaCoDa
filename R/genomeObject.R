#' Genome Initialization
#' 
#' \code{initializeGenomeObject} initializes the Rcpp Genome object
#' 
#' @param file A file of coding sequences in fasta or RFPData format
#' 
#' @param genome A genome object can be passed in to concatenate the input file to it (optional).
#' 
#' @param observed.expression.file A string containing the location of a file containing
#'  empirical expression rates (optional).
#' 
#' @param fasta A boolean value which decides whether to initialize with a
#'  fasta file or an RFPData file. (TRUE for fasta, FALSE for RFPData)
#' 
#' @param match.expression.by.id If TRUE (default), observed expression values will be assigned by matching sequence identifier.
#' If FALSE, observed expression values will be assigned by order.
#' 
#' @param append If TRUE (FALSE is default), function will read in additional genome data to append to an existing genome.
#' If FALSE, genome data is cleared before reading in data (no preexisting data). 
#' 
#' @return This function returns the initialized Genome object.
#' 
initializeGenomeObject <- function(file, genome=NULL, observed.expression.file=NULL, fasta=TRUE, match.expression.by.id=TRUE, append=FALSE) {
  if (is.null(genome)){ 
    genome <- new(Genome)
  }

  if (fasta == TRUE) {
    genome$readFasta(file, append)
  } else {
    genome$readRFPData(file, append)
  }
  if(!is.null(observed.expression.file)) {
    genome$readObservedPhiValues(observed.expression.file, match.expression.by.id)
  }
  return(genome)
}


#' Get Codon Counts For Each Amino Acid 
#' 
#' @param genome A genome object from which the counts of each
#' codon can be obtained.
#'  
#' @param aa A one character representation of an amino acid.
#' 
#' @return Returns a matrix storing the codonCounts for the given amino acid. 
#' 
#' @description provides the codon counts for a fiven amino acid across all genes
#' 
#' @details The returned matrix containes a row for each gene and a coloumn for each codon.
#' 
getCodonCountsForAA <- function(genome, aa){
  # get codon count for aa
  codons <- AAToCodon(aa, F)
  codonCounts <- lapply(codons, function(codon){
    codonCounts <- genome$getCodonCountsPerGene(codon)
  })
  codonCounts <- do.call("cbind", codonCounts)
  return(codonCounts)
}



#' Length of Genome
#' 
#' \code{length} gives the length of a genome
#' 
#' @param x A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @return returns the number of genes in a genome
length.Rcpp_Genome <- function(x) {
  return(x$getGenomeSize(F))
}

#' Summary of Genome
#' 
#' \code{summary} summarizes the description of a genome, such as number of genes and average gene length.
#' 
#' @param object A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @param ... Optional, additional arguments to be passed to the main summary function 
#' that affect the summary produced.
#'
#' @return This function returns by default an object of class c("summaryDefault", table").
summary.Rcpp_Genome <- function(object, ...) {
  # TODO output stuff like:
  # - no. of genes
  # - avg. gene length
  # - avg. A,C,G,T content
  # - avg. AA composition
  # - ...
  summary(object, ...)
}

#' Gene Names of Genome
#' 
#' 
#' @param genome A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @param simulated A logical value denoting if the gene names to be listed are simulated or not.
#' The default value is FALSE.
#' 
#' @description returns the identifiers of the genes within the genome specified.
#' 
#' @return gene.names Returns the names of the genes as a vector of strings.
#' 
getNames <- function(genome, simulated = FALSE)
{
  genes <- genome$getGenes(simulated)
  gene.names <- unlist(lapply(1:length(genes), function(i){return(genes[[i]]$id)}))
  return(gene.names)
}


#' Add gene observed synthesis rates
#' 
#' \code{addObservedSynthesisRateSet} returns the observed synthesis rates of the genes within the genome specified.
#' 
#' @param genome A genome object initialized with \code{\link{initializeGenomeObject}} to add observed expression data.
#' 
#' @param observed.expression.file A string containing the location of a file containing
#'  empirical expression rates (optional).
#' 
#' @param match.expression.by.id If TRUE (default) observed expression values will be assigned by matching sequence identifier.
#' If FALSE observed expression values will be assigned by order
#' 
#' @return gene.names Returns the names of the genes as a vector of strings.
#' 
addObservedSynthesisRateSet <- function(genome, observed.expression.file, match.expression.by.id=TRUE)
{
  genome$readObservedPhiValues(observed.expression.file, match.expression.by.id)
  return(genome)
}


#' Get gene observed synthesis rates
#' 
#' \code{getObservedSynthesisRateSet} returns the observed synthesis rates of the genes within the genome specified.
#' 
#' @param genome A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @param simulated A logical value denoting if the synthesis rates to be listed are simulated or not.
#' The default value is FALSE.
#' 
#' @return gene.names Returns the names of the genes as a vector of strings.
#' 
getObservedSynthesisRateSet <- function(genome, simulated = FALSE)
{
  genes <- genome$getGenes(simulated)
  expression <- lapply(1:length(genes), function(i){return(genes[[i]]$getObservedSynthesisRateValues())})
  
  return(do.call(rbind, expression))
}


