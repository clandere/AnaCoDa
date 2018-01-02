#' Initialize Parameter 
#' 
#' @param genome An object of type Genome necessary for the initialization of the Parameter object.
#' The default value is NULL.
#' 
#' @param sphi Initial values for sphi. Expected is a vector of length numMixtures.
#' The default value is NULL.
#' 
#' @param num.mixtures The number of mixtures elements for the underlying mixture distribution (numMixtures > 0).
#' The default value is 1.
#' 
#' @param gene.assignment A vector holding the initial mixture assignment for each gene. 
#' The vector length has to equal the number of genes in the genome.
#' Valid values for the vector range from 1 to numMixtures. 
#' It is possible but not advised to leave a mixture element empty.
#' The default Value is NULL.
#' 
#' @param initial.expression.values (Optional) A vector with intial phi values.
#' The length of the vector has to equal the number of genes in the Genome object.
#' The default value is NULL.
#' 
#' @param model Specifies the model used. Valid options are "ROC", "PA", "PANSE", or "FONSE".
#' The default model is "ROC".
#' ROC is described in Gilchrist et al. 2015.
#' PA, PANSE and FONSE are currently unpublished.
#' 
#' @param split.serine Whether serine should be considered as 
#' one or two amino acids when running the model.
#' TRUE and FALSE are the only valid values.
#' The default value for split.serine is TRUE.
#' 
#' @param mixture.definition A string describing how each mixture should
#' be treated with respect to mutation and selection.
#' Valid values consist of "allUnique", "mutationShared", and "selectionShared".
#' The default value for mixture.definition is "allUnique".
#' See details for more information.
#' 
#' @param mixture.definition.matrix A matrix representation of how
#' the mutation and selection categories correspond to the mixtures.
#' The default value for mixture.definition.matrix is NULL.
#' If provided, the model will use the matrix to initialize the mutation and selection
#' categories instead of the definition listed directly above.
#' See details for more information.
#' 
#' @param init.with.restart.file File name containing information to reinitialize a 
#' previous Parameter object.
#' If given, all other arguments will be ignored.
#' The default value for init.with.restart.file is NULL.
#' 
#' @param mutation.prior.sd Controlling the standard deviation of the normal 
#' prior on the mutation parameters
#' 
#' @param init.csp.variance specifies the initial proposal width for codon specific parameter (default is 0.0025). 
#' The proposal width adapts during the runtime to reach a taget acceptance rate of ~0.25
#' 
#' @param init.sepsilon specifies the initial value for sepsilon. default is 0.1
#' 
#' @return parameter Returns an initialized Parameter object.
#' 
#' @description \code{initializeParameterObject} initializes a new parameter object or reconstructs one from a restart file
#' 
#' @details \code{initializeParameterObject} checks the values of the arguments 
#' given to insure the values are valid.
#' 
#' The mixture definition and mixture definition matrix describes how the mutation
#' and selection categories are set up with respect to the number of mixtures. For
#' example, if mixture.definition = "allUnique" and numMixtures = 3, a matrix
#' representation would be \code{matrix(c(1,2,3,1,2,3), ncol=2)}
#' where each row represents a mixture, the first column represents the mutation
#' category, and the second column represents the selection category.
#' Another example would be mixture.definition = "selectionShared" and numMixtures = 4 (
#' \code{matrix(c(1,2,3, 4,1,1,1,1), ncol=2)}).
#' In this case, the selection category is the same for every mixture. If a matrix
#' is given, and it is valid, then the mutation/selection relationship will be
#' defined by the given matrix and the keyword will be ignored. A matrix should only
#' be given in cases where the keywords would not create the desired matrix.

initializeParameterObject <- function(genome = NULL, sphi = NULL, num.mixtures = 1, 
                                    gene.assignment = NULL, initial.expression.values = NULL,
                                    model = "ROC", split.serine = TRUE, 
                                    mixture.definition = "allUnique", 
                                    mixture.definition.matrix = NULL,
                                    init.with.restart.file = NULL, mutation.prior.sd = 0.35, 
				    init.csp.variance = 0.0025, init.sepsilon = 0.1){
  # check input integrity
  if(is.null(init.with.restart.file)){
    if(length(sphi) != num.mixtures){
      stop("Not all mixtures have an Sphi value assigned!\n")
    }
    if(length(genome) != length(gene.assignment)){
      stop("Not all Genes have a mixture assignment!\n")
    }
    if(max(gene.assignment) > num.mixtures){
      stop("Gene is assigned to non existing mixture!\n")
    }
    if(num.mixtures < 1){
      stop("num. mixture has to be a positive non-zero value!\n")
    }    
    #TODO: should we check integrity of other values, such as numMixtures being
    #positive?
  }

  
  
  
  if(model == "ROC"){
    if(is.null(init.with.restart.file)){
      parameter <- initializeROCParameterObject(genome, sphi, num.mixtures, 
                                                gene.assignment, initial.expression.values, split.serine, 
                            mixture.definition, mixture.definition.matrix, 
                            mutation.prior.sd, init.csp.variance, init.sepsilon)    
    }else{
      parameter <- new(ROCParameter, init.with.restart.file)
    }
  }else if(model == "FONSE"){
    if(is.null(init.with.restart.file)){
      parameter <- initializeFONSEParameterObject(genome, sphi, num.mixtures, 
                                                  gene.assignment, initial.expression.values, split.serine, 
                            mixture.definition, mixture.definition.matrix, init.csp.variance)
    }else{
      parameter <- new(FONSEParameter, init.with.restart.file)
    }
  }else if(model == "PA"){
    if(is.null(init.with.restart.file)){
      parameter <- initializePAParameterObject(genome, sphi, num.mixtures, 
                                                gene.assignment, initial.expression.values, split.serine, 
                            mixture.definition, mixture.definition.matrix, init.csp.variance) 
    }else{
      parameter <- new(PAParameter, init.with.restart.file)
    }
  }else if(model == "PANSE"){
    if(is.null(init.with.restart.file)){
      parameter <- initializePANSEParameterObject(genome, sphi, num.mixtures, 
                                                gene.assignment, initial.expression.values, split.serine, 
                            mixture.definition, mixture.definition.matrix, init.csp.variance) 
    }else{
      parameter <- new(PANSEParameter, init.with.restart.file)
    }
  }else{
    stop("Unknown model.")
  }
  
  return(parameter)
}


#Called from initializeParameterObject. 
initializeROCParameterObject <- function(genome, sphi, numMixtures, geneAssignment,
                      expressionValues = NULL, split.serine = TRUE,
                      mixture.definition = "allUnique", 
                      mixture.definition.matrix = NULL, mutation_prior_sd = 0.35, init.csp.variance = 0.0025, init.sepsilon = 0.1){

  if(is.null(mixture.definition.matrix)){ 
    # keyword constructor
    parameter <- new(ROCParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(ROCParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeSynthesisRateByGenome(genome, mean(sphi))
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  n.obs.phi.sets <- ncol(getObservedSynthesisRateSet(genome))
  parameter$setNumObservedSynthesisRateSets(n.obs.phi.sets)

  parameter$mutation_prior_sd <- mutation_prior_sd

  if (n.obs.phi.sets != 0){
   parameter$setInitialValuesForSepsilon(as.vector(init.sepsilon))
  }

  parameter <- initializeCovarianceMatrices(parameter, genome, numMixtures, geneAssignment, init.csp.variance)
   
  return(parameter)
}


#Called from initializeParameterObject.
initializePAParameterObject <- function(genome, sphi, numMixtures, geneAssignment, 
                          expressionValues = NULL, split.serine = TRUE, 
                          mixture.definition = "allUnique", 
                          mixture.definition.matrix = NULL, init.csp.variance){

  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(PAParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(PAParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeSynthesisRateByGenome(genome, sphi)
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  ## TODO (Cedric): use init.csp.variance to set initial proposal width for CSP parameters
  
  return (parameter)
}

#Called from initializeParameterObject.
initializePANSEParameterObject <- function(genome, sphi, numMixtures, geneAssignment, 
                          expressionValues = NULL, split.serine = TRUE, 
                          mixture.definition = "allUnique", 
                          mixture.definition.matrix = NULL, init.csp.variance){

  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(PANSEParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(PANSEParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeSynthesisRateByGenome(genome, sphi)
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  return (parameter)
}

#Called from initializeParameterObject.
initializeFONSEParameterObject <- function(genome, sphi, numMixtures, 
                        geneAssignment, expressionValues = NULL, split.serine = TRUE,
                        mixture.definition = "allUnique", 
                        mixture.definition.matrix = NULL, init.csp.variance){

  # create Parameter object
  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(FONSEParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(FONSEParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeSynthesisRateByGenome(genome, sphi)
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  parameter <- initializeCovarianceMatrices(parameter, genome, numMixtures, geneAssignment, init.csp.variance)
  
  return(parameter)
}



#' Return Codon Specific Paramters (or write to csv) estimates as data.frame
#' 
#' @param parameter parameter on object created by \code{initializeParameterObject}.
#' 
#' @param filename Posterior estimates will be written to file instead of returned if specified (format: csv).
#' 
#' @param CSP which type of codon specific parameter should be returned (mutation (default) or selection)
#' 
#' @param mixture estimates for which mixture should be returned
#' 
#' @param samples The number of samples used for the posterior estimates.
#' 
#' @return returns a data.frame with the posterior estimates of the models 
#' codon specific parameters or writes it directly to a csv file if \code{filename} is specified
#' 
#' @description \code{getCSPEstimates} returns the codon specific
#' parameter estimates for a given parameter and mixture or write it to a csv file.
#' 
#' 
getCSPEstimates <- function(parameter, filename=NULL, CSP="Mutation", mixture = 1, samples = 10){
  Amino_Acid <- c()
  Value <- c()
  Codon <- c()
  quantile_list <- vector("list")
  
  if (class(parameter) == "Rcpp_ROCParameter" || class(parameter) == "Rcpp_FONSEParameter"){
    names.aa <- aminoAcids()
    
    for(aa in names.aa){
      if(aa == "M" || aa == "W" || aa == "X") next
      codons <- AAToCodon(aa, T)
    
      for(i in 1:length(codons)){
        Amino_Acid <- c(Amino_Acid, aa)
        Codon <- c(Codon, codons[i])
      
        if(CSP == "Mutation"){
          Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, TRUE))
          quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, codons[i], 0, c(0.025, 0.975), TRUE))
        }
        else if(CSP == "Selection"){
          Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, TRUE))
          quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, codons[i], 1, c(0.025, 0.975), TRUE))
        }
        else {
          stop("Unknown parameter type given with argument: CSP")
        }
      }
    }
  }
  else if (class(parameter) == "Rcpp_PAParameter"){
    groupList <- parameter$getGroupList()

    for(i in 1:length(groupList)){
      aa <- codonToAA(groupList[i])
      Codon <- c(Codon, groupList[i])
      Amino_Acid <- c(Amino_Acid, aa)
       
      if(CSP == "Alpha"){
        Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, FALSE))
        quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, codons[i], 0, c(0.025, 0.975), FALSE))
        }
      else if(CSP == "Lambda Prime"){
        Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, FALSE))
        quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, codons[i], 1, c(0.025, 0.975), FALSE))
      }
      else {
        stop("Unknown parameter type given with argument: CSP")
      }
    }
  }
  else if (class(parameter) == "Rcpp_PANSEParameter"){
    groupList <- parameter$getGroupList()

    for(i in 1:length(groupList)){
      aa <- codonToAA(groupList[i])
      Codon <- c(Codon, groupList[i])
      Amino_Acid <- c(Amino_Acid, aa)
       
      if(CSP == "Alpha"){
        Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, FALSE))
        quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, codons[i], 0, c(0.025, 0.975), FALSE))
        }
      else if(CSP == "Lambda Prime"){
        Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, FALSE))
        quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, codons[i], 1, c(0.025, 0.975), FALSE))
      }
      else {
        stop("Unknown parameter type given with argument: CSP")
      }
    }
  }
  else{
    stop("Unknown object provided with argument: parameter")
  }

  quantile_list <- matrix(unlist(quantile_list), nrow = 2)
  data <- data.frame(Amino_Acid, Codon, Value, Lower=quantile_list[1,], Upper=quantile_list[2,])
  colnames(data) <- c("AA", "Codon", "Posterior", "0.025%", "0.975%")
  if(is.null(filename))
  {
    return(data)
  }else {
    write.csv(data, file = filename, row.names = FALSE, quote=FALSE)
  }
}



# Uses a multinomial logistic regression to estimate the codon specific parameters for every category.
# Delta M is the intercept - and Delta eta is the slope of the regression.
# The package VGAM is used to perform the regression.
getCSPbyLogit <- function(codonCounts, phi, coefstart = NULL, x.arg = FALSE, 
                          y.arg = FALSE, qr.arg = FALSE){
  #avoid cases with 0 aa count
  idx <- rowSums(codonCounts) != 0
  
  # performs the regression and returns Delta M and Delta eta as well as other information no used here
  ret <- VGAM::vglm(codonCounts[idx, ] ~ phi[idx],
                    VGAM::multinomial, coefstart = coefstart,
                    x.arg = x.arg, y.arg = y.arg, qr.arg = qr.arg)
  coefficients <- ret@coefficients
  
  ## convert delta.t to delta.eta
  coefficients <- -coefficients
  
  ret <- list(coefficients = coefficients,
              coef.mat = matrix(coefficients, nrow = 2, byrow = TRUE),
              R = ret@R)
  return(ret)
}


#TODO: Need comments explaining what is going on
subMatrices <- function(M, r, c){
  rg <- (row(M) - 1) %/% r + 1
  cg <- (col(M) - 1) %/% c + 1
  rci <- (rg - 1) * max(cg) + cg
  return(rci)
}



#TODO: Need comments explaining what is going on
splitMatrix <- function(M, r, c){
  rci <- subMatrices(M, r, c)
  N <- prod(dim(M)) / r / c
  cv <- lapply(1:N, function(x) M[rci==x])
  
  return(lapply(1:N, function(i) matrix(cv[[i]], nrow = r)))
} 



#######
### CURRENTLY NOT EXPOSED
#######
#' Initialize Covariance Matrices
#' 
#' @param parameter A Parameter object that corresponds to one of the model types. 
#' Valid values are "ROC", "PA", and "FONSE".
#' 
#' @param genome An object of type Genome necessary for the initialization of the Parameter object. 
#' 
#' @param numMixtures The number of mixture elements for the underlying mixture distribution (numMixtures > 0).
#' 
#' @param geneAssignment A vector holding the initial mixture assignment for each gene.
#' The vector length has to equal the number of genes in the genome.
#' Valid values for the vector range from 1 to numMixtures.
#' It is possible but not advised to leave a mixture element empty.
#' 
#' @param init.csp.variance initial proposal variance for codon specific parameter, default is 0.0025.
#' 
#' @return parameter Returns the Parameter argument, now modified with initialized mutation, selection, and covariance matrices.
#' 

# Also initializes the mutaiton and selection parameter
initializeCovarianceMatrices <- function(parameter, genome, numMixtures, geneAssignment, init.csp.variance = 0.0025) {
  numMutationCategory <- parameter$numMutationCategories
  numSelectionCategory <- parameter$numSelectionCategories
  
  phi <- parameter$getCurrentSynthesisRateForMixture(1) # phi values are all the same initially

  names.aa <- aminoAcids()
 # ct <- getInstance()
#  names.aa <- ct$getGroupList()
  
  for(aa in names.aa){
    if(aa == "M" || aa == "W" || aa == "X") next
    #should go away when CT is up and running
    
    codonCounts <- getCodonCountsForAA(genome, aa)
    numCodons <- dim(codonCounts)[2] - 1
    #-----------------------------------------
    # TODO WORKS CURRENTLY ONLY FOR ALLUNIQUE!
    #-----------------------------------------
    covmat <- vector("list", numMixtures)
    for(mixElement in 1:numMixtures){    
      idx <- geneAssignment == mixElement
      csp <- getCSPbyLogit(codonCounts[idx, ], phi[idx])
      
      parameter$initMutation(csp$coef.mat[1,], mixElement, aa)
      parameter$initSelection(csp$coef.mat[2,], mixElement, aa)
    }
    
    # One covariance matrix for all mixtures.
    # Currently only variances used.
    compl.covMat <- diag((numMutationCategory + numSelectionCategory) * numCodons) * init.csp.variance
    parameter$initCovarianceMatrix(compl.covMat, aa)
  }
  
  
  
  #for(aa in names.aa){
  # if(aa == "M" || aa == "W" || aa == "X") next
  #should go away when CT is up and running
  
  #codonCounts <- getCodonCountsForAA(genome, aa)
  #numCodons <- dim(codonCounts)[2] - 1
  #-----------------------------------------
  # TODO WORKS CURRENTLY ONLY FOR ALLUNIQUE!
  #-----------------------------------------
  # covmat <- vector("list", numMixtures)
  #for(mixElement in 1:numMixtures){    
  # idx <- geneAssignment == mixElement
  #csp <- getCSPbyLogit(codonCounts[idx, ], phi[idx])
  
  #   parameter$initMutation(csp$coef.mat[1,], mixElement, aa)
  #  parameter$initSelection(csp$coef.mat[2,], mixElement, aa)
  # split matrix into sup matrices (dM and dEta)
  # covmat[[mixElement]] <- splitMatrix(t(csp$R) %*% csp$R, numCodons, numCodons)  # we expect the covariance matrix, but get the decomposition.
  #  }
  # compl.covMat <- matrix(0, ncol = numMixtures * numCodons * 2, nrow = numMixtures * numCodons * 2)
  #matrix.positions <- subMatrices(compl.covMat, numCodons, numCodons)
  
  #compl.seq <- seq(1, dim(compl.covMat)[1], numCodons)
  #mut.seq <- compl.seq[1:(length(compl.seq)/2)]
  #i <- 1
  #for(pos in mut.seq){ 
  # compl.covMat[matrix.positions == matrix.positions[pos, pos]] <- unlist(covmat[[i]][1])
  #  i <- i + 1
  # i <- ifelse(i > numMutationCategory, 1, i)
  #  }
  # sel.seq <- compl.seq[(length(compl.seq)/2 + 1):length(compl.seq)]
  #  i <- 1
  # for(pos in sel.seq){ 
  #  compl.covMat[matrix.positions == matrix.positions[pos, pos]] <- unlist(covmat[[i]][4])
  # i <- i + 1
  #i <- ifelse(i > numMutationCategory, 1, i)
  #}
  
  #ofdiag.seq <- mut.seq + numCodons*numMutationCategory
  #for(i in 1:length(mut.seq)){
  #  compl.covMat[matrix.positions == matrix.positions[mut.seq[i], ofdiag.seq[i]]] <- unlist(covmat[[i]][2])
  # compl.covMat[matrix.positions == matrix.positions[ofdiag.seq[i], mut.seq[i]]] <- unlist(covmat[[i]][3])
  #}
  #for testing - in actuality this is used, it is currently overwriting 
  #previous steps.
  #compl.covMat <- diag((numMutationCategory + numSelectionCategory) * numCodons) * 0.05
  #compl.covMat / max(compl.covMat)
  #parameter$initCovarianceMatrix(compl.covMat, aa)
  #}
    
  return(parameter)
}


#' Returns mixture assignment estimates for each gene
#' 
#' @param parameter on object created by \code{initializeParameterObject}
#' 
#' @param gene.index a integer or vector of integers representing the gene(s) of interesst.
#' 
#' @param samples number of samples for the posterior estimate
#' 
#' @return returns a vector with the mixture assignment of each gene corresbonding to \code{gene.index} in the same order as the genome. 
#'
#' @description Posterior estimates for the mixture assignment of specified genes
#' 
#' @details The returned vector is unnamed as gene ids are only stored in the \code{genome} object, 
#' but the \code{gene.index} vector can be used to match the assignment to the genome.
#' 
#' @examples 
#' \dontrun{
#' # get the mixture assignment for all genes
#' getMixtureAssignmentEstimate(parameter, 1:length(genome), 1000)
#' 
#' # get the mixture assignment for a subsample
#' getMixtureAssignmentEstimate(parameter, 5:100, 1000)
#' # or
#' getMixtureAssignmentEstimate(parameter, c(10, 30:50, 3, 90), 1000)
#' }

getMixtureAssignmentEstimate <- function(parameter, gene.index, samples)
{
  mixtureAssignment <- unlist(lapply(gene.index,  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
  return(mixtureAssignment)
}


#' Returns the estimated phi posterior for a gene
#' 
#' @param parameter on object created by \code{initializeParameterObject}.
#' 
#' @param gene.index a integer or vector of integers representing the gene(s) of interesst.
#' 
#' @param mixtureAssignment a vector with the same length as \code{gene.index} decsribing the mixture assignment of the gene
#' 
#' @param samples number of samples for the posterior estimate
#' 
#' @return returns a vector with the mixture assignment of each gene corresbonding to \code{gene.index} in the same order as the genome. 
#'
#' @description Posterior estimates for the phi value of specified genes
#' 
#' @details The returned vector is unnamed as gene ids are only stored in the \code{genome} object, 
#' but the \code{gene.index} vector can be used to match the assignment to the genome.
#' 
#' @examples 
#' \dontrun{
#' # get the estimated phi values for all genes assuming they are in mixture 1
#' getExpressionEstimatesForMixture(parameter, 1:length(genome), rep(1, length(genome), 1000)
#' 
#' # get the estimated phi values for all genes for the estimated mixture element
#' mix.assign <- getMixtureAssignmentEstimate(parameter, 1:length(genome), 1000)
#' getExpressionEstimatesForMixture(parameter, 1:length(genome), mix.assign, 1000)
#' }

getExpressionEstimatesForMixture <- function(parameter, gene.index, mixtureAssignment, samples)
{
  expressionValues <- unlist(lapply(gene.index, function(geneIndex){ 
    expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex]) 
    parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples, geneIndex, expressionCategory) 
  }))
  
  expressionStdErr <- sqrt(unlist(lapply(gene.index, function(geneIndex){ 
    expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex]) 
    parameter$getSynthesisRateVarianceByMixtureElementForGene(samples, geneIndex, expressionCategory, T) 
  }))) / samples
  
  expr.mat <- cbind(expressionValues, expressionStdErr)
  colnames(expr.mat) <- c("PHI", "Std_Error")
  return(expr.mat)
}

#' Write Parameter Object to a File
#' 
#' @param parameter parameter on object created by \code{initializeParameterObject}.
#' 
#' @param file A filename that where the data will be stored.
#' 
#' @return This function has no return value.
#' 
#' @description \code{writeParameterObject} will write the parameter object as binary to the filesystem
#' 
#' @details As Rcpp object are not serializable with the default R \code{save} function, 
#' therefore this custom save function is provided (see \link{loadParameterObject}).
#' 
#' @examples 
#' \dontrun{
#' writeParameterObject(parameter, "file.Rda")
#' }
#' 
writeParameterObject <- function(parameter, file)
{
  UseMethod("writeParameterObject", parameter)
}


# extracts traces and parameter information from the base class Parameter
extractBaseInfo <- function(parameter){
  trace <- parameter$getTraceObject()
  stdDevSynthesisRateTraces <- trace$getStdDevSynthesisRateTraces()
  stdDevSynthesisRateAcceptRatTrace <- trace$getStdDevSynthesisRateAcceptanceRateTrace()
  synthRateTrace <- trace$getSynthesisRateTrace()
  synthAcceptRatTrace <- trace$getSynthesisRateAcceptanceRateTrace()
  mixAssignTrace <- trace$getMixtureAssignmentTrace()
  mixProbTrace <- trace$getMixtureProbabilitiesTrace()
  codonSpecificAcceptRatTrace <- trace$getCodonSpecificAcceptanceRateTrace()
  numMix <- parameter$numMixtures
  numMut <- parameter$numMutationCategories
  numSel <- parameter$numSelectionCategories
  categories <- parameter$getCategories()
  curMixAssignment <- parameter$getMixtureAssignment()
  lastIteration <- parameter$getLastIteration()
  
  varList <- list(stdDevSynthesisRateTraces = stdDevSynthesisRateTraces, 
                    stdDevSynthesisRateAcceptRatTrace = stdDevSynthesisRateAcceptRatTrace,
                    synthRateTrace = synthRateTrace,
                    synthAcceptRatTrace = synthAcceptRatTrace,
                    mixAssignTrace = mixAssignTrace,
                    mixProbTrace = mixProbTrace,
                    codonSpecificAcceptRatTrace = codonSpecificAcceptRatTrace,
                    numMix = numMix,
                    numMut = numMut,
                    numSel = numSel,
                    categories = categories,
                    curMixAssignment = curMixAssignment,
                    lastIteration = lastIteration
                    )
  return(varList)
}


#called from "writeParameterObject."
writeParameterObject.Rcpp_ROCParameter <- function(parameter, file){
  paramBase <- extractBaseInfo(parameter)
  
  currentMutation <- parameter$currentMutationParameter
  currentSelection <- parameter$currentSelectionParameter
  proposedMutation <- parameter$proposedMutationParameter
  proposedSelection <- parameter$proposedSelectionParameter
  model = "ROC"
  mutationPrior <- parameter$getMutationPriorStandardDeviation()
  
  trace <- parameter$getTraceObject()
  
  mutationTrace <- trace$getCodonSpecificParameterTrace(0)
  selectionTrace <- trace$getCodonSpecificParameterTrace(1)
  synthesisOffsetAcceptRatTrace <- trace$getSynthesisOffsetAcceptanceRateTrace()
  synthesisOffsetTrace <- trace$getSynthesisOffsetTrace()
  observedSynthesisNoiseTrace <- trace$getObservedSynthesisNoiseTrace()
  if (length(synthesisOffsetTrace) == 0){
    withPhi = FALSE
  }else{
    withPhi = TRUE
  }
  
  save(list = c("paramBase", "currentMutation", "currentSelection",
                "proposedMutation", "proposedSelection", "model",  
                "mutationPrior", "mutationTrace", "selectionTrace", 
                "synthesisOffsetAcceptRatTrace", "synthesisOffsetTrace", 
                "observedSynthesisNoiseTrace", "withPhi"),
       file=file)
}


#called from "writeParameterObject."
writeParameterObject.Rcpp_PAParameter <- function(parameter, file){
  paramBase <- extractBaseInfo(parameter)
  
  currentAlpha <- parameter$currentAlphaParameter
  currentLambdaPrime <- parameter$currentLambdaPrimeParameter
  proposedAlpha <- parameter$proposedAlphaParameter
  proposedLambdaPrime <- parameter$proposedLambdaPrimeParameter
  model = "PA"
  
  
  trace <- parameter$getTraceObject()
  alphaTrace <- trace$getCodonSpecificParameterTrace(0)
  lambdaPrimeTrace <- trace$getCodonSpecificParameterTrace(1)

  save(list = c("paramBase", "currentAlpha", "currentLambdaPrime", "proposedAlpha",
                "proposedLambdaPrime", "model", "alphaTrace", "lambdaPrimeTrace"),
       file=file)
}

#called from "writeParameterObject."
writeParameterObject.Rcpp_PANSEParameter <- function(parameter, file){
  paramBase <- extractBaseInfo(parameter)
  
  currentAlpha <- parameter$currentAlphaParameter
  currentLambdaPrime <- parameter$currentLambdaPrimeParameter
  proposedAlpha <- parameter$proposedAlphaParameter
  proposedLambdaPrime <- parameter$proposedLambdaPrimeParameter
  model = "PANSE"
  
  
  trace <- parameter$getTraceObject()
  alphaTrace <- trace$getCodonSpecificParameterTrace(0)
  lambdaPrimeTrace <- trace$getCodonSpecificParameterTrace(1)

  save(list = c("paramBase", "currentAlpha", "currentLambdaPrime", "proposedAlpha",
                "proposedLambdaPrime", "model", "alphaTrace", "lambdaPrimeTrace"),
       file=file)
}


#called from "writeParameterObject."
writeParameterObject.Rcpp_FONSEParameter <- function(parameter, file)
{
  paramBase <- extractBaseInfo(parameter)
  
  currentMutation <- parameter$currentMutationParameter
  currentSelection <- parameter$currentSelectionParameter

  model = "FONSE"
  mutationPrior <- parameter$getMutationPriorStandardDeviation()
  
  trace <- parameter$getTraceObject()
  
  mutationTrace <- trace$getCodonSpecificParameterTrace(0)
  selectionTrace <- trace$getCodonSpecificParameterTrace(1)
  
  save(list = c("paramBase", "currentMutation", "currentSelection",
                "model","mutationPrior", "mutationTrace", "selectionTrace"),
       file=file)
}




#' Load Parameter Object
#'  
#' @param files A list of parameter filenames to be loaded. If multiple files are given,
#' the parameter objects will be concatenated in the order provided
#' 
#' @return Returns an initialized Parameter object.
#' 
#' @description \code{loadParameterObject} will load a parameter object from the filesystem
#' 
#' @details The function loads one or multiple files. In the case of multiple file, e.g. due to the use of check pointing, the files will
#' be concatenated to one parameter object. See \link{writeParameterObject} for the writing of parameter objects
#' 
#' @examples 
#' \dontrun{
#' # load a single parameter object
#' parameter <- loadParameterObject("path/to/rda/file.Rda")
#' 
#' # load and concatenate multiple parameter object
#' parameter <- loadParameterObject(c("path/to/rda/file.Rda", "path/to/rda/file2.Rda"))
#' }
#' 

loadParameterObject <- function(files)
{
  #A temporary env is set up to stop R errors.
  firstModel <- "Invalid model"
  for (i in 1:length(files)){
    tempEnv <- new.env();
    load(file = files[i], envir = tempEnv)
    if (i == 1){
      firstModel <- tempEnv$model
    }else{
      if (firstModel != tempEnv$model){
        stop("The models do not match between files")
      }#end of error check
    }#end of if-else
  }#end of for
  
#  browser()
  if (firstModel == "ROC"){
    parameter <- new(ROCParameter)
    parameter <- loadROCParameterObject(parameter, files)
  }else if (firstModel == "PA") {
    parameter <- new(PAParameter)
    parameter <- loadPAParameterObject(parameter, files)
  }else if (firstModel == "PANSE") {
    parameter <- new(PANSEParameter)
    parameter <- loadPANSEParameterObject(parameter, files)
  }else if (firstModel == "FONSE") {
    parameter <- new(FONSEParameter)
    parameter <- loadFONSEParameterObject(parameter, files)
  }else{
    stop("File data corrupted")
  }
  return(parameter)
}

#Sets all the common variables in the Parameter objects.
setBaseInfo <- function(parameter, files)
{
  for (i in 1:length(files)) {
    tempEnv <- new.env();
    load(file = files[i], envir = tempEnv)
    if (i == 1) {
      categories <- tempEnv$paramBase$categories
      categories.matrix <- do.call("rbind", tempEnv$paramBase$categories)
      numMixtures <- tempEnv$paramBase$numMix
      numMutationCategories <- tempEnv$paramBase$numMut
      numSelectionCategories <- tempEnv$paramBase$numSel
      mixtureAssignment <- tempEnv$paramBase$curMixAssignment
      lastIteration <- tempEnv$paramBase$lastIteration
      max <- tempEnv$paramBase$lastIteration + 1
      
      stdDevSynthesisRateTraces <- vector("list", length = numSelectionCategories)
      for (j in 1:numSelectionCategories) {
        stdDevSynthesisRateTraces[[j]] <- tempEnv$paramBase$stdDevSynthesisRateTraces[[j]][1:max]
      }
      stdDevSynthesisRateAcceptanceRateTrace <- tempEnv$paramBase$stdDevSynthesisRateAcceptRatTrace
      synthesisRateTrace <- vector("list", length = numSelectionCategories)
      for (j in 1:numSelectionCategories) {
        for (k in 1:length(tempEnv$paramBase$synthRateTrace[[j]])){
          synthesisRateTrace[[j]][[k]] <- tempEnv$paramBase$synthRateTrace[[j]][[k]][1:max]
        }
      }
      synthesisRateAcceptanceRateTrace <- tempEnv$paramBase$synthAcceptRatTrace
      mixtureAssignmentTrace <- vector("list", length = length(tempEnv$paramBase$mixAssignTrace))
      for (j in 1:length(tempEnv$paramBase$mixAssignTrace)){
        mixtureAssignmentTrace[[j]] <- tempEnv$paramBase$mixAssignTrace[[j]][1:max]
      }
      mixtureProbabilitiesTrace <- c()
      for (j in 1:numMixtures) {
        mixtureProbabilitiesTrace[[j]] <- tempEnv$paramBase$mixProbTrace[[j]][1:max]
      }
      codonSpecificAcceptanceRateTrace <- tempEnv$paramBase$codonSpecificAcceptRatTrace
    } else {
      if (sum(categories.matrix != do.call("rbind", tempEnv$paramBase$categories)) != 0){
          stop("categories is not the same between all files")
      }#end of error check

      if (numMixtures != tempEnv$paramBase$numMix){
        stop("The number of mixtures is not the same between files")
      }
      
      if (numMutationCategories != tempEnv$paramBase$numMut){
        stop("The number of mutation categories is not the same between files")
      }
      
      if (numSelectionCategories != tempEnv$paramBase$numSel){
        stop("The number of selection categories is not the same between files")
      }
      
      if (length(mixtureAssignment) != length(tempEnv$paramBase$curMixAssignment)){
        stop("The length of the mixture assignment is not the same between files. 
             Make sure the same genome is used on each run.")
      }
      
      curStdDevSynthesisRateTraces <- tempEnv$paramBase$stdDevSynthesisRateTraces
      curStdDevSynthesisRateAcceptanceRateTrace <- tempEnv$paramBase$stdDevSynthesisRateAcceptRatTrace
      curSynthesisRateTrace <- tempEnv$paramBase$synthRateTrace
      curSynthesisRateAcceptanceRateTrace <- tempEnv$paramBase$synthAcceptRatTrace
      curMixtureAssignmentTrace <- tempEnv$paramBase$mixAssignTrace
      curMixtureProbabilitiesTrace <- tempEnv$paramBase$mixProbTrace
      curCodonSpecificAcceptanceRateTrace <- tempEnv$paramBase$codonSpecificAcceptRatTrace
      
      lastIteration <- lastIteration + tempEnv$paramBase$lastIteration
      
      
      #assuming all checks have passed, time to concatenate traces
      max <- tempEnv$paramBase$lastIteration + 1
      combineTwoDimensionalTrace(stdDevSynthesisRateTraces, curStdDevSynthesisRateTraces, max)

      size <- length(curStdDevSynthesisRateAcceptanceRateTrace)
      stdDevSynthesisRateAcceptanceRateTrace <- c(stdDevSynthesisRateAcceptanceRateTrace, 
                                      curStdDevSynthesisRateAcceptanceRateTrace[2:size])

      
      combineThreeDimensionalTrace(synthesisRateTrace, curSynthesisRateTrace, max)
      size <- length(curSynthesisRateAcceptanceRateTrace)
      combineThreeDimensionalTrace(synthesisRateAcceptanceRateTrace, curSynthesisRateAcceptanceRateTrace, size)
      
      combineTwoDimensionalTrace(mixtureAssignmentTrace, curMixtureAssignmentTrace, max)
      combineTwoDimensionalTrace(mixtureProbabilitiesTrace, curMixtureProbabilitiesTrace, max)
      size <- length(curCodonSpecificAcceptanceRateTrace)
      combineTwoDimensionalTrace(codonSpecificAcceptanceRateTrace, curCodonSpecificAcceptanceRateTrace, size)
    }
  }

  parameter$setCategories(categories)
  parameter$setCategoriesForTrace()  
  parameter$numMixtures <- numMixtures
  parameter$numMutationCategories <- numMutationCategories
  parameter$numSelectionCategories <- numSelectionCategories
  parameter$setMixtureAssignment(tempEnv$paramBase$curMixAssignment) #want the last in the file sequence
  parameter$setLastIteration(lastIteration)
  
  trace <- parameter$getTraceObject()
  trace$setStdDevSynthesisRateTraces(stdDevSynthesisRateTraces)
  trace$setStdDevSynthesisRateAcceptanceRateTrace(stdDevSynthesisRateAcceptanceRateTrace)
  trace$setSynthesisRateTrace(synthesisRateTrace)
  trace$setSynthesisRateAcceptanceRateTrace(synthesisRateAcceptanceRateTrace)
  trace$setMixtureAssignmentTrace(mixtureAssignmentTrace)
  trace$setMixtureProbabilitiesTrace(mixtureProbabilitiesTrace)
  trace$setCodonSpecificAcceptanceRateTrace(codonSpecificAcceptanceRateTrace)
  
  parameter$setTraceObject(trace)
  return(parameter)
}


#Called from "loadParameterObject."
loadROCParameterObject <- function(parameter, files)
{
  parameter <- setBaseInfo(parameter, files)
  for (i in 1:length(files)){
    tempEnv <- new.env();
    load(file = files[i], envir = tempEnv)

    numMutationCategories <- tempEnv$paramBase$numMut
    numSelectionCategories <- tempEnv$paramBase$numSel
    max <- tempEnv$paramBase$lastIteration + 1
  
    if (i == 1){
      withPhi <- tempEnv$withPhi
      if (withPhi){
        phiGroups <- length(tempEnv$synthesisOffsetTrace)
        synthesisOffsetTrace <- c()
        for (j in 1:phiGroups) {
          synthesisOffsetTrace[[j]] <- tempEnv$synthesisOffsetTrace[[j]][1:max]
        }
        
        
        synthesisOffsetAcceptanceRateTrace <- tempEnv$synthesisOffsetAcceptRatTrace
        
        
        observedSynthesisNoiseTrace <- c()
        for (j in 1:phiGroups) {
          observedSynthesisNoiseTrace[[j]] <- tempEnv$observedSynthesisNoiseTrace[[j]][1:max]
        }
        #need number of phi groups, not the number of mixtures apparently.
      }else {
        synthesisOffsetTrace <- c()
        synthesisOffsetAcceptanceRateTrace <- c()
        observedSynthesisNoiseTrace <- c()
      }
      
      codonSpecificParameterTraceMut <- vector("list", length=numMutationCategories)
      for (j in 1:numMutationCategories) {
	codonSpecificParameterTraceMut[[j]] <- vector("list", length=length(tempEnv$mutationTrace[[j]]))
        for (k in 1:length(tempEnv$mutationTrace[[j]])){
          codonSpecificParameterTraceMut[[j]][[k]] <- tempEnv$mutationTrace[[j]][[k]][1:max]
          #codonSpecificParameterTraceSel[[j]][[k]] <- tempEnv$selectionTrace[[j]][[k]][1:max]
        }
      }

      codonSpecificParameterTraceSel <- vector("list", length=numSelectionCategories)
      for (j in 1:numSelectionCategories) {
	codonSpecificParameterTraceSel[[j]] <- vector("list", length=length(tempEnv$selectionTrace[[j]]))
        for (k in 1:length(tempEnv$selectionTrace[[j]])){
          #codonSpecificParameterTraceMut[[j]][[k]] <- tempEnv$mutationTrace[[j]][[k]][1:max]
          codonSpecificParameterTraceSel[[j]][[k]] <- tempEnv$selectionTrace[[j]][[k]][1:max]
        }
      }
    }else{
      curSynthesisOffsetTrace <- tempEnv$synthesisOffsetTrace
      curSynthesisOffsetAcceptanceRateTrace <- tempEnv$synthesisOffsetAcceptRatTrace
      curObservedSynthesisNoiseTrace <- tempEnv$observedSynthesisNoiseTrace
      curCodonSpecificParameterTraceMut <- tempEnv$mutationTrace
      curCodonSpecificParameterTraceSel <- tempEnv$selectionTrace
      if (withPhi != tempEnv$withPhi){
        stop("Runs do not match in concern in with.phi")
      }
      
      
      if (withPhi){
        combineTwoDimensionalTrace(synthesisOffsetTrace, curSynthesisOffsetTrace, max)
        size <- length(curSynthesisOffsetAcceptanceRateTrace)
        combineTwoDimensionalTrace(synthesisOffsetAcceptanceRateTrace, curSynthesisOffsetAcceptanceRateTrace, size)
        combineTwoDimensionalTrace(observedSynthesisNoiseTrace, curObservedSynthesisNoiseTrace, max)
      }
      
      combineThreeDimensionalTrace(codonSpecificParameterTraceMut, curCodonSpecificParameterTraceMut, max)
      combineThreeDimensionalTrace(codonSpecificParameterTraceSel, curCodonSpecificParameterTraceSel, max)
    }#end of if-else
  }#end of for loop (files)
  
  trace <- parameter$getTraceObject()
  trace$setSynthesisOffsetTrace(synthesisOffsetTrace)
  trace$setSynthesisOffsetAcceptanceRateTrace(synthesisOffsetAcceptanceRateTrace)
  trace$setObservedSynthesisNoiseTrace(observedSynthesisNoiseTrace)
  trace$setCodonSpecificParameterTrace(codonSpecificParameterTraceMut, 0)
  trace$setCodonSpecificParameterTrace(codonSpecificParameterTraceSel, 1)
  
  parameter$currentMutationParameter <- tempEnv$currentMutation
  parameter$currentSelectionParameter <- tempEnv$currentSelection
  parameter$proposedMutationParameter <- tempEnv$proposedMutation
  parameter$proposedSelectionParameter <- tempEnv$proposedSelection
  parameter$setTraceObject(trace)
  return(parameter) 
}


#Called from "loadParameterObject."
loadPAParameterObject <- function(parameter, files)
{
  parameter <- setBaseInfo(parameter, files)
  
  for (i in 1:length(files)){
    tempEnv <- new.env();
    load(file = files[i], envir = tempEnv)
  
    max <- tempEnv$paramBase$lastIteration + 1
    numMixtures <- tempEnv$paramBase$numMix
    numMutationCategories <- tempEnv$paramBase$numMut
    numSelectionCategories <- tempEnv$paramBase$numSel

    if (i == 1){
      #for future use: This may break if PA is ran with more than
      #one mixture, in this case just follow the format of the 
      #ROC CSP parameters.
      alphaTrace <- vector("list", length=numMutationCategories)
      for (j in 1:numMutationCategories) {
        for (k in 1:length(tempEnv$alphaTrace[[j]])){
          alphaTrace[[j]][[k]] <- tempEnv$alphaTrace[[j]][[k]][1:max]
        }
      }
      lambdaPrimeTrace <- vector("list", length=numSelectionCategories)
      for (j in 1:numSelectionCategories) {
        for (k in 1:length(tempEnv$lambdaPrimeTrace[[j]])){
          lambdaPrimeTrace[[j]][[k]] <- tempEnv$lambdaPrimeTrace[[j]][[k]][1:max]
        }
      }
    }else{
      
      curAlphaTrace <- tempEnv$alphaTrace
      curLambdaPrimeTrace <- tempEnv$lambdaPrimeTrace
      
      combineThreeDimensionalTrace(alphaTrace, curAlphaTrace, max)
      combineThreeDimensionalTrace(lambdaPrimeTrace, curLambdaPrimeTrace, max)
    }
  }#end of for loop (files)
  
  
  parameter$currentAlphaParameter <- tempEnv$currentAlpha
  parameter$proposedAlphaParameter <- tempEnv$proposedAlpha
  parameter$currentLambdaPrimeParameter <- tempEnv$currentLambdaPrime
  parameter$proposedLambdaPrimeParameter <- tempEnv$proposedLambdaPrime
  trace <- parameter$getTraceObject()
  trace$setCodonSpecificParameterTrace(alphaTrace, 0)
  trace$setCodonSpecificParameterTrace(lambdaPrimeTrace, 1)
  
  parameter$setTraceObject(trace)
  return(parameter) 
}

loadPANSEParameterObject <- function(parameter, files)
{
  parameter <- setBaseInfo(parameter, files)
  
  for (i in 1:length(files)){
    tempEnv <- new.env();
    load(file = files[i], envir = tempEnv)
  
    max <- tempEnv$paramBase$lastIteration + 1
    numMixtures <- tempEnv$paramBase$numMix
    numMutationCategories <- tempEnv$paramBase$numMut
    numSelectionCategories <- tempEnv$paramBase$numSel

    if (i == 1){
      #for future use: This may break if PANSE is ran with more than
      #one mixture, in this case just follow the format of the 
      #ROC CSP parameters.
      alphaTrace <- vector("list", length=numMutationCategories)
      for (j in 1:numMutationCategories) {
        for (k in 1:length(tempEnv$alphaTrace[[j]])){
          alphaTrace[[j]][[k]] <- tempEnv$alphaTrace[[j]][[k]][1:max]
        }
      }
      lambdaPrimeTrace <- vector("list", length=numSelectionCategories)
      for (j in 1:numSelectionCategories) {
        for (k in 1:length(tempEnv$lambdaPrimeTrace[[j]])){
          lambdaPrimeTrace[[j]][[k]] <- tempEnv$lambdaPrimeTrace[[j]][[k]][1:max]
        }
      }
    }else{
      
      curAlphaTrace <- tempEnv$alphaTrace
      curLambdaPrimeTrace <- tempEnv$lambdaPrimeTrace
      
      combineThreeDimensionalTrace(alphaTrace, curAlphaTrace, max)
      combineThreeDimensionalTrace(lambdaPrimeTrace, curLambdaPrimeTrace, max)
    }
  }#end of for loop (files)
  
  
  parameter$currentAlphaParameter <- tempEnv$currentAlpha
  parameter$proposedAlphaParameter <- tempEnv$proposedAlpha
  parameter$currentLambdaPrimeParameter <- tempEnv$currentLambdaPrime
  parameter$proposedLambdaPrimeParameter <- tempEnv$proposedLambdaPrime
  trace <- parameter$getTraceObject()
  trace$setCodonSpecificParameterTrace(alphaTrace, 0)
  trace$setCodonSpecificParameterTrace(lambdaPrimeTrace, 1)
  
  parameter$setTraceObject(trace)
  return(parameter) 
}


#Called from "loadParameterObject."
loadFONSEParameterObject <- function(parameter, files)
{
  parameter <- setBaseInfo(parameter, files)
  for (i in 1:length(files)){
    tempEnv <- new.env();
    load(file = files[i], envir = tempEnv)
    
    numMutationCategories <- tempEnv$paramBase$numMut
    numSelectionCategories <- tempEnv$paramBase$numSel
    max <- tempEnv$paramBase$lastIteration + 1

    if (i == 1){
      
      codonSpecificParameterTraceMut <- vector("list", length=numMutationCategories)
      for (j in 1:numMutationCategories) {
	codonSpecificParameterTraceMut[[j]] <- vector("list", length=length(tempEnv$mutationTrace[[j]]))
        for (k in 1:length(tempEnv$mutationTrace[[j]])){
          codonSpecificParameterTraceMut[[j]][[k]] <- tempEnv$mutationTrace[[j]][[k]][1:max]
          #codonSpecificParameterTraceSel[[j]][[k]] <- tempEnv$selectionTrace[[j]][[k]][1:max]
        }
      }

      codonSpecificParameterTraceSel <- vector("list", length=numSelectionCategories)
      for (j in 1:numSelectionCategories) {
	codonSpecificParameterTraceSel[[j]] <- vector("list", length=length(tempEnv$selectionTrace[[j]]))
        for (k in 1:length(tempEnv$selectionTrace[[j]])){
          #codonSpecificParameterTraceMut[[j]][[k]] <- tempEnv$mutationTrace[[j]][[k]][1:max]
          codonSpecificParameterTraceSel[[j]][[k]] <- tempEnv$selectionTrace[[j]][[k]][1:max]
        }
      }

    }else{
      curCodonSpecificParameterTraceMut <- tempEnv$mutationTrace
      curCodonSpecificParameterTraceSel <- tempEnv$selectionTrace

      
      combineThreeDimensionalTrace(codonSpecificParameterTraceMut, curCodonSpecificParameterTraceMut, max)
      combineThreeDimensionalTrace(codonSpecificParameterTraceSel, curCodonSpecificParameterTraceSel, max)
    }#end of if-else
  }#end of for loop (files)
  
  trace <- parameter$getTraceObject()
  trace$setCodonSpecificParameterTrace(codonSpecificParameterTraceMut, 0)
  trace$setCodonSpecificParameterTrace(codonSpecificParameterTraceSel, 1)
  
  parameter$currentMutationParameter <- tempEnv$currentMutation
  parameter$currentSelectionParameter <- tempEnv$currentSelection
  parameter$setTraceObject(trace)
  return(parameter)  
}


#Intended to combine 2D traces (vector of vectors) read in from C++. The first
#element of the second trace is omited since it should be the same as the 
#last value of the first trace.
combineTwoDimensionalTrace <- function(trace1, trace2, max){
  for (size in 1:length(trace1))
  {
    trace1[[size]]<- c(trace1[[size]], trace2[[size]][2:max])
  }
}


#Intended to combine 3D traces (vector of vectors of vectors) read in from C++. The first
#element of the second trace is omited since it should be the same as the 
#last value of the first trace.
combineThreeDimensionalTrace <- function(trace1, trace2, max){
  
  for (size in 1:length(trace1)){
    for (sizeTwo in 1:length(trace1[[size]])){
      trace1[[size]][[sizeTwo]] <- c(trace1[[size]][[sizeTwo]], 
                          trace2[[size]][[sizeTwo]][2:max])
    }
  }
}
