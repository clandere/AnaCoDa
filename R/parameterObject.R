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
#' @param init.w.obs.phi TRUE: initialize phi values with observed phi values 
#' (data from RNAseq, mass spectrometry, ribosome footprinting) Default is FALSE. 
#' If multiple observed phi values exist for a gene, the geometric mean of these values is used as initial phi.
#' When using this function, one should remove any genes with 
#' missing phi values, as these genes will not have an initial phi value.
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
#' \code{matrix(c(1,2,3,4,1,1,1,1), ncol=2)}).
#' In this case, the selection category is the same for every mixture. If a matrix
#' is given, and it is valid, then the mutation/selection relationship will be
#' defined by the given matrix and the keyword will be ignored. A matrix should only
#' be given in cases where the keywords would not create the desired matrix.
#' 
#' @examples 
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#' restart_file <- system.file("extdata", "restart_file.rst", package = "AnaCoDa")
#' 
#' genome <- initializeGenomeObject(file = genome_file)
#' 
#' ## initialize a new parameter object
#' sphi_init <- 1
#' numMixtures <- 1
#' geneAssignment <- rep(1, length(genome))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' 
#' ## re-initialize a parameter object from a restart file. Useful for checkpointing
#' parameter <- initializeParameterObject(init.with.restart.file = restart_file)
#' 
#' ## initialize a parameter object with a custon mixture definition matrix
#' def.matrix <- matrix(c(1,1,1,2), ncol=2)
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = c(0.5, 2), num.mixtures = 2,
#'                                        gene.assignment = geneAssignment,
#'                                        mixture.definition.matrix = def.matrix)
#' 

initializeParameterObject <- function(genome = NULL, sphi = NULL, num.mixtures = 1, 
                                      gene.assignment = NULL, initial.expression.values = NULL,
                                      model = "ROC", split.serine = TRUE, 
                                      mixture.definition = "allUnique", 
                                      mixture.definition.matrix = NULL,
                                      init.with.restart.file = NULL, mutation.prior.sd = 0.35, 
                                      init.csp.variance = 0.0025, init.sepsilon = 0.1, 
                                      init.w.obs.phi=FALSE){
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
    if (!is.null(sphi)) {
      if (length(sphi) != num.mixtures) {
        stop("sphi must be a vector of length numMixtures\n")
      }
    }
    if (!is.null(initial.expression.values)) {
      if (length(initial.expression.values) != length.Rcpp_Genome(genome)) {
        stop("initial.expression.values must have length equal to the number of genes in the Genome object\n")
      }
    }
    if (!identical(split.serine, TRUE) && !identical(split.serine, FALSE)) {
      stop("split.serine must be a boolean value\n")
    }
    if (mixture.definition != "allUnique" && mixture.definition != "mutationShared" &&
        mixture.definition != "selectionShared") {
      stop("mixture.definition must be \"allUnique\", \"mutationShared\", or \"selectionShared\". Default is \"allUnique\"\n")
    }
    if (mutation.prior.sd < 0) {
      stop("mutation.prior.sd should be positive\n")
    }
    if (init.csp.variance < 0) {
      stop("init.csp.variance should be positive\n")
    } 
    if (any(init.sepsilon < 0)) {
      stop("init.sepsilon should be positive\n")
    }
  } else {
    if (!file.exists(init.with.restart.file)) {
      stop("init.with.restart.file provided does not exist\n")
    }
  }
  
  
  
  
  if(model == "ROC"){
    if(is.null(init.with.restart.file)){
      parameter <- initializeROCParameterObject(genome, sphi, num.mixtures, 
                                                gene.assignment, initial.expression.values, split.serine, 
                                                mixture.definition, mixture.definition.matrix, 
                                                mutation.prior.sd, init.csp.variance, init.sepsilon,init.w.obs.phi)    
    }else{
      parameter <- new(ROCParameter, init.with.restart.file)
    }
  }else if(model == "FONSE"){
    if(is.null(init.with.restart.file)){
      parameter <- initializeFONSEParameterObject(genome, sphi, num.mixtures, 
                                                  gene.assignment, initial.expression.values, split.serine, 
                                                  mixture.definition, mixture.definition.matrix, init.csp.variance,init.w.obs.phi)
    }else{
      parameter <- new(FONSEParameter, init.with.restart.file)
    }
  }else if(model == "PA"){
    if(is.null(init.with.restart.file)){
      parameter <- initializePAParameterObject(genome, sphi, num.mixtures, 
                                               gene.assignment, initial.expression.values, split.serine, 
                                               mixture.definition, mixture.definition.matrix, init.csp.variance,init.w.obs.phi) 
    }else{
      parameter <- new(PAParameter, init.with.restart.file)
    }
  }else if(model == "PANSE"){
    if(is.null(init.with.restart.file)){
      parameter <- initializePANSEParameterObject(genome, sphi, num.mixtures, 
                                                  gene.assignment, initial.expression.values, split.serine, 
                                                  mixture.definition, mixture.definition.matrix, init.csp.variance,init.w.obs.phi) 
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
                                         mixture.definition.matrix = NULL, mutation_prior_sd = 0.35, init.csp.variance = 0.0025, init.sepsilon = 0.1,init.w.obs.phi=FALSE){
  
  if(is.null(mixture.definition.matrix)){ 
    # keyword constructor
    parameter <- new(ROCParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(ROCParameter, as.vector(sphi), geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  if(is.null(expressionValues) && init.w.obs.phi == F)
  {
    parameter$initializeSynthesisRateByGenome(genome, mean(sphi))
    
  } 
  else if(init.w.obs.phi == T && is.null(expressionValues))
  {
    observed.phi <- getObservedSynthesisRateSet(genome)
    if (ncol(observed.phi)-1 > 1)
    {
      observed.phi <- apply(observed.phi[,2:ncol(observed.phi)],geomMean,MARGIN = 1)
    }
    else
    {
      observed.phi <- observed.phi[,2]
    }
    parameter$initializeSynthesisRateByList(observed.phi)
  }
  else if (!is.null(expressionValues) && init.w.obs.phi == F)
  {
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  else
  {
    stop("expressionValues is not NULL and init.w.obs.phi == TRUE. Please choose only one of these options.")
  }
  
  n.obs.phi.sets <- ncol(getObservedSynthesisRateSet(genome)) - 1
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
                                        mixture.definition.matrix = NULL, init.csp.variance,init.w.obs.phi=FALSE){
  
  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(PAParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(PAParameter, as.vector(sphi), geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  # initialize expression values
  if(is.null(expressionValues) && init.w.obs.phi == F)
  {
    parameter$initializeSynthesisRateByGenome(genome, mean(sphi))
    
  } 
  else if(init.w.obs.phi == T && is.null(expressionValues))
  {
    observed.phi <- getObservedSynthesisRateSet(genome)
    if (ncol(observed.phi)-1 > 1)
    {
      observed.phi <- apply(observed.phi[,2:ncol(observed.phi)],geomMean,MARGIN = 1)
    }
    else
    {
      observed.phi <- observed.phi[,2]
    }
    parameter$initializeSynthesisRateByList(observed.phi)
  }
  else if (!is.null(expressionValues) && init.w.obs.phi == F)
  {
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  else
  {
    stop("expressionValues is not NULL and init.w.obs.phi == TRUE. Please choose only one of these options.")
  }
  
  ## TODO (Cedric): use init.csp.variance to set initial proposal width for CSP parameters
  
  return (parameter)
}

#Called from initializeParameterObject.
initializePANSEParameterObject <- function(genome, sphi, numMixtures, geneAssignment, 
                                           expressionValues = NULL, split.serine = TRUE, 
                                           mixture.definition = "allUnique", 
                                           mixture.definition.matrix = NULL, init.csp.variance,init.w.obs.phi=FALSE){
  
  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(PANSEParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(PANSEParameter, as.vector(sphi), geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  # initialize expression values
  if(is.null(expressionValues) && init.w.obs.phi == F)
  {
    parameter$initializeSynthesisRateByGenome(genome, mean(sphi))
    
  } 
  else if(init.w.obs.phi == T && is.null(expressionValues))
  {
    observed.phi <- getObservedSynthesisRateSet(genome)
    if (ncol(observed.phi)-1 > 1)
    {
      observed.phi <- apply(observed.phi[,2:ncol(observed.phi)],geomMean,MARGIN = 1)
    }
    else
    {
      observed.phi <- observed.phi[,2]
    }
    parameter$initializeSynthesisRateByList(observed.phi)
  }
  else if (!is.null(expressionValues) && init.w.obs.phi == F)
  {
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  else
  {
    stop("expressionValues is not NULL and init.w.obs.phi == TRUE. Please choose only one of these options.")
  }
  
  return (parameter)
}

#Called from initializeParameterObject.
initializeFONSEParameterObject <- function(genome, sphi, numMixtures, 
                                           geneAssignment, expressionValues = NULL, split.serine = TRUE,
                                           mixture.definition = "allUnique", 
                                           mixture.definition.matrix = NULL, init.csp.variance,init.w.obs.phi=FALSE){
  
  # create Parameter object
  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(FONSEParameter, as.vector(sphi), numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(FONSEParameter, as.vector(sphi), geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  # initialize expression values
  if(is.null(expressionValues) && init.w.obs.phi == F)
  {
    parameter$initializeSynthesisRateByGenome(genome, mean(sphi))
    
  } 
  else if(init.w.obs.phi == T && is.null(expressionValues))
  {
    observed.phi <- getObservedSynthesisRateSet(genome)
    if (ncol(observed.phi)-1 > 1)
    {
      observed.phi <- apply(observed.phi[,2:ncol(observed.phi)],geomMean,MARGIN = 1)
    }
    else
    {
      observed.phi <- observed.phi[,2]
    }
    parameter$initializeSynthesisRateByList(observed.phi)
  }
  else if (!is.null(expressionValues) && init.w.obs.phi == F)
  {
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  else
  {
    stop("expressionValues is not NULL and init.w.obs.phi == TRUE. Please choose only one of these options.")
  }
  parameter <- initializeCovarianceMatrices(parameter, genome, numMixtures, geneAssignment, init.csp.variance)
  
  return(parameter)
}

 
# getCSPEstimates <- function(parameter, filename=NULL, CSP="Mutation", mixture = 1, samples = 10, optref = FALSE){
#   Amino_Acid <- c()
#   Value <- c()
#   Codon <- c()
#   quantile_list <- vector("list")
#   
#   if (class(parameter) == "Rcpp_ROCParameter" || class(parameter) == "Rcpp_FONSEParameter"){
#     names.aa <- aminoAcids()
#     
#     for(aa in names.aa){
#       if(aa == "M" || aa == "W" || aa == "X") next
#       codons <- AAToCodon(aa, T)
#       
#       for(i in 1:length(codons)){
#         Amino_Acid <- c(Amino_Acid, aa)
#         Codon <- c(Codon, codons[i])
#         
#         if(CSP == "Mutation"){
#           Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, TRUE))
#           quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, codons[i], 0, c(0.025, 0.975), TRUE))
#         }
#         else if(CSP == "Selection"){
#           Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, TRUE))
#           quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, codons[i], 1, c(0.025, 0.975), TRUE))
#         }
#         else {
#           stop("Unknown parameter type given with argument: CSP")
#         }
#       }
#     }
#   }
#   else if (class(parameter) == "Rcpp_PAParameter"){
#     groupList <- parameter$getGroupList()
#     
#     for(i in 1:length(groupList)){
#       aa <- codonToAA(groupList[i])
#       Codon <- c(Codon, groupList[i])
#       Amino_Acid <- c(Amino_Acid, aa)
#       
#       if(CSP == "Alpha"){
#         Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, groupList[i], 0, FALSE))
#         quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, groupList[i], 0, c(0.025, 0.975), FALSE))
#       }
#       else if(CSP == "Lambda Prime"){
#         Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, groupList[i], 1, FALSE))
#         quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, groupList[i], 1, c(0.025, 0.975), FALSE))
#       }
#       else {
#         stop("Unknown parameter type given with argument: CSP")
#       }
#     }
#   }
#   else if (class(parameter) == "Rcpp_PANSEParameter"){
#     groupList <- parameter$getGroupList()
#     
#     for(i in 1:length(groupList)){
#       aa <- codonToAA(groupList[i])
#       Codon <- c(Codon, groupList[i])
#       Amino_Acid <- c(Amino_Acid, aa)
#       
#       if(CSP == "Alpha"){
#         Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, groupList[i], 0, FALSE))
#         quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, groupList[i], 0, c(0.025, 0.975), FALSE))
#       }
#       else if(CSP == "Lambda Prime"){
#         Value <- c(Value, parameter$getCodonSpecificPosteriorMean(mixture, samples, groupList[i], 1, FALSE))
#         quantile_list <- c(quantile_list, parameter$getCodonSpecificQuantile(mixture, samples, groupList[i], 1, c(0.025, 0.975), FALSE))
#       }
#       else {
#         stop("Unknown parameter type given with argument: CSP")
#       }
#     }
#   }
#   else{
#     stop("Unknown object provided with argument: parameter")
#   }
#   
#   quantile_list <- matrix(unlist(quantile_list), nrow = 2)
#   data <- data.frame(Amino_Acid, Codon, Value, Lower=quantile_list[1,], Upper=quantile_list[2,])
#   colnames(data) <- c("AA", "Codon", "Posterior", "0.025%", "0.975%")
#   data <- add.ref.codon(data, optref)
#   
#   if(is.null(filename))
#   {
#     return(data)
#   }else {
#     write.csv(data, file = filename, row.names = FALSE, quote=FALSE)
#   }
# }
# 
# add.ref.codon <- function(csp, opt.ref = FALSE)
# {
#   fill.in <- data.frame(Codon=as.character(AnaCoDa::codons()))
#   csp <- merge(x = csp, y = fill.in, by = "Codon", all = T)
#   for(i in 1:64)
#     csp$AA[i] <- AnaCoDa::codonToAA(as.character(csp$Codon[i]))
#   csp[is.na(csp)] <- 0
#   csp <- csp[order(csp$AA), ]
#   
#   if(opt.ref)
#   {
#     aas <- aminoAcids()
#     n.aa <- length(aas)
#     for(j in 1:n.aa)
#     {
#       aa <- aas[j]
#       if(aa == "W" || aa == "M" || aa == "X") next
#       aa.pos <- which(csp$AA == aa)
#       csp[aa.pos, 3] <- csp[aa.pos, 3] - min(csp[aa.pos, 3])
#       csp[aa.pos, 4] <- csp[aa.pos, 4] - min(csp[aa.pos, 4])
#       csp[aa.pos, 5] <- csp[aa.pos, 5] - min(csp[aa.pos, 5])
#     }
#   }
#   return(csp)
# }



#' Calculates the marginal log-likelihood for a set of parameters
#'
#' @param parameter An object created with \code{initializeParameterObject}. 
#'  
#' @param mcmc An object created with \code{initializeMCMCObject} 
#'  
#' @param mixture determines for which mixture the marginal log-likelihood should be calculated
#'    
#' @param n.samples How many samples should be used for the calculation 
#'  
#' @param divisor A value > 1 in order to scale down the tails of the importance distribution
#'  
#' @return This function returns the model object created. 
#'  
#' @description initializes the model object. 
#' 
#' @details calculateMarginalLogLikelihood Calculate marginal log-likelihood for
#' calculation of the Bayes factor using a generalized harmonic mean estimator of the 
#' marginal likelihood. See Gronau et al. (2017) for details
#'
#' @examples 
#' \dontrun{
#' # Calculate the log-marginal likelihood
#' parameter <- loadParameterObject("parameter.Rda")
#' mcmc <- loadMCMCObject("mcmc.Rda")
#' calculateMarginalLogLikelihood(parameter, mcmc, mixture = 1, 
#' samples = 500, scaling = 1.5)
#'
#' # Calculate the Bayes factor for two models
#' parameter1 <- loadParameterObject("parameter1.Rda")
#' parameter2 <- loadParameterObject("parameter2.Rda")
#' mcmc1 <- loadMCMCObject("mcmc1.Rda")
#' mcmc2 <- loadMCMCObject("mcmc2.Rda")
#' mll1 <- calculateMarginalLogLikelihood(parameter1, mcmc1, mixture = 1, 
#' samples = 500, scaling = 1.5)
#' mll2 <- calculateMarginalLogLikelihood(parameter2, mcmc2, mixture = 1, 
#' samples = 500, scaling = 1.5)
#' cat("Bayes factor: ", mll1 - mll2, "\n")
#' }
#'  
calculateMarginalLogLikelihood <- function(parameter, mcmc, mixture, n.samples, divisor)
{  
  if(divisor < 1) stop("Generalized Harmonic Mean Estimation of Marginal Likelihood requires importance sampling distribution variance divisor be greater than 1")

  ## Collect information from AnaCoDa objects
  trace <- parameter$getTraceObject()
  ## This should be the posterior instead of the log_posterior but this causes an overflow, find fix!!!
  log_posterior <- mcmc$getLogPosteriorTrace()
  log_posterior <- log_posterior[(length(log_posterior) - n.samples+1):(length(log_posterior))]
  
  ### HANDLE CODON SPECIFIC PARAMETERS
  log_imp_dens_sample <- rep(0, n.samples)
  for(ptype in 0:1) # for all parameter types (mutation/selection parameters)
  {
    for(aa in AnaCoDa::aminoAcids()) # for all amino acids
    {
      if(aa == "M" || aa == "W" || aa == "X") next # ignore amino acids with only one codon or stop codons 
      codons <- AnaCoDa::AAToCodon(aa, focal = T)
      ## get covariance matrix and mean of importance distribution
      sample_mat <- matrix(NA, ncol = length(codons), nrow = n.samples)
      mean_vals <- rep(NA, length(codons))
      for(i in 1:length(codons)) # for all codons
      {
        vec <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], ptype, TRUE)
        vec <- vec[(length(vec) - n.samples+1):(length(vec))]
        sample_mat[,i] <- vec
        mean_vals[i] <- mean(vec)
      }
      ## scale/shrinked covariance matrix
      cov_mat <- cov(sample_mat) / divisor
      
      for(i in 1:n.samples)
      {
        ## calculate importance density for collected samples
        ## mikeg: We should double check with Russ that it is okay to use the full covariance matrix of the sample (even though we have no prior on the cov structure) when constructing the importance density function.
        ## It seems logical to do so
        log_imp_dens_aa <- dmvnorm(x = sample_mat[i,], mean = mean_vals, sigma = cov_mat, log = TRUE)
        log_imp_dens_sample[i] = log_imp_dens_sample[i] + log_imp_dens_aa
      }
    }
  }
  
  ## HANDLE GENE SPECIFIC PARAMETERS
  
  # phi values are stored on natural scale.
  synt_trace <- trace$getSynthesisRateTrace()[[mixture]]
  n_genes = length(synt_trace);
  
  sd_vals <- rep(NA, n_genes)
  mean_vals <- rep(NA, n_genes)
  for(i in 1:n_genes) ## i is indexing across genes
  {
    vec <- synt_trace[[i]]
    
    vec <- vec[(length(vec) - n.samples+1):(length(vec))]
    
    sd_vals[i] <- sd(vec)
    mean_vals[i] <- mean(vec)
    log_mean_vals = log(mean_vals) - 0.5 * log(1+(sd_vals^2/mean_vals^2))
    log_sd_vals = sqrt(log(1+(sd_vals^2/mean_vals^2)))
    
    ## Calculate vector of importance density for entire \phi trace of gene.
    log_imp_dens_phi <- dlnorm(x = vec, meanlog = log_mean_vals[i], sdlog = log_sd_vals[i]/divisor, log = TRUE)
    
    
    ## update importance density function vector of sample with current gene;
    log_imp_dens_sample = log_imp_dens_sample + log_imp_dens_phi
    
  } ## end synth_trace loop
  
  ## Scale importance density for each sample by its posterior probability (on log scale)
  log_imp_dens_over_posterior = log_imp_dens_sample - log_posterior
  
  ## now scale by max term to facilitate summation
  max_log_term = max(log_imp_dens_over_posterior)
  ## Y = X - max_X
  offset_log_imp_dens_over_posterior <- log_imp_dens_over_posterior - max_log_term
  ## Z = sum(exp(vec(Y)))
  offset_sum_imp_dens_over_posterior <- sum(exp(offset_log_imp_dens_over_posterior))/n.samples
  log_sum_imp_dens_over_posterior <- log(offset_sum_imp_dens_over_posterior) + max_log_term
  ## ln(ML) = ln(n) -(Z + max_X)
  log_marg_lik <- log(n.samples) - log_sum_imp_dens_over_posterior
  ##marg_lik = 1.0/(log_inv_marg_lik/n.samples) # equation 9
  return(log_marg_lik)
}



#' Find and return list of optimal codons
#' 
#' \code{findOptimalCodon} extracrs the optimal codon for each amino acid.
#' 
#' @param csp a \code{data.frame} as returned by \code{getCSPEstimates}.
#'
#' @return A named list with with optimal codons for each amino acid.
#'
#' @examples 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- 1
#' numMixtures <- 1
#' geneAssignment <- rep(1, length(genome))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' model <- initializeModelObject(parameter = parameter, model = "ROC")
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, 
#'                              adaptive.width=adaptiveWidth, est.expression=TRUE, 
#'                              est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#' divergence.iteration <- 10
#' \dontrun{
#' runMCMC(mcmc = mcmc, genome = genome, model = model, 
#'         ncores = 4, divergence.iteration = divergence.iteration)
#' 
#' csp_mat <- getCSPEstimates(parameter, CSP="Selection")
#' opt_codons <- findOptimalCodon(csp_mat)
#' }

findOptimalCodon <- function(csp)
{
  aas <- aminoAcids()
  n.aa <- length(aas)
  result <- vector("list", length(aas))
  names(result) <- aas
  for(j in 1:n.aa)
  {
    aa <- aas[j]
    if(aa == "W" || aa == "M" || aa == "X") next
    aa.pos <- which(csp$AA == aa)
    opt.codon.pos <- which(csp[aa.pos, 3] == min(csp[aa.pos, 3]))
    result[[j]] <- csp$Codon[aa.pos[opt.codon.pos]]
  }
  return(result)
}


#' Return Codon Specific Paramters (or write to csv) estimates as data.frame
#' 
#' @param parameter parameter an object created by \code{initializeParameterObject}.
#' 
#' @param filename Posterior estimates will be written to file (format: csv). Filename will be in the format <parameter_name>_<filename>.csv.
#' 
#' @param mixture estimates for which mixture should be returned
#' 
#' @param samples The number of samples used for the posterior estimates.
#'
#' @param relative.to.optimal.codon Boolean determining if parameters should be relative to the preferred codon or the alphabetically last codon (Default=TRUE). Only applies to ROC and FONSE models 
#'
#' @param report.original.ref Include the original reference codon (Default = TRUE). Note this is only included for the purposes of simulations, which expect the input parameter file to be in a specific format. Later version of AnaCoDa will remove this. 

#' @return returns a list data.frame with the posterior estimates of the models 
#' codon specific parameters or writes it directly to a csv file if \code{filename} is specified
#' 
#' @description \code{getCSPEstimates} returns the codon specific
#' parameter estimates for a given parameter and mixture or write it to a csv file.
#'
#' @examples  
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- c(1,1)
#' numMixtures <- 2
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' model <- initializeModelObject(parameter = parameter, model = "ROC")
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, 
#'                              adaptive.width=adaptiveWidth, est.expression=TRUE, 
#'                              est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#' divergence.iteration <- 10
#' \dontrun{
#' runMCMC(mcmc = mcmc, genome = genome, model = model, 
#'         ncores = 4, divergence.iteration = divergence.iteration)
#' 
#' ## return estimates for codon specific parameters
#' csp_mat <- getCSPEstimates(parameter)
#' 
#' # write the result directly to the filesystem as a csv file. No values are returned
#' getCSPEstimates(parameter, filename=file.path(tempdir(), "test.csv"))
#' 
#' }

getCSPEstimates <- function(parameter, filename=NULL, mixture = 1, samples = 10, relative.to.optimal.codon=T, report.original.ref = T)
{
  model.conditions <- checkModel(parameter)
  model.uses.ref.codon <- model.conditions$model.uses.ref.codon
  names.aa <- model.conditions$aa
  codons <- model.conditions$codons
  parameter.names <- model.conditions$parameter.names
  
  ## Creates empty vector of 0 for initial dataframes
  init <- rep(0.0,length(codons))
  
  param.1<- data.frame(Codon=codons,AA=names.aa,Posterior=init,Lower.quant=init, Upper.quant=init,stringsAsFactors = F,row.names = codons)
  param.2 <- data.frame(Codon=codons,AA=names.aa,Posterior=init,Lower.quant=init, Upper.quant=init,stringsAsFactors = F,row.names=codons)
  if (model.uses.ref.codon)
  {
    codons <- codons[which(codons %in% unlist(lapply(X = names.aa,FUN = AAToCodon,T)))]
  }
  ## Get parameter estimate for each codon
  for (codon in codons)
  {
    param.1[codon,"Posterior"] <- parameter$getCodonSpecificPosteriorMean(mixtureElement=mixture,samples=samples,codon=codon,paramType=0,withoutReference=model.uses.ref.codon)
    param.2[codon,"Posterior"] <- parameter$getCodonSpecificPosteriorMean(mixtureElement=mixture,samples=samples,codon=codon,paramType=1,withoutReference=model.uses.ref.codon)
    param.1[codon,c("Lower.quant","Upper.quant")] <- parameter$getCodonSpecificQuantile(mixtureElement=mixture, samples=samples,codon=codon,paramType=0, probs=c(0.025, 0.975),withoutReference=model.uses.ref.codon)
    param.2[codon,c("Lower.quant","Upper.quant")]  <- parameter$getCodonSpecificQuantile(mixtureElement=mixture, samples=samples,codon=codon,paramType=1, probs=c(0.025, 0.975),withoutReference=model.uses.ref.codon)
  }
  colnames(param.1) <- c("Codon", "AA", "Posterior", "0.025%", "0.975%")
  colnames(param.2) <- c("Codon", "AA", "Posterior", "0.025%", "0.975%")
  
  ## Only called if model actually uses reference codon
  if(relative.to.optimal.codon && model.uses.ref.codon)
  {
    csp.param <- optimalAsReference(param.1,param.2,parameter.names,report.original.ref)
  } else if (relative.to.optimal.codon == F || model.uses.ref.codon == F ){
    ## This is just in case the user wants to exclude the original reference codon
    ## TO DO: update C++ function which might expect certain format for the input CSP file parameters
    if (model.uses.ref.codon && !report.original.ref)
    {
      param.1 <- param.1[-which(param.1[,"Posterior"]==0),]
      param.2 <- param.2[-which(param.2[,"Posterior"]==0),]
    }
    csp.param <- vector(mode="list",length=2)
    names(csp.param) <- parameter.names
    csp.param[[parameter.names[1]]] <- param.1[,c("AA", "Codon", "Posterior", "0.025%", "0.975%")]
    csp.param[[parameter.names[2]]] <- param.2[,c("AA", "Codon", "Posterior", "0.025%", "0.975%")]
  }
  if(is.null(filename))
  {
    return(csp.param)
  }else {
    write.csv(csp.param[[parameter.names[1]]], file = paste0(filename,"_",parameter.names[1],".csv"), row.names = FALSE, quote=FALSE)
    write.csv(csp.param[[parameter.names[2]]], file = paste0(filename,"_",parameter.names[2],".csv"), row.names = FALSE, quote=FALSE)
  }
}

## NOT EXPOSED
optimalAsReference <- function(param.1,param.2,parameter.names,report.original.ref)
{
  updated.param.1 <- data.frame()
  updated.param.2 <- data.frame()
  aa <- unique(param.2[,"AA"])
  for (a in aa)
  {
    codons <- AAToCodon(a)
    ## Create temporary data frames for modifying values
    tmp.1 <- param.1[codons,] ## "Mutation" parameter
    tmp.2 <- param.2[codons,] ## "Selection" parameter
    current.reference.row <- which(tmp.2[,"Posterior"]==0)
    optimal.parameter.value <- min(tmp.2[,"Posterior"])
    ## No reason to do anything if optimal value is 0
    if (optimal.parameter.value != 0.0)
    {
      tmp.2[,c("Posterior","0.025%","0.975%")] <- tmp.2[,c("Posterior","0.025%","0.975%")] - optimal.parameter.value
      ##Get row of the optimal codon, which should be 0
      optimal.codon.row <- which(tmp.2[,"Posterior"]==0.0)
      tmp.2[current.reference.row,c("0.025%","0.975%")] <- tmp.2[optimal.codon.row,c("0.025%","0.975%")] + tmp.2[current.reference.row,"Posterior"]
      ## Can now change optimal codon values to 0.0
      tmp.2[optimal.codon.row,c("Posterior","0.025%","0.975%")] <- 0.0
      
      ## Find corresponding reference value for other parameter
      optimal.parameter.value <- tmp.1[optimal.codon.row,"Posterior"]
      tmp.1[,c("Posterior","0.025%","0.975%")] <- tmp.1[,c("Posterior","0.025%","0.975%")] - optimal.parameter.value
      tmp.1[current.reference.row,c("0.025%","0.975%")] <- tmp.1[optimal.codon.row,c("0.025%","0.975%")] + tmp.1[current.reference.row,"Posterior"]
      tmp.1[optimal.codon.row,c("Posterior","0.025%","0.975%")] <- 0.0
    }
    if (!report.original.ref)
    {
      tmp.1 <- tmp.1[-current.reference.row,]
      tmp.2 <- tmp.2[-current.reference.row,]
    }
    updated.param.1 <- rbind(updated.param.1,tmp.1)
    updated.param.2 <- rbind(updated.param.2,tmp.2)
  }
  csp.param <- vector(mode="list",length=2)
  names(csp.param) <- parameter.names
  csp.param[[parameter.names[1]]] <- updated.param.1[,c("AA", "Codon", "Posterior", "0.025%", "0.975%")]
  csp.param[[parameter.names[2]]] <- updated.param.2[,c("AA", "Codon", "Posterior", "0.025%", "0.975%")]
  return(csp.param)
}

## NOT EXPOSED
checkModel <- function(parameter)
{
  class.type = class(parameter)
  if(class(parameter)=="Rcpp_ROCParameter" || class(parameter)=="Rcpp_FONSEParameter")
  {
    model.uses.ref.codon <- TRUE
    names.aa <- parameter$getGroupList()
    codons <- unlist(lapply(names.aa,AAToCodon))
    aa <- unlist(lapply(codons,codonToAA))
    parameter.names <- c("Mutation","Selection")
    
  } else
  {
    model.uses.ref.codon <- FALSE
    codons <- parameter$getGroupList()
    aa <- unlist(lapply(codons,codonToAA))
    parameter.names <- c("Alpha","Lambda Prime")
  }
  return(list(aa=aa,codons=codons,model.uses.ref.codon=model.uses.ref.codon,parameter.names=parameter.names))
}


#' Calculate Selection coefficients
#' 
#' \code{getSelectionCoefficients} calculates the selection coefficient of each codon in each gene.
#' 
#' @param genome A genome object initialized with 
#' \code{\link{initializeGenomeObject}} to add observed expression data.
#' 
#' @param parameter an object created by \code{initializeParameterObject}.
#' 
#' @param samples The number of samples used for the posterior estimates.
#' 
#' @return A matrix with selection coefficients.
#' 
#' @examples 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- 1
#' numMixtures <- 1
#' geneAssignment <- rep(1, length(genome))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' model <- initializeModelObject(parameter = parameter, model = "ROC")
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, 
#'                              adaptive.width=adaptiveWidth, est.expression=TRUE, 
#'                              est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#' divergence.iteration <- 10
#' \dontrun{
#' runMCMC(mcmc = mcmc, genome = genome, model = model, 
#'         ncores = 4, divergence.iteration = divergence.iteration)
#' 
#' ## return estimates for selection coefficients s for each codon in each gene
#' selection.coefficients <- getSelectionCoefficients(genome = genome, 
#'                                                    parameter = parameter, samples = 1000)
#' }
#' 
getSelectionCoefficients <- function(genome, parameter, samples = 100)
{
  sel.coef <- parameter$calculateSelectionCoefficients(samples)
  grouplist <- parameter$getGroupList()
  codon.names <- NULL
  if(class(parameter) == "Rcpp_ROCParameter" || class(parameter) == "Rcpp_FONSEParameter")
  {
    for(aa in grouplist)
      codon.names <- c(codon.names, AAToCodon(aa))
    
    sel.coef <- sel.coef[, -c(60, 61)] # The matrix is to large as it could store M and W which is not used here.
  }else{
    codon.names <- grouplist
  }
  
  gene.names <- getNames(genome)
  
  colnames(sel.coef) <- codon.names
  rownames(sel.coef) <- gene.names
  return(sel.coef)  
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

#' extracts an object of traces from a parameter object.
#'
#' @param parameter A Parameter object that corresponds to one of the model types.
#'
#' @return trace Returns an object of type Trace extracted from the given parameter object
#'
#' @examples  
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- c(1,1)
#' numMixtures <- 2
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' 
#' trace <- getTrace(parameter) # empty trace object since no MCMC was perfomed
#' 
getTrace <- function(parameter){
  return(parameter$getTraceObject())
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
    
    codonCounts <- getCodonCountsForAA(aa, genome) # ignore column with gene ids
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
  
  #codonCounts <- getCodonCountsForAA(aa, genome)
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
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- c(1,1)
#' numMixtures <- 2
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' model <- initializeModelObject(parameter = parameter, model = "ROC")
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width=adaptiveWidth, 
#'                              est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#' divergence.iteration <- 10
#' \dontrun{
#' runMCMC(mcmc = mcmc, genome = genome, model = model, 
#'         ncores = 4, divergence.iteration = divergence.iteration)
#' 
#' # get the mixture assignment for all genes
#' mixAssign <- getMixtureAssignmentEstimate(parameter = parameter, 
#'                                           gene.index = 1:length(genome), samples = 1000)
#' 
#' # get the mixture assignment for a subsample
#' mixAssign <- getMixtureAssignmentEstimate(parameter = parameter, 
#'                                           gene.index = 5:100, samples = 1000)
#' # or
#' mixAssign <- getMixtureAssignmentEstimate(parameter = parameter, 
#'                                           gene.index = c(10, 30:50, 3, 90), samples = 1000)
#' }
#' 
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
#' @param samples number of samples for the posterior estimate
#'
#' @param quantiles vector of quantiles, (default: c(0.025, 0.975))
#' 
#' @return returns a vector with the mixture assignment of each gene corresbonding to \code{gene.index} in the same order as the genome. 
#'
#' @description Posterior estimates for the phi value of specified genes
#' 
#' @details The returned vector is unnamed as gene ids are only stored in the \code{genome} object, 
#' but the \code{gene.index} vector can be used to match the assignment to the genome.
#' 
#' @examples  
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- c(1,1)
#' numMixtures <- 2
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' model <- initializeModelObject(parameter = parameter, model = "ROC")
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, 
#'                              adaptive.width=adaptiveWidth, est.expression=TRUE,
#'                              est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#' divergence.iteration <- 10
#' \dontrun{
#' runMCMC(mcmc = mcmc, genome = genome, model = model, 
#'         ncores = 4, divergence.iteration = divergence.iteration)
#' 
#' # get the estimated expression values for all genes based on the mixture 
#' # they are assigned to at each step
#' estimatedExpression <- getExpressionEstimates(parameter, 1:length(genome), 1000)
#' }
#' 
getExpressionEstimates <- function(parameter, gene.index, samples, quantiles=c(0.025, 0.975))
{
  expressionValues <- unlist(lapply(gene.index, function(geneIndex){ 
    parameter$getSynthesisRatePosteriorMeanForGene(samples, geneIndex, FALSE) 
  }))
  
  expressionValuesLog <- unlist(lapply(gene.index, function(geneIndex){ 
    parameter$getSynthesisRatePosteriorMeanForGene(samples, geneIndex, TRUE) 
  }))
  
  expressionStdErr <- sqrt(unlist(lapply(gene.index, function(geneIndex){ 
    parameter$getSynthesisRateVarianceForGene(samples, geneIndex, TRUE, FALSE) 
  }))) / samples
  
  expressionStdErrLog <- sqrt(unlist(lapply(gene.index, function(geneIndex){ 
    parameter$getSynthesisRateVarianceForGene(samples, geneIndex, TRUE, TRUE) 
  }))) / samples
  
  expressionQuantile <- lapply(gene.index, function(geneIndex){ 
    parameter$getExpressionQuantile(samples, geneIndex, quantiles, FALSE) 
  })
  expressionQuantile <- do.call(rbind, expressionQuantile)
  
  expressionQuantileLog <- lapply(gene.index, function(geneIndex){ 
    parameter$getExpressionQuantile(samples, geneIndex, quantiles, TRUE) 
  })
  expressionQuantileLog <- do.call(rbind, expressionQuantileLog)
  
  expr.mat <- cbind(expressionValues, expressionValuesLog, expressionStdErr, expressionStdErrLog, expressionQuantile, expressionQuantileLog)
  colnames(expr.mat) <- c("PHI", "log10.PHI", "Std.Error", "log10.Std.Error", quantiles, paste("log10.", quantiles, sep=""))
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
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- c(1,1)
#' numMixtures <- 2
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' 
#' ## writing an empty parameter object as the runMCMC routine was not called yet
#' writeParameterObject(parameter = parameter, file = file.path(tempdir(), "file.Rda"))
#' 
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
  grouplist <- parameter$getGroupList()
  
  
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
                  lastIteration = lastIteration,
                  grouplist = grouplist
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
#' parameter <- loadParameterObject("parameter.Rda")
#' 
#' # load and concatenate multiple parameter object
#' parameter <- loadParameterObject(c("parameter1.Rda", "parameter2.Rda"))
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
      grouplist <- tempEnv$paramBase$grouplist
      
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
      
      if(length(grouplist) != length(tempEnv$paramBase$grouplist)){
        stop("Number of Amino Acids/Codons is not the same between files.")	
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
  parameter$setGroupList(grouplist)
  
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



#' Take the geometric mean of a vector
#'  
#' @param x A vector of numerical .
#' 
#' @param rm.invalid Boolean value for handling 0, negative, or NA values in the vector. Default is TRUE and will not
#' include these values in the calculation. If FALSE, these values will be replaced by the value give to \code{default} and will
#' be included in the calculation.
#' 
#' @param default Numerical value that serves as the value to replace 0, negative, or NA values in the calculation when rm.invalid is FALSE.
#' Default is 1e-5.
#' 
#' @return Returns the geometric mean of a vector.
#' 
#' @description \code{geomMean} will calculate the geometric mean of a list of numerical values.
#' 
#' @details This function is a special version of the geometric mean specifically for AnaCoda.
#' Most models in Anacoda assume a log normal distribution for phi values, thus all values in \code{x} are expectd to be positive.
#' geomMean returns the geometric mean of a vector and can handle 0, negative, or NA values. 
#' 
#' @examples 
#' x <- c(1, 2, 3, 4)
#' geomMean(x)
#' 
#' y<- c(1, NA, 3, 4, 0, -1)
#' # Only take the mean of non-Na values greater than 0
#' geomMean(y)
#' 
#' # Replace values <= 0 or NAs with a default value 0.001 and then take the mean
#' geomMean(y, rm.invalid = FALSE, default = 0.001)
#' 
geomMean <- function(x, rm.invalid = TRUE, default = 1e-5)
{
  if(!rm.invalid)
  {
    x[x <= 0 | is.na(x)] <- default
  } else{
    x <- x[which(x > 0 & !is.na(x))]
  }
  total <- prod(x) ^ (1/length(x))
  return(total)
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




