#'  Model Initialization 
#'
#' @param parameter An object created with \code{initializeParameterObject}. 
#'  
#' @param model A string containing the model to run (ROC, FONSE, or PA), has to match parameter object. 
#'  
#' @param with.phi (ROC only) A boolean that determines whether or not to include empirical
#'    phi values (expression rates) for the calculations. 
#'    
#' @param fix.observation.noise (ROC only) Allows to fix the noise in the observed expression dataset to the initial condition.
#'	The initial condition for the observed expression noise can be set in the parameter object. 
#'  
#' @param rfp.count.column (PA and PANSE only) A number representing the RFP count column to use.
#'  
#' @return This function returns the model object created. 
#'  
#' @description initializes the model object. 
#' 
#' @details initializeModelObject initializes a model. The type of model is determined based on the string passed to the \code{model} argument.
#'  The Parameter object has to match the model that is initialized. E.g. to initialize a ROC model, 
#'  it is required that a ROC parameter object is passed to the function.
#'
#' @examples 
#' 
#' #initializing a model object
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#' expression_file <- system.file("extdata", "expression.csv", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file, 
#'                                  observed.expression.file = expression_file)
#' sphi_init <- c(1,1)
#' numMixtures <- 2
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' 
#' # initializing a model object assuming we have observed expression (phi) 
#' # values stored in the genome object.
#' initializeModelObject(parameter = parameter, model = "ROC", with.phi = TRUE)
#'
#' # initializing a model object ignoring observed expression (phi) 
#' # values stored in the genome object.
#' initializeModelObject(parameter = parameter, model = "ROC", with.phi = FALSE)
#' 
initializeModelObject <- function(parameter, model = "ROC", with.phi = FALSE, fix.observation.noise = FALSE, rfp.count.column = 1) {
  if(model == "ROC") {
    c.model <- new(ROCModel, with.phi, fix.observation.noise)
  } else if (model == "FONSE") {
    c.model = new(FONSEModel)
  } else if (model == "PA") {
    c.model <- new(PAModel, rfp.count.column)
  } else if (model == "PANSE") {
    c.model <- new(PANSEModel, rfp.count.column)
  } else {
    stop("Unknown model.")
  }
  c.model$setParameter(parameter)
  return(c.model)
}
