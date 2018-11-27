library(testthat)
library(AnaCoDa)

context("MCMC with ROC")

# This file currently checks the logPosterior recorded at iteration 100, between an old, hard-coded test and a current Unit Test.
# Two tests are run: one with Phi, one without Phi. The existence of the relevant input files is also checked.

# Possible implementation change: take the logPosterior value and hard code it here, and compare via
# mcmc$getLogPosteriorTrace(), which returns a vector. Get the average of these values
# and compare it with the hard-coded average of logPosteriorTrace.

# In R, file.path is faster than paste
fileName = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
expressionFile = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR_phi_withPhiSet.csv")

# Ensure the input files exist.
test_that("file exists: simulatedAllUniqueR.fasta", {
  expect_equal(file.exists(fileName), T)
})
test_that("file exists: simulatedAllUniqueR_phi_withPhiSet.csv", {
  expect_equal(file.exists(expressionFile), T)
})

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "allUnique"

samples <- 10
thinning <- 10
adaptiveWidth <- 10
divergence.iteration <- 0

mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)

### With Phi
set.seed(446141)

genome <- initializeGenomeObject(file = fileName, observed.expression.file = expressionFile, match.expression.by.id=FALSE)

geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

model <- initializeModelObject(parameter, "ROC", with.phi = TRUE) 

outFile = file.path("UnitTestingOut", "testMCMCROCLogPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()

test_that("identical MCMC-ROC input with Phi, same log posterior", {
  knownLogPosterior <- -828047
  testLogPosterior <- round(mcmc$getLogPosteriorTrace()[10])
  expect_equal(knownLogPosterior, testLogPosterior)
})

### Without Phi
set.seed(446141)

genome <- initializeGenomeObject(file = fileName) 

geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

model <- initializeModelObject(parameter, "ROC", with.phi = FALSE) 

outFile = file.path("UnitTestingOut", "testMCMCROCLogWithoutPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()

test_that("identical MCMC-ROC input without Phi, same log posterior", {
  knownLogPosterior <- -831329
  testLogPosterior <- round(mcmc$getLogPosteriorTrace()[10])
  expect_equal(knownLogPosterior, testLogPosterior)
})

