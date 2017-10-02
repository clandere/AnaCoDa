[![Build Status](https://travis-ci.org/clandere/RibModelFramework.svg)](https://travis-ci.org/clandere/RibModelFramework)

# AnaCoDa

* AnaCoDa is a collection of codon models. 
* the release version can be obtained from ...


# Examples: Running models
## Example 1: Using codon data in the form of CDS in fasta format with one mixture (ROC)
The following example illustrates how you would estimates parameters under the ROC model of a given set of protein coding genes, assuming the same mutation and selection regime for all genes.

```{r, echo = FALSE}
genome <- initializeGenomeObject(file = "genome.fasta")
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, geneAssignment = rep(1, length(genome)))
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)
model <- initializeModelObject(parameter = parameter, model = "ROC")
runMCMC(mcmc = mcmc, genome = genome, model = model)
```

## Example 2: Using codon data in the form of CDS in fasta format with one mixture (FONSE)
The following example illustrates how you would estimates parameters under the FONSE model of a given set of protein coding genes, assuming the same mutation and selection regime for all genes.

```{r, echo = FALSE}
genome <- initializeGenomeObject(file = "genome.fasta")
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, geneAssignment = rep(1, length(genome)))
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)
model <- initializeModelObject(parameter = parameter, model = "FONSE")
runMCMC(mcmc = mcmc, genome = genome, model = model)
```

## Example 3: Using codon data in the form of Ribosome footprints with one mixture (PA)
The following example illustrates how you would estimates parameters under the PA model of a given set of protein coding genes, assuming the same mutation and selection regime for all genes.

```{r, echo = FALSE}
genome <- initializeGenomeObject(file = "rfpcounts.tsv", fasta = FALSE)
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)))
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)
model <- initializeModelObject(parameter = parameter, model = "PA")
runMCMC(mcmc = mcmc, genome = genome, model = model)
```

# Examples: Advanced examples
* As the above examples illustrated the commonalities in the way all models are called. The following example will use the default ROC model for illustration purposes
## Example 4
* multiple mixture distributions with genes being initially randomly assigned to a mixture distribution. The mixture assignment of each gene will be estimated. As the below example shows, only arguments passed to the parameter object have to be adjusted to reflect a change in the number of assumed mixture distributions.

```{r, echo = FALSE}
genome <- initializeGenomeObject(file = "genome.fasta")
parameter <- initializeParameterObject(genome = genome, sphi = c(1,2,3), num.mixtures = 3, geneAssignment = sample(1:3, length(genome), replace=TRUE))
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)
model <- initializeModelObject(parameter = parameter, model = "ROC")
runMCMC(mcmc = mcmc, genome = genome, model = model)
```

## Example 5
* This example is based on the previous one, but instead of estimating the assignemnt of each gene to one of the three mixture distributions, we will fix the mixture assignemt to the initial assignment

```{r, echo = FALSE}
genome <- initializeGenomeObject(file = "genome.fasta")
parameter <- initializeParameterObject(genome = genome, sphi = c(1,2,3), num.mixtures = 3, geneAssignment = sample(1:3, length(genome), replace=TRUE))
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50, est.mix = FALSE)
model <- initializeModelObject(parameter = parameter, model = "ROC")
runMCMC(mcmc = mcmc, genome = genome, model = model)
```



