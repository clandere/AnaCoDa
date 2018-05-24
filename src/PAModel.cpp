#include "include/PA/PAModel.h"

//R runs only
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//--------------------------------------------------//
//----------- Constructors & Destructors ---------- //
//--------------------------------------------------//

PAModel::PAModel(unsigned _RFPCountColumn) : Model()
{
	parameter = NULL;
	RFPCountColumn = _RFPCountColumn;
	//ctor
}


PAModel::~PAModel()
{
	//dtor
	//TODO: call Parent's deconstructor
	//delete parameter;
}


double PAModel::calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime, unsigned currRFPValue,
                                                      unsigned currNumCodonsInMRNA, double phiValue)
{
	double logLikelihood = ((std::lgamma((currNumCodonsInMRNA * currAlpha) + currRFPValue)) - (std::lgamma(currNumCodonsInMRNA * currAlpha)))
						   + (currRFPValue * (std::log(phiValue) - std::log(currLambdaPrime + phiValue)))
						   + ((currNumCodonsInMRNA * currAlpha) * (std::log(currLambdaPrime) - std::log(currLambdaPrime + phiValue)));
    //my_print("alpha = % \n lambda = % \n phivalue = % \n loglikelihood % \n", currAlpha, currLambdaPrime, phiValue, logLikelihood);

	return logLikelihood;
}





//------------------------------------------------//
//---------- Likelihood Ratio Functions ----------//
//------------------------------------------------//


void PAModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
{
	double logLikelihood = 0.0;
	double logLikelihood_proposed = 0.0;

	unsigned alphaCategory = parameter->getMutationCategory(k);
	unsigned lambdaPrimeCategory = parameter->getSelectionCategory(k);
	unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(k);


	double phiValue = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, false);
	double phiValue_proposed = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, true);

#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for reduction(+:logLikelihood,logLikelihood_proposed)
#endif
	for (unsigned index = 0; index < getGroupListSize(); index++) //number of codons, without the stop codons
	{
		std::string codon = getGrouping(index);

		double currAlpha = getParameterForCategory(alphaCategory, PAParameter::alp, codon, false);
		double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PAParameter::lmPri, codon, false);
		unsigned currRFPValue = gene.geneData.getRFPValue(index, RFPCountColumn);

		unsigned currNumCodonsInMRNA = gene.geneData.getCodonCountForCodon(index);
		if (currNumCodonsInMRNA == 0) continue;

		logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPValue, currNumCodonsInMRNA, phiValue);
		logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPValue, currNumCodonsInMRNA, phiValue_proposed);
	}

	double stdDevSynthesisRate = parameter->getStdDevSynthesisRate(lambdaPrimeCategory, false);
	double logPhiProbability = Parameter::densityLogNorm(phiValue, (-(stdDevSynthesisRate * stdDevSynthesisRate) / 2), stdDevSynthesisRate, true);
	double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, (-(stdDevSynthesisRate * stdDevSynthesisRate) / 2), stdDevSynthesisRate, true);
	double currentLogPosterior = (logLikelihood + logPhiProbability);
	double proposedLogPosterior = (logLikelihood_proposed + logPhiProbability_proposed);

	CSHyperParameters[0] = logProbabilityRatio[0] = (proposedLogPosterior - currentLogPosterior) - (std::log(phiValue) - std::log(phiValue_proposed));//Is recalulcated in MCMC
	CSHyperParameters[1] = logProbabilityRatio[1] = currentLogPosterior - std::log(phiValue_proposed);
	CSHyperParameters[2] = logProbabilityRatio[2] = proposedLogPosterior - std::log(phiValue);
	CSHyperParameters[3] = logProbabilityRatio[3] = currentLogPosterior;
	CSHyperParameters[4] = logProbabilityRatio[4] = proposedLogPosterior;
	logProbabilityRatio[5] = logLikelihood;
	logProbabilityRatio[6] = logLikelihood_proposed;

}


void PAModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome,
                                                                std::vector<double> &logAcceptanceRatioForAllMixtures)
{
	double logLikelihood = 0.0;
	double logLikelihood_proposed = 0.0;
    ExpectedZ =  calculateExpectedZ(genome);
    Y = calculateY(genome);
	Gene *gene;
	unsigned index = SequenceSummary::codonToIndex(grouping);

#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(gene) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		gene = &genome.getGene(i);
		// which mixture element does this gene belong to
		unsigned mixtureElement = parameter->getMixtureAssignment(i);
		// how is the mixture element defined. Which categories make it up
		unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
		unsigned lambdaPrimeCategory = parameter->getSelectionCategory(mixtureElement);
		unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
		// get non codon specific values, calculate likelihood conditional on these
		double phiValue = parameter->getSynthesisRate(i, /*synthesisRateCategory*/mixtureElement, false);
		unsigned currRFPValue = gene->geneData.getRFPValue(index, RFPCountColumn);
		unsigned currNumCodonsInMRNA = gene->geneData.getCodonCountForCodon(index);
		if (currNumCodonsInMRNA == 0) continue;


		double currAlpha = getParameterForCategory(alphaCategory, PAParameter::alp, grouping, false);
		double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PAParameter::lmPri, grouping, false);
        currLambdaPrime *= (ExpectedZ/Y);

		double propAlpha = getParameterForCategory(alphaCategory, PAParameter::alp, grouping, true);
		double propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PAParameter::lmPri, grouping, true);


		logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPValue, currNumCodonsInMRNA, phiValue);
		logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(propAlpha, propLambdaPrime, currRFPValue, currNumCodonsInMRNA, phiValue);
	}
	logAcceptanceRatioForAllMixtures[0] = logLikelihood_proposed - logLikelihood;
	logAcceptanceRatioForAllMixtures[1] = logLikelihood;
	logAcceptanceRatioForAllMixtures[2] = logLikelihood_proposed;
	logAcceptanceRatioForAllMixtures[3] = logLikelihood;
	logAcceptanceRatioForAllMixtures[4] = logLikelihood_proposed;
}


void PAModel::calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> & logProbabilityRatio)
{

	double lpr = 0.0; // this variable is only needed because OpenMP doesn't allow variables in reduction clause to be reference

	unsigned selectionCategory = getNumSynthesisRateCategories();
	std::vector<double> currentStdDevSynthesisRate(selectionCategory, 0.0);
	std::vector<double> currentMphi(selectionCategory, 0.0);
	std::vector<double> proposedStdDevSynthesisRate(selectionCategory, 0.0);
	std::vector<double> proposedMphi(selectionCategory, 0.0);
	for (unsigned i = 0u; i < selectionCategory; i++)
	{
		currentStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, false);
		currentMphi[i] = -((currentStdDevSynthesisRate[i] * currentStdDevSynthesisRate[i]) / 2);
		proposedStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, true);
		proposedMphi[i] = -((proposedStdDevSynthesisRate[i] * proposedStdDevSynthesisRate[i]) / 2);
		// take the Jacobian into account for the non-linear transformation from logN to N distribution
		lpr -= (std::log(currentStdDevSynthesisRate[i]) - std::log(proposedStdDevSynthesisRate[i]));
		// take prior into account
		//TODO(Cedric): make sure you can control that prior from R
		lpr -= Parameter::densityNorm(currentStdDevSynthesisRate[i], 1.0, 0.1, true) - Parameter::densityNorm(proposedStdDevSynthesisRate[i], 1.0, 0.1, true);
	}

	logProbabilityRatio.resize(1);
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for reduction(+:lpr)
#endif
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		unsigned mixture = getMixtureAssignment(i);
		mixture = getSynthesisRateCategory(mixture);
		double phi = getSynthesisRate(i, mixture, false);
		if (i == 0)
		{
			//my_print("proposed: %\n", Parameter::densityLogNorm(phi, proposedMPhi[mixture], proposedStdDevSynthesisRate[mixture], false));
			//my_print("current: %\n", Parameter::densityLogNorm(phi, currentMPhi, currentStdDevSynthesisRate, false));
		}
		lpr += Parameter::densityLogNorm(phi, proposedMphi[mixture], proposedStdDevSynthesisRate[mixture], true) -
				Parameter::densityLogNorm(phi, currentMphi[mixture], currentStdDevSynthesisRate[mixture], true);
		//my_print("LPR: %\n", lpr);
	}

	logProbabilityRatio[0] = lpr;
}





//----------------------------------------------------------//
//---------- Initialization and Restart Functions ----------//
//----------------------------------------------------------//


void PAModel::initTraces(unsigned samples, unsigned num_genes)
{
	parameter->initAllTraces(samples, num_genes);
}


void PAModel::writeRestartFile(std::string filename)
{
	return parameter->writeEntireRestartFile(filename);
}





//----------------------------------------//
//---------- Category Functions ----------//
//----------------------------------------//


double PAModel::getCategoryProbability(unsigned i)
{
	return parameter->getCategoryProbability(i);
}


unsigned PAModel::getMutationCategory(unsigned mixture)
{
	return parameter->getMutationCategory(mixture);
}


unsigned PAModel::getSelectionCategory(unsigned mixture)
{
	return parameter->getSelectionCategory(mixture);
}


unsigned PAModel::getSynthesisRateCategory(unsigned mixture)
{
	return parameter->getSynthesisRateCategory(mixture);
}


std::vector<unsigned> PAModel::getMixtureElementsOfSelectionCategory(unsigned k)
{
	return parameter->getMixtureElementsOfSelectionCategory(k);
}





//------------------------------------------//
//---------- Group List Functions ----------//
//------------------------------------------//


unsigned PAModel::getGroupListSize()
{
	return parameter->getGroupListSize();
}


std::string PAModel::getGrouping(unsigned index)
{
	return parameter->getGrouping(index);
}





//---------------------------------------------------//
//---------- stdDevSynthesisRate Functions ----------//
//---------------------------------------------------//


double PAModel::getStdDevSynthesisRate(unsigned selectionCategory, bool proposed)
{
	return parameter->getStdDevSynthesisRate(selectionCategory, proposed);
}


double PAModel::getCurrentStdDevSynthesisRateProposalWidth()
{
	return parameter->getCurrentStdDevSynthesisRateProposalWidth();
}


void PAModel::updateStdDevSynthesisRate()
{
	parameter->updateStdDevSynthesisRate();
}





//----------------------------------------------//
//---------- Synthesis Rate Functions ----------//
//----------------------------------------------//


double PAModel::getSynthesisRate(unsigned index, unsigned mixture, bool proposed)
{
	return parameter->getSynthesisRate(index, mixture, proposed);
}


void PAModel::updateSynthesisRate(unsigned i, unsigned k)
{
	parameter->updateSynthesisRate(i, k);
}





//-----------------------------------------//
//---------- Iteration Functions ----------//
//-----------------------------------------//


unsigned PAModel::getLastIteration()
{
	return parameter->getLastIteration();
}


void PAModel::setLastIteration(unsigned iteration)
{
	parameter->setLastIteration(iteration);
}





//-------------------------------------//
//---------- Trace Functions ----------//
//-------------------------------------//


void PAModel::updateStdDevSynthesisRateTrace(unsigned sample)
{
	parameter->updateStdDevSynthesisRateTrace(sample);
}


void PAModel::updateSynthesisRateTrace(unsigned sample, unsigned i)
{
	parameter->updateSynthesisRateTrace(sample, i);
}


void PAModel::updateMixtureAssignmentTrace(unsigned sample, unsigned i)
{
	parameter->updateMixtureAssignmentTrace(sample, i);
}


void PAModel::updateMixtureProbabilitiesTrace(unsigned sample)
{
	parameter->updateMixtureProbabilitiesTrace(sample);
}


void PAModel::updateCodonSpecificParameterTrace(unsigned sample, std::string codon)
{
	parameter->updateCodonSpecificParameterTrace(sample, codon);
}

void PAModel::updateCodonSpecificHyperParameter(std::string aa, double randomNumber)
{
	parameter->updateCodonSpecificHyperParameter(aa, randomNumber, CSHyperParameters[0], CSHyperParameters[3], CSHyperParameters[4],
        CSHyperParameters[1],CSHyperParameters[2]);
}


void PAModel::updateHyperParameterTraces(unsigned sample)
{
	updateStdDevSynthesisRateTrace(sample);
}


void PAModel::updateTracesWithInitialValues(Genome & genome)
{
	std::vector <std::string> groupList = parameter->getGroupList();

	for (unsigned i = 0; i < genome.getGenomeSize(); i++)
	{
		parameter->updateSynthesisRateTrace(0, i);
		parameter->updateMixtureAssignmentTrace(0, i);
	}

	for (unsigned i = 0; i < groupList.size(); i++)
	{
		parameter->updateCodonSpecificParameterTrace(0, getGrouping(i));
	}
}





//----------------------------------------------//
//---------- Adaptive Width Functions ----------//
//----------------------------------------------//


void PAModel::adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void PAModel::adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void PAModel::adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt)
{
	parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth, lastIteration, adapt);
}


void PAModel::adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt)
{
	adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


void PAModel::proposeCodonSpecificParameter()
{
	parameter->proposeCodonSpecificParameter();
}


void PAModel::proposeHyperParameters()
{
	parameter->proposeStdDevSynthesisRate();
}


void PAModel::proposeSynthesisRateLevels()
{
	parameter->proposeSynthesisRateLevels();
}


unsigned PAModel::getNumPhiGroupings()
{
	return parameter->getNumObservedPhiSets();
}


unsigned PAModel::getMixtureAssignment(unsigned index)
{
	return parameter->getMixtureAssignment(index);
}


unsigned PAModel::getNumMixtureElements()
{
	return parameter->getNumMixtureElements();
}


unsigned PAModel::getNumSynthesisRateCategories()
{
	return parameter->getNumSynthesisRateCategories();
}


void PAModel::setNumPhiGroupings(unsigned value)
{
	parameter->setNumObservedPhiSets(value);
}


void PAModel::setMixtureAssignment(unsigned i, unsigned catOfGene)
{
	parameter->setMixtureAssignment(i, catOfGene);
}


void PAModel::setCategoryProbability(unsigned mixture, double value)
{
	parameter->setCategoryProbability(mixture, value);
}


void PAModel::updateCodonSpecificParameter(std::string aa)
{
	parameter->updateCodonSpecificParameter(aa);
}

void PAModel::updateGibbsSampledHyperParameters(Genome &genome)
{
	//TODO: fill in
}


void PAModel::updateAllHyperParameter()
{
	updateStdDevSynthesisRate();
}


void PAModel::updateHyperParameter(unsigned hp)
{
	// NOTE: when adding additional hyper parameter, also add to updateAllHyperParameter()
	switch (hp)
	{
		case 0:
			updateStdDevSynthesisRate();
			break;
		default:
			updateStdDevSynthesisRate();
			break;
	}
}


void PAModel::simulateGenome(Genome &genome)
{
	for (unsigned geneIndex = 0u; geneIndex < genome.getGenomeSize(); geneIndex++)
	{
		unsigned mixtureElement = getMixtureAssignment(geneIndex);
		Gene gene = genome.getGene(geneIndex);
		double phi = parameter->getSynthesisRate(geneIndex, mixtureElement, false);
		Gene tmpGene = gene;
		for (unsigned codonIndex = 0u; codonIndex < 61; codonIndex++)
		{
			std::string codon = SequenceSummary::codonArray[codonIndex];
			unsigned alphaCat = parameter->getMutationCategory(mixtureElement);
			unsigned lambdaPrimeCat = parameter->getSelectionCategory(mixtureElement);

			double alpha = getParameterForCategory(alphaCat, PAParameter::alp, codon, false);
			double lambdaPrime = getParameterForCategory(lambdaPrimeCat, PAParameter::lmPri, codon, false);

			double alphaPrime = alpha * gene.geneData.getCodonCountForCodon(codon);

#ifndef STANDALONE
			RNGScope scope;
			NumericVector xx(1);
			xx = rgamma(1, alphaPrime, 1.0/lambdaPrime);
			xx = rpois(1, xx[0] * phi);
			tmpGene.geneData.setRFPValue(codonIndex, xx[0], RFPCountColumn);
#else
			std::gamma_distribution<double> GDistribution(alphaPrime, 1.0/lambdaPrime);
			double tmp = GDistribution(Parameter::generator);
			std::poisson_distribution<unsigned> PDistribution(phi * tmp);
			unsigned simulatedValue = PDistribution(Parameter::generator);
			tmpGene.geneData.setRFPValue(codonIndex, simulatedValue, RFPCountColumn);
#endif
		}
		genome.addGene(tmpGene, true);
	}
}


void PAModel::printHyperParameters()
{
	for (unsigned i = 0u; i < getNumSynthesisRateCategories(); i++)
	{
		//my_print("stdDevSynthesisRate posterior estimate for selection category %: %\n", i, parameter->getStdDevSynthesisRate(i));
	}
	//my_print("\t current stdDevSynthesisRate proposal width: %\n", getCurrentStdDevSynthesisRateProposalWidth());
}


/* getParameter (RCPP EXPOSED)
 * Arguments: None
 *
 * Returns the PAParameter of the model.
*/
PAParameter* PAModel::getParameter()
{
	return parameter;
}


/* setParameter (RCPP EXPOSED)
 * Arguments: The PAParameter object to be referenced and stored.
 * Sets the private parameter object to the parameter object specified.
*/
void PAModel::setParameter(PAParameter &_parameter)
{
	parameter = &_parameter;
}


double PAModel::calculateAllPriors()
{
	return 0.0; //TODO(Cedric): implement me, see ROCModel
}

//TODO: Debug Function
double PAModel::calculateExpectedZ(Genome &genome)
{
    return parameter->calculateExpectedZ(genome);
}

double PAModel::getParameterForCategory(unsigned category, unsigned param, std::string codon, bool proposal)
{
	return parameter->getParameterForCategory(category, param, codon, proposal);
}

//DEBUG function
double PAModel::calculateY(Genome &genome)
{
    std::vector <Gene> tmp = genome.getGenes();
    int sum = 0;
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++){
        sum += tmp[i].getTotalRFPCount();    
    }
    return sum;
}
