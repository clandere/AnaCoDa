#include "include/FONSE/FONSEModel.h"


//--------------------------------------------------//
//----------- Constructors & Destructors ---------- //
//--------------------------------------------------//


FONSEModel::FONSEModel() : Model()
{
	parameter = 0;
}


FONSEModel::~FONSEModel()
{
	//dtor
	//delete parameter;
}


double FONSEModel::calculateLogLikelihoodRatioPerAA(Gene& gene, std::string grouping, double *mutation, double *selection, double phiValue)
{
	unsigned numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
	double logLikelihood = 0.0;

	std::vector <unsigned> *positions;
	std::vector <double> codonProb(6, 0);

	//Find the maximum index
	unsigned minIndexVal = 0u;
	for (unsigned i = 1; i < (numCodons - 1); i++)
	{
		if (selection[minIndexVal] > selection[i])
		{
			minIndexVal = i;
		}
	}

	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(grouping, aaStart, aaEnd, false);
	for (unsigned i = aaStart, k = 0; i < aaEnd; i++, k++)
	{
		positions = gene.geneData.getCodonPositions(i); 
		for (unsigned j = 0; j < positions->size(); j++)
		{
			calculateLogCodonProbabilityVector(numCodons, positions->at(j), minIndexVal, mutation, selection, phiValue, codonProb); 
			if (codonProb[k] == 0) continue;
			logLikelihood += codonProb[k];
		} 
		//positions->clear();
	}
 	return logLikelihood;
}


double FONSEModel::calculateMutationPrior(std::string grouping, bool proposed)
{
	unsigned numCodons = SequenceSummary::GetNumCodonsForAA(grouping, true);
	double mutation[5];

	double priorValue = 0.0;

	unsigned numMutCat = parameter->getNumMutationCategories();
	double mutation_prior_sd = parameter->getMutationPriorStandardDeviation();
	for (unsigned i = 0u; i < numMutCat; i++)
	{
		parameter->getParameterForCategory(i, FONSEParameter::dM, grouping, proposed, mutation);
		for (unsigned k = 0u; k < numCodons; k++)
		{
			priorValue += Parameter::densityNorm(mutation[k], 0.0, mutation_prior_sd, true);
		}
	}
	return priorValue;
}





//------------------------------------------------//
//---------- Likelihood Ratio Functions ----------//
//------------------------------------------------//


void FONSEModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
{
	double likelihood = 0.0;
	double likelihood_proposed = 0.0;
	std::string curAA;
	std::vector <unsigned> positions;
	double mutation[5];
	double selection[5];

	//SequenceSummary *sequenceSummary = gene.getSequenceSummary(); //currently unused

	// get correct index for everything
	unsigned mutationCategory = parameter->getMutationCategory(k);
	unsigned selectionCategory = parameter->getSelectionCategory(k);
	unsigned expressionCategory = parameter->getSynthesisRateCategory(k);

	double phiValue = parameter->getSynthesisRate(geneIndex, expressionCategory, false);
	double phiValue_proposed = parameter->getSynthesisRate(geneIndex, expressionCategory, true);


	/* TODO: This loop causes a compiler warning because i is an int, but openMP won't compile if I change i to unsigned.
		Maybe worth looking into? */
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(mutation, selection, positions, curAA) reduction(+:likelihood,likelihood_proposed)
#endif
	for (unsigned i = 0u; i < getGroupListSize(); i++)
	{
		curAA = getGrouping(i);

		parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, curAA, false, mutation);
		parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, curAA, false, selection);

		likelihood += calculateLogLikelihoodRatioPerAA(gene, curAA, mutation, selection, phiValue);
		likelihood_proposed += calculateLogLikelihoodRatioPerAA(gene, curAA, mutation, selection, phiValue_proposed);
	}

	//my_print("% %\n", logLikelihood, logLikelihood_proposed);

	double stdDevSynthesisRate = parameter->getStdDevSynthesisRate(selectionCategory, false);
	double logPhiProbability = Parameter::densityLogNorm(phiValue, (-(stdDevSynthesisRate * stdDevSynthesisRate) / 2), stdDevSynthesisRate, true);
	double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, (-(stdDevSynthesisRate * stdDevSynthesisRate) / 2), stdDevSynthesisRate, true);
	double currentLogLikelihood = (likelihood + logPhiProbability);
	double proposedLogLikelihood = (likelihood_proposed + logPhiProbability_proposed);
	if (phiValue == 0) {
		my_print("phiValue is 0\n");
	}
	if (phiValue_proposed == 0) {
		my_print("phiValue_prop is 0\n");
	}
	logProbabilityRatio[0] = (proposedLogLikelihood - currentLogLikelihood) - (std::log(phiValue) - std::log(phiValue_proposed));
	logProbabilityRatio[1] = currentLogLikelihood - std::log(phiValue_proposed);
	if (std::isinf(logProbabilityRatio[1])) {
		my_print("logProb1 inf\n");
	}
	logProbabilityRatio[2] = proposedLogLikelihood - std::log(phiValue);
	if (std::isinf(logProbabilityRatio[2])) {
		my_print("logProb2 inf\n");
	}

	//------------NOTE: Jeremy, Cedric changed the reverse jump to where we DON'T include it. I had my PA
	//LogLikelihood go to 0 because of now missing terms. I have added the code underneath to where we calculate it---/

	logProbabilityRatio[3] = currentLogLikelihood;
	logProbabilityRatio[4] = proposedLogLikelihood;
}


void FONSEModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, std::vector<double> &logAcceptanceRatioForAllMixtures)
{
	int numGenes = genome.getGenomeSize();
	//int numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
	double likelihood = 0.0;
	double likelihood_proposed = 0.0;

	double mutation[5];
	double selection[5];
	double mutation_proposed[5];
	double selection_proposed[5];

	std::string curAA;

	Gene *gene;
	SequenceSummary *sequenceSummary;
	unsigned aaIndex = SequenceSummary::AAToAAIndex(grouping);

#ifdef _OPENMP
//#ifndef __APPLE__
	#pragma omp parallel for private(mutation, selection, mutation_proposed, selection_proposed, curAA, gene, sequenceSummary) reduction(+:likelihood,likelihood_proposed)
#endif
	for (unsigned i = 0u; i < numGenes; i++)
	{
		gene = &genome.getGene(i);
		sequenceSummary = gene->getSequenceSummary();
		if (sequenceSummary->getAACountForAA(aaIndex) == 0) continue;

		// which mixture element does this gene belong to
		unsigned mixtureElement = parameter->getMixtureAssignment(i);
		// how is the mixture element defined. Which categories make it up
		unsigned mutationCategory = parameter->getMutationCategory(mixtureElement);
		unsigned selectionCategory = parameter->getSelectionCategory(mixtureElement);
		unsigned expressionCategory = parameter->getSynthesisRateCategory(mixtureElement);
		// get phi value, calculate likelihood conditional on phi
		double phiValue = parameter->getSynthesisRate(i, expressionCategory, false);


		// get current mutation and selection parameter
		parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, grouping, false, mutation);
		parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, grouping, false, selection);

		// get proposed mutation and selection parameter
		parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, grouping, true, mutation_proposed);
		parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, grouping, true, selection_proposed);
		
		likelihood += calculateLogLikelihoodRatioPerAA(*gene, grouping, mutation, selection, phiValue);
		likelihood_proposed += calculateLogLikelihoodRatioPerAA(*gene, grouping, mutation_proposed, selection_proposed, phiValue);
	}
	//likelihood_proposed = likelihood_proposed + calculateMutationPrior(grouping, true);
	//likelihood = likelihood + calculateMutationPrior(grouping, false);

	logAcceptanceRatioForAllMixtures[0] = (likelihood_proposed - likelihood);
}


void FONSEModel::calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> & logProbabilityRatio)
{
	double lpr = 0.0;
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
		lpr += Parameter::densityLogNorm(phi, proposedMphi[mixture], proposedStdDevSynthesisRate[mixture], true)
			   - Parameter::densityLogNorm(phi, currentMphi[mixture], currentStdDevSynthesisRate[mixture], true);
	}
	logProbabilityRatio[0] = lpr;
}





//----------------------------------------------------------//
//---------- Initialization and Restart Functions ----------//
//----------------------------------------------------------//


void FONSEModel::initTraces(unsigned samples, unsigned num_genes)
{
	parameter->initAllTraces(samples, num_genes);
}


void FONSEModel::writeRestartFile(std::string filename)
{
	return parameter->writeEntireRestartFile(filename);
}





//----------------------------------------//
//---------- Category Functions ----------//
//----------------------------------------//


double FONSEModel::getCategoryProbability(unsigned i)
{
	return parameter->getCategoryProbability(i);
}


unsigned FONSEModel::getMutationCategory(unsigned mixture)
{
	return parameter->getMutationCategory(mixture);
}


unsigned FONSEModel::getSelectionCategory(unsigned mixture)
{
	return parameter->getSelectionCategory(mixture);
}


unsigned FONSEModel::getSynthesisRateCategory(unsigned mixture)
{
	return parameter->getSynthesisRateCategory(mixture);
}


std::vector<unsigned> FONSEModel::getMixtureElementsOfSelectionCategory(unsigned k)
{
	return parameter->getMixtureElementsOfSelectionCategory(k);
}





//------------------------------------------//
//---------- Group List Functions ----------//
//------------------------------------------//


unsigned FONSEModel::getGroupListSize()
{
	return parameter->getGroupListSize();
}


std::string FONSEModel::getGrouping(unsigned index)
{
	return parameter->getGrouping(index);
}





//---------------------------------------------------//
//---------- stdDevSynthesisRate Functions ----------//
//---------------------------------------------------//


double FONSEModel::getStdDevSynthesisRate(unsigned selectionCategory, bool proposed)
{
	return parameter->getStdDevSynthesisRate(selectionCategory, proposed);
}


double FONSEModel::getCurrentStdDevSynthesisRateProposalWidth()
{
	return parameter->getCurrentStdDevSynthesisRateProposalWidth();
}


void FONSEModel::updateStdDevSynthesisRate()
{
	parameter->updateStdDevSynthesisRate();
}





//----------------------------------------------//
//---------- Synthesis Rate Functions ----------//
//----------------------------------------------//


double FONSEModel::getSynthesisRate(unsigned index, unsigned mixture, bool proposed)
{
	return parameter->getSynthesisRate(index, mixture, proposed);
}


void FONSEModel::updateSynthesisRate(unsigned i, unsigned k)
{
	parameter->updateSynthesisRate(i, k);
}





//-----------------------------------------//
//---------- Iteration Functions ----------//
//-----------------------------------------//


unsigned FONSEModel::getLastIteration()
{
	return parameter->getLastIteration();
}


void FONSEModel::setLastIteration(unsigned iteration)
{
	parameter->setLastIteration(iteration);
}





//-------------------------------------//
//---------- Trace Functions ----------//
//-------------------------------------//


void FONSEModel::updateStdDevSynthesisRateTrace(unsigned sample)
{
	parameter->updateStdDevSynthesisRateTrace(sample);
}


void FONSEModel::updateSynthesisRateTrace(unsigned sample, unsigned i)
{
	parameter->updateSynthesisRateTrace(sample, i);
}


void FONSEModel::updateMixtureAssignmentTrace(unsigned sample, unsigned i)
{
	parameter->updateMixtureAssignmentTrace(sample, i);
}


void FONSEModel::updateMixtureProbabilitiesTrace(unsigned sample)
{
	parameter->updateMixtureProbabilitiesTrace(sample);
}


void FONSEModel::updateCodonSpecificParameterTrace(unsigned sample, std::string grouping)
{
	parameter->updateCodonSpecificParameterTrace(sample, grouping);
}


void FONSEModel::updateHyperParameterTraces(unsigned sample)
{
	updateStdDevSynthesisRateTrace(sample);
}


void FONSEModel::updateTracesWithInitialValues(Genome & genome)
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


void FONSEModel::adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void FONSEModel::adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void FONSEModel::adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt)
{
	parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth, lastIteration, adapt);
}


void FONSEModel::adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt)
{
	adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


void FONSEModel::proposeCodonSpecificParameter()
{
	parameter->proposeCodonSpecificParameter();
}


void FONSEModel::proposeHyperParameters()
{
	parameter->proposeHyperParameters();
}


void FONSEModel::proposeSynthesisRateLevels()
{
	parameter->proposeSynthesisRateLevels();
}


unsigned FONSEModel::getNumPhiGroupings()
{
	return parameter->getNumObservedPhiSets();
}


unsigned FONSEModel::getMixtureAssignment(unsigned index)
{
	return parameter->getMixtureAssignment(index);
}


unsigned FONSEModel::getNumMixtureElements()
{
	return parameter->getNumMixtureElements();
}


unsigned FONSEModel::getNumSynthesisRateCategories()
{
	return parameter->getNumSynthesisRateCategories();
}


void FONSEModel::setNumPhiGroupings(unsigned value)
{
	parameter->setNumObservedPhiSets(value);
}


void FONSEModel::setMixtureAssignment(unsigned i, unsigned catOfGene)
{
	parameter->setMixtureAssignment(i, catOfGene);
}


void FONSEModel::setCategoryProbability(unsigned mixture, double value)
{
	parameter->setCategoryProbability(mixture, value);
}


void FONSEModel::updateCodonSpecificParameter(std::string grouping)
{
	parameter->updateCodonSpecificParameter(grouping);
}


void FONSEModel::updateGibbsSampledHyperParameters(Genome &genome)
{
	//TODO: fill in
}


void FONSEModel::updateAllHyperParameter()
{
	updateStdDevSynthesisRate();
}

void FONSEModel::updateHyperParameter(unsigned hp)
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


void FONSEModel::simulateGenome(Genome & genome)
{
	unsigned codonIndex;
	std::string curAA;

	std::string tmpDesc = "Simulated Gene";

	for (unsigned geneIndex = 0; geneIndex < genome.getGenomeSize(); geneIndex++) //loop over all genes in the genome
	{
		if (geneIndex % 100 == 0) my_print("Simulating Gene %\n", geneIndex);
		Gene gene = genome.getGene(geneIndex);
		SequenceSummary sequenceSummary = gene.geneData;
		std::string tmpSeq = "ATG"; //Always will have the start amino acid


		unsigned mixtureElement = getMixtureAssignment(geneIndex);
		unsigned mutationCategory = getMutationCategory(mixtureElement);
		unsigned selectionCategory = getSelectionCategory(mixtureElement);
		unsigned synthesisRateCategory = getSynthesisRateCategory(mixtureElement);
		double phi = getSynthesisRate(geneIndex, synthesisRateCategory, false);

		std::string geneSeq = gene.getSequence();
		for (unsigned position = 1; position < (geneSeq.size() / 3); position++)
		{
			std::string codon = geneSeq.substr((position * 3), 3);
			curAA = SequenceSummary::codonToAA(codon);

			//TODO: Throw an error here instead
			if (curAA == "X") {
				if (position < (geneSeq.size() / 3) - 1) my_print("Warning: Internal stop codon found in gene % at position %. Ignoring and moving on.\n", gene.getId(), position);
				continue;
			}

			unsigned numCodons = SequenceSummary::GetNumCodonsForAA(curAA);

			double* codonProb = new double[numCodons](); //size the arrays to the proper size based on # of codons.
			double* mutation = new double[numCodons - 1]();
			double* selection = new double[numCodons - 1]();


			if (curAA == "M" || curAA == "W")
			{
				codonProb[0] = 1;
			}
			else
			{
				getParameterForCategory(mutationCategory, FONSEParameter::dM, curAA, false, mutation);
				getParameterForCategory(selectionCategory, FONSEParameter::dOmega, curAA, false, selection);
				calculateCodonProbabilityVector(numCodons, position, mutation, selection, phi, codonProb);
			}


			codonIndex = Parameter::randMultinom(codonProb, numCodons);
			unsigned aaStart, aaEnd;
			SequenceSummary::AAToCodonRange(curAA, aaStart, aaEnd, false);  //need the first spot in the array where the codons for curAA are
			codon = sequenceSummary.indexToCodon(aaStart + codonIndex);//get the correct codon based off codonIndex
			tmpSeq += codon;
		}
		std::string codon = sequenceSummary.indexToCodon((unsigned)Parameter::randUnif(61.0, 64.0)); //randomly choose a stop codon, from range 61-63
		tmpSeq += codon;
		Gene simulatedGene(tmpSeq, tmpDesc, gene.getId());
		genome.addGene(simulatedGene, true);
	}
}


void FONSEModel::printHyperParameters()
{
	for (unsigned i = 0u; i < getNumSynthesisRateCategories(); i++)
	{
		my_print("stdDevSynthesisRate posterior estimate for selection category %: %\n", i, getStdDevSynthesisRate(i));
	}
	my_print("\t current stdDevSynthesisRate proposal width: %\n", getCurrentStdDevSynthesisRateProposalWidth());
}


/* getParameter (RCPP EXPOSED)
* Arguments: None
*
* Returns the FONSEParameter of the model.
*/
FONSEParameter* FONSEModel::getParameter()
{
	return parameter;
}

void FONSEModel::setParameter(FONSEParameter &_parameter)
{
	parameter = &_parameter;
}


double FONSEModel::calculateAllPriors()
{
	double priorRatio = 0.0;
	unsigned size = getGroupListSize();

	for (unsigned i = 0; i < size; i++)
	{
		std::string grouping = getGrouping(i);
		priorRatio += calculateMutationPrior(grouping, false);
	}

	// add more priors if necessary.

	return priorRatio;
}

//Calculates the log probability of each codon for an amino acid and puts them in a vector.
void FONSEModel::calculateLogCodonProbabilityVector(unsigned numCodons, unsigned position, unsigned minIndexValue,
												 double *mutation, double *selection, double phi, std::vector <double> &codonProb)
{
	double denominator;

	/* log(c_i) = \Delta M - (\phi * \beta(i) * \Delta \omega),                 *
	 * where \beta(i) = a_1 + (i * a_2)                                         *
	 *                                                                          *
	 * Right now a_1 and a_2 are set to 4.0. However, we are planning on making *
	 * them hyperparameters in the future, since they are constant for the      *
	 * entire genome.                                                           */

	 // If the min(selection) is less than zero than we have to adjust the reference codon.
	 // If the reference codon is the min value (0) then we do not have to adjust the reference codon.
	 // This is necessary to deal with very large phi values (> 10^4) and avoid producing Inf which then
	 // causes the denominator to be Inf (Inf / Inf = NaN).
	if (selection[minIndexValue] < 0.0) {
		denominator = 0.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = -(mutation[i] - mutation[minIndexValue]) - (phi * (4.0 + (4.0 * position)) * (selection[i] - selection[minIndexValue]));
			denominator += std::exp(codonProb[i]);
		}
		//Alphabetically, the last codon is the reference codon.
		codonProb[numCodons - 1] = (mutation[minIndexValue]) + (phi * (4.0 + (4.0 * position)) * selection[minIndexValue]);
		denominator += std::exp(codonProb[numCodons - 1]);
	}
	else
	{
		denominator = 1.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = -(mutation[i]) - (phi * (4.0 + (4.0 * position)) * selection[i]);
			denominator += std::exp(codonProb[i]);
		}
		//Again, the last codon is the reference codon
		codonProb[numCodons - 1] = 0.0;
	}

	//Here we take the log of the denominator (the summation term) so that we can finish calculating 
	//the log probabilities simple by subtracting the log of the denominator from each element.
	denominator = std::log(denominator);
	for (unsigned i = 0; i < numCodons; i++)
	{
		codonProb[i] = codonProb[i] - denominator;
	}
}


//Since the simulateGenome function utilizes the codon probability vector, but doesn't deal with the log values,
//this function simply returns the vector with each codon's probability.
void FONSEModel::calculateCodonProbabilityVector(unsigned numCodons, unsigned position,
													double *mutation, double *selection, double phi, double codonProb[])
{
	double denominator;
	unsigned minIndexValue = 0u;
	for (unsigned i = 1; i < (numCodons - 1); i++)
	{
		if (selection[minIndexValue] < selection[i])
		{
			minIndexValue = i;
		}
	}

	/* c_i = exp[\Delta M - (\phi * \beta(i) * \Delta \omega)],                 *
	* where \beta(i) = a_1 + (i * a_2)                                         *
	*                                                                          *
	* Right now a_1 and a_2 are set to 4.0. However, we are planning on making *
	* them hyperparameters in the future, since they are constant for the      *
	* entire genome.                                                           */


	// If the min(selection) is less than zero than we have to adjust the reference codon.
	// If the reference codon is the min value (0) then we do not have to adjust the reference codon.
	// This is necessary to deal with very large phi values (> 10^4) and avoid producing Inf which then
	// causes the denominator to be Inf (Inf / Inf = NaN).
	if (selection[minIndexValue] < 0.0) {
		denominator = 0.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = std::exp(-(mutation[i] - mutation[minIndexValue]) - (phi * (4.0 + (4.0 * position)) * (selection[i] - selection[minIndexValue])));
			denominator += codonProb[i];
		}
		//Alphabetically, the last codon is the reference codon.
		codonProb[numCodons - 1] = std::exp((mutation[minIndexValue]) + (phi * (4.0 + (4.0 * position)) * selection[minIndexValue]));
		denominator += codonProb[numCodons - 1];
	}
	else
	{
		denominator = 1.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = std::exp(-(mutation[i]) - (phi * (4.0 + (4.0 * position)) * selection[i]));
			denominator += codonProb[i];
		}
		//Again, the last codon is the reference codon
		codonProb[numCodons - 1] = 1.0;
	}

	//As is found in ROCModel.cpp, multiplication is a faster operation than division so we 
	//save time here by dividing once and then multiplying numCodons times instead of dividing
	//numCodons times.
	denominator = 1 / denominator;
	for (unsigned i = 0; i < numCodons; i++)
	{
		codonProb[i] *= denominator;
	}
}

void FONSEModel::getParameterForCategory(unsigned category, unsigned param, std::string aa, bool proposal, double* returnValue)
{
	parameter->getParameterForCategory(category, param, aa, proposal, returnValue);
}








// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//


#ifndef STANDALONE


/*std::vector<double> FONSEModel::CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection, double phi)
{
	unsigned numCodons = mutation.size() + 1;
	double* _mutation = &mutation[0];
	double* _selection = &selection[0];
	double* codonProb = new double[numCodons]();
	calculateCodonProbabilityVector(numCodons, _mutation, _selection, phi, codonProb);
	std::vector<double> returnVector(codonProb, codonProb + numCodons);
	return returnVector;
}*/


#endif
