#include "include/base/Model.h"

Model::Model()
{
//ctor
}

// TODO: Rule of Three dictates we may need a copy assignment operator as well (operator=)

Model::~Model()
{
//dtor
}

//Cedric: This functions will replace calculateMutationPrior in ROC/FONSE model and allows us to more generally use priors on codon specific parameters.
//			We have to first change how current and proposed csp values are stored to move the function getParameterForCategory up into the base parameter class.

double Model::calculatePriorForCodonSpecificParam(Parameter *parameter, std::string grouping, unsigned paramType, bool proposed)
{
	unsigned numCodons = SequenceSummary::GetNumCodonsForAA(grouping, true); // TODO(Cedric): rename getNumCodonsForGrouping and have it return 1 if grouping is a codon to make it applicable for RFP
	double parameterValues[5];

	double priorValue = 0.0;

	unsigned numCat = parameter->getNumMutationCategories(); // TODO(Cedric): rename this function getNumCSPone or something like that
	double prior_sd = parameter->getCodonSpecificPriorStdDev(paramType);
	for (unsigned i = 0u; i < numCat; i++)
	{
		//parameter->getParameterForCategory(i, paramType, grouping, proposed, mutation); // TODO(Cedric): we have to change how csp are stored first!
		for (unsigned k = 0u; k < numCodons; k++)
		{
			priorValue += Parameter::densityNorm(parameterValues[k], 0.0, prior_sd, true);
		}
	}
	return priorValue;
}

//Temporary for debugging
void Model::updateCodonSpecificHyperParameter(std::string grouping, double randomNumber){
    return;
}
