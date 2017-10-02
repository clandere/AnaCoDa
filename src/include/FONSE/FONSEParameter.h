#ifndef FONSEPARAMETER_H
#define FONSEPARAMETER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>


#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "../base/Parameter.h"
#include "../base/Trace.h"

class FONSEParameter : public Parameter
{
	private:

		double bias_csp;
		double mutation_prior_sd;

		std::vector <double> propose(std::vector <double> currentParam, double(*proposal)(double a, double b), double A, std::vector <double> B);


	public:

		//Constructors & Destructors:
		FONSEParameter();
		explicit FONSEParameter(std::string filename);
		FONSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector <unsigned> geneAssignment,
				   std::vector <std::vector <unsigned> > thetaKmatrix, bool splitSer = true,
				   std::string _mutationSelectionState = "allUnique");
		FONSEParameter& operator=(const FONSEParameter& rhs);
		FONSEParameter(const FONSEParameter &other); //TODO: No longer needed?
		virtual ~FONSEParameter();



		//Initialization, Restart, Index Checking:
		void initFONSEParameterSet();
		void initFONSEValuesFromFile(std::string filename);
		void writeEntireRestartFile(std::string filename);
		void writeFONSERestartFile(std::string filename);
		void initFromRestartFile(std::string filename);

		void initAllTraces(unsigned samples, unsigned num_genes);
		void initMutationCategories(std::vector<std::string> files, unsigned numCategories);
		void initSelectionCategories(std::vector<std::string> files, unsigned numCategories);



		//Trace Functions:
		void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping);


		//Covariance Functions:
		CovarianceMatrix& getCovarianceMatrixForAA(std::string aa);

		//CSP Functions:
		double getCurrentCodonSpecificProposalWidth(unsigned aa);
		void proposeCodonSpecificParameter();
		void updateCodonSpecificParameter(std::string grouping);


		//Prior Functions:
		double getMutationPriorStandardDeviation();
		void setMutationPriorStandardDeviation(double _mutation_prior_sd);


		//Other functions:
		void getParameterForCategory(unsigned category, unsigned paramType, std::string aa, bool proposal, double *returnSet);
		void proposeHyperParameters();





		//R Section:

#ifndef STANDALONE
		//Constructors & Destructors:
		FONSEParameter(std::vector<double> stdDevSynthesisRate, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix,
		 				bool splitSer = true);
		FONSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
						bool splitSer = true, std::string _mutationSelectionState = "allUnique");



		//Initialization, Restart, Index Checking:
		void initCovarianceMatrix(SEXP matrix, std::string aa);
		void initMutation(std::vector <double> mutationValues, unsigned mixtureElement, std::string aa);
		void initSelection(std::vector <double> selectionValues, unsigned mixtureElement, std::string aa);
		void initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories,
						std::string paramType);


		//CSP Functions:
		std::vector< std::vector <double> > getCurrentMutationParameter();
		std::vector< std::vector <double> > getCurrentSelectionParameter();
		void setCurrentMutationParameter(std::vector<std::vector<double>> _currentMutationParameter);
		void setCurrentSelectionParameter(std::vector<std::vector<double>> _currentSelectionParameter);


#endif //STANDALONE

	protected:
};
#endif // FONSEPARAMETER_H
