#ifndef PANSEMODEL_H
#define PANSEMODEL_H


#include "../base/Model.h"
#include "PANSEParameter.h"
#include "../SequenceSummary.h"

#include <sstream>

class PANSEModel: public Model
{
	private:
		PANSEParameter *parameter;
		unsigned RFPCountColumn;

		double calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime,
				unsigned currRFPObserved, unsigned currNumCodonsInMRNA, double phiValue);


	public:
		//Constructors & Destructors:
		explicit PANSEModel(unsigned RFPCountColumn = 0u);
		virtual ~PANSEModel();


		//Likelihood Ratio Functions:
		virtual void calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k,
				double* logProbabilityRatio); // Depends on RFPCountColumn
		virtual void calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome,
				std::vector<double> &logAcceptanceRatioForAllMixtures); // Depends on RFPCountColumn
		virtual void calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration,
				std::vector <double> &logProbabilityRatio);


		//Initialization and Restart Functions:
		virtual void initTraces(unsigned samples, unsigned num_genes);
		virtual void writeRestartFile(std::string filename);


		//Category Functions:
		virtual double getCategoryProbability(unsigned i);
		virtual unsigned getMutationCategory(unsigned mixture);
		virtual unsigned getSelectionCategory(unsigned mixture);
		virtual unsigned getSynthesisRateCategory(unsigned mixture);
		virtual std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned k);


		//Group List Functions:
		virtual unsigned getGroupListSize();
		virtual std::string getGrouping(unsigned index);


		//stdDevSynthesisRate Functions:
		virtual double getStdDevSynthesisRate(unsigned selectionCategory, bool proposed = false);
		virtual double getCurrentStdDevSynthesisRateProposalWidth();
		virtual void updateStdDevSynthesisRate();


		//Synthesis Rate Functions:
		virtual double getSynthesisRate(unsigned index, unsigned mixture, bool proposed = false);
		virtual void updateSynthesisRate(unsigned i, unsigned k);


		//Iteration Functions:
		virtual unsigned getLastIteration();
		virtual void setLastIteration(unsigned iteration);


		//Trace Functions:
		virtual void updateStdDevSynthesisRateTrace(unsigned sample);
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned i);
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned i);
		virtual void updateMixtureProbabilitiesTrace(unsigned sample);
		virtual void updateCodonSpecificParameterTrace(unsigned sample, std::string codon);
		virtual void updateHyperParameterTraces(unsigned sample);
		virtual void updateTracesWithInitialValues(Genome &genome);


		//Adaptive Width Functions:
		virtual void adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt = true);
		virtual void adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt = true);


		//Other Functions:
		virtual void proposeCodonSpecificParameter();
		virtual void proposeHyperParameters();
		virtual void proposeSynthesisRateLevels();

		virtual unsigned getNumPhiGroupings();
		virtual unsigned getMixtureAssignment(unsigned index);
		virtual unsigned getNumMixtureElements();
		virtual unsigned getNumSynthesisRateCategories();

		virtual void setNumPhiGroupings(unsigned value);
		virtual void setMixtureAssignment(unsigned i, unsigned catOfGene);
		virtual void setCategoryProbability(unsigned mixture, double value);

		virtual void updateCodonSpecificParameter(std::string aa);
		virtual void updateGibbsSampledHyperParameters(Genome &genome);
		virtual void updateAllHyperParameter();
		virtual void updateHyperParameter(unsigned hp);

		virtual void simulateGenome(Genome &genome); // Depends on RFPCountColumn
		virtual void printHyperParameters();
		PANSEParameter* getParameter();

		void setParameter(PANSEParameter &_parameter);
		virtual double calculateAllPriors();
		virtual double getParameterForCategory(unsigned category, unsigned param, std::string codon, bool proposal);
		double u_gamma_helper(double s, double x);

		double u_gamma(double s, double x);
		double u_gamma_log(double s, double x);
        double generalized_integral(double p, double z);

        double generalized_integral_log(double p, double z);
        double prob_elongation(double curralpha, double currlambda, double currv);
        double prob_elongation_log(double curralpha, double currlambda, double currv);
        double delta_g(int i, int g, double *lambda, double *v_g, double *alpha);
        
        double delta_g_log(int i, int g, double *lambda, double *v_g, double *alpha);
        double prob_Y_g(double curralpha, int sample_size, double lambda_prime, double psi, double prevdelta);
        double prob_Y_g_log(double curralpha, int sample_size, double lambda_prime, double psi, double prevdelta);
        
        //Psi-Phi Conversion Functions
        double psi2phi(double psi, double sigma);
        double phi2psi(double phi, double sigma);

	protected:
};

#endif // PANSEMODEL_H
