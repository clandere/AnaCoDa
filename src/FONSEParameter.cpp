#include "include/FONSE/FONSEParameter.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


FONSEParameter::FONSEParameter() : Parameter()
{
	//CTOR
	bias_csp = 0;
	mutation_prior_sd = 0.35;
	currentCodonSpecificParameter.resize(2);
	proposedCodonSpecificParameter.resize(2);
}


FONSEParameter::FONSEParameter(std::string filename) : Parameter(22)
{
	currentCodonSpecificParameter.resize(2);
	proposedCodonSpecificParameter.resize(2);
	initFromRestartFile(filename);
}


FONSEParameter::FONSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
	std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer, std::string _mutationSelectionState) :
	Parameter(22)
{
	initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initFONSEParameterSet();
}


FONSEParameter& FONSEParameter::operator=(const FONSEParameter& rhs)
{
	if (this == &rhs)
		return *this; // handle self assignment

	Parameter::operator=(rhs);

	// proposal bias and std for codon specific parameter
	bias_csp = rhs.bias_csp;
	std_csp = rhs.std_csp;

	mutation_prior_sd = rhs.mutation_prior_sd;

	return *this;
}


FONSEParameter::FONSEParameter(const FONSEParameter &other) : Parameter(other)
{
	bias_csp = other.bias_csp;
	std_csp = other.std_csp;

	mutation_prior_sd = other.mutation_prior_sd;

}


FONSEParameter::~FONSEParameter()
{
	// destructor
}





//---------------------------------------------------------------//
// ---------- Initialization, Restart, Index Checking -----------//
//---------------------------------------------------------------//


void FONSEParameter::initFONSEParameterSet()
{
	mutation_prior_sd = 0.35;
	groupList = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y", "Z" };
	// proposal bias and std for codon specific parameter
	bias_csp = 0;
	std_csp.resize(numParam, 0.1);
	//may need getter fcts

	currentCodonSpecificParameter.resize(2);
	proposedCodonSpecificParameter.resize(2);

	currentCodonSpecificParameter[dM].resize(numMutationCategories);
	proposedCodonSpecificParameter[dM].resize(numMutationCategories);

	currentCodonSpecificParameter[dOmega].resize(numSelectionCategories);
	proposedCodonSpecificParameter[dOmega].resize(numSelectionCategories);

	unsigned biggestCat = std::max(std::max(numMutationCategories, numSelectionCategories), maxGrouping);
	for(unsigned i = 0u; i < biggestCat; i++)
	{
		std::vector<double> tmp(numParam, 0.0);
		if(i < numMutationCategories)
		{
			currentCodonSpecificParameter[dM][i] = tmp;
			proposedCodonSpecificParameter[dM][i] = tmp;
		}
		if(i < numSelectionCategories)
		{
			proposedCodonSpecificParameter[dOmega][i] = tmp;
			currentCodonSpecificParameter[dOmega][i] = tmp;
		}
		if(i < maxGrouping)
		{
			std::string aa = SequenceSummary::AminoAcidArray[i];
			//TODO: Explain this
			unsigned numCodons = SequenceSummary::GetNumCodonsForAA(aa, true);
			CovarianceMatrix m((numMutationCategories + numSelectionCategories) * numCodons);
			m.choleskyDecomposition();
			covarianceMatrix.push_back(m);
		}
	}
}


void FONSEParameter::initFONSEValuesFromFile(std::string filename)
{
	std::ifstream input;
	std::vector <double> mat;
	covarianceMatrix.resize(maxGrouping);
	input.open(filename.c_str());
	if (input.fail())
		my_printError("ERROR: Could not open RestartFile.txt to initialize FONSE values\n");
	else
	{
		std::string tmp, variableName;
		unsigned cat = 0;
		while (getline(input, tmp))
		{
			int flag;
			if (tmp[0] == '>')
				flag = 1;
			else if (input.eof() || tmp == "\n")
				flag = 2;
			else if (tmp[0] == '#')
				flag = 3;
			else
				flag = 4;

			if (flag == 1)
			{
				mat.clear();
				cat = 0;
				variableName = tmp.substr(1, tmp.size() - 2);
				if (variableName == "covarianceMatrix")
				{
					getline(input, tmp);
					//char aa = tmp[0];
					cat = SequenceSummary::AAToAAIndex(tmp); //getting the index value of the amino acid
				}
			}
			else if (flag == 2)
			{
				my_print("here\n");
			}
			else if (flag == 3) //user comment, continue
			{
				continue;
			}
			else
			{
				std::istringstream iss;
				if (variableName == "currentMutationParameter")
				{
					if (tmp == "***")
					{
						currentCodonSpecificParameter[dM].resize(currentCodonSpecificParameter[dM].size() + 1);
						cat++;
					}
					else if (tmp == "\n")
						continue;
					else
					{
						double val;
						iss.str(tmp);
						while (iss >> val)
						{
							currentCodonSpecificParameter[dM][cat - 1].push_back(val);
						}
					}
				}
				else if (variableName == "currentSelectionParameter")
				{
					if (tmp == "***")
					{
						currentCodonSpecificParameter[dOmega].resize(currentCodonSpecificParameter[dOmega].size() + 1);
						cat++;
					}
					else if (tmp == "\n")
						continue;
					else
					{
						double val;
						iss.str(tmp);
						while (iss >> val)
						{
							currentCodonSpecificParameter[dOmega][cat - 1].push_back(val);
						}
					}
				}
				else if (variableName == "covarianceMatrix")
				{
					if (tmp == "***") //end of matrix
					{
						CovarianceMatrix CM(mat);
						CM.choleskyDecomposition();//Solving a system of equations
						covarianceMatrix[cat] = CM;
					}
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						mat.push_back(val);
					}
				}
				else if (variableName == "std_csp")
				{
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						std_csp.push_back(val);
					}
				}
				else if (variableName == "mutation_prior_sd")
				{
					iss.str(tmp);
					iss >> mutation_prior_sd;
				}
			}
		}
	}
	input.close();

	//init other values
	bias_csp = 0;
	proposedCodonSpecificParameter[dM].resize(numMutationCategories);
	//looping through the bigger of the two categories
	unsigned biggerCat = std::max(numMutationCategories, numSelectionCategories);
	for (unsigned i = 0; i < biggerCat; i++)
	{
		//making sure not going out of bounds on either of the vectors
		if(i < numMutationCategories)
			proposedCodonSpecificParameter[dM][i] = currentCodonSpecificParameter[dM][i];
		if(i < numSelectionCategories)
			proposedCodonSpecificParameter[dOmega][i] = currentCodonSpecificParameter[dOmega][i];
	}

	groupList = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y", "Z" };
	//groupList = { "C", "D", "E", "F", "H", "K", "M", "N", "Q", "W", "Y" };
}


void FONSEParameter::writeEntireRestartFile(std::string filename)
{
	writeBasicRestartFile(filename);
	writeFONSERestartFile(filename);
}


void FONSEParameter::writeFONSERestartFile(std::string filename)
{
  std::ofstream out(filename.c_str(), std::ofstream::app);
  //out.open(filename.c_str(), std::ofstream::app);
  if (out.fail())
      my_printError("ERROR: Could not open RestartFile.txt to append\n");
	else
	{
		std::ostringstream oss;
		unsigned j;
		oss << ">mutation_prior_sd:\n" << mutation_prior_sd << "\n";
		oss << ">std_csp:\n";
		for (unsigned i = 0; i < std_csp.size(); i++)
		{
			oss << std_csp[i];
			if ((i + 1) % 10 == 0)
				oss << "\n";
			else
				oss << " ";
		}
		oss << ">currentMutationParameter:\n";
		for (unsigned i = 0; i < currentCodonSpecificParameter[dM].size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < currentCodonSpecificParameter[dM][i].size(); j++)
			{
				oss << currentCodonSpecificParameter[dM][i][j];
				if ((j + 1) % 10 == 0)
					oss << "\n";
				else
					oss << " ";
			}
			if (j % 10 != 0)
				oss << "\n";
		}

		oss << ">currentSelectionParameter:\n";
		for (unsigned i = 0; i < currentCodonSpecificParameter[dOmega].size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < currentCodonSpecificParameter[dOmega][i].size(); j++)
			{
				oss << currentCodonSpecificParameter[dOmega][i][j];
				if ((j + 1) % 10 == 0)
					oss << "\n";
				else
					oss << " ";
			}
			if (j % 10 != 0)
				oss << "\n";
		}
		for (unsigned i = 0; i < groupList.size(); i++)
		{
			std::string aa = groupList[i];
			oss << ">covarianceMatrix:\n" << aa << "\n";
			CovarianceMatrix m = covarianceMatrix[SequenceSummary::AAToAAIndex(aa)];
			std::vector<double>* tmp = m.getCovMatrix();
			int size = m.getNumVariates();
			for (unsigned k = 0; k < size * size; k++)
			{
				if (k % size == 0 && k != 0) { oss << "\n"; }
				oss << tmp->at(k) << "\t";
			}
			oss << "\n***\n";
		}
		out << oss.str();
	}
    out.close();
}


void FONSEParameter::initFromRestartFile(std::string filename)
{
    initBaseValuesFromFile(filename);
    initFONSEValuesFromFile(filename);
}


void FONSEParameter::initAllTraces(unsigned samples, unsigned num_genes)
{
    traces.initializeFONSETrace(samples, num_genes, numMutationCategories, numSelectionCategories, numParam,
                         numMixtures, categories, maxGrouping,currentSynthesisRateLevel[0],mixtureAssignment);
}

void FONSEParameter::initMutationCategories(std::vector<std::string> files, unsigned numCategories)
{
    for (unsigned category = 0; category < numCategories; category++)
    {
      //Open the file for the category
      std::ifstream currentFile;
      currentFile.open(files[category].c_str());
      if (currentFile.fail())
        my_printError("Error opening file % to initialize mutation values.\n", category);
			else
			{
				std::string tmp;
				currentFile >> tmp; //The first line is a header (Amino Acid, Codon, Value, Std_deviation)

				while (currentFile >> tmp)
				{
					//Get the Codon and Index
					std::size_t pos = tmp.find(',', 2); //Amino Acid and a comma will always be the first 2 characters
					std::string codon = tmp.substr(2, pos - 2);
					unsigned codonIndex = SequenceSummary::codonToIndex(codon, true);

					//get the value to store
					std::size_t pos2 = tmp.find(',', pos + 1);
					//my_print("%\n", tmp.substr(pos + 1, pos2 - pos - 1 ));
					double value = std::atof(tmp.substr(pos + 1, pos2 - pos - 1).c_str());

					currentCodonSpecificParameter[dM][category][codonIndex] = value;
					proposedCodonSpecificParameter[dM][category][codonIndex] = value;
				}
			}
      currentFile.close();
    } //END OF A CATEGORY/FILE
}


void FONSEParameter::initSelectionCategories(std::vector<std::string> files, unsigned numCategories)
{
  for (unsigned category = 0; category < numCategories; category++)
  {
	  //Open the file for the category
	  std::ifstream currentFile;
	  currentFile.open(files[category].c_str());
	  if (currentFile.fail())
	      my_printError("Error opening file % to initialize mutation values.\n", category);
		else
		{
			std::string tmp;
			currentFile >> tmp; //The first line is a header (Amino Acid, Codon, Value, Std_deviation)

			while (currentFile >> tmp)
			{
				//Get the Codon and Index
				std::size_t pos = tmp.find(',', 2); //Amino Acid and a comma will always be the first 2 characters
				std::string codon = tmp.substr(2, pos - 2);
				unsigned codonIndex = SequenceSummary::codonToIndex(codon, true);

				//get the value to store
				std::size_t pos2 = tmp.find(',', pos + 1);
				//	my_print("%\n", tmp.substr(pos + 1, pos2 - pos - 1 ));
				double value = std::atof(tmp.substr(pos + 1, pos2 - pos - 1).c_str());

				currentCodonSpecificParameter[dOmega][category][codonIndex] = value;
				proposedCodonSpecificParameter[dOmega][category][codonIndex] = value;
			}
		}
    currentFile.close();
  } //END OF A CATEGORY/FILE
}






// --------------------------------------//
// ---------- Trace Functions -----------//
// --------------------------------------//


void FONSEParameter::updateCodonSpecificParameterTrace(unsigned sample, std::string grouping)
{
	traces.updateCodonSpecificParameterTraceForAA(sample, grouping, currentCodonSpecificParameter[dM], dM);
	traces.updateCodonSpecificParameterTraceForAA(sample, grouping, currentCodonSpecificParameter[dOmega], dOmega);
}


// ------------------------------------------//
// ---------- Covariance Functions ----------//
// ------------------------------------------//


CovarianceMatrix& FONSEParameter::getCovarianceMatrixForAA(std::string aa)
{
  aa[0] = (char)std::toupper(aa[0]);
  unsigned aaIndex = SequenceSummary::aaToIndex.find(aa)->second;
  return covarianceMatrix[aaIndex];
}


// -----------------------------------//
// ---------- CSP Functions ----------//
// -----------------------------------//


double FONSEParameter::getCurrentCodonSpecificProposalWidth(unsigned aa)
{
    unsigned aaStart, aaEnd;
    //Gets the codon range based on the Amino Acid
    SequenceSummary::AAIndexToCodonRange(aa, aaStart, aaEnd, false);
    return std_csp[aaStart];
}


void FONSEParameter::proposeCodonSpecificParameter()
{
  for (unsigned k = 0; k < getGroupListSize(); k++)
  {
    std::vector<double> iidProposed;
    std::string aa = getGrouping(k);
	
    unsigned aaStart, aaEnd;
    SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
    unsigned numCodons = aaEnd - aaStart;
    for (unsigned i = 0u; i < (numCodons * (numMutationCategories + numSelectionCategories)); i++)
    {
      iidProposed.push_back(randNorm(0.0, 1.0)); //Random distribution between 0 and 1
    }

    std::vector<double> covaryingNums;
	//TODO: Explain the following line
    covaryingNums = covarianceMatrix[SequenceSummary::AAToAAIndex(aa)].transformIidNumbersIntoCovaryingNumbers(iidProposed);
		unsigned biggestCat = std::max(numMutationCategories, numSelectionCategories);

		for (unsigned i = 0; i < biggestCat; i++)
    {
      for (unsigned j = i * numCodons, l = aaStart; j < (i * numCodons) + numCodons; j++, l++)
      {
				if(i < numMutationCategories)
        	proposedCodonSpecificParameter[dM][i][l] = currentCodonSpecificParameter[dM][i][l] + covaryingNums[j];

				if(i < numSelectionCategories)
					proposedCodonSpecificParameter[dOmega][i][l] = currentCodonSpecificParameter[dOmega][i][l]
					+ covaryingNums[(numMutationCategories * numCodons) + j];
      }
    }
  }
}


void FONSEParameter::updateCodonSpecificParameter(std::string grouping)
{
	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(grouping, aaStart, aaEnd, true);
    unsigned aaIndex = SequenceSummary::aaToIndex.find(grouping)->second;
	
    numAcceptForCodonSpecificParameters[aaIndex]++;

    unsigned biggestCat = std::max(numMutationCategories, numSelectionCategories);
    for (unsigned k = 0u; k < biggestCat; k++)
    {
        for (unsigned i = aaStart; i < aaEnd; i++)
        {
		    if (i < numMutationCategories)
      	        currentCodonSpecificParameter[dM][k][i] = proposedCodonSpecificParameter[dM][k][i];
			if (i < numSelectionCategories)
			    currentCodonSpecificParameter[dOmega][k][i] = proposedCodonSpecificParameter[dOmega][k][i];
        }
    }
}





// -------------------------------------//
// ---------- Prior Functions ----------//
// -------------------------------------//
double FONSEParameter::getMutationPriorStandardDeviation()
{
	return mutation_prior_sd;
}


void FONSEParameter::setMutationPriorStandardDeviation(double _mutation_prior_sd)
{
	mutation_prior_sd = _mutation_prior_sd;
}




// -------------------------------------//
// ---------- Other Functions ----------//
// -------------------------------------//

void FONSEParameter::getParameterForCategory(unsigned category, unsigned paramType, std::string aa, bool proposal,
                                             double *returnSet)
{
    std::vector<double> *tempSet;
	tempSet = proposal ? &proposedCodonSpecificParameter[paramType][category] : &currentCodonSpecificParameter[paramType][category];
	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);

    unsigned j = 0u;
    for (unsigned i = aaStart; i < aaEnd; i++, j++)
    {
        returnSet[j] = tempSet->at(i);
    }
}


void FONSEParameter::proposeHyperParameters()
{
    for (unsigned i = 0u; i < numSelectionCategories; i++)
    {
		stdDevSynthesisRate_proposed[i] = std::exp(randNorm(std::log(stdDevSynthesisRate[i]), std_stdDevSynthesisRate));
    }
}








// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//


#ifndef STANDALONE


//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


FONSEParameter::FONSEParameter(std::vector<double> stdDevSynthesisRate, std::vector<unsigned> geneAssignment,
                               std::vector<unsigned> _matrix, bool splitSer) : Parameter(22)
{
    unsigned _numMixtures = _matrix.size() / 2;
    std::vector<std::vector<unsigned>> thetaKMatrix;
    thetaKMatrix.resize(_numMixtures);

    unsigned index = 0;
    for (unsigned i = 0; i < _numMixtures; i++)
    {
        for (unsigned j = 0; j < 2; j++, index++)
        {
            thetaKMatrix[i].push_back(_matrix[index]);
        }
    }
    initParameterSet(stdDevSynthesisRate, _matrix.size() / 2, geneAssignment, thetaKMatrix, splitSer);
    initFONSEParameterSet();

}


FONSEParameter::FONSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
                               bool splitSer, std::string _mutationSelectionState) : Parameter(22)
{
    std::vector<std::vector<unsigned>> thetaKMatrix;
    initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
    initFONSEParameterSet();
}





//---------------------------------------------------------------//
// ---------- Initialization, Restart, Index Checking ---------- //
//---------------------------------------------------------------//


void FONSEParameter::initCovarianceMatrix(SEXP _matrix, std::string aa)
{
    std::vector<double> tmp;
    NumericMatrix matrix(_matrix);

    for (unsigned i = 0u; i < aa.length(); i++)	aa[i] = (char)std::toupper(aa[i]);

    unsigned aaIndex = SequenceSummary::aaToIndex.find(aa)->second;
    unsigned numRows = matrix.nrow();
    std::vector<double> covMatrix(numRows * numRows);

    //NumericMatrix stores the matrix by column, not by row. The loop
    //below transposes the matrix when it stores it.
    unsigned index = 0;
    for (unsigned i = 0; i < numRows; i++)
    {
        for (unsigned j = i; j < numRows * numRows; j += numRows, index++)
        {
            covMatrix[index] = matrix[j];
        }
    }
    CovarianceMatrix m(covMatrix);
    m.choleskyDecomposition();
    covarianceMatrix[aaIndex] = m;
}

void FONSEParameter::initMutation(std::vector<double> mutationValues, unsigned mixtureElement, std::string aa)
{
    //TODO: seperate out the R wrapper functionality and make the wrapper
    //currentMutationParameter
    bool check = checkIndex(mixtureElement, 1, numMixtures);
    if (check)
    {
        mixtureElement--;

        unsigned category = getMutationCategory(mixtureElement);
        aa[0] = (char)std::toupper(aa[0]);
	    
        unsigned aaStart, aaEnd;
	    SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
        for (unsigned i = aaStart, j = 0; i < aaEnd; i++, j++)
        {
            currentCodonSpecificParameter[dM][category][i] = mutationValues[j];
        }
    }
}


void FONSEParameter::initSelection(std::vector<double> selectionValues, unsigned mixtureElement, std::string aa)
{
    //TODO: seperate out the R wrapper functionality and make the wrapper
    //currentSelectionParameter
    bool check = checkIndex(mixtureElement, 1, numMixtures);
    if (check)
    {
        mixtureElement--;

        int category = getSelectionCategory(mixtureElement);

        aa[0] = (char)std::toupper(aa[0]);
	    
        unsigned aaStart, aaEnd;
        SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
        for (unsigned i = aaStart, j = 0; i < aaEnd; i++, j++)
        {
            currentCodonSpecificParameter[dOmega][category][i] = selectionValues[j];
        }
    }
}


// -----------------------------------//
// ---------- CSP Functions ----------//
// -----------------------------------//


std::vector< std::vector <double> > FONSEParameter::getCurrentMutationParameter()
{
    return currentCodonSpecificParameter[dM];
}


std::vector< std::vector <double> > FONSEParameter::getCurrentSelectionParameter()
{
    return currentCodonSpecificParameter[dOmega];
}


void FONSEParameter::setCurrentMutationParameter(std::vector<std::vector<double>> _currentMutationParameter)
{
	currentCodonSpecificParameter[dM] = _currentMutationParameter;
}


void FONSEParameter::setCurrentSelectionParameter(std::vector<std::vector<double>> _currentSelectionParameter)
{
	currentCodonSpecificParameter[dOmega] = _currentSelectionParameter;
}
#endif
