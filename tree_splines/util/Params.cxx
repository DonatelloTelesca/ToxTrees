/*
 * observation.cpp
 *
 *  Created on: 13 juil. 2011
 *      Author
 */


#include "Params.h"
#include "Utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

using namespace std;
using namespace nsTreeSplines;
using namespace nsUtils;

#define PAR nsTreeSplines::Params
#define BNU boost::numeric::ublas

PAR::Params (void) throw ()
{
	nb_indiv = 0; 
	nb_dose = 0; 
	nb_repeat = 0; 
	nb_assay = 0; 
	nb_time = 0;
	nb_predictors = 0;
	nb_test_indiv = 0;
	nb_test_step = 0;
	max_depth = 0;
	delta = 0;
	xi = 0;
	nb_loop = 0;
	phi_d = 0;
	phi_t = 0;
	penalty_matrix_det = 0;
	predictors = NULL;
	alpha = 0.95;
	beta2 = 2;
	pgrow = 0.50;
	pprune = 0.25;
	pchange = 0.0;
	a1 = 0.0;
	b1 = 0.0;
	lambda1 = 0.0;
	lambda2 = 0.0;
	lambda3 = 0.0;
}

PAR::~Params() throw() 
{
	for (int i = 0; i < nb_indiv; i++)
		delete [] predictors[i];
	delete [] predictors;
	delete [] tau;
}

void PAR::ChooseSplitter (int & index, double & value) throw ()
{
	// Choose a split dimension
	double u = rand() / double(RAND_MAX);
	index =  (int)floor(u*nb_predictors);
	  
	// Choose split value
	u = rand() / double(RAND_MAX);
	int index2 = (int)floor(u*nb_indiv);
	value = predictors[index2][index];
	
	//index = 0;
	//value = 10;
}

// TODO this function speed can be improved, by using map structure for example
std::vector <int> PAR::GetAvailablePredictors (const std::vector <int> & indiv) throw ()
{
	std::vector <int> result;
	/*for (unsigned int i = 0; i < indiv.size(); i++)
	{
		vector <double> uniquesort;
		for (int j = 0; j < nb_predictors; j++)
			uniquesort.push_back (predictors[indiv[i]][j]);
		sort(uniquesort.begin(), uniquesort.end());
		vector<double>::iterator it = unique(uniquesort.begin(), uniquesort.end());
		uniquesort.resize( it - uniquesort.begin() );
		uniquesort.erase(uniquesort.begin()); // Erase the minimal value to be shure there is at leat one split

		if (uniquesort.size() > 1)
			result.push_back(i);
	}*/
	for (int i = 0; i < nb_predictors; i++)
	{
		vector <double> uniquesort;
		for (unsigned int j = 0; j < indiv.size(); j++)
			uniquesort.push_back (predictors[indiv[j]][i]);
		sort(uniquesort.begin(), uniquesort.end());
		vector<double>::iterator it = unique(uniquesort.begin(), uniquesort.end());
		uniquesort.resize( it - uniquesort.begin() );
		uniquesort.erase(uniquesort.begin()); // Erase the minimal value to be shure there is at leat one split

		if (uniquesort.size() > 0)
			result.push_back(i);
	}
	return result;
}

std::vector<double> PAR::GetAvailableValues (const int & index) throw ()
{
	std::vector<double> result;
	for (int j = 0; j < nb_indiv; j++)
		result.push_back (predictors[j][index]);
	sort (result.begin(), result.end());
	vector<double>::iterator it = unique (result.begin(), result.end());
	result.resize(it - result.begin());
	result.erase(result.begin());
	return result;
}

/*BNU::matrix <double> PAR::GetSigma_d_Inv () throw ()
{
	BNU::matrix <double> result (sigma_d.size1(), sigma_d.size2());
	for (unsigned int i = 0; i < result.size1(); i++)
		for (unsigned int j = 0; j < result.size2(); j++)
			result (i,j) = 0.0;

	double col0 = 1.0 / (1.0 - (phi_d * phi_d));
	double col1 = (phi_d * -1) / (1.0 - (phi_d * phi_d));
	double col2 = (1.0 + (phi_d * phi_d)) / (1.0 - (phi_d * phi_d));
	for (unsigned int i = 1; i < sigma_d.size1()-1; i++)
	{
		result (i, i-1) = col1;
		result (i, i) = col2;
		result (i, i+1) = col1;
	}
	result (0,0) = col0;
	result (0,1) = col1;
	result (sigma_d.size1()-1, sigma_d.size2()-1) = col0;
	result (sigma_d.size1()-1, sigma_d.size2()-2) = col1;
	return result;
}*/

std::ostream& nsTreeSplines::operator<< (std::ostream & os, nsTreeSplines::Params * params) throw ()
{   
	os << "nb_indiv = " << params->nb_indiv << endl;
	os << "nb_dose = " << params->nb_dose << endl;
	os << "nb_repeat = " << params->nb_repeat << endl;
	os << "nb_assay = " << params->nb_assay << endl;
	os << "nb_Time = " << params->nb_time << endl;
	os << "nb_predictors = " << params->nb_predictors << endl;
	os << "max depth = " << params->max_depth << endl;
	os << "delta = " << params->delta << endl;
	os << "xi = " << params->xi << endl;
	os << "nb_loop = " << params->nb_loop << endl;
	os << "phi_d = " << params->phi_d << endl;
	os << "a1 = " << params->a1 << endl;
	os << "b1 = " << params->b1 << endl;
	os << "aplha = " << params->alpha << endl;
	os << "beta = " << params->beta2 << endl;
	os << "w = " << params->w << endl;
	os << "pgrow pprune pchange = " << params->pgrow << " " << params->pprune << " " << params->pchange << endl;
	os << "lambda1 = " << params->lambda1 << endl;
	os << "lambda2 = " << params->lambda2 << endl;
	os << "lambda3 = " << params->lambda3 << endl;

	os << "\n------------- X -------------- " << endl;
	Display2DArray (params->predictors, params->nb_indiv, params->nb_predictors);

	os << "\n--------- penalty_matrix ----------- " << endl;
	os << params->penalty_matrix << endl;
	os << "Determinant " << params->penalty_matrix_det << endl;

	os << "\n--------- sigma_j ----------- " << endl;
	os << params->sigma_j << endl;
	os << "\n--------- sigma_d ----------- " << endl;
	os << params->sigma_d << endl;
	os << "\n--------- spline base matrix ----------- " << endl;
	os << params->spline_base_matrix << endl;
	os << "\n--------- lambda matrix ----------- " << endl;
	os << params->lambda << endl;
	return os;
}

#undef PAR
