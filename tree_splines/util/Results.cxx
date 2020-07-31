/*
 * observation.cpp
 *
 *  Created on: 13 juil. 2011
 *      Author
 */

#include "Results.h"
#include "Utils.h"

#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace nsTreeSplines;
using namespace nsUtils;

#define RES nsTreeSplines::Results
#define BNU boost::numeric::ublas

RES::Results (const int & assay, const int & _predictors, const char * _fileout) throw ()
{
	nb_assay = assay; 
	nb_predictors = _predictors;

	predictors = new int * [assay];
	current_predictors = new int * [assay];
	predictor_values = new map<double, int> * [assay];
	current_predictor_values = new map<double, int> * [assay];
	for (int i = 0; i < assay; i++)
	{
		predictors[i] = new int [_predictors];
		current_predictors[i] = new int [_predictors];
		for (int j = 0; j < _predictors; j++)
		{
			current_predictors[i][j] = 0;
			predictors[i][j] = 0;
		}

		predictor_values [i] = new map<double, int>  [_predictors];
		current_predictor_values [i] = new map<double, int>  [_predictors];
	}

	fileout.open(_fileout);
}

RES::~Results() throw() 
{
	FlushMatrixes ();
	beta_matrixes.clear();

	// Write predictor count
	for (int i = 0; i < nb_assay; i++)
	{
		for (int j = 0; j < nb_predictors; j++)
			fileout << (j == 0 ? "" : " ") << predictors[i][j];
		fileout << endl;
	}

	// Write values count
	for (int i = 0; i < nb_assay; i++)
	{
		for (int j = 0; j < nb_predictors; j++)
		{
			map <double, int>::iterator it;
			for (it = predictor_values[i][j].begin(); it != predictor_values[i][j].end(); it++)
				fileout << "(" << (*it).first << "," << (*it).second << ")";
			fileout << endl;
		}
	}

	// Write tau
	for (unsigned int i = 0; i < taus.size(); i++)
	{
		for (unsigned int j = 0; j < taus[i].size(); j++)
			fileout << (j == 0 ? "" : " ") << taus[i][j];
		fileout << endl;
	}

	// Write phi_d
	for (unsigned int i = 0; i < phi_d.size(); i++)
		fileout << (i == 0 ? "" : " ") << phi_d[i];
	fileout << endl;
	
	// Write phi_t
	for (unsigned int i = 0; i < phi_t.size(); i++)
		fileout << (i == 0 ? "" : " ") << phi_t[i];
	fileout << endl;

	// Write sigma_j
	for (unsigned int i = 0; i < sigmas.size (); i++)
		fileout << sigmas[i] << endl;
	
	// Write s_res
	for (unsigned int i = 0; i < s_res.size (); i++)
		fileout << s_res[i] << endl;

	// Write nb_sons
	for (unsigned int i = 0; i < nb_sons.size(); i++)
	{
		for (unsigned int j = 0; j < nb_sons[i].size(); j++)
			fileout << (j == 0 ? "" : " ") << nb_sons[i][j];
		fileout << endl;
	}
	for (int i = 0; i < nb_assay; i++)
	{
		delete [] predictors[i];
		delete [] current_predictors[i];
		delete [] predictor_values [i];
		delete [] current_predictor_values [i];
	}

	delete [] predictors;
	delete [] current_predictors;
	delete [] predictor_values;
	delete [] current_predictor_values;
	fileout.close();
}

void RES::FlushMatrixes () throw ()
{
	for (unsigned int i = 0; i < beta_matrixes.size(); i++)
	{
		for (unsigned int j = 0; j < beta_matrixes[i].size(); j++)
			fileout << (j == 0 ? "" : " " ) << beta_matrixes[i][j];
		fileout << endl;
	}
	beta_matrixes.clear ();
}

void RES::IncreaseAll (const int & assay) throw ()
{
	for (int i = 0; i < nb_predictors; i++)
	{
		if (predictors[assay][i] > 0)
			predictors[assay][i]++;
		map<double,int>::iterator it;
		for ( it= predictor_values[assay][i].begin() ; it != predictor_values[assay][i].end(); it++ )
			if ((*it).second > 0)
				(*it).second++;
	}
}

void RES::AddAll (const int & assay) throw ()
{
	for (int i = 0; i < nb_predictors; i++)
	{
		predictors[assay][i] += current_predictors[assay][i];
		map<double,int>::iterator it;
		for (it = current_predictor_values[assay][i].begin() ; it != current_predictor_values[assay][i].end(); it++ )
		{
			if (predictor_values[assay][i].find ((*it).first) == predictor_values[assay][i].end())
				predictor_values[assay][i][(*it).first] = 0;
			predictor_values[assay][i][(*it).first] += (*it).second;
		}
	}
}

void RES::CheckResults () throw ()
{
	for (int i = 0; i < nb_assay; i++)
	{
		for (int j = 0; j < nb_predictors; j++)
		{
			int value_total = 0;
			map <double, int>::iterator it;
			for (it = predictor_values[i][j].begin(); it != predictor_values[i][j].end(); it++)
				value_total += (*it).second;
			if (value_total != predictors[i][j])
				cout << "Incorrect value assay " << i << " predictor " << j << endl;
		}
	}
}

std::ostream& nsTreeSplines::operator<< (std::ostream & os, nsTreeSplines::Results * results) throw ()
{   
	os << "\n------------- Predictors -------------- " << endl;
	Display2DArray (results->predictors, results->nb_assay, results->nb_predictors);

	os << "\n------------- Predictors Values ------------- " << endl;
	for (int i = 0; i < results->nb_assay; i++)
	{
		for (int j = 0; j < results->nb_predictors; j++)
		{
			cout << "PV [" << i << "," << j << "] -> ";
			map <double, int>::iterator it;
			for (it = results->predictor_values[i][j].begin(); it != results->predictor_values[i][j].end(); it++)
				cout << ", (" << (*it).first << "," << (*it).second << ")";
			cout << endl;
		}
	}

	return os;
}

#undef RES
