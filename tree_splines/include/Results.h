/**
 * @File     : Result,h
 * @Authors  : C. Low-Kam, L. Di-Jorio
 * @Date     : 21/06/2012
 * @Version  : V1.0
 * @Synopsis : Structure storing and managing results
**/

#if !defined __RESULTS_H__
#define      __RESULTS_H__

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include <stdio.h>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#define BNU boost::numeric::ublas

namespace nsTreeSplines
{
    class Results
    {
		public :
			int nb_assay; // Number of assays (J)
			int nb_predictors; // Number of predictors
			int ** predictors; // Count of used predictors
			int ** current_predictors; // Count of used predictors
			std::map <double, int> ** predictor_values; // Count of 
			std::map <double, int> ** current_predictor_values; // Count of 
			std::vector < BNU::vector<double> > beta_matrixes; // Beta matrix
			std::vector < std::vector <double> > taus; // tau
			std::ofstream fileout; // file for storing result
			std::vector <double> phi_d;
			std::vector <double> phi_t;
			std::vector < BNU::matrix <double> > sigmas;
			std::vector < BNU::matrix <double> > s_res;
			std::vector < std::vector <int> > nb_sons;

			Results (const int & assay, const int & _predictors, const char * _fileout) throw ();
			~Results (void) throw ();

			// -------------------------------------- //
			//                METHODS                 //
			// -------------------------------------- //
			void FlushMatrixes () throw ();
			void IncreaseAll (const int & assay) throw ();
			void AddAll (const int & assay) throw ();
			void CheckResults () throw ();
	};

	std::ostream& operator<< (std::ostream & os, nsTreeSplines::Results * result) throw ();
}

#undef BNU
#endif /* __PARAMS_H__ */
