/*
 * observation.cpp
 *
 *  Created on: 13 juil. 2011
 *      Author
 */

#include "Utils.h"
#include "Cholesky.hpp"

#include <string.h>
#include <iostream>
#include <iomanip>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/bindings/lapack/gesvd.hpp>
//#include "boost/numeric/ublas/matrix_proxy.hpp"
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

using namespace std;
using namespace nsUtils;
using namespace boost::numeric::ublas;

#define NU nsUtils
#define BNU boost::numeric::ublas;

#ifndef BOOST_UBLAS_TYPENAME
#define BOOST_UBLAS_TYPENAME typename
#endif

static boost::mt19937 gen (time(0));
static boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator_norm (boost::mt19937(time(0)),boost::normal_distribution<>());

std::vector<string> NU::Split (string splitter, string str) throw ()
{
	std::vector<string> result;
	result.push_back ("");
	char s = splitter[0];

	for (unsigned int i = 0; i < str.length(); i++)
	{
		char c = str[i];
		if (c == s)
			result.push_back("");
		else
			result[result.size()-1] += str[i];
	}
	return result;
}

void NU::Split (string splitter, string str, int * result) throw ()
{
	string tmp = "";
	char s = splitter[0];
	int index = 0;
	for (unsigned int i = 0; i < str.length(); i++)
	{
		char c = str[i];
		if (c == s)
		{
			result[index] = atoi(tmp.c_str());
			tmp = "";
			index++;
		}
		else
			tmp += str[i];
	}
	result[index] = atoi(tmp.c_str());
}

void NU::Split (string splitter, string str, double * result) throw ()
{
	string tmp = "";
	char s = splitter[0];
	int index = 0;
	for (unsigned int i = 0; i < str.length(); i++)
	{
		char c = str[i];
		if (c == s)
		{
			result[index] = atof(tmp.c_str());
			tmp = "";
			index++;
		}
		else
			tmp += str[i];
	}
	result[index] = atof(tmp.c_str());
}

void NU::Split (std::string splitter, std::string str, mapped_matrix<double> & result, int & line) throw ()
{
	string tmp = "";
	char s = splitter[0];
	int index = 0;
	for (unsigned int i = 0; i < str.length(); i++)
	{
		char c = str[i];
		if (c == s)
		{
			result(line, index) = atof(tmp.c_str());
			tmp = "";
			index++;
		}
		else
			tmp += str[i];
	}
	result(line, index) = atof(tmp.c_str());
}

void NU::Split (std::string splitter, std::string str, symmetric_matrix<double, lower> & result, int & line) throw ()
{
	string tmp = "";
	char s = splitter[0];
	int index = 0;
	for (unsigned int i = 0; i < str.length(); i++)
	{
		char c = str[i];
		if (c == s)
		{
			result (line, index) = atof (tmp.c_str());
			tmp = "";
			index++;
			if (index > line)
				return;
		}
		else
			tmp += str[i];
	}
	result(line, index) = atof(tmp.c_str());
}

/*void NU::Split (std::string splitter, std::string str, diagonal_matrix<double, lower> & result, int & line) throw ()
{
	string tmp = "";
	char s = splitter[0];
	int index = 0;
	for (unsigned int i = 0; i < str.length(); i++)
	{
		char c = str[i];
		if (c == s)
		{
			result (line, index) = atof (tmp.c_str());
			tmp = "";
			index++;
			if (index > line)
				return;
		}
		else
			tmp += str[i];
	}
	result(line, index) = atof(tmp.c_str());
}*/

void NU::Display2DArray (int ** array, int _i, int _j) throw ()
{
	for (int i = 0; i < _i; i++)
	{
		for (int j = 0; j < _j; j++)
			cout << setw(10) << array[i][j];
		cout << endl;
	}
}

void NU::Display2DArray (double ** array, int _i, int _j) throw ()
{
	for (int i = 0; i < _i; i++)
	{
		for (int j = 0; j < _j; j++)
			cout << setw(10) << array[i][j];
		cout << endl;
	}
}

int NU::DeterminantSign(const permutation_matrix<std ::size_t>& pm) throw ()
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}
 
double NU::Determinant(mapped_matrix<double>& m) throw ()
{
	mapped_matrix<double> copied_m (m);
    permutation_matrix<std ::size_t> pm (copied_m.size1());
    double det = 1.0;
    if(lu_factorize (copied_m, pm)) 
        det = 0.0;
    else 
	{
        for(unsigned int i = 0; i < m.size1(); i++)
            det *= copied_m(i,i); // multiply by elements on diagonal
        det = det * DeterminantSign(pm);
    }
    return det;
}

double NU::Determinant(matrix<double>& m) throw ()
{
	mapped_matrix<double> copied_m (m);
    permutation_matrix<std ::size_t> pm (copied_m.size1());
    double det = 1.0;
    if(lu_factorize (copied_m, pm)) 
        det = 0.0;
    else 
	{
        for(unsigned int i = 0; i < m.size1(); i++)
            det *= copied_m(i,i); // multiply by elements on diagonal
        det = det * DeterminantSign(pm);
    }
    return det;
}

bool NU::InvertMatrix(const matrix<double> & input, matrix<double>& inverse) throw ()
{
	typedef permutation_matrix<std::size_t> pmatrix;

	matrix<double> A(input); // create a working copy of the input
	pmatrix pm(A.size1()); // create a permutation matrix for the LU-factorization

	int res = lu_factorize(A, pm); // perform LU-factorization
	if (res != 0)
		return false;
	
	inverse.assign(identity_matrix<double> (A.size1())); // create identity matrix of "inverse"
	lu_substitute(A, pm, inverse); // backsubstitute to get the inverse

	return true;
}

bool NU::InvertMatrix(symmetric_matrix<double, lower>& input, symmetric_matrix<double, lower>& inverse)
{
	matrix <double> input2 (input);
	matrix <double> inverse2 (inverse);
	bool test = InvertMatrix (input2, inverse2);
	inverse = symmetric_matrix<double, lower> (inverse2);
	return test;
}

bool NU::InvertMatrix(matrix<double>& input, symmetric_matrix<double, lower>& inverse)
{
	matrix <double> inverse2 (inverse);
	bool test = InvertMatrix (input, inverse2);
	inverse = symmetric_matrix<double, lower> (inverse2);
	return test;
}

void NU::InvertMatrixChol (const matrix<double>& input, symmetric_matrix<double, lower>& inverse)
{
	matrix <double> inverse2 (inverse);
	InvertMatrixChol (input, inverse2);
	inverse = symmetric_matrix<double, lower> (inverse2);
}

void NU::InvertMatrixChol (const matrix<double>& input, matrix<double>& inverse) throw ()
{
	matrix <double> tmp (input.size1(), input.size2());
	cholesky_decompose(input, tmp);
	inverse.resize (input.size1(), input.size2(), false);
	for (unsigned int i = 0; i < input.size1(); i++)
		for (unsigned int j = 0; j < input.size2(); j++)
			inverse (i,j) = 0.0;
	for (unsigned int i = 0; i < input.size1(); i++)
		inverse (i, i) = 1.0;
	cholesky_solve(tmp, inverse, ublas::lower());
}

/** \brief make a immutable symmetric adaptor from a matrix
  *
  * \usage: 
	<code>
	 A = symmetric< lower >(B);
	 A = symmetric(B, lower());
	</code>
  */
template < class TYPE, class MATRIX > ublas::symmetric_adaptor<const MATRIX, TYPE> symmetric(const MATRIX & A, const TYPE& uplo = TYPE())
{
	  return ublas::symmetric_adaptor<const MATRIX, TYPE>(A);
}

void NU::RInvertWishart (symmetric_matrix<double> &wishart, symmetric_matrix <double> &wishart_inv, const matrix<double> & lambda, const int & m) throw ()
{
	matrix <double> tmp (lambda);
	matrix <double> tmp_inv (wishart_inv);
	//InvertMatrix ( lambda, tmp_inv);
	InvertMatrixChol ( lambda, tmp_inv);
	Wishart (tmp, tmp_inv, m);
	wishart_inv = symmetric <ublas::lower> (tmp_inv);
	//InvertMatrix (tmp, wishart); // should be here, but removed for optimisation (we need the inverse of the result)
	InvertMatrixChol (tmp, wishart); // should be here, but removed for optimisation (we need the inverse of the result)
}

void NU::Wishart (matrix<double> &res, const matrix<double> & lambda, const int & m) throw ()
{
	matrix <double> X (m, lambda.size2());
	for (unsigned int i = 0; i < X.size1(); i++)
		for (unsigned int j = 0; j < X.size2(); j++)
			X(i,j) = RNorm();

	matrix <double> L (lambda.size1(), lambda.size1());                                        
	cholesky_decompose(lambda, L);
	for (unsigned int i = 0; i < L.size1(); i++)
		for (unsigned int j = i+1; j < L.size2(); j++)
			L(i,j) = 0;
	//L = trans(L);
	//cout << L << endl;
	res = prod (trans(X), X);
	res = prod (L, res);
	res = prod (res, trans(L));
}

const symmetric_matrix<double, ublas::lower> NU::kron (const symmetric_matrix <double, ublas::lower> & lhs, const symmetric_matrix<double, ublas::lower> & rhs) throw ()
{
    typedef typename symmetric_matrix<double, ublas::lower>::value_type scalar_t;
	typedef symmetric_matrix<double, ublas::lower>::const_iterator1 t_it1;
	typedef symmetric_matrix<double, ublas::lower>::const_iterator2 t_it2;

    const unsigned int lhs_nrows = lhs.size1();
    const unsigned int lhs_ncols = lhs.size2();
    const unsigned int rhs_nrows = rhs.size1();
    const unsigned int rhs_ncols = rhs.size2();

    symmetric_matrix <double, ublas::lower> result (lhs_nrows*rhs_nrows, lhs_ncols*rhs_ncols);

    for (t_it1 i = lhs.begin1(); i != lhs.end1(); ++i)
		for (t_it2 j = i.begin(); j != i.end(); ++j)
		{
			const unsigned int rj1 = j.index1()*rhs_nrows;
			const unsigned int cj2 = j.index2()*rhs_ncols;
			const scalar_t lhs_ij = (*j);
			
			for (t_it1 k = rhs.begin1(); k != rhs.end1(); ++k)
				for (t_it2 l = k.begin(); l != k.end(); ++l)
					result(rj1 + l.index1(), cj2 + l.index2()) = lhs_ij * (*l);
		}
    return result;
}

matrix <double> NU::InvARMatrix (const unsigned int & size1, const unsigned int & size2, const double & phi) throw ()
{
    matrix <double> result (size1, size2);
	if (size1 == 1 && size2 == 1)
	{
		result (0,0) = 1;
		return result;
	}

    for (unsigned int i = 0; i < result.size1(); i++)
        for (unsigned int j = 0; j < result.size2(); j++)
            result (i,j) = 0.0;

    double col0 = 1.0 / (1.0 - (phi * phi));
    double col1 = (phi * -1) / (1.0 - (phi * phi));
    double col2 = (1.0 + (phi * phi)) / (1.0 - (phi * phi));
    for (unsigned int i = 1; i < size1-1; i++)
    {   
        result (i, i-1) = col1;
        result (i, i) = col2;
        result (i, i+1) = col1;
    }   
    result (0,0) = col0;
    result (0,1) = col1;

    result (size1-1, size2-1) = col0;
    result (size1-1, size2-2) = col1;
    return result;
}

double NU::RNorm () throw ()
{
	return generator_norm();
}

double NU::RUnif (const double & min, const double & max) throw ()
{
	boost::random::uniform_real_distribution <> d8 (min,max);
	boost::variate_generator<boost::mt19937&,boost::random::uniform_real_distribution<> > var_runif( gen, d8 );
	return var_runif();
}

double NU::RGamma (const double & shape, const double & scale) throw ()
{
	boost::gamma_distribution<> gd( shape );
	boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma( gen, gd );

	return 1.0 / (scale*var_gamma());
}

void NU::AddMatrix (matrix <double> &a, matrix <double> b) throw ()
{
	for (unsigned int i = 0; i < a.size1(); i++)
		for (unsigned int j = 0; j < a.size2(); j++)
			a (i, j) += b (i,j);
}

#undef BNU
#undef NU
