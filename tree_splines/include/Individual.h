/**
 * @file Individual.h
 * @synopsis Individual description for handling dataset
 **/

#ifndef __INDIVIDUAL_H__
#define __INDIVIDUAL_H__

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "Params.h"

namespace nsTreeSplines
{
	class Individual
	{
		private:
			nsTreeSplines::Params * m_Params; // access to the program parameters (contains the size needed for the array)
			double * m_Array;  // table of observations
			int		 m_Id;	   // id of the observation in the database

		public:
			Individual (nsTreeSplines::Params * params) throw ();
			Individual (nsTreeSplines::Params * params, const int & id) throw ();
			Individual (nsTreeSplines::Individual * individual) throw ();
			virtual ~Individual ()			throw ();

			int    GetId      ()		  throw ();
			double * GetValues ()		  throw ();
			double GetValueAt (const int & i) throw ();
			void   SetId	  (const int & id) throw ();
			nsTreeSplines::Params * GetParams () throw ();

			double& operator() (unsigned dose, unsigned time, unsigned repeat);
			double  operator() (unsigned dose, unsigned time, unsigned repeat) const;
	};
	std::ostream& operator<< (std::ostream & os, nsTreeSplines::Individual * individual) throw ();
}

#include "Individual.hxx"

#endif /* __INDIVIDUAL_H__ */
