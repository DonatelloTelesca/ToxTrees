/*
 * observation.h
 *
 *  Created on: 13 juil. 2011
 *      Author
 */

#ifndef __OBSERVATION_H__
#define __OBSERVATION_H__

#include <stdio.h>
#include <stdlib.h>
#include <string>

namespace nsLBart
{
	class Observation 
	{
		private:
			int      m_Length; // length of observation (number of columns)
			double * m_Array;  // table of observations (always the same length)
			int		 m_Id;	   // id of the observation in the database

		public:
			//Observation (std::string line, const bool & y) throw ();
			Observation (const int & length) throw ();
			Observation (const int & length, const int & id) throw ();
			Observation (Observation * o)    throw ();
			virtual ~Observation()			 throw ();

			void   Print      ()          throw ();
			int    GetLength  ()          throw ();
			double GetValueAt (int index) throw ();
			int    GetId      ()		  throw ();
			double * GetValues ()		  throw ();
			void   SetLength  (const int & l) throw ();
			void   SetId	  (const int & id) throw ();
			void   SetValueAt (const int & i, const double & val) throw ();
			//void   Substract  (const int & i, const double & val) throw ();
	};
}

#include "Observation.hxx"

#endif /* __OBSERVATION_H__ */
