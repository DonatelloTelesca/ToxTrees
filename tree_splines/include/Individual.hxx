/**
 * @File     : Individual.hxx
 * @Authors  : C. Low-Kam, L. Di-Jorio
 * @Date     : 29/10/2011
 * @Version  : V1.0
 * @Synopsis : 
**/

#if !defined __INDIVIDUAL_HXX__
#define      __INDIVIDUAL_HXX__

#include "Individual.h"

#define IND nsTreeSplines::Individual

inline int IND::GetId () throw ()
{
	return m_Id;
}

inline double * IND::GetValues () throw ()
{
	return m_Array;
}

inline void IND::SetId (const int & id) throw ()
{
	m_Id = id;
}

inline nsTreeSplines::Params * IND::GetParams (void) throw ()
{
	return m_Params;
}

inline double IND::GetValueAt (const int & i) throw ()
{
	return m_Array[i];
}

inline double& IND::operator() (unsigned dose, unsigned time, unsigned repeat)
{
	return m_Array [ ((dose*m_Params->nb_time)+time) * m_Params->nb_repeat + repeat];
}

inline double IND::operator() (unsigned dose, unsigned time, unsigned repeat) const
{
	return m_Array [ ((dose*m_Params->nb_time)+time) * m_Params->nb_repeat + repeat];
}

#undef IND
#endif
