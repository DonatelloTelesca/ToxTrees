/**
 * @File     : Observation.hxx
 * @Authors  : C. Low-Kam, L. Di-Jorio
 * @Date     : 29/10/2011
 * @Version  : V1.0
 * @Synopsis : 
**/

#if !defined __OBSERVATION_HXX__
#define      __OBSERVATION_HXX__

#include "Observation.h"

#define OBS nsLBart::Observation

inline int OBS::GetLength() throw ()
{
    return m_Length;
}
 
inline double OBS::GetValueAt (int index) throw ()
{
	return m_Array[index];
}

inline int OBS::GetId () throw ()
{
	return m_Id;
}

inline double * OBS::GetValues () throw ()
{
	return m_Array;
}

inline void OBS::SetValueAt (const int & i, const double & val) throw ()
{
	m_Array[i] = val;
}

inline void OBS::SetId (const int & id) throw ()
{
	m_Id = id;
}

/*inline void OBS::Substract (const int & i, const double & val) throw ()
{
	m_Array[i] -= val;
}*/

#undef OBS
#endif
