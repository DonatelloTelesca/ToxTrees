/**
 * @File     : YNode.hxx
 * @Authors  : C. Low-Kam, L. Di-Jorio
 * @Date     : 29/10/2011
 * @Version  : V1.0
 * @Synopsis : 
**/

#if !defined __YNODE_HXX__
#define      __YNODE_HXX__

#include "YNode.h"

#define YNODE nsLBart::YNode

inline YNODE::~YNode () throw ()
{}

inline int YNODE::GetObject () throw ()
{
    return m_Object;
}

inline void YNODE::SetObject (int i) throw ()
{
    m_Object = i;
}

#undef YNODE

#endif
