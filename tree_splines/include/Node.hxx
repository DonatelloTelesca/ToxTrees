/**
 * @File     : Node.hxx
 * @Authors  : C. Low-Kam, L. Di-Jorio
 * @Date     : 29/10/2011
 * @Version  : V1.0
 * @Synopsis : 
**/

#if !defined __NODE_HXX__
#define      __NODE_HXX__

#include "Node.h"

#define NODE nsTreeSplines::Node
#define IND  nsTreeSplines::Individual
#define BNU boost::numeric::ublas

inline IND * NODE::GetIndividualAt(int i) throw ()
{
    return m_AllIndividuals[i];
}

inline void NODE::AddIndividual (Individual* o) throw ()
{
    m_AllIndividuals.push_back(o);
}

inline int NODE::GetRule() throw ()
{
    return m_RuleIndex;
}

inline void NODE::SetRule( int i) throw ()
{
    m_RuleIndex = i; 
}

inline double NODE::GetValue() throw ()
{
    return m_RuleSplit;
}

inline void NODE::SetValue(double j) throw ()
{
    m_RuleSplit=j;
}

/*inline vector<int> NODE::GetAvailableRules() throw ()
{
    return availablepredictor();
}*/

inline void NODE::SetFather(NODE * n) throw ()
{
    m_Father = n;
}

inline NODE * NODE::GetFather() throw ()
{
    return m_Father;
}

inline void NODE::SetLChild(NODE * n) throw ()
{
    m_Lchild = n;
}

inline NODE * NODE::GetLChild() throw ()
{
    return m_Lchild;
}

inline void NODE::SetRChild(NODE * n) throw ()
{
    m_Rchild = n;
}

inline NODE * NODE::GetRChild() throw ()
{
    return m_Rchild;
}

inline void NODE::SetIsLeftChild() throw ()
{
    m_IsLeftChild = true;
}

inline bool NODE::GetIsLeftChild() throw ()
{
    return m_IsLeftChild;
}

inline int NODE::GetNinds() throw ()
{
    return m_AllIndividuals.size();
}

inline void NODE::SetDepth (int d) throw ()
{
    m_Depth = d;
}

inline int NODE::GetDepth() throw ()
{
    return m_Depth;
}

/*inline int NODE::GetObject() throw ()
{
	return -1;
}*/

inline int NODE::GetIndividualSize () throw ()
{
	return m_AllIndividuals.size ();
}

inline BNU::vector <double> & NODE::GetMu () throw ()
{
	return m_Mu;
}

inline BNU::vector <double> & NODE::GetBeta () throw ()
{
	return m_Beta;
}

inline BNU::matrix <double> & NODE::GetSigmaInvert () throw ()
{
	return m_SigmaInvert;
}

#undef NODE
#undef IND
#undef BNU

#endif
