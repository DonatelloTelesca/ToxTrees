/*
 * node.cpp
 *
 *  Created on: 13 july 2011
 *      Author
 */

#include "YNode.h"
#include <iostream>

using namespace std;
using namespace nsLBart;

#define YNODE nsLBart::YNode

YNODE::YNode() throw ()
	: Node (),
	m_Object (-1)
{ }

YNODE::YNode(const int & obj) throw ()
	: Node (),
	m_Object (obj)
{ }

void YNODE::PrintTree() throw ()
{
	// Print node and all descendants
	cout << "--------" << endl;
	cout << "depth: " << m_Depth << " index " << m_RuleIndex << " rule " << m_RuleSplit << " obj " << m_Object << endl;
	this->Print();
	if (m_Lchild != NULL)
		m_Lchild->PrintTree();
	if (m_Rchild != NULL)
		m_Rchild->PrintTree();
}

double YNODE::GetSum () throw ()
{
	double sum = 0.0;
	for (unsigned int i = 0; i < m_AllObservations.size(); i++)
		sum += m_AllObservations[i]->GetValueAt(2);
	return sum;
}
