/*
 * node.cpp
 *
 *  Created on: 13 july 2011
 *      Author
 */

#include "Node.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <map>

#include "Cholesky.hpp"
#include "Utils.h"

using namespace std;
using namespace nsTreeSplines;
using namespace nsUtils;

#define NODE nsTreeSplines::Node
#define BNU boost::numeric::ublas

// --------------------------------- //
// ---------   METHODS    ---------- //
// --------------------------------- //

double unifRand()
{
    return rand() / double(RAND_MAX);
}

boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator(boost::mt19937(time(0)),boost::normal_distribution<>());
// --------------------------------- //

NODE::Node() throw ()
	:m_Father(NULL),
	m_Lchild (NULL),
	m_Rchild (NULL),
	m_RuleIndex (-1),
	m_Depth (0),
	m_RuleSplit (100000000.0),
	m_IsLeftChild (false)
{ }

NODE::Node (NODE * n, const bool & fill_obs) throw ()
{
	m_Father = n->GetFather ();
	m_Lchild = n->GetLChild ();
	m_Rchild = n->GetRChild ();
	m_RuleIndex = n->GetRule ();
	m_Depth = n->GetDepth ();
	m_RuleSplit = n->GetValue ();
	m_IsLeftChild = n->GetIsLeftChild ();
	m_SigmaInvert = n->GetSigmaInvert();
	m_Mu = n->GetMu();
	m_Beta = n->GetBeta();

	if (fill_obs)
		for (int i = 0; i < n->GetNinds(); i++)
			m_AllIndividuals.push_back(n->GetIndividualAt(i));//new Individual(n->GetIndividualAt(i)));
}

NODE::Node(NODE * father, const int & depth, const bool & islc) throw ()
	:m_Father(father),
	m_Lchild (NULL),
	m_Rchild (NULL),
	m_RuleIndex (-1),
	m_Depth (depth),
	m_RuleSplit (100000000.0),
	m_IsLeftChild (islc)
{}

NODE::~Node() throw ()
{
	m_AllIndividuals.clear ();

	if (m_Lchild != NULL)
		delete m_Lchild;
	if (m_Rchild != NULL)
		delete m_Rchild;
}

void NODE::Print() throw ()
{
	for (unsigned int i=0; i< m_AllIndividuals.size(); i++)
		cout << m_AllIndividuals[i] << endl;
}

void NODE::PrintTree() throw ()
{
	// Print node and all descendants
	cout << "--------" << endl;
	cout << "depth: " << m_Depth << " index " << m_RuleIndex << " rule " << m_RuleSplit << endl;
	this->Print();
	if (m_Lchild != NULL)
		m_Lchild->PrintTree();
	if (m_Rchild != NULL)
		m_Rchild->PrintTree();
}

void NODE::GetInternalNodes (vector <Node *> & all_nodes) throw ()
{
	if (m_Rchild != NULL || m_Lchild != NULL)
		all_nodes.push_back(this);
	if (m_Lchild != NULL)
		m_Lchild->GetInternalNodes (all_nodes);
	if (m_Rchild != NULL)
		m_Rchild->GetInternalNodes (all_nodes);
}

void NODE::GetLeaves (std::vector <Node *> & leaves) throw ()
{
	if (m_Lchild == NULL && m_Rchild == NULL)
		leaves.push_back (this);
	if (m_Lchild != NULL)
		m_Lchild->GetLeaves (leaves);
	if (m_Rchild != NULL)
		m_Rchild->GetLeaves (leaves);
}

vector<int> NODE::AvailablePredictor() throw ()
{
	vector <int> ind;
	for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
		ind.push_back(m_AllIndividuals[i]->GetId());
	// Return vector of indexes of available variables for splits
	return m_AllIndividuals[0]->GetParams()->GetAvailablePredictors (ind);
}

bool NODE::Grow() throw ()
{
    vector <int> indiv;
	for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
		indiv.push_back(m_AllIndividuals[i]->GetId());

	vector<int> tmp = m_AllIndividuals[0]->GetParams()->GetAvailablePredictors(indiv);
	
	//cout << "number of observations: " << m_AllIndividuals.size() << endl;
	//cout << "number of available predictors: " << tmp.size() << endl;

    if (m_AllIndividuals.size () < 2 || tmp.size () < 1)
	{
		//cout << "can not grow" << endl;
		return false;
	}

	m_AllIndividuals[0]->GetParams()->ChooseSplitter (m_RuleIndex, m_RuleSplit);
	//cout << "Chosen splitters " << m_RuleIndex << " " << m_RuleSplit << endl;
	
	// Create left child and right child
	m_Lchild = new Node (this, m_Depth+1, true);
	m_Rchild = new Node (this, m_Depth+1, false);
	for (unsigned int i=0; i < m_AllIndividuals.size(); i++)
	{
		if (m_AllIndividuals[0]->GetParams()->predictors[m_AllIndividuals[i]->GetId()][m_RuleIndex] < m_RuleSplit)
			m_Lchild->AddIndividual (m_AllIndividuals[i]);//new Individual (m_AllIndividuals[i]));
		else m_Rchild->AddIndividual(m_AllIndividuals[i]);//new Individual(m_AllIndividuals[i]));
	}

	//m_Rchild->Print();
	//cout << "---------------" << endl;
	//m_Rchild->Print();
	//cout << "---------------" << endl;
	//cout << "depth of parent: "      << m_Depth << endl;
	//cout << "depth of left child: "  << m_Lchild->GetDepth() << endl;
	//cout << "depth of right child: " << m_Rchild->GetDepth()<< endl;

	return true;
}

void NODE::Grow (const int & index, const double & value) throw ()
{
	cout << "Appel Grow" << endl;
	m_RuleIndex = index;
	m_RuleSplit = value;

	m_Lchild = new Node (this, m_Depth+1, false);
	m_Rchild = new Node (this, m_Depth+1, true);
	for (unsigned int i=0; i < m_AllIndividuals.size(); i++)
	{
		if (m_AllIndividuals[0]->GetParams()->predictors[m_AllIndividuals[i]->GetId()][m_RuleIndex] < m_RuleSplit)
			m_Lchild->AddIndividual (m_AllIndividuals[i]);//new Individual (m_AllIndividuals[i]));
		else m_Rchild->AddIndividual(m_AllIndividuals[i]);//new Individual(m_AllIndividuals[i]));
	}
}

void NODE::Prune(const bool & keep) throw ()
{
	if (m_Lchild->GetRule() == -1 && m_Rchild->GetRule()== -1) // Prune a pair of terminal nodes only
	{
		if (keep)
		{
			for (int i = 0; i < m_Lchild->GetNinds (); i++)
				m_AllIndividuals.push_back(m_Lchild->GetIndividualAt(i));//new Individual(m_Lchild->GetIndividualAt(i)));
			for (int i = 0; i < m_Rchild->GetNinds (); i++)
				m_AllIndividuals.push_back(m_Rchild->GetIndividualAt(i));//new Individual(m_Rchild->GetIndividualAt(i)));
		}

		delete m_Lchild;
		m_Lchild = NULL;

		delete m_Rchild;
		m_Rchild = NULL;
		
		m_RuleIndex = -1;
	    m_RuleSplit = 100000000.0;
	}
}

void NODE::ClearIndividuals () throw ()
{
	/*for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
		delete m_AllIndividuals[i];*/
	m_AllIndividuals.clear ();
}

void NODE::DeleteIndividuals () throw ()
{
	for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
		delete m_AllIndividuals[i];
}

NODE * NODE::CopyBranch (NODE * np, NODE * nc, const bool & fill_obs) throw ()
{
	Node * branch = new Node (nc, fill_obs);
	branch->SetFather (np);
	if (branch->GetLChild() != NULL)
		branch->SetLChild(CopyBranch (branch, nc->GetLChild(), fill_obs));
	if (branch->GetRChild() != NULL)
		branch->SetRChild(CopyBranch (branch, nc->GetRChild(), fill_obs));

	return branch;
}

void NODE::Empty() throw ()
{
	//for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
	//	delete m_AllIndividuals[i];
	m_AllIndividuals.clear();
	if (m_Lchild != NULL)
		m_Lchild->Empty();
	if (m_Rchild != NULL)
		m_Rchild->Empty();
}

void NODE::ReassembleIndividuals (const bool & keep = false, const bool & original = true) throw ()
{
	if (m_Lchild == NULL && m_Rchild == NULL) // terminal node
	{
		for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
			m_Father->AddIndividual(m_AllIndividuals[i]);//new Individual (m_AllIndividuals[i]));
		if (!keep)
			m_AllIndividuals.clear ();
	}
	else
	{
		m_Lchild->ReassembleIndividuals (keep, false);
		m_Rchild->ReassembleIndividuals (keep, false);
		if (m_Father != NULL && !original)
		{
			for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
				m_Father->AddIndividual(m_AllIndividuals[i]);//new Individual (m_AllIndividuals[i]));
			if (!keep)
				m_AllIndividuals.clear();
		}
	}
}

void NODE::DisassembleIndividuals () throw ()
{
	if (m_Lchild == NULL && m_Rchild == NULL)
		return;

	ClearIndividuals ();
	
	if (m_Lchild != NULL)
		m_Lchild->DisassembleIndividuals ();
	if (m_Rchild != NULL)
		m_Rchild->DisassembleIndividuals ();
}

NODE * NODE::Change () throw ()
{
	Node * branch = CopyBranch (GetFather(), this, true);
	branch->ReassembleIndividuals(false, true);

	//cout << "  ---------- BRANCH ------------  " << endl;
	//branch->PrintTree ();
	//cout << "  ---------- END BRANCH ------------  " << endl;

	// Draw a new rule and splitting value
	m_AllIndividuals[0]->GetParams()->ChooseSplitter (branch->m_RuleIndex, branch->m_RuleSplit);

	//cout << "new rule index: " << branch->GetRule() << endl;

	vector <double> values = m_AllIndividuals[0]->GetParams()->GetAvailableValues (branch->GetRule());
	//for (unsigned int i = 0; i < values.size(); i++)
	//	cout << " " << values[i];
	//cout << endl;

	while (true)
	{
		//cout << "VALUES {";
		//for (unsigned int i = 0; i < values.size(); i++)
		//	cout << (i == 0 ? "" : ",") << values[i];
		//cout << "}" << endl;
		// Randomly choose a splitting value
		double u = unifRand();
		int index2 = (int)floor(u*values.size());
		branch->SetValue(values[index2]);

		//cout << "--- choix valeur " << index2 << " = " << values[index2] << endl;
		if (branch->CanChange ())
		{
			//cout << "changement ok" << endl;
			return branch;
		}
		else
		{
			//cout << "--------- On ne peut changer " << endl;
			if (values.size() == 1)
			{
				//cout << "#################### DELETE BRANCH " << endl;
				delete branch;
				return NULL;
			}
			values.erase(values.begin()+index2);
			//cout << "  ---------- BEFORE ------------  " << endl;
			//branch->PrintTree ();
			//cout << "  ---------- END BEFORE ------------  " << endl;
			branch->GetLChild()->Empty();
			branch->GetRChild()->Empty();
			//branch->PrintTree();
		}
	}
	delete branch;
	return NULL;
}

bool NODE::CanChange () throw ()
{
	//cout << " ---- NODE ---- " << endl;
	//Print();
	//cout << " ---- END NODE ---- " << endl;
	if (m_Lchild == NULL && m_Rchild == NULL)
		return true;
	if (m_AllIndividuals.size() == 1)
		return false;

	//cout << "ICI" << endl;
	
	/*vector <double> values = m_AllIndividuals[0]->GetParams()->GetAvailableValues (m_RuleIndex);
	if (!((findRule(AvailablePredictor(), m_RuleIndex) != -1) && (values[0] < m_RuleSplit) && (values[values.size()-1] >= m_RuleSplit)))
		return false;*/

	//m_Lchild->Empty();
	//m_Rchild->Empty();
	// Fill children with values
	for (unsigned int i=0; i < m_AllIndividuals.size(); i++)
	{
		//if (m_AllIndividuals[0]->GetParams()->predictors[i][m_RuleIndex] < m_RuleSplit)
		if (m_AllIndividuals[0]->GetParams()->predictors[m_AllIndividuals[i]->GetId()][m_RuleIndex] < m_RuleSplit)
			m_Lchild->AddIndividual (m_AllIndividuals[i]);//new Individual (m_AllIndividuals[i]));
		else m_Rchild->AddIndividual(m_AllIndividuals[i]);//new Individual(m_AllIndividuals[i]));
	}

	if (m_Lchild->GetNinds() == 0 || m_Rchild->GetNinds() == 0)
		return false;
	return (m_Lchild->CanChange () && m_Rchild->CanChange ());
}

bool NODE::CanGrow () throw ()
{
	return (m_AllIndividuals[0] != NULL && m_AllIndividuals.size() > 4 && (AvailablePredictor()).size() > 1);
}

// Change rule with father
// If before swap the brother has the same rule, swap with the two brothers
// Then check if swap is possible
NODE * NODE::Swap (const bool & swap_with_left) throw ()
{
	Node * branch = CopyBranch (GetFather(), this, true);
	branch->ReassembleIndividuals(false, true);
	for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
		branch->AddIndividual(new Individual(m_AllIndividuals[i]));

	if (m_Rchild->GetRule() == m_Lchild->GetRule() && m_Rchild->GetValue() == m_Lchild->GetValue()) // two sons has same rules
	{
		branch->GetRChild()->SetRule(m_RuleIndex); 
		branch->GetRChild()->SetValue(m_RuleSplit); 
		branch->GetLChild()->SetRule(m_RuleIndex); 
		branch->GetLChild()->SetValue(m_RuleSplit); 
	}
	else if (swap_with_left)
	{
		branch->GetLChild()->SetRule(m_RuleIndex); 
		branch->GetLChild()->SetValue(m_RuleSplit); 
	}
	else
	{
		branch->GetRChild()->SetRule(m_RuleIndex); 
		branch->GetRChild()->SetValue(m_RuleSplit); 
	}
	
	branch->SetRule (swap_with_left ? m_Lchild->GetRule() : m_Rchild->GetRule());
	branch->SetValue (swap_with_left ? m_Lchild->GetValue() : m_Rchild->GetValue());
	
	if (branch->CanChange ())
		return branch;

	delete branch;
	return NULL;
}

double NODE::GetLogILik(BNU::symmetric_matrix <double, BNU::lower> & sigma_hat_j_inv, const double & tau) throw ()
{
	double logilik = 0.0;
	ComputeSigmaMatrix (sigma_hat_j_inv, tau);
	ComputeMuMatrix (sigma_hat_j_inv);
// 	cout << "------------------ NOEUD INDIV [ ";
	for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
// 		cout << (i == 0 ? "" : ",") << m_AllIndividuals[i]->GetId();
// 	cout << "] m_beta " << m_Beta << endl;
	logilik += (0.5 * log (m_AllIndividuals[0]->GetParams()->penalty_matrix_det)) - (0.5 * log (Determinant (m_SigmaInvert)))
		- ((double)((m_AllIndividuals[0]->GetParams()->nb_dose + m_AllIndividuals[0]->GetParams()->delta) * 
					(m_AllIndividuals[0]->GetParams()->nb_time > 1 ? (m_AllIndividuals[0]->GetParams()->nb_time + m_AllIndividuals[0]->GetParams()->delta) : 1)/ 2.0)) * log (tau);
	ublas::vector<double> b_copy (m_Mu);
	b_copy = prod (trans(m_Mu), m_SigmaInvert);
	
	//b_copy = prod (b_copy, b); TODO does not work, making it by hand
	double x = 0.0;
	for (unsigned int i = 0; i < m_Mu.size(); i++)
		x += m_Mu(i) * b_copy(i);
	
	return (logilik + (0.5 * x));
}

double NODE::LogPriT(const double & alpha, const double & beta) throw ()
{
	//cout << " ------------------- TREATING ------------------" << endl;
	//Print ();
	//cout << " -----------------------------------------------" << endl;
	double pgrow = (AvailablePredictor().size () <= 1 ? 0 : ((m_AllIndividuals.size () < 2 ? .001 : 1) * alpha/ (pow(1.0 + m_Depth, beta))));
	//cout << "pgrow " << pgrow << " alpha " << alpha << " beta " << beta << " " << AvailablePredictor().size () << endl;
	// Prior probability of branch
	if(m_Lchild == NULL)
		return log (1.0 - pgrow);
	
	double retval = log( pgrow );
	//cout << "nbr predic " << AvailablePredictor().size() << endl
	retval -= log(AvailablePredictor().size());
	retval -= log((double)(m_AllIndividuals[0]->GetParams()->GetAvailableValues(m_RuleIndex)).size());
	retval += m_Lchild->LogPriT(alpha,beta) + m_Rchild->LogPriT(alpha,beta);
	//cout << "logPriT: " << retval << endl;

	return retval;
}

void NODE::ComputeSigmaMatrix (BNU::symmetric_matrix <double, BNU::lower> & sigma_hat_j_inv, const double & tau) throw ()
{
	//cout << "$$$$$$$$$$$$$$$$$$$$$ COMP " << endl;
	//Print ();
	int size = (m_AllIndividuals[0]->GetParams()->nb_dose + m_AllIndividuals[0]->GetParams()->delta) * (m_AllIndividuals[0]->GetParams()->nb_time + m_AllIndividuals[0]->GetParams()->delta);
	m_SigmaInvert.resize (size, size, false);
	//cout << "test 2" << endl;
	//matrix <double> trans_spline_base_matrix = trans (params->spline_base_matrix);
	//cout << m_AllIndividuals[0]->GetParams()->nb_repeat << " * " << m_AllIndividuals.size() << endl << " --------TTT ------ ";
	//Print ();
	m_SigmaInvert = m_AllIndividuals[0]->GetParams()->nb_repeat * m_AllIndividuals.size() * 
		prod (sigma_hat_j_inv, m_AllIndividuals[0]->GetParams()->spline_base_matrix);
	//cout << "test 3" << endl;
	m_SigmaInvert = prod (trans (m_AllIndividuals[0]->GetParams()->spline_base_matrix), m_SigmaInvert);
	//cout << "test 1" << endl;
	//cout << "TAU4" <<  m_AllIndividuals[0]->GetParams()->tau[j] << endl;
	m_SigmaInvert += (1 / tau) * m_AllIndividuals[0]->GetParams()->penalty_matrix;
	//cout << "test 5" << endl;
}

void NODE::ComputeMuMatrix (BNU::symmetric_matrix <double, BNU::lower> & sigma_hat_j_inv) throw ()
{
	// Cholesky decompose    
	BNU::matrix <double> L (m_SigmaInvert.size1(), m_SigmaInvert.size2());
	size_t res = cholesky_decompose(m_SigmaInvert, L);

	if (res != 0) // TODO throw an exeption here
	{
		cerr << "Cholesky decomposition failure" << endl;
		return;
	}


	// Compute beta here (optimisation)
	int size = 0;
	if (m_AllIndividuals[0]->GetParams()->nb_time == 1)
		size = m_AllIndividuals[0]->GetParams()->nb_dose + m_AllIndividuals[0]->GetParams()->delta;
	else
		size = (m_AllIndividuals[0]->GetParams()->nb_dose + m_AllIndividuals[0]->GetParams()->delta) * (m_AllIndividuals[0]->GetParams()->nb_time + m_AllIndividuals[0]->GetParams()->delta);
	m_Beta.resize(size, false);
	for (unsigned int i = 0; i < m_Beta.size(); i++)
		m_Beta(i) = RNorm();

	solve(trans(L), m_Beta, ublas::lower() );
	
	m_Mu.resize (m_AllIndividuals[0]->GetParams()->nb_dose * m_AllIndividuals[0]->GetParams()->nb_time, false);
	for (int d = 0; d < m_AllIndividuals[0]->GetParams()->nb_dose; d++)
	{                        
		for (int t = 0; t < m_AllIndividuals[0]->GetParams()->nb_time; t++)
		{
			double sum = 0;
			for (unsigned int i = 0; i < m_AllIndividuals.size (); i++)
				for (int k = 0; k < m_AllIndividuals[0]->GetParams()->nb_repeat; k++)
					sum += (*m_AllIndividuals[i])(d, t, k);
			//cout << "size " << m_Mu.size() << " nb dose " << m_AllIndividuals[0]->GetParams()->nb_dose << " " <<
			//	(d * (m_AllIndividuals[0]->GetParams()->nb_dose -1) + t) << " d " << d << " t " << t << endl;
			m_Mu((d * m_AllIndividuals[0]->GetParams()->nb_time) + t) = sum; 
		}
	}

	BNU::matrix <double> tmp = prod(trans (m_AllIndividuals[0]->GetParams()->spline_base_matrix), sigma_hat_j_inv); 
	m_Mu = prod (tmp, m_Mu);
	 
	cholesky_solve(L, m_Mu, ublas::lower());
	//cout << "preb " << m_Beta;
	m_Beta += m_Mu;
	//cout << " - m_Beta --" << m_Beta << " m_Mu " << m_Mu << endl;
}

void NODE::ComputeResiduals (std::vector < BNU::matrix <double> > & s_res_vect, const int & index, const int & repeat) throw ()
{
	BNU::vector <double> tmp = prod (m_AllIndividuals[0]->GetParams()->spline_base_matrix, m_Beta);
	for (int d = 0; d < m_AllIndividuals[0]->GetParams()->nb_dose; d++)
	{
		for (int t = 0; t < m_AllIndividuals[0]->GetParams()->nb_time; t++)
		{
			//double sum = 0.0;
			for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
			{
				//for (int k = 0; k < m_AllIndividuals[0]->GetParams()->nb_repeat; k++)
				//{
// 				cout << "BEFORE "  << "index " << index << " ind " << i  << " repet " << repeat << " dose " << d  << " " << (*m_AllIndividuals[i])(d, t, repeat) << " - " << tmp((d * m_AllIndividuals[0]->GetParams()->nb_time) + t);
				(*m_AllIndividuals[i])(d, t, repeat) = (*m_AllIndividuals[i])(d, t, repeat) - tmp((d * m_AllIndividuals[0]->GetParams()->nb_time) + t); 
// 				cout << " = " << (*m_AllIndividuals[i])(d, t, repeat) << endl;
				s_res_vect [(m_AllIndividuals[i]->GetId() * m_AllIndividuals[0]->GetParams()->nb_repeat) + repeat] ((d * m_AllIndividuals[0]->GetParams()->nb_time) + t, index) = (*m_AllIndividuals[i])(d, t, repeat);
					//sum += (*m_AllIndividuals[i])(d, t, k);
				//}
			}
			//s_res ((d * m_AllIndividuals[0]->GetParams()->nb_time) + t, index) += sum;
		}
	}
}

void NODE::ComputeNu (std::vector<Individual *> & new_database, const double & current_sigma) throw ()
{
	for (unsigned int i = 0; i < m_AllIndividuals.size(); i++)
	{
		for (int d = 0; d < m_AllIndividuals[0]->GetParams()->nb_dose; d++)
		{
			for (int t = 0; t < m_AllIndividuals[0]->GetParams()->nb_time; t++)
			{
				for (int k = 0; k < m_AllIndividuals[0]->GetParams()->nb_repeat; k++)
				{
					//cout << "[" << m_AllIndividuals[i]->GetId() << "] (" << d << "," << t << "," << k << ") * " << current_sigma << endl;
					(*new_database[ m_AllIndividuals[i]->GetId() ])(d, t, k) += (current_sigma * (*m_AllIndividuals[i])(d, t, k)); // TODO incorporate time
				}
			}
		}
	}
	//cout << " ---------------------------- " << endl;
}

void NODE::ComputeR (std::vector< Individual *> & new_database, std::vector< std::vector < Individual *> > & database, const int & tree) throw ()
{
	for (unsigned int i = 0; i < m_AllIndividuals.size (); i++)
	{
		for (int d = 0; d < m_AllIndividuals[0]->GetParams()->nb_dose; d++)
		{
			for (int t = 0; t < m_AllIndividuals[0]->GetParams()->nb_time; t++)
			{
				for (int k = 0; k < m_AllIndividuals[0]->GetParams()->nb_repeat; k++)
				{
					(*m_AllIndividuals[i])(d, t, k) = (*database[ m_AllIndividuals[i]->GetId() ][tree])(d, t, k)  //TODO Time to change here
						- (*new_database[ m_AllIndividuals[i]->GetId() ])(d, t, k);
				}
			}
		}
	}
}

void NODE::InitResults (nsTreeSplines::Results * results, const int & current) throw ()
{
	if (m_Father == NULL && m_Lchild == NULL && m_Rchild == NULL)
		return;
	if (m_Lchild == NULL && m_Rchild == NULL)
		return;

	results->current_predictors [current][m_RuleIndex]++;
	map <double, int>::iterator it = results->current_predictor_values [current][m_RuleIndex].find (m_RuleSplit);
	if (it == results->current_predictor_values [current][m_RuleIndex].end())
		results->current_predictor_values [current][m_RuleIndex][m_RuleSplit] = 1;
	else
		results->current_predictor_values [current][m_RuleIndex][m_RuleSplit]++;
	if (m_Lchild != NULL)
		m_Lchild->InitResults (results, current);
	if (m_Rchild != NULL)
		m_Rchild->InitResults (results, current);
}

void NODE::TestBetaMatrix (nsTreeSplines::Results * results, double ** test_predictors, const int & current) throw ()
{
	if (m_Lchild == NULL && m_Rchild == NULL)
	{
		results->beta_matrixes.push_back (m_Beta);
		return;
	}

	if (test_predictors[current][m_RuleIndex] < m_RuleSplit)
		m_Lchild->TestBetaMatrix (results, test_predictors, current);
	else m_Rchild->TestBetaMatrix(results, test_predictors, current);
}

#undef BNU
#undef NODE
