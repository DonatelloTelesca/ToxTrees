/*
 * Node.h
 *
 *  Created on: 13 juil. 2011
 *      Author
 */

#ifndef __NODE_H__
#define __NODE_H__

#include <string>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <boost/numeric/ublas/matrix.hpp>

#include "Individual.h"
#include "Results.h"

#define BNU boost::numeric::ublas

namespace nsTreeSplines
{
	class Node 
	{
		protected:

			Node * m_Father;      // Father of Node
			Node * m_Lchild;      // Left child of Node
			Node * m_Rchild;      // Right child of Node
			int    m_RuleIndex;   // Index of splitting rule
			int    m_Depth;       // Depth of Node, root has 0
			double m_RuleSplit;   // Value of splitting rule
			bool   m_IsLeftChild; // True if is a left child
			
			std::vector<Individual*> m_AllIndividuals; // All Individuals contained in Node, size changes
			BNU::matrix <double> m_SigmaInvert;
			BNU::vector <double> m_Mu;
			BNU::vector <double> m_Beta;

		public:
			Node () throw ();
			Node (Node * n, const bool & fill_obs) throw ();
			Node (Node * father, const int & depth, const bool & islc) throw ();
			virtual ~Node() throw ();

			// -------------------------------------- //
			//                GETTERS                 //
			// -------------------------------------- //

			int    GetRule () throw ();
			int    GetNinds () throw ();
			int    GetDepth() throw ();
			double GetValue() throw ();
			bool   GetIsLeftChild() throw ();
			Node * GetFather() throw ();
			Node * GetLChild() throw ();
			Node * GetRChild() throw ();
			Individual* GetIndividualAt (int i)   throw ();
			std::vector<int> GetAvailableRules () throw ();
			int GetIndividualSize () throw ();
			BNU::vector <double> & GetMu () throw ();
			BNU::vector <double> & GetBeta () throw ();
			BNU::matrix <double> & GetSigmaInvert () throw ();

			// -------------------------------------- //
			//                SETTERS                 //
			// -------------------------------------- //

			void AddIndividual (Individual*) throw ();
			void SetRule   (int i)    throw ();
			void SetNobs   (int i)    throw ();
			void SetDepth  (int d)    throw ();
			void SetValue  (double j) throw ();
			void SetFather (Node * n) throw ();
			void SetLChild (Node * n) throw ();
			void SetRChild (Node * n) throw ();
			void SetIsLeftChild ()    throw ();

			// -------------------------------------- //
			//                METHODS                 //
			// -------------------------------------- //

			void Print() throw ();
			virtual void PrintTree() throw ();
			void GetInternalNodes (std::vector <Node *> & all_nodes) throw ();
			void GetLeaves (std::vector <Node *> & leaves) throw ();
			std::vector<int> AvailablePredictor() throw ();
			bool Grow() throw ();
			void Grow (const int & index, const double & value) throw ();
			void Prune(const bool & keep) throw ();
			Node * CopyBranch(Node * np, Node * nc, const bool & fill_obs) throw ();
			void ClearIndividuals () throw ();
			void DeleteIndividuals () throw ();
			void Empty() throw ();
			void ReassembleIndividuals (const bool & keep, const bool & original) throw ();
			void DisassembleIndividuals () throw ();
			Node * Change() throw ();
			bool CanChange() throw ();
			bool CanGrow() throw ();
			Node * Swap(const bool & swap_with_left) throw ();
			double GetLogILik(BNU::symmetric_matrix <double, BNU::lower> & sigma_hat_j_inv, const double & tau) throw ();
			double LogPriT   (const double & alpha, const double & beta) throw ();
			void ComputeSigmaMatrix (BNU::symmetric_matrix <double, BNU::lower> & sigma_hat_j_inv, const double & tau) throw ();
			void ComputeMuMatrix (BNU::symmetric_matrix <double, BNU::lower> & sigma_hat_j_inv) throw ();
			void ComputeResiduals (std::vector < BNU::matrix <double> > & s_res, const int & index, const int & repeat) throw ();
			void ComputeNu (std::vector<Individual *> & new_database, const double & current_sigma) throw ();
			void ComputeR (std::vector< Individual *> & new_database, std::vector< std::vector < Individual *> > & database, const int & tree) throw ();
			void InitResults (nsTreeSplines::Results * results, const int & current) throw ();
			void TestBetaMatrix (nsTreeSplines::Results * results, double ** test_predictors, const int & current) throw ();
	};
}

#include "Node.hxx"

#undef BNU
#endif /* __NODE_H__ */
