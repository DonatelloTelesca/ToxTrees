/**
 * @file tree_splines.cpp
 * @author Cecile Low-Kam, Lisa Di Jorio, based on the code of Hugh Chipman, Robert McCulloch
 * @date 24/05/2012
 * Main program, initialize param
 **/

#include "Utils.h"
#include "Params.h"
#include "Individual.h"
#include "Node.h"
#include "Cholesky.hpp"
#include "Results.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>

#include <boost/random/variate_generator.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#define MATHLIB_STANDALONE 1
#define BNU boost::numeric::ublas

using namespace std;
using namespace nsUtils;
using namespace nsTreeSplines;
using namespace BNU;

Params * params = new Params ();
Node ** currentTrees;
std::vector < std::vector <Node*> > treesLeaves;
std::vector < std::vector <Individual * > >database;

std::vector < matrix <double> > s_res_vect;
symmetric_matrix<double, BNU::lower> sigma_hat_j_inv; //(params->nb_dose);
matrix<double> sigma_tmp_1;
matrix<double> sigma_tmp_2_inv;

Results * results;
double ** test_predictors;

int current_loop;

bool ReadInitFile (string filename, string fileout) throw ()
{
	ifstream f (filename.c_str());
	string line = "";

	for (;getline (f, line);)
	{
		if (strstr (line.c_str(), "-- I J K D T P D"))
		{
			getline (f, line);
			std::vector<string> result = Split (" ", line);
			params->nb_indiv = atoi(result[0].c_str());
			params->nb_assay = atoi(result[1].c_str());
			params->nb_repeat = atoi(result[2].c_str());
			params->nb_dose = atoi(result[3].c_str());
			params->nb_time = atoi(result[4].c_str());
			params->nb_predictors = atoi(result[5].c_str());
			params->max_depth = atoi(result[6].c_str());
			if (params->nb_indiv == 0 || params->nb_dose == 0 || params->nb_repeat == 0 
					|| params->nb_assay == 0 || params->nb_predictors == 0 || params->nb_time == 0
					|| params->max_depth == 0)
			{
				cerr << "None of the -- I J K D T P D can be set to zero" << endl;
				return false;
			}

			currentTrees = new Node * [params->nb_assay];
			for (int i = 0; i < params->nb_assay; i++)
			{
				currentTrees[i] = new Node ();
				std::vector<Node *> tmp;
				treesLeaves.push_back (tmp);
				treesLeaves[treesLeaves.size()-1].push_back(currentTrees[i]);
			}

			for (unsigned int i = 0; i < (params->nb_indiv * params->nb_repeat); i++)
				s_res_vect.push_back(matrix <double> (params->nb_dose * params->nb_time, params->nb_assay));
			sigma_hat_j_inv.resize (params->nb_dose);
		}
		if (strstr (line.c_str(), "-- pgrow pprune pchange"))
		{
			getline (f, line);
			std::vector<string> result = Split (" ", line);
			params->pgrow = atof (result[0].c_str());
			params->pprune = atof (result[1].c_str());
			params->pchange = atof (result[2].c_str());
		}
		else if (strstr (line.c_str(), "-- X")) // predictor matrix
		{
			if (params->nb_indiv == 0 || params->nb_dose == 0 || params->nb_repeat == 0 
					|| params->nb_assay == 0 || params->nb_predictors == 0 || params->nb_time == 0)
			{
				cerr << "Please set the -- I J K D T P option before setting X" << endl;
				return false;
			}

			params->predictors = new double * [params->nb_indiv];
			test_predictors = new double * [params->nb_test_indiv];
			
			for (int i = 0; i < params->nb_indiv; i++)
			{
				getline (f, line);
				params->predictors[i] = new double [params->nb_predictors];
				Split (" ", line, params->predictors[i]);
			}
			for (int i = 0; i < params->nb_test_indiv; i++)
			{
				getline (f, line);
				test_predictors[i] = new double [params->nb_predictors];
				Split (" ", line, test_predictors[i]);
			}
		}
		else if (strstr (line.c_str(), "-- delta"))
		{
			getline (f, line);
			params->delta = atoi(line.c_str());
		}
		else if (strstr (line.c_str(), "-- w"))
		{
			getline (f, line);
			params->w = atof (line.c_str());
		}
		else if (strstr (line.c_str(), "-- phi_d"))
		{
			getline (f, line);
			params->phi_d = atof(line.c_str());
		}
		else if (strstr (line.c_str(), "-- phi_t"))
		{
			getline (f, line);
			params->phi_t = atof(line.c_str());
		}
		else if (strstr (line.c_str(), "-- nb_loop"))
		{
			getline (f, line);
			params->nb_loop = atoi (line.c_str());
		}
		else if (strstr (line.c_str(), "-- train_nb_loop"))
		{
			getline (f, line);
			params->nb_loop_train = atoi (line.c_str());
		}
		else if (strstr (line.c_str(), "-- alpha"))
		{
			getline (f, line);
			params->alpha = atof (line.c_str());
		}
		else if (strstr (line.c_str(), "-- beta"))
		{
			getline (f, line);
			params->beta2 = atof (line.c_str());
		}
		else if (strstr (line.c_str(), "-- nb_test_indiv"))
		{
			getline (f, line);
			params->nb_test_indiv = atoi (line.c_str());
		}
		else if (strstr (line.c_str(), "-- nb_step_test"))
		{
			getline (f, line);
			params->nb_test_step = atoi (line.c_str());
		}
		else if (strstr(line.c_str(), "-- penalty"))
		{
			int size = 0;
			if (params->nb_time == 1)
				size = (params->nb_dose + params->delta);
			else
				size = (params->nb_dose + params->delta) * (params->nb_time + params->delta);
			params->penalty_matrix.resize(size, size, false);
			for (int i = 0; i < size; i++) // TODO test delta
			{
				getline (f, line);
				Split (" ", line, params->penalty_matrix, i);
			}
			params->penalty_matrix_det = Determinant (params->penalty_matrix);
		}
		else if (strstr(line.c_str(), "-- sigma_j"))
		{
			params->sigma_j.resize (params->nb_assay);
			for (int i = 0; i < params->nb_assay; i++)
			{
				getline (f, line);
				Split (" ", line, params->sigma_j, i);
			}
		}
		else if (strstr (line.c_str(), "-- B"))
		{
			if (params->nb_indiv == 0 || params->nb_dose == 0 || params->nb_repeat == 0 
					|| params->nb_assay == 0 || params->nb_predictors == 0 || params->nb_time == 0)
			{
				cerr << "Please set the -- I J K D T P option before setting B" << endl;
				return false;
			}
			int size = 0;
			if (params->nb_time == 1)
				size = (params->nb_dose + params->delta);
			else
				size = (params->nb_dose + params->delta) * (params->nb_time + params->delta);
			params->spline_base_matrix.resize (params->nb_dose * params->nb_time, size, false);
			for (int i = 0; i < (params->nb_dose * params->nb_time); i++)
			{
				getline (f, line);
				Split (" ", line, params->spline_base_matrix, i);
			}
		}
		else if (strstr(line.c_str(), "-- a1 b1"))
		{
			if (params->nb_indiv == 0 || params->nb_dose == 0 || params->nb_repeat == 0 
					|| params->nb_assay == 0 || params->nb_predictors == 0 || params->nb_time == 0)
			{
				cerr << "Please set the -- I J K D T P option before setting B" << endl;
				return false;
			}

			getline (f, line);
			std::vector<string> result = Split (" ", line);
			params->a1 = atof(result[0].c_str());
			params->b1 = atof(result[1].c_str());
			
			params->tau = new double [params->nb_assay];
			for (int i = 0; i < params->nb_assay; i++)
			{
				//cout << "gene " << i << " " << params->a1 << " " << (1.0 / params->b1) << endl;
				params->tau[i] = RGamma (params->a1, 1.0 / params->b1);
			}
		}
		else if (strstr (line.c_str(), "-- xi"))
		{
			getline (f, line);
			params->xi = atof (line.c_str());
		}
		else if (strstr (line.c_str(), "-- lambda1 lambda2 lambda3"))
		{
			getline (f, line);
			std::vector<string> result = Split (" ", line);
			params->lambda1 = atoi(result[0].c_str());
			params->lambda2 = atoi(result[1].c_str());
			params->lambda3 = atoi(result[2].c_str());
		}
		else if (strstr(line.c_str(), "-- lambda"))
		{
			if (params->nb_indiv == 0 || params->nb_dose == 0 || params->nb_repeat == 0 
					|| params->nb_assay == 0 || params->nb_predictors == 0 || params->nb_time == 0)
			{
				cerr << "Please set the -- I J K D T P option before setting lambda" << endl;
				return false;
			}

			params->lambda.resize(params->nb_assay, params->nb_assay, false);
			for (int i = 0; i < params->nb_assay; i++)
			{
				getline (f, line);
				Split (" ", line, params->lambda, i);
			}
		}
		else if (line == ""); // Empty line : do noting
	}
	f.close ();

	results = new Results (params->nb_assay, params->nb_predictors, fileout.c_str());
	return true;
}

void ReadDataFile (string filename) throw ()
{
	ifstream f (filename.c_str());
	string line = "";

	database.resize (params->nb_indiv);

	for (int ind = 0; ind < params->nb_indiv; ind++)
	{
		database[ind].resize(params->nb_assay);
		for (int assay = 0; assay < params->nb_assay; assay++)
		{
			//currentTrees[assay]->AddIndividual (new Individual (params, ind));
			database[ind][assay] = new Individual (params, ind);
			currentTrees[assay]->AddIndividual (new Individual (params, ind));
		}

		for (int dose = 0; dose < params->nb_dose; dose++)
		{
			for (int time = 0; time < params->nb_time; time++)
			{
				getline (f, line);
				std::vector <string> result = Split (" ", line);
				for (int assay = 0; assay < params->nb_assay; assay++)
				{
					for (int repeat = 0; repeat < params->nb_repeat; repeat++)
					{
						(*currentTrees[assay]->GetIndividualAt(ind))(dose, time, repeat) = atof(result[assay*params->nb_repeat + repeat].c_str());
						(*database[ind][assay])(dose, time, repeat) = atof(result[assay*params->nb_repeat + repeat].c_str());
					}
				}
			}
		}
	}

	f.close ();
}

void InitSResVect (void) throw ()
{
	for (int a = 0; a < params->nb_assay; a++)
		for (int d = 0; d < params->nb_dose; d++)
			for (int t = 0; t < params->nb_time; t++)
				for (int i = 0; i < params->nb_indiv; i++)
					for (int r = 0; r < params->nb_repeat; r++)
						s_res_vect [(i * params->nb_repeat) + r] ((d * params->nb_time) + t, a) = (*database[i][a])(d, t, r);
}

void DisplayDatabase () throw ()
{
	for (unsigned int i = 0; i < database.size(); i++)
		for (unsigned int j = 0; j < database[i].size (); j++)
			cout << "j=" << j << " " << database[i][j];
}

bool InitParam (char ** argv) throw ()
{
	string init_file, data_file, fileout;
	for (int i = 1; i < 7; i+=2)
	{
		string tmp = "";
		if (argv[i])
			tmp = string (argv[i]);
		if (tmp == "-i") // init file following
			init_file = string (argv[i+1]);
		else if (tmp == "-f")
			data_file = string (argv[i+1]);
		else if (tmp == "-o")
			fileout = string (argv[i+1]);
		else if (tmp == "")
		{
			cerr << "Unknown option " << tmp << endl;
			return false;
		}
	}
	
	if (!ReadInitFile (init_file, fileout))
		return false;

	ReadDataFile (data_file);
	return true;
}

void InitSigmaD (void) throw ()
{
	params->sigma_d.resize (params->nb_dose, false);
	for (int i = 0; i < params->nb_dose; i++)
		for (int j = 0; j <= i; j++)
			params->sigma_d (i, j) = pow (params->phi_d, (i-j));
	
	params->sigma_t.resize (params->nb_time, false);
	for (int i = 0; i < params->nb_time; i++)
		for (int j = 0; j <= i; j++)
			params->sigma_t (i, j) = pow (params->phi_t, (i-j));
}

// routine DT
void Mh (const double & l1, const double & l2, const double & l3, double & phi, const double & special) throw ()
{
	int j = (params->nb_indiv * params->nb_assay * params->nb_repeat * special) + (params->nb_assay * special);
	double R = (params->phi_d + params->w > 1 ? 1 : params->phi_d + params->w);
	double L = (params->phi_d - params->w < -1 ? -1 : params->phi_d - params->w);

	double xStar = RUnif (L, R); 
	double R1 = (xStar + params->w > 1 ? 1 : xStar + params->w);
	double L1 = (xStar - params->w < -1 ? -1 : xStar - params->w);
	double f1 = - (j * (params->nb_dose - 1) / 2)*log(1 - params->phi_d * params->phi_d)
		+(l1 - params->phi_d * l2 + params->phi_d * params->phi_d * l3) / (2* (1 - params->phi_d * params->phi_d));
	f1 += (j * (params->nb_dose - 1)/2)*log(1 - xStar * xStar) -(l1 - xStar * l2 + xStar * xStar * l3) / (2 * ( 1 - xStar * xStar));
	f1 += -log(R1 - L1) + log(R - L);
	if (log (RUnif (0, 1)) < f1)
		phi = xStar;
}

// Returns the number of current leaves that can grow
int GetToGrow (const int & index) throw ()
{
	int size = 0;
	for (unsigned int i = 0; i < treesLeaves[index].size (); i++)
		if (treesLeaves[index][i]->GetDepth() < params->max_depth && treesLeaves[index][i]->CanGrow())
			size++;
	return size;
}

/**
 * Randomly choose a leave to grow. The leave has to be growable
 **/
int GetRandomGrow (const int & index) throw ()
{
	std::vector <int> tmp;
	for (unsigned int i = 0; i < treesLeaves[index].size (); i++)
		if (treesLeaves[index][i]->GetDepth() < params->max_depth && treesLeaves[index][i]->GetNinds() >= 2 
				&& treesLeaves[index][i]->AvailablePredictor().size() >= 1)
			tmp.push_back(i);
	if (tmp.size() == 0)
		return -1;
	return (tmp.size() == 0 ? -1 : tmp[((int)floor((rand() / double(RAND_MAX))*tmp.size()))]);
}

int GetRandomPrune (const int & index) throw ()
{
	std::vector <int> tmp;
	for (unsigned int i = 0; i < treesLeaves[index].size ();i++)
	{
		if (/*treesLeaves[index][i]->GetDepth() < params->max_depth &&*/ treesLeaves[index][i]->GetFather()->GetLChild() != NULL 
				&& treesLeaves[index][i]->GetFather()->GetRChild() != NULL // father has two sons
				&& treesLeaves[index][i]->GetFather()->GetLChild()->GetLChild() == NULL && treesLeaves[index][i]->GetFather()->GetLChild()->GetRChild() == NULL
				&& treesLeaves[index][i]->GetFather()->GetRChild()->GetLChild() == NULL && treesLeaves[index][i]->GetFather()->GetRChild()->GetRChild() == NULL)
			tmp.push_back(i);
	}
	if (tmp.size() == 0)
		return -1;
	return (tmp.size() == 0 ? -1 : tmp[((int)floor((rand() / double(RAND_MAX))*tmp.size()))]);
}

Node * GetRandomChange (const int & index) throw ()
{
	std::vector <Node *> internal;
	currentTrees[index]->GetInternalNodes (internal);
	return (internal.size() == 0 ? NULL : internal[((int)floor((rand() / double(RAND_MAX))*internal.size()))]);
}

Node * GetRandomSwap (const int & index, Node ** OldNode) throw ()
{
	std::vector <Node *> internal;
	if (currentTrees[index]->GetLChild() != NULL)
		currentTrees[index]->GetLChild()->GetInternalNodes (internal);
	if (currentTrees[index]->GetRChild() != NULL)
		currentTrees[index]->GetRChild()->GetInternalNodes (internal);

	//cout << "Nb internal " << internal.size() << endl;

	while (internal.size() > 0)
	{
		int random = (int)floor((rand() / double(RAND_MAX))*internal.size());
		Node * n = internal[random]->GetFather()->Swap (internal[random]->GetIsLeftChild());
		if (n != NULL)
		{
			(*OldNode) = internal[random]->GetFather();
			return n;
		}
		else
			internal.erase(internal.begin() + random);
	}
	return NULL;
}

/**
 * Return the number of parents of all leave (size of level n-1)
 **/
int BeforeLeaves (const int & index) throw ()
{
    int size = 0;
    if (treesLeaves[index][0]->GetFather() == NULL)
        return 0;
    for (unsigned int i = 0; i < treesLeaves[index].size ();)
    {   
        if (treesLeaves[index][i]->GetFather()->GetLChild() != NULL && treesLeaves[index][i]->GetFather()->GetRChild() != NULL
                && treesLeaves[index][i]->GetFather()->GetLChild()->GetLChild() == NULL && treesLeaves[index][i]->GetFather()->GetLChild()->GetRChild() == NULL
                && treesLeaves[index][i]->GetFather()->GetRChild()->GetLChild() == NULL && treesLeaves[index][i]->GetFather()->GetRChild()->GetRChild() == NULL)
        {   
            size++;
            i++;
        }   
        i++;
    }   
    return size;
}

/**
 * @param index l'arbre qui est traite
 **/
double GetLogILik (const int & index) throw ()
{
	double Lx = 0.0;
	for (unsigned int i = 0; i < treesLeaves[index].size (); i++)
		Lx += treesLeaves[index][i]->GetLogILik(sigma_hat_j_inv, params->tau[index]);
	return Lx;
}

bool Grow (const int & index, const int & node) throw ()
{
	Node * n = treesLeaves[index][node];
    bool grown = n->Grow (); 
    if (grown && n->GetLChild()->GetNinds() >= 1 && n->GetRChild()->GetNinds () >= 1) // Grow could be accepted
    {   
		double PGn = params->alpha * pow ((1 + n->GetDepth()), params->beta2 * -1); // Prior probability of node to grow
        //cout << "PGn: " << PGn << endl;
        double Lx = GetLogILik(index); // Log likelihood of tree (before growing)
        //cout << "Lx: " << Lx << endl;
    
        double Pbot = 1.0 / (double) GetToGrow(index); // Probability for this node to grow
		//cout << "Pbot " << Pbot << endl;
        double PBx = (n->GetFather () == NULL ? 1 : params->pgrow); // Probability of grow move
		double PGl = params->alpha * pow ((1 + n->GetLChild()->GetDepth()), params->beta2 * -1);
		double PGr = params->alpha * pow ((1 + n->GetRChild()->GetDepth()), params->beta2 * -1);
        //cout << "PGl: " << PGl << " " << endl;
        //cout << "PGr: " << PGr << endl;

        // Log-likelihood of tree after growing
		double Ly = Lx - n->GetLogILik (sigma_hat_j_inv, params->tau[index]) + n->GetLChild()->GetLogILik (sigma_hat_j_inv, params->tau[index])
			+ n->GetRChild()->GetLogILik(sigma_hat_j_inv, params->tau[index]);
        //cout << "Ly: " << Ly << endl;

        double Pnog = 1.0 / (double) (BeforeLeaves(index) + 1); // Probability of node already grown to be pruned
        //cout << "Pnog: " << Pnog << endl;

        double PDy = params->pprune; // Probability of step prune
        //cout << "PDy: " << PDy << endl;

        // Probability of move
        double alpha1 = (PGn * (1.0 - PGl) * (1.0 - PGr) * PDy * Pnog) / ((1.0 - PGn) * PBx * Pbot);
        //cout << "alpha1: " << alpha1 << endl;
        double alpha2 = alpha1 * exp(Ly - Lx);
        //cout << "alpha2: " << alpha2 << endl;
        double alpha = min(1.0, alpha2);
        //cout << "alpha: " << alpha << endl;

        double v = rand() / double(RAND_MAX);
        if (v < alpha)
        {   
            //cout << "=> GROW MOVE IS ACCEPTED" << endl;
			if (current_loop > params->nb_loop_train) // store this result
			{
				results->current_predictors[index][n->GetRule()]++;
				if (results->current_predictor_values[index][n->GetRule()].find(n->GetValue()) == results->current_predictor_values[index][n->GetRule()].end())
					results->current_predictor_values[index][n->GetRule()][n->GetValue()] = 1;
				else
					results->current_predictor_values[index][n->GetRule()][n->GetValue()]++;

				results->AddAll(index);
			}
			n->ClearIndividuals ();
            return true;
        }   
        else
        {   
            //cout << "=> GROW MOVE IS REFUSED v=" << v << " alpha= " << alpha << endl;
            n->Prune (false); 
			if (current_loop > params->nb_loop_train) // store this result
				results->AddAll(index);
        }   
    }   
    else
    {
		if (current_loop > params->nb_loop_train) // store this result
			results->AddAll(index);
		if (grown)
			n->Prune(false);
	}
	return false;
}

bool Prune (const int & index, const int & node) throw ()
{
	Node * n = treesLeaves[index][node]->GetFather(); // We will prune this node

	int rindex = n->GetRule ();
	double rvalue = n->GetValue ();

	double PDx = params->pprune; // Probability of step prune
	//cout << "PDx: " << PDx << endl;
	
	double Pnog = 1 / (double) (BeforeLeaves(index) + 1); // Probability of node to prune to be chosen
	//cout << "Pnog: " << Pnog << endl;
	
	double PGl = params->alpha * pow ((1 + n->GetLChild()->GetDepth()), params->beta2 * -1);
	double PGr = params->alpha * pow ((1 + n->GetRChild()->GetDepth()), params->beta2 * -1);
	//cout << "PGl: " << PGl << endl;
	//cout << "PGr: " << PGr << endl;
	
	double Lx = GetLogILik(index); // Log likelihood of tree
	//cout << "Lx: " << Lx << endl;
		
	Node * node_left = new Node (n->GetLChild(), true);
	Node * node_right = new Node (n->GetRChild(), true);

	std::vector <int> tmp;
	for (unsigned int i = 0; i < treesLeaves[index].size(); i++)
	{
		if (treesLeaves[index][i]->GetIndividualAt(0) == node_left->GetIndividualAt(0))
			tmp.push_back(i);
		else if (treesLeaves[index][i]->GetIndividualAt(0) == node_right->GetIndividualAt(0))
			tmp.push_back(i);
	}
	n->Prune (true);

	treesLeaves[index].erase (treesLeaves[index].begin()+(tmp[0] > tmp[1] ? tmp[0] : tmp[1]));//treesLeaves[index].begin() + (isLeft ? node+1 : node));
	treesLeaves[index].erase (treesLeaves[index].begin()+(tmp[0] > tmp[1] ? tmp[1] : tmp[0]));//treesLeaves[index].begin() + (isLeft ? node+1 : node));

	// Log-likelihood of tree
	double Ly = Lx + n->GetLogILik (sigma_hat_j_inv, params->tau[index]) - 
		node_right->GetLogILik(sigma_hat_j_inv, params->tau[index]) - node_left->GetLogILik(sigma_hat_j_inv, params->tau[index]);
	//cout << "Ly: " << Ly << endl;
	
	double PBy = params->pgrow;
    //cout << "PBy: " << PBy << endl;
     
    double PGn = params->alpha * pow ((1 + n->GetDepth()), params->beta2 * -1); // Prior probability of node to grow
    //cout << "PGn: " << PGn << endl;
    
	double Pbot = 1.0 / (double)(GetToGrow(index) + 1); // Probability for this node to grow
	//cout << "Pbot: " << Pbot << endl;//" // " << GetToGrow(index) << " " << canRGrow << " " << canLGrow << " " << treesLeaves[index].size() << endl;

	// Probability of move
	double alpha1 =((1.0-PGn)*PBy*Pbot)/(PGn*(1.0-PGl)*(1.0-PGr)*PDx*Pnog);
	double alpha2 = alpha1*exp(Ly-Lx);
	//cout << "alpha2: " << alpha2 << endl;
	double alpha = min(1.0,alpha2);
	//cout << "alpha: " << alpha << endl;

	double v = rand() / double(RAND_MAX);
	if (v < alpha)
	{
		//cout << "=> PRUNE MOVE IS ACCEPTED" << endl;
		if (current_loop > params->nb_loop_train) // store this result
		{
			results->current_predictors[index][rindex]--;
			results->current_predictor_values[index][rindex][rvalue]--;
			results->AddAll(index);
		}
		delete node_left;
		delete node_right;
		treesLeaves[index].push_back(n);
		return true;
	}
	else
	{
		//cout << "=> PRUNE MOVE IS REFUSED v=" << v << " alpha= " << alpha << endl;
		n->SetLChild (node_left);
		n->SetRChild (node_right);
		n->SetRule (rindex);
		n->SetValue (rvalue);
		treesLeaves[index].push_back(n->GetLChild ());
		treesLeaves[index].push_back(n->GetRChild ());
		n->ClearIndividuals ();
		if (current_loop > params->nb_loop_train) // store this result
			results->AddAll(index);
	}
	return false;
}

bool Change (Node * n, const int & index) throw ()
{
	Node * changed = n->Change ();
	/*cout << " ---------------- BEFORE ------------------" <<endl;
	n->PrintTree();
	cout << " ------------------------------------------" <<endl;
	changed->PrintTree();
	cout << " ------------------------------------------" <<endl;*/

	if (changed == NULL)
	{
		//cout << "CHANGED RETURNED FALSE" << endl;
		return false;
	}

	n->ReassembleIndividuals(true, true);
	double XLogPi = n->LogPriT(params->alpha, params->beta2); // Get log-likelihood and log prior of old tree
	n->DisassembleIndividuals ();
	//cout << "XLogPi: " << XLogPi << endl;
	double XLogL = GetLogILik (index);
	//cout << "XLogL: " << XLogL << endl;

	// Make the change effective
	if (n->GetFather() != NULL && n->GetIsLeftChild ())
		n->GetFather()->SetLChild (changed);
	else if (n->GetFather () != NULL)
		n->GetFather()->SetRChild (changed);

	double YLogPi = changed->LogPriT(params->alpha, params->beta2); // Get log-likelihood and log prior of new tree
	changed->DisassembleIndividuals();
	//cout << "YLogPi: " << YLogPi << endl;
    
	std::vector <Node *> new_leaves;
	if (n->GetFather () == NULL)
		changed->GetLeaves(new_leaves);
	else
		currentTrees[index]->GetLeaves (new_leaves);
	double YLogL = 0.0;
    for (unsigned int i = 0; i < new_leaves.size(); i++)
		YLogL += new_leaves[i]->GetLogILik(sigma_hat_j_inv, params->tau[index]);
	//cout << "YLogL: " << YLogL << endl;

	// Compute probability of move
	//cout << "alpha2: " << exp(YLogPi+YLogL-XLogPi-XLogL) << endl;
	double alpha = min(1.0,exp(YLogPi+YLogL-XLogPi-XLogL));
	//cout << "alpha: " << alpha << endl;
	double v = rand() / double(RAND_MAX);
	if (v < alpha)
	{
		//cout << "=> CHANGE MOVE IS ACCEPTED" << endl;
		if (n->GetFather() == NULL)
			currentTrees[index] = changed;

		treesLeaves[index].clear ();
		for (unsigned int i = 0; i < new_leaves.size(); i++)
			treesLeaves[index].push_back(new_leaves[i]);
		if (current_loop > params->nb_loop_train) // store this result
		{
			results->current_predictors[index][n->GetRule()]--;
			results->current_predictors[index][changed->GetRule()]++;
			results->current_predictor_values[index][n->GetRule()][n->GetValue()]--;
			results->current_predictor_values[index][changed->GetRule()][changed->GetValue()]++;
			results->AddAll(index);
		}
	
		delete n;
		return true;
	}
	else
	{
		//cout << "=> CHANGE MOVE IS REFUSED" << endl;
		delete changed;
		
		if (n->GetFather() != NULL && n->GetIsLeftChild ())
			n->GetFather()->SetLChild (n);
		else if (n->GetFather() != NULL)
			n->GetFather()->SetRChild (n);
		if (current_loop > params->nb_loop_train) // store this result
			results->AddAll(index);
	}
	return false;
}

bool Swap (Node * n, Node * OldNode, const int & index) throw ()
{
	// Special case : the node has already been swaped,
	// so only compute probabilities
	OldNode->ReassembleIndividuals(true, true);
	double XLogPi = OldNode->LogPriT( params->alpha, params->beta2); // Compute log-likelihood and log of prior of old tree
	OldNode->DisassembleIndividuals();
	//cout << "XLogPi: " << XLogPi << endl;
	double XLogL = GetLogILik (index);
	//cout << "XLogL: " << XLogL << endl;
	
	// Make the change effective
	if (OldNode->GetFather() != NULL && OldNode->GetIsLeftChild ())
		OldNode->GetFather()->SetLChild (n);
	else if (OldNode->GetFather () != NULL)
		OldNode->GetFather()->SetRChild (n);

	double YLogPi = n->LogPriT(params->alpha, params->beta2); // Get log-likelihood and log prior of new tree
	n->DisassembleIndividuals();
	//cout << "YLogPi: " << YLogPi << endl;
	
	std::vector <Node *> new_leaves;
	if (n->GetFather () == NULL)
		n->GetLeaves(new_leaves);
	else
		currentTrees[index]->GetLeaves (new_leaves);

	double YLogL = 0.0;
    for (unsigned int i = 0; i < new_leaves.size(); i++)
		YLogL += new_leaves[i]->GetLogILik(sigma_hat_j_inv, params->tau[index]);
	//cout << "YLogL: " << YLogL << endl;

	//cout << "alpha2: " << exp(YLogPi+YLogL-XLogPi-XLogL) << endl; // Compute probability of move
	double alpha = min(1.0,exp(YLogPi+YLogL-XLogPi-XLogL));
	//cout << "alpha: " << alpha << endl;
	double v = rand() / double(RAND_MAX);
	if (v < alpha)
	{
		//cout << "=> SWAP MOVE IS ACCEPTED" << endl;
		if (OldNode->GetFather() == NULL)
			currentTrees[index] = n;

		treesLeaves[index].clear ();
		for (unsigned int i = 0; i < new_leaves.size(); i++)
			treesLeaves[index].push_back(new_leaves[i]);

		if (current_loop > params->nb_loop_train) // store this result
		{
			results->current_predictors[index][OldNode->GetRule()]--;
			results->current_predictors[index][n->GetRule()]++;
			results->current_predictor_values[index][OldNode->GetRule()][OldNode->GetValue()]--;
			results->current_predictor_values[index][n->GetRule()][n->GetValue()]++;
	
			if (n->GetRChild()->GetRule() != -1)
			{	
				results->current_predictors[index][OldNode->GetRChild()->GetRule()]--;
				results->current_predictors[index][n->GetRChild()->GetRule()]++;
				results->current_predictor_values[index][OldNode->GetRChild()->GetRule()][OldNode->GetRChild()->GetValue()]--;
				results->current_predictor_values[index][n->GetRChild()->GetRule()][n->GetRChild()->GetValue()]++;
			}
			
			if (n->GetLChild()->GetRule() != -1)
			{	
				results->current_predictors[index][OldNode->GetLChild()->GetRule()]--;
				results->current_predictors[index][n->GetLChild()->GetRule()]++;
				results->current_predictor_values[index][OldNode->GetLChild()->GetRule()][OldNode->GetLChild()->GetValue()]--;
				results->current_predictor_values[index][n->GetLChild()->GetRule()][n->GetLChild()->GetValue()]++;
			}

			results->AddAll(index);
			//currentTrees[index]->InitResults(results, index);
		}
	
		delete OldNode;
	}
	else
	{
		//cout << "=> SWAP MOVE IS REFUSED" << endl;
		delete n;
		
		if (OldNode->GetFather() != NULL && OldNode->GetIsLeftChild ())
			OldNode->GetFather()->SetLChild (OldNode);
		else if (OldNode->GetFather() != NULL)
			OldNode->GetFather()->SetRChild (OldNode);
			
		if (current_loop > params->nb_loop_train) // store this result
			//currentTrees[index]->InitResults(results, index);
			results->AddAll (index);
	}
	return false;
}

void Metrop (const int & index) throw ()
{
	// jth column without the jth element
	sigma_tmp_1.resize (params->sigma_j.size1()-1, 1);
	//cout << "sigma1 " << sigma_tmp_1 << endl;
	for (unsigned int i = 0; i < params->sigma_j.size1(); i++)
		if (i != (unsigned)index)
			sigma_tmp_1(i < (unsigned)index ? i : i-1, 0) = params->sigma_j(i, index);
	//cout << "sigma1 " << sigma_tmp_1 << endl;

	// without jth column and jth line
	matrix <double> sigma2 (params->sigma_j.size1()-1, params->sigma_j.size2()-1);
	sigma_tmp_2_inv.resize (params->sigma_j.size1()-1, params->sigma_j.size2()-1);
	for (unsigned int i = 0; i < params->sigma_j.size1(); i++)
		for (unsigned int j = 0; j < params->sigma_j.size2(); j++)
			if (i != (unsigned)index && j != (unsigned)index)
				sigma2 (i < (unsigned)index ? i : i-1, j < (unsigned)index ? j : j-1) = params->sigma_j(i,j);

	//cout << "sigma2 " << sigma2 << endl;
	InvertMatrixChol (sigma2, sigma_tmp_2_inv);
	//cout << "sigma2_inv " << sigma_tmp_2_inv << endl;
	//cout << "tau " << params->tau[index] << endl;

	sigma2 = prod (sigma_tmp_2_inv, sigma_tmp_1);
	sigma2 = prod (trans(sigma_tmp_1), sigma2);
	
	BNU::symmetric_matrix<double, BNU::lower> sigma_hat_j (params->nb_dose);
	sigma_hat_j = (params->sigma_j (index, index) - sigma2(0,0)) * kron (params->sigma_d, params->sigma_t);
	//cout << "SIGMA HAT J " << sigma_hat_j << endl;
	
	InvertMatrixChol (sigma_hat_j, sigma_hat_j_inv);
	//cout << "sigma_hat_j_inv " << sigma_hat_j_inv << endl;

	//cout << treesLeaves[index].size() << " " << treesLeaves[index][0]->GetFather () << endl;
	if (treesLeaves[index].size() == 1 && treesLeaves[index][0]->GetFather () == NULL) // current node is a root, this is the first move
	{
		//cout << "#### GROW RACINE " << treesLeaves[index].size() << endl;

		if (Grow (index, 0))
		{
			treesLeaves[index].pop_back(); // The only node present in the vector is the root of the tree
			//cout << "GROW ACCEPTED" << endl;
			treesLeaves[index].push_back(currentTrees[index]->GetLChild());
			treesLeaves[index].push_back(currentTrees[index]->GetRChild());
		}
		else
			GetLogILik (index); // Just to set beta
	}
	else // randomly choose 
	{
		double u = rand() / double(RAND_MAX);
		//cout << "RANDOM CHOISE " << u << endl;
			
		if(u < params->pgrow)
		{
			int i = GetRandomGrow (index);
			//cout << "#### GROW RANDOM " << i << endl;
			if (i != -1 && Grow (index, i))
			{
				//cout << "GROW ACCEPTED " << endl;
				treesLeaves[index].push_back(treesLeaves[index][i]->GetLChild());
				treesLeaves[index].push_back(treesLeaves[index][i]->GetRChild());
				treesLeaves[index].erase(treesLeaves[index].begin() + i); // The only node present in the vector is the root of the tree
			}
		}
		else if (u < (params->pgrow + params->pprune)) // prune movement
		{
			int i = GetRandomPrune (index);
			//cout << "#### PRUNE RANDOM " << i << " size before " << treesLeaves[index].size() << endl;
			Prune (index, i);
		}
		else if (u < (params->pgrow + params->pprune + params->pchange)) // change movement
		{
			//cout << "#### CHANGE" << endl;
			//currentTrees[index]->PrintTree();
			Node * n = GetRandomChange(index);
			//n->Print ();
			Change (n, index);
			//currentTrees[index]->PrintTree();
		}
		else
		{
			//cout << "#### SWAP" << endl;
			//currentTrees[index]->PrintTree();
			Node * OldNode;
			Node * n = GetRandomSwap (index, &OldNode);
			
			if (n != NULL)
			{
				//cout << "SWAP ACCEPTED" << endl;
				Swap (n, OldNode, index);
			}
			else if (current_loop > params->nb_loop_train)
				results->AddAll(index);
		}
	}
}

void ComputeTreeResiduals (const int & index) throw ()
{
	for (unsigned int i = 0; i < treesLeaves[index].size (); i++)
	{
		for (int r = 0; r < params->nb_repeat; r++)
			treesLeaves[index][i]->ComputeResiduals (s_res_vect, index, r);
	}
}

/**
 * m_kron = kroneker, ne change pas
 * Correspond a page 5 du old, derniere formule
 **/
void ComputePsi (const matrix <double> & m_s_res, double & phi, const double & special) throw ()
{
	//matrix <double> psi = prod (m_s_res, m_kron);
	//psi = prod (psi, trans (m_s_res)); // Ici il va falloir une boucle. Il faut remplacer m_s_res

	double lambda1 = 0.0, lambda2 = 0.0 , lambda3 = 0.0;

	for (unsigned int i = 0; i < m_s_res.size1(); i++)
	{
		lambda1 += m_s_res(i,i);
		if (i < m_s_res.size1()-1)
			lambda2 += m_s_res(i,i+1) + m_s_res(i+1,i);
	}
	lambda3 = lambda1 - m_s_res(0,0) - m_s_res(m_s_res.size1()-1, m_s_res.size2()-1) + params->lambda3;
	lambda1 += params->lambda1;
	lambda2 += params->lambda2;

	Mh (lambda1, lambda2, lambda3, phi, special);
}

void PrintSResVect (void) throw ()
{
	for (int i = 0; i < params->nb_indiv; i++)
		for (int j = 0; j < params->nb_repeat; j++)
			cout << "(" << i << "," << j << ") -- " << s_res_vect [(i * params->nb_repeat) + j] << endl;
}

/**
 * Compute m_kron = inv_sigma_d kron inv_sigma_t. m_kron ne doit pas bouger
 **/
void ComputeResiduals () throw ()
{
	//cout << "----------------------- RESIDUS ----------------------- " << endl;
	//PrintSResVect ();
	//cout << "------------------------------------------------------- " << endl;
	// Compute new sigma_j
	matrix <double> m_s_res (params->nb_assay, params->nb_assay);
	for (unsigned int i = 0; i < m_s_res.size1(); i++)
		for (unsigned int j = 0; j < m_s_res.size2(); j++)
			m_s_res (i, j) = 0.0;
	matrix <double> m_kron = kron (InvARMatrix (params->sigma_d.size1(), params->sigma_d.size2(), params->phi_d), InvARMatrix (params->sigma_t.size1(), params->sigma_t.size2(), params->phi_t));

	for (int i = 0; i < params->nb_indiv; i++)
	{
		for (int k = 0; k < params->nb_repeat; k++)
		{
			matrix <double> tmp = prod (trans (s_res_vect[(i * params->nb_repeat) + k]), m_kron);
			tmp = prod (tmp, s_res_vect[(i * params->nb_repeat) + k]);
			AddMatrix (m_s_res, tmp);
		}
	}

	//cout << "SOMME 1 " << m_s_res << endl;
	m_s_res += params->lambda;
	int param_xi = params->xi + (params->nb_indiv * params->nb_dose * params->nb_repeat * params->nb_time);
	symmetric_matrix <double> tmp_sigma_j_inv (params->nb_assay, params->nb_assay);
	RInvertWishart (params->sigma_j, tmp_sigma_j_inv, m_s_res, param_xi);
	InvertMatrixChol (params->sigma_j, tmp_sigma_j_inv);
	//cout << "SIGMA J " << params->sigma_j << endl;
	//cout << "SIGMA J INV " << tmp_sigma_j_inv << endl;

	if (current_loop > params->nb_loop_train)
		results->sigmas.push_back (params->sigma_j);
	//cout << "Nouveau sigma_j " << params->sigma_j << endl;
	//InvertMatrix (tmp_sigma_j_inv, params->sigma_j);

	//matrix <double> psi = prod (tmp_sigma_j_inv, trans(s_res));
	//psi = prod (s_res, psi);

	// Compute phi_d
	m_s_res.resize (params->nb_dose, params->nb_dose);
	for (unsigned int i = 0; i < m_s_res.size1(); i++)
		for (unsigned int j = 0; j < m_s_res.size2(); j++)
			m_s_res (i, j) = 0.0;
	m_kron = kron (tmp_sigma_j_inv, InvARMatrix (params->sigma_t.size1(), params->sigma_t.size2(), params->phi_t));
	//cout << "NOUVEAU KRON " << m_kron << endl;
	for (int i = 0; i < params->nb_indiv; i++)
	{
		for (int k = 0; k < params->nb_repeat; k++)
		{
			matrix < double > tmp (params->nb_dose, params->nb_assay * params->nb_time);
			for (int j = 0; j < params->nb_assay; j++)
				for (int t = 0; t < params->nb_time; t++)
					for (int d = 0; d < params->nb_dose; d++)
						tmp (d, (j * params->nb_time) + t) = s_res_vect[(i * params->nb_repeat) + k] ((d * params->nb_time) + t, j);
			//cout << "Mat reorganisee pour i=" << i << ", k=" << k << " -> " << tmp << endl;
			matrix < double > tmp2 = matrix <double> (tmp);
			tmp2 = prod (tmp, m_kron);
			tmp2 = prod (tmp2, trans (tmp));
			AddMatrix (m_s_res, tmp2);
		}
	}
	//cout << "SOMME 2 " << m_s_res << endl;
	ComputePsi (m_s_res, params->phi_d, params->nb_time);
	results->phi_d.push_back(params->phi_d);

	if (current_loop > params->nb_loop_train)
		results->s_res.push_back (m_s_res);
	//cout << "PSI_D = " << params->phi_d << endl;

	// compute phi_t
	if (params->nb_time >1)
	{
		m_s_res.resize (params->nb_time, params->nb_time, false);
		m_kron = kron (tmp_sigma_j_inv, InvARMatrix (params->sigma_d.size1(), params->sigma_d.size2(), params->phi_d));
		for (unsigned int i = 0; i < m_s_res.size1(); i++)
			for (unsigned int j = 0; j < m_s_res.size2(); j++)
				m_s_res (i, j) = 0.0;
		for (int i = 0; i < params->nb_indiv; i++)
		{
			for (int k = 0; k < params->nb_repeat; k++)
			{
				matrix < double > tmp (params->nb_time, (params->nb_assay * params->nb_dose));
				for (int j = 0; j < params->nb_assay; j++)
					for (int d = 0; d < params->nb_dose; d++)
						for (int t = 0; t < params->nb_time; t++)
							tmp (t, (j * params->nb_dose) + d) = s_res_vect[(i * params->nb_repeat) + k] ((d * params->nb_time) + t, j);
				//cout << "Mat reorganisee pour i=" << i << ", k=" << k << " -> " << tmp << endl;
				matrix < double > tmp2 = matrix <double> (tmp);
				tmp2 = prod (tmp, m_kron);
				tmp2 = prod (tmp2, trans (tmp));
				AddMatrix (m_s_res, tmp2);
			}
		}
		//cout << "SOMME 3 " << m_s_res << endl;
	
		ComputePsi (m_s_res, params->phi_t, params->nb_dose);
		results->phi_t.push_back(params->phi_t);
		//cout << "PSI_T = " << params->phi_t << endl;
	}
	//cout << params->phi_d << endl;
	//cout << " sigma_j " << params->sigma_j << endl;
}

void ComputeR (const int & index) throw ()
{
	sigma_tmp_1.resize (params->sigma_j.size1()-1, 1);
	for (unsigned int i = 0; i < params->sigma_j.size1(); i++)
		if (i != (unsigned)index)
			sigma_tmp_1(i < (unsigned)index ? i : i-1, 0) = params->sigma_j(i, index);
	//cout << "sigma1 " << sigma_tmp_1 << endl;

	// without jth column and jth line
	matrix <double> sigma2 (params->sigma_j.size1()-1, params->sigma_j.size2()-1);
	sigma_tmp_2_inv.resize (params->sigma_j.size1()-1, params->sigma_j.size2()-1);
	for (unsigned int i = 0; i < params->sigma_j.size1(); i++)
		for (unsigned int j = 0; j < params->sigma_j.size2(); j++)
			if (i != (unsigned)index && j != (unsigned)index)
				sigma2 (i < (unsigned)index ? i : i-1, j < (unsigned)index ? j : j-1) = params->sigma_j(i,j);
	//cout << "SIGMA2 " << sigma2 << endl;
	InvertMatrixChol (sigma2, sigma_tmp_2_inv);
	//cout << "sigma2_inv " << sigma_tmp_2_inv << endl;
	//cout << "LAAAAA " << sigma_tmp_1 << " " << sigma_tmp_2_inv << endl;
	sigma2 = prod (trans(sigma_tmp_1), sigma_tmp_2_inv);
	//cout << "kronecker " << sigma2 << endl;
	//cout << endl;

	std::vector < Individual *> new_database;
	new_database.resize (params->nb_indiv);
	for (unsigned int i = 0; i < new_database.size(); i++)
		new_database[i] = new Individual (params, i);

	for (int i = 0; i < params->nb_assay; i++)
	{
		if (i != index)
		{
			for (unsigned int j = 0; j < treesLeaves[i].size(); j++)
			{
				//cout << "(0, " << (i < index ? i : i-1) << ") i=" << i << ", j=" << j << endl;
				treesLeaves[i][j]->ComputeNu (new_database, sigma2(0, i < index ? i : i-1));
			}
		}
	}

	//cout << "------ DISPLAY NEW DB FOR " << index << endl;
	//for (unsigned int i = 0; i < new_database.size(); i++)
	//	cout << new_database[i] << endl;
	//cout << "------------------ " << endl;
	//DisplayDatabase ();

	for (unsigned int i = 0; i < treesLeaves[index].size(); i++)
		treesLeaves[index][i]->ComputeR (new_database, database, index);

	for (unsigned int i = 0; i < new_database.size (); i++)
		delete new_database[i];
}

void SmoothParameters (const int & index) throw ()
{
	//cout << "---------------------------------------- PRINT TREE -------------------------------------" << endl;
	//currentTrees[index]->PrintTree();
	//cout << "-----------------------------------------------------------------------------------------" << endl;
	double a = params->a1 + ((double)(treesLeaves[index].size() * (params->nb_dose + params->delta) * (params->nb_time > 1 ? (params->nb_time + params->delta) : 1)/2.0));
	//cout << "a " << a << endl;
	double b = 0.0;
	for (unsigned int i = 0; i < treesLeaves[index].size(); i++) // For each terminal leaves
	{
		//treesLeaves[index][i]->Print();
		BNU::vector <double> tmp = prod (trans(treesLeaves[index][i]->GetBeta()), params->penalty_matrix);
		double x = 0.0;
		for (unsigned int j = 0; j < treesLeaves[index][i]->GetBeta().size(); j++)
			x += (treesLeaves[index][i]->GetBeta())(j) * tmp (j);
		b += 0.5 * x;
	}
	b += params->b1;
	//cout << "b " << b << endl;
	params->tau[index] = RGamma (a, 1.0/b); // TODO change by rgamma here
	//cout << "tau[" << index << "]" << params->tau[index] << endl;

	if (current_loop > params->nb_loop_train) // Store tau
	{
		if (index == 0) // new loop
		{
			std::vector <double> t;
			results->taus.push_back (t);
		}
		results->taus[current_loop - params->nb_loop_train - 1].push_back(params->tau[index]);
	}
	//cout << "tau [" << index << "]" << params->tau[index] << endl;
}

void Run (void) throw ()
{
	for (current_loop = 0; current_loop < params->nb_loop; current_loop++)
	{
		if (current_loop%10 == 0)
			cout << "\n/********************************** ITER " << current_loop << " **********************************/" << endl;

		//cout << "\n/********************************** ITER " << i << " **********************************/" << endl;
		//for (unsigned int _i = 0; _i < s_res.size1(); _i++)
		//	for (unsigned int _j = 0; _j < s_res.size2(); _j++)
		//		s_res (_i, _j) = 0.0;

		for (int j = 0; j < params->nb_assay; j++)
		{
			if (current_loop > 0)
			{
				//currentTrees[j]->PrintTree();
				//cout << "Compute R " << i << " " << j << endl;
				ComputeR (j);
				//currentTrees[j]->PrintTree();
			}
			//cout << " ------------------- METROP ARBRE " << j << endl;
			Metrop (j);
			SmoothParameters (j);
			ComputeTreeResiduals(j);
			//currentTrees[j]->PrintTree();
			//cout << " -------------------------FIN ABR " << j << " --------------- " << endl << endl;
			if (current_loop == params->nb_loop_train) // First init
			{
				currentTrees[j]->InitResults (results, j);
				results->AddAll (j);
			}
			if (current_loop > params->nb_loop_train && params->nb_test_step > 0 && (current_loop-params->nb_loop_train)%params->nb_test_step == 0)
			{
				for (int i = 0; i < params->nb_test_indiv; i++)
					currentTrees[j]->TestBetaMatrix (results, test_predictors, i);

				if (results->beta_matrixes.size() > 100)
					results->FlushMatrixes();
			}
			if (current_loop > params->nb_loop_train)
			{
				if (j == 0)
				{
					std::vector <int> t;
					results->nb_sons.push_back (t);
				}
				results->nb_sons[current_loop - params->nb_loop_train -1].push_back (treesLeaves[j].size());
			}
		}
		ComputeResiduals ();
	}
}

void BuiltTestTree (void) throw ()
{
	currentTrees[0]->Grow (2, 2);
	currentTrees[0]->ClearIndividuals ();
	currentTrees[0]->GetRChild()->Grow (0,2);
	//currentTrees[0]->GetRChild()->GetLChild()->Grow (0,6);
	currentTrees[0]->GetRChild()->ClearIndividuals ();
	//currentTrees[0]->GetRChild()->GetLChild()->ClearIndividuals ();
	currentTrees[0]->GetLChild()->Grow (1,2);
	currentTrees[0]->GetLChild()->ClearIndividuals ();

	treesLeaves[0].push_back(currentTrees[0]->GetRChild()->GetRChild());
	//treesLeaves[0].push_back(currentTrees[0]->GetRChild()->GetLChild());
	//treesLeaves[0].push_back(currentTrees[0]->GetLChild()->GetRChild());
	//treesLeaves[0].push_back(currentTrees[0]->GetLChild()->GetRChild());
}

int main( int argc, char* argv [])
{
	if (argc < 5)
	{
		cerr << "Usage: " << argv[0] << " <-f filename><-i init_file><-o file_out>" << endl;
		//DisplayHelp ();
		return -1;
	}

	if (!InitParam (argv))
		return -1;

	InitSigmaD ();
	InitSResVect ();
	cout << params << endl;
	//DisplayDatabase ();
	//cout << "------------------------------------------------" << endl;

	srand(time(0));

	Run ();

	/*for (int i = 0; i < params->nb_assay; i++)
	{
		cout << " -------------------- ARBRE " << i << " ------------------------ " << endl;
		currentTrees[i]->PrintTree();
		cout << " -------------------------------------------------- " << endl;
	}*/

	cout << results << endl;

	//treesLeaves[0].pop_back();
	//BuiltTestTree ();
	//currentTrees[0]->PrintTree();

	cout << "CHECKING RESULTS COHERENCY" << endl;
	results->CheckResults();
	cout << "DONE" << endl;

	for (int i = 0; i < params->nb_assay; i++)
	{
		for (unsigned int j = 0; j < treesLeaves[i].size(); j++) // clear leaves individuals
			treesLeaves[i][j]->DeleteIndividuals();
		delete currentTrees[i];
	}
	delete [] currentTrees;

	for (unsigned int i = 0; i < database.size(); i++)
		for (unsigned int j = 0; j < database[i].size(); j++)
			delete database[i][j];

	for (int i = 0; i < params->nb_test_indiv; i++)
		delete [] test_predictors[i];
	delete [] test_predictors;

	delete params;
	delete results;

	return 0;
}

#undef BNU
