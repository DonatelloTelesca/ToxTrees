/*
 * observation.cpp
 *
 *  Created on: 13 juil. 2011
 *      Author
 */


#include "Observation.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace nsLBart;

#define OBS nsLBart::Observation

OBS::Observation (const int & length) throw ()
{
	m_Length = length;
	m_Array = new double [m_Length];
	m_Id = -1;
}

OBS::Observation (const int & length, const int & id) throw ()
{
	m_Length = length;
	m_Array = new double [m_Length];
	m_Id = id;
}

OBS::Observation (Observation * o) throw ()
{
	m_Length = o->GetLength ();
	m_Id = o->GetId ();
	m_Array = new double [m_Length];
	for (int i = 0; i < m_Length; i++)
		m_Array[i] = o->GetValueAt(i);
}

OBS::~Observation() throw()
{
	if (m_Array != NULL)
	{
		delete [] m_Array;
		m_Array = NULL;
	}
}

void OBS::SetLength (const int & l) throw ()
{
	delete [] m_Array;
	m_Length = l;
	m_Array = new double [m_Length];
}

void OBS::Print() throw ()
{
	cout << "id(" << m_Id << ") ";
	for (int i=0; i< m_Length; i++)
	{
		cout << m_Array[i]<< " ";
	}
	cout << endl;
}

#undef OBS
