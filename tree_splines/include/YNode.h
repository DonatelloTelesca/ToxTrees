/*
 * YNode.h
 *
 *  Created on: 13 juil. 2011
 *      Author
 */

#ifndef __YNODE_H__
#define __YNODE_H__

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "Node.h"

namespace nsLBart
{
	class YNode : public Node 
	{
		private:
			int m_Object;

		public:
			YNode () throw ();
			YNode (const int & obj) throw ();
			//Node (Node * n, const bool & fill_obs) throw ();
			//Node (Node * father, const int & col, const int & depth, const bool & islc) throw ();
			virtual ~YNode() throw ();

			// -------------------------------------- //
			//                GETTERS                 //
			// -------------------------------------- //

			int    GetObject() throw ();

			// -------------------------------------- //
			//                SETTERS                 //
			// -------------------------------------- //

			void SetObject   (int i)    throw ();

			// -------------------------------------- //
			//                METHODS                 //
			// -------------------------------------- //

			//virtual void ReadFile(char * filename) throw ();
			virtual void PrintTree() throw ();
			double GetSum () throw ();
	};
}

#include "YNode.hxx"

#endif /* __YNODE_H__ */
