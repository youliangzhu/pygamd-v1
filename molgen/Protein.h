/*
PYGAMD - Python GPU-Accelerated Molecular Dynamics Software
VERSION 1
COPYRIGHT
	PYGAMD Copyright (c) (2021) You-Liang Zhu, Zhong-Yuan Lu
LICENSE
	This program is a free software: you can redistribute it and/or 
	modify it under the terms of the GNU General Public License. 
	This program is distributed in the hope that it will be useful, 
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
	See the General Public License v3 for more details.
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
DISCLAIMER
	The authors of PYGAMD do not guarantee that this program and its 
	derivatives are free from error. In no event shall the copyright 
	holder or contributors be liable for any indirect, incidental, 
	special, exemplary, or consequential loss or damage that results 
	from its use. We also have no responsibility for providing the 
	service of functional extension of this program to general users.
USER OBLIGATION 
	If any results obtained with PYGAMD are published in the scientific 
	literature, the users have an obligation to distribute this program 
	and acknowledge our efforts by citing the paper "Y.-L. Zhu, H. Liu, 
	Z.-W. Li, H.-J. Qian, G. Milano, and Z.-Y. Lu, J. Comput. Chem. 2013,
	34, 2197-2211" in their article.
CORRESPONDENCE
	Dr. You-Liang Zhu
	Email: ylzhu@pygamd.com
*/

#include "Molecule.h"

using namespace std;

#ifndef __PROTEIN_H__
#define __PROTEIN_H__

class Protein : public Molecule 
	{
	public:

		Protein(const std::string& fname);

		 virtual ~Protein(){}
		void allocate_amino_acid_data();
		unsigned int readSequences(std::string fname);
		void generateType();
		void generateTopology();
		void generateSites();
		void readPos(std::string fname);
		unsigned int getIndex(std::string name);
		void assignTypes();
		virtual void generate();
	protected:   	
		unsigned int m_N_amino_acid;
		unsigned int m_cg_nmax;
		unsigned int m_n_amino_acid_types;
		unsigned int m_atom_nmax;		

		std::vector<std::string> m_name;
		std::vector<unsigned int> m_cg_np;
		std::vector<std::string> m_cg_ptype;
		std::vector<std::string> m_sequence;
		std::vector<unsigned int> m_natom;
		std::vector<std::string> m_read_atom_type;
		std::vector<unsigned int> m_read_atom_molid;
		std::vector<std::string> m_read_atom_molname;
		std::vector<vec>m_read_atom_pos;		
		};

void export_Protein(pybind11::module &m);
#endif

  

