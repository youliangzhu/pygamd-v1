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
#include "Object.h"
#include "DNAchain.h"
using namespace std;

#ifndef __GENERATORS_H__
#define __GENERATORS_H__

class Generators
	{
	public:
		//! Constructs the compute
		Generators(double Lx, double Ly, double Lz);	
		//! Destructor
		virtual ~Generators();
		 
        void addMolecule(Molecule* Mol, unsigned int Nm);
		void setMinimumDistance(double mini_dis);
		void setMinimumDistance(const std::string& name1, const std::string& name2, double mini_dis);
		void outPutXml(std::string fname);
		void outPutMol2(std::string fname);	
		void outPutMST(std::string fname);		
        void updatePos(std::vector<vec>& pos, std::vector<std::string>& type, std::vector<vec>& ori, std::vector<vec4>& quat,
						std::vector<unsigned int>& body, std::vector<unsigned int>& molecule);
		void initiateList();
		unsigned int cellid(int i, int j, int k);
		void setPrecision(unsigned int npre)
			{
			m_nprecision = npre;
			}	

		void setDimension(unsigned int ndn)
			{
			m_dimension = ndn;
			}			
			
		void setHead(unsigned int nhead)
			{
			m_nhead = nhead;
			}
		void setParam(const string& name1, const string& name2, double epsilon, double sigma, double rcut);
		unsigned int switchNametoType(const string& name);
		void generate();
    protected:	
        std::vector< Molecule* > m_molecules; 
        std::vector< unsigned int > m_Nmol;	
        std::vector<vec> m_pos_all;
        std::vector<unsigned int> m_body_all;
        std::vector<vec> m_ori_all;
        std::vector<vec4> m_quat_all;		
        std::vector<unsigned int> m_molecule_all;
		unsigned int m_N;
		unsigned int m_Num_bonds;
		unsigned int m_Num_mol;
		double m_Lx;
		double m_Ly;
		double m_Lz;
        std::vector<double> m_min_dis;
        std::vector<vec> m_params;		
        std::vector< string > m_type_mapping_all;
		double m_rcut_max;
		unsigned int m_NBtype;	
		
		vec_uint m_dim;
		vec m_width;
		unsigned int* m_list;
		unsigned int* m_head;
		
		double m_edge;
		double R2S();
		unsigned int m_dimension;
		unsigned int m_nprecision;
		unsigned int m_nhead;
		bool m_generated;
	};

void export_Generators(pybind11::module &m);		

#endif

