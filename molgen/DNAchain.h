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

#ifndef __DNACHAIN_H__
#define __DNACHAIN_H__

struct nucleotide
	{
		//! The enum
		enum Enum
			{
			A = 0,
			G,
			C,
			T,
			};
	}; 
			
struct sites
	{
		//! The enum
		enum Enum
			{
			Su = 0,
			Ph,
			Ab,
			Gb,
			Cb,
			Tb,
			};
	}; 	

class DNAchain : public Molecule 
	{
	public:
    enum Strand
        {
        ss=0,
		ds,
        };
		DNAchain(unsigned int Nbp, DNAchain::Strand s);
		DNAchain(const std::string& fname, unsigned int NatomPerMole, DNAchain::Strand s);
		DNAchain(unsigned int Nbp, DNAchain::Strand s, std::string form);
		DNAchain(const std::string& fname, unsigned int NatomPerMole, DNAchain::Strand s, std::string form);

		 virtual ~DNAchain(){}
		void allocate_DNAdata();
		void setSequences(std::string types, std::string fname);
		void setSequences(std::string fname);
		void readSequences(std::string fname);
		void generateType();
		void generateTopology();
		void setScale(double scale)
			{
			m_scale = scale;
			}
		void setStartPoint(double x, double y, double z)
			{
			m_start_point = vec(x, y, z);
			m_set_start_point=true;
			}	
		void setDirection(double x, double y, double z)
			{
			m_direction = vec(x, y, z);
			m_set_direction=true;
			}				
		void generateSites();
		void placeSites();
		virtual void generate();		
	protected:   	
		unsigned int m_N_bps;
		unsigned int m_N_sites;	
		vector<std::string> m_sequence;
		vector<vec4> m_site_data1;					// the data of six sites
		vector<vec4> m_site_data2;					// the data of six sites
		vector<std::string> m_site_name;
		vector<vec> m_site_xyz;	
		vector<int> m_site_id;	
		double m_delt_degree;
		double m_delt_z;
		bool m_circle;
		double m_radius;
		double m_scale;
		Strand m_s;
		vector<vec> m_xyz_temp;
		vec m_start_point;
		vec m_direction;
		bool m_set_start_point;
		bool m_set_direction;
		std::string m_form;
		};
		

void export_DNAchain(pybind11::module &m);
#endif

  

