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

#include "DCDBuilder.h"
#include <cuda_runtime.h>

#ifndef __MOLINFO_H__
#define __MOLINFO_H__

class MolInfo
    {
    public:
		MolInfo(XMLBuilder& build);
		
		void initialize();
		void buildMol();
		void buildmap();
		void bondExclude();
		void angleExclude();
		void dihedralExclude();
		void buildExclusionList();
		void buildExclusionGPUList();
		unsigned int switchNameToIndex(const std::string &name);
		
        const std::vector<unsigned int >& getMolIdPerParticle() const 
			{
			return m_mol_id_particle;
			}
        unsigned int* getMolIdPerParticleGPU() 
			{
			if(m_mol_id_particle_device)
				{
				cudaMalloc(&d_mol_id_particle, sizeof(unsigned int)*m_N);
				cudaMemcpy(d_mol_id_particle, &m_mol_id_particle[0], sizeof(unsigned int)*m_N, cudaMemcpyHostToDevice);
				m_mol_id_particle_device=false;
				}
			return d_mol_id_particle;
			}
        const std::vector<unsigned int >& getMolTypeId() const 
			{
			return m_mol_type_id;
			}
		unsigned int getKindNumofMol() const 
			{
			return m_mol_type_exchmap.size();
			}
		unsigned int getNumofMol() const 
			{
			return m_n_mol;
			}
        const std::vector<vec >& getPos0() const 
			{
			return m_pos0;
			}
        const std::vector<vec >& getPos() const 
			{
			return m_pos;
			}			
        const std::vector<string >& getMoltypeMap() const 
			{
			return m_mol_type_exchmap;
			}
        const std::vector<unsigned int>& getNmolperKind() const 
			{
			return m_n_mol_per_kind;
			}
       const std::vector<unsigned int>& getMolsize() const 
			{
			return m_mol_size;
			}
       const std::vector<uint_2>& getMolStartEnd() const 
			{
			return m_mol_stat_end;
			}
       const std::vector<unsigned int>& getList() const 
			{
			return m_list;
			}
       const std::vector<unsigned int>& getHead() const 
			{
			return m_head;
			}
       const vec& getWith() const 
			{
			return m_width;
			}
       const vec_uint& getDim() const 
			{
			return m_dim;
			}	
       const std::vector<vec_int>& getMap() const 
			{
			return m_map;
			}	

       const std::vector<unsigned int>& getNExclusion() 
			{
			buildExclusionList();
			return m_n_exclusion;
			}		
       const std::vector<unsigned int>& getExclusionList()
			{
			buildExclusionList();
			return m_exclusion_list;
			}	

		unsigned int* getNExclusionGPU()
			{
			buildExclusionGPUList();
			return d_n_exclusion;
			}		
		unsigned int* getExclusionListGPU()
			{
			buildExclusionGPUList();
			return d_exclusion_list;
			}
			
		void addFreeParticleType(const std::string& name) 
			{
			// search for the type map
			for (unsigned int i = 0; i < free_particle_types.size(); i++)
				{
				if (free_particle_types[i] == name)
					return;
				}
			// add a new one if it is not found
			free_particle_types.push_back(name);
			return;
			}
			
		unsigned int getFreeParticleTypeId(const std::string& name) 
			{
			// search for the type map
			for (unsigned int i = 0; i < free_particle_types.size(); i++)
				{
				if (free_particle_types[i] == name)
					return i;
				}
			return NO_INDEX;	
			}
			
		std::vector< std::string>& getFreeParticleTypes()
			{
			return free_particle_types;
			}
	
       bool ifexclude(unsigned int a, unsigned int b);
	   unsigned int cellid(int i, int j, int k);
       void computeList(double r_cut);	
	   void updatePosition0();
       void outPutInfo();	   
    private:
		XMLBuilder& m_build;
		std::vector<unsigned int > m_mol_id_particle;
		unsigned int* d_mol_id_particle;
		bool m_mol_id_particle_device;
		unsigned int m_n_mol;
		std::vector< std::string> free_particle_types;
		std::vector<vec> m_pos;
		BoxSize m_box;
		std::vector<std::string> m_mol_type;
		std::vector<unsigned int > m_mol_type_id;
		std::vector<std::string> m_mol_type_exchmap;
		std::vector<vec> m_pos0;
		std::vector<unsigned int > m_n_mol_per_kind;
		std::vector<unsigned int> m_mol_size;	
		std::vector<uint_2> m_mol_stat_end;
		std::vector<vec_int > m_map;		
		std::vector<unsigned int> m_exclusion_list;	
		std::vector<unsigned int> m_n_exclusion;
		unsigned int* d_exclusion_list;
		unsigned int* d_n_exclusion;
		unsigned int m_n_max_bond;
		unsigned int m_n_max_angle;
		unsigned int m_n_max_dihedral;
		unsigned int m_n_max_ex;
		std::vector<uint_2> m_bonds;
		std::vector<unsigned int> m_n_bond;
		unsigned int m_N;
		unsigned int m_dimension;
		bool m_bond_ex;
		bool m_angle_ex;
		bool m_dihedral_ex;
		bool m_build_ex_host_table;
		bool m_build_ex_device_table;
		
		vec_uint m_dim;
		vec m_width;
		std::vector<unsigned int> m_list;
		std::vector<unsigned int> m_head;
		
	};

#endif
	