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

#include <memory>
#include "XMLBuilder.h"

using namespace std;

#ifndef __MSTREADER_H__
#define __MSTREADER_H__

class MSTReader : public XMLBuilder
    {
    public:
        MSTReader();
		virtual bool readDataFromMST(const string &fname);
        virtual void outPutInfo();

    protected:
        void reset_params();
		void clear_data();

		bool m_mst_read;
		bool m_invariant_data;
		bool m_variant_data;
		bool m_num_particles_read;
		bool m_timestep_read;
		bool m_dimension_read;
		bool m_bond_read;
		bool m_angle_read;
		bool m_dihedral_read;
		bool m_vsite_read;		
		bool m_box_read;		
		bool m_position_read;
		bool m_type_read;
		bool m_image_read;
		bool m_mass_read;
		bool m_velocity_read;
		bool m_charge_read;
		bool m_body_read;
		bool m_diameter_read;		
		bool m_rotangle_read;
		bool m_force_read;
		bool m_virial_read;			
		bool m_molecule_read;
		bool m_init_read;
		bool m_cris_read;
		bool m_orientation_read;
		bool m_quaternion_read;	
		bool m_rotation_read;
		bool m_inert_read;
		bool m_asphere_read;
		bool m_patch_read;		
		std::map< std::string, bool > m_read_indicator;
		ifstream::pos_type m_sp;
    };

#endif



