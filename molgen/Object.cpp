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

#include "Object.h"
#include<time.h> 
#include<stdlib.h> 
using namespace std;

Object::Object(unsigned int NatomPerMole, Object::Shape s): Molecule(NatomPerMole), m_s(s)
    {
	m_radius = 1.0;
	m_xyz_temp.resize(NatomPerMole);
    }

Object::Object(const std::string& fname, unsigned int NatomPerMole, Object::Shape s): Molecule(fname, NatomPerMole), m_s(s)
    {
	m_radius = 1.0;
	m_xyz_temp.resize(NatomPerMole);
    }

void Object::removeCM()
	{
	if(m_Nread_particle!=m_NatomPerMole)
		{
		cerr << endl << "***Error! The number of read particles " << m_Nread_particle <<" is different from initialized set value "<< m_NatomPerMole<<" !"<< endl << endl;
		throw runtime_error("Object::removeCM error");	
		}

	vec center_mass;
	center_mass.x=0.0;
	center_mass.y=0.0;
	center_mass.z=0.0;
	for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
		{
		center_mass.x += m_xyz_temp[taga].x;
		center_mass.y += m_xyz_temp[taga].y; 
		center_mass.z += m_xyz_temp[taga].z;
		}
	center_mass.x /= double(m_NatomPerMole);
	center_mass.y /= double(m_NatomPerMole);
	center_mass.z /= double(m_NatomPerMole);
	for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
		{
		m_xyz_temp[taga].x -= center_mass.x;
		m_xyz_temp[taga].y -= center_mass.y; 
		m_xyz_temp[taga].z -= center_mass.z;
		}
	}

void Object::generateSphere()
	{
	for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
		{
		if(!m_be_generated[taga])
			{
			bool fitpoint = false;
			bool success = true;
			unsigned int turns =0;
			vec posa;
			while(!fitpoint)
				{
				vector<vec> testpos;
				for(unsigned int i=0; i<m_testnum; i++)
					{
					double theta = R2S()*M_PI;
					double ran = R2S();
					double phi = 2.0*asin(sqrt(ran));
					ran=R2S()-0.5;
					if(ran>0)
						phi = -phi;
					posa.x = sin(phi)*cos(theta)*m_radius;
					posa.y = sin(phi)*sin(theta)*m_radius;
					posa.z = cos(phi)*m_radius;								
					testpos.push_back(posa);
					}
				vector<unsigned int> vd;							
				fitpoint = checkdistance(taga, taga, vd, testpos, posa, 1);
				if(turns>10000)
					{
					success = false;
					fitpoint = true;
					}
				turns +=1;								
				}
			if(!success)
				{
				cerr << endl << "***Error! Can not generate this object!"<<endl << endl;
				throw runtime_error("Object::generate error");	
				}
			m_xyz[taga].x = posa.x;
			m_xyz[taga].y = posa.y; 
			m_xyz[taga].z = posa.z;
			
			m_xyz_temp[taga].x = posa.x;
			m_xyz_temp[taga].y = posa.y; 
			m_xyz_temp[taga].z = posa.z;
			
			m_be_generated[taga]= true;
			}
		}
	vec center_mass;
	center_mass.x=0.0;
	center_mass.y=0.0;
	center_mass.z=0.0;
	for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
		{
		center_mass.x += m_xyz_temp[taga].x;
		center_mass.y += m_xyz_temp[taga].y; 
		center_mass.z += m_xyz_temp[taga].z;
		}
	center_mass.x /= double(m_NatomPerMole);
	center_mass.y /= double(m_NatomPerMole);
	center_mass.z /= double(m_NatomPerMole);
	for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
		{
		m_xyz_temp[taga].x -= center_mass.x;
		m_xyz_temp[taga].y -= center_mass.y; 
		m_xyz_temp[taga].z -= center_mass.z;
		}
	}

void Object::placeObject()
	{
	bool fitpoint = false;
	bool success = true;
	unsigned int turns =0;
	vec posa, centermass, angle;
	vec4 q, ar0, ar1, ar2;
	while(!fitpoint)
		{
		centermass.x = ( R2S() - 0.5 )*m_mol_Lx + m_shift_Lx;
		centermass.y = ( R2S() - 0.5 )*m_mol_Ly + m_shift_Ly;
		centermass.z = ( R2S() - 0.5 )*m_mol_Lz + m_shift_Lz;
		angle.x = R2S()*M_PI;
		angle.y = R2S()*M_PI;
		angle.z = R2S()*M_PI;
		if(m_mol_Lz<1.0e-6)
			angle.x = 0.0;
		q.x = cos(angle.x/2.0)*cos((angle.y+angle.z)/2.0);
		q.y = sin(angle.x/2.0)*cos((angle.y-angle.z)/2.0);
		q.z = sin(angle.x/2.0)*sin((angle.y-angle.z)/2.0);
		q.w = cos(angle.x/2.0)*sin((angle.y+angle.z)/2.0);
		double rsq = q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w;
		double r = sqrt(rsq);
		q.x /= r;
		q.y /= r;
		q.z /= r;
		q.w /= r;
		exyzFromQuaternion(q, ar0, ar1, ar2);
		for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
			{
			vec particle_pos = m_xyz_temp[taga];
			vec ri, posa;
			ri.x = ar0.x * particle_pos.x + ar1.x * particle_pos.y + ar2.x * particle_pos.z;
			ri.y = ar0.y * particle_pos.x + ar1.y * particle_pos.y + ar2.y * particle_pos.z;
			ri.z = ar0.z * particle_pos.x + ar1.z * particle_pos.y + ar2.z * particle_pos.z;
			posa.x = ri.x + centermass.x;
			posa.y = ri.y + centermass.y; 
			posa.z = ri.z + centermass.z;
			vector<vec> testpos;
			for(unsigned int i=0; i<m_testnum; i++)
				{	
				testpos.push_back(posa);
				}
			vector<unsigned int> vd;							
			fitpoint = checkdistance(taga, taga, vd, testpos, posa, 0);	
			if(!fitpoint)
				break;
			}
		if(turns>10000)
			{
			success = false;
			fitpoint = true;
			}
		turns +=1;
		}
	if(!success)
		{
		cerr << endl << "***Error! Can not generate this object!"<<endl << endl;
		throw runtime_error("Object::generate error");	
		}
	for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
		{
		vec particle_pos = m_xyz_temp[taga];
		vec ri, posa;
		ri.x = ar0.x * particle_pos.x + ar1.x * particle_pos.y + ar2.x * particle_pos.z;
		ri.y = ar0.y * particle_pos.x + ar1.y * particle_pos.y + ar2.y * particle_pos.z;
		ri.z = ar0.z * particle_pos.x + ar1.z * particle_pos.y + ar2.z * particle_pos.z;
		posa.x = ri.x + centermass.x;
		posa.y = ri.y + centermass.y; 
		posa.z = ri.z + centermass.z;
		m_xyz[taga] = posa;
		
		if (m_orientation[taga]==2)
			{
			vec ori = m_ori_vec_read[taga];

			ori.x = ar0.x * ori.x + ar1.x * ori.y + ar2.x * ori.z;
			ori.y = ar0.y * ori.x + ar1.y * ori.y + ar2.y * ori.z;
			ori.z = ar0.z * ori.x + ar1.z * ori.y + ar2.z * ori.z;

			m_ori_vec[taga] = ori;
			}		
		
		if (m_quaternion[taga]==2)
			{
			vec4 ars0, ars1, ars2, ars_new0, ars_new1, ars_new2;
			vec4 quat = m_quat_vec_read[taga];
			if(quat.x!=0.0||quat.y!=0.0||quat.z!=0.0||quat.w!=0.0)
				{
				exyzFromQuaternion(quat, ars0, ars1, ars2);
	//		cout<<"initial "<<quat.x<<" "<<quat.y<<" "<<quat.z<<" "<<quat.w<<endl;			
				ars_new0.x = ar0.x * ars0.x + ar1.x * ars0.y + ar2.x * ars0.z;
				ars_new0.y = ar0.y * ars0.x + ar1.y * ars0.y + ar2.y * ars0.z;
				ars_new0.z = ar0.z * ars0.x + ar1.z * ars0.y + ar2.z * ars0.z;
	
				ars_new1.x = ar0.x * ars1.x + ar1.x * ars1.y + ar2.x * ars1.z;
				ars_new1.y = ar0.y * ars1.x + ar1.y * ars1.y + ar2.y * ars1.z;
				ars_new1.z = ar0.z * ars1.x + ar1.z * ars1.y + ar2.z * ars1.z;	
	
				ars_new2.x = ar0.x * ars2.x + ar1.x * ars2.y + ar2.x * ars2.z;
				ars_new2.y = ar0.y * ars2.x + ar1.y * ars2.y + ar2.y * ars2.z;
				ars_new2.z = ar0.z * ars2.x + ar1.z * ars2.y + ar2.z * ars2.z; 
				vec4 quat2;
				quaternionFromEXYZ(quat2, ars_new0, ars_new1, ars_new2);
	//			cout<<"after "<<quat2.x<<" "<<quat2.y<<" "<<quat2.z<<" "<<quat2.w<<endl;
				m_quat_vec[taga] = quat2;				
				}
			else
				m_quat_vec[taga] = quat;
			}
		}
	}

void Object::generate()
	{
	m_xyz.clear();
	m_xyz.resize(m_NatomPerMole);

	if (m_firststep)
		{
		initData();
		genName();
		cout<<"Molecule: "<<m_mol_name<<endl;
		cout<<"-- statistics --"<<endl;
		cout<<"The number of particles: "<<m_NatomPerMole<<endl;
		cout<<"The number of types: "<<m_Ntypes<<endl;
		for(unsigned int i=0; i<m_Ntypes; i++)
			cout<<m_type_mapping[i]<<endl;	
		cout<<"The number of bonds in a molecule: "<<m_bond.size()<<endl;
		generateAngle();
		generateDihedral();	
		cout<<"generating ..."<<endl;		
		if(!m_set_mol_box)
			{
			m_mol_Lx=m_Lx;
			m_mol_Ly=m_Ly;
			m_mol_Lz=m_Lz;
			}
		for(unsigned int i=0; i<m_NatomPerMole; i++)
			{
			if(m_be_generated_read[i])
				{
				m_be_generated[i]= true;
				m_xyz[i]=m_xyz_read[i];
				m_xyz_temp[i]=m_xyz_read[i];
				}
			else
				{
				m_be_generated[i]= false;
				}
			}
		if (m_s==none)
			{
			removeCM();
			}
		else if(m_s==sphere)
			{
			generateSphere();
			}
		m_firststep = false;
		}
	placeObject();
	}

void export_Object(pybind11::module &m)
	{
	pybind11::class_<Object, Molecule, std::shared_ptr<Object> >(m, "Object")
		.def(pybind11::init<unsigned int, Object::Shape >())
		.def(pybind11::init<const std::string&, unsigned int, Object::Shape >())
		.def("setRadius", &Object::setRadius)	
		;
    pybind11::enum_<Object::Shape>(m, "Shape")
    .value("none",Object::Shape::none)
    .value("sphere",Object::Shape::sphere)
	.export_values()
	;
	} 

