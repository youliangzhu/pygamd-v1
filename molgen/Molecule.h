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
#include "mst_reader.h"
#include<pybind11/pybind11.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>

#include <iomanip> 
using namespace std;

#ifndef __MOLECULE_H__
#define __MOLECULE_H__

struct Cylinder
    {
    Cylinder(double ox=0.0, double oy=0.0, double oz=0.0, double nx=1.0, double ny=0.0, double nz=0.0, double r_min=0.0, double r_max=0.0)
            : origin_x(ox), origin_y(oy), origin_z(oz), radius_min(r_min), radius_max(r_max)
        {
        // Directionize nx, ny, nz
        double len = sqrt(nx*nx + ny*ny + nz*nz);
        direction_x = nx / len;
        direction_y = ny / len;
        direction_z = nz / len;
        }
        
    double origin_x;   
    double origin_y;   
    double origin_z;   
    
    double direction_x;   
    double direction_y;   
    double direction_z; 
	double radius_min;
	double radius_max;
    };	
	
struct Sphere
    {
    Sphere(double ox=0.0, double oy=0.0, double oz=0.0, double r_min=0.0, double r_max=0.0)
            : origin_x(ox), origin_y(oy), origin_z(oz), radius_min(r_min), radius_max(r_max)
        {
        }
        
    double origin_x;   
    double origin_y;   
    double origin_z;
	double radius_min;
	double radius_max;
    };		
		
	
void exyzFromQuaternion(vec4 &quat, vec4 &ex_space, vec4 &ey_space, vec4 &ez_space);

void quaternionFromEXYZ(vec4 &quat, vec4 &ex_space, vec4 &ey_space, vec4 &ez_space);


class Molecule
	{
	public:
		//! Constructs the compute
		Molecule(unsigned int natomPerMole);
		//! Constructs the compute
		Molecule(const std::string& fname, unsigned int NatomPerMole);		
		//! Destructor
		 virtual ~Molecule(){}

		void setParticleTypes(std::string types);
		void setTopology(std::string topology);
		void setPosition(unsigned int i, double px,double py, double pz);		
		void setIsotactic(bool iso);
		void setBondLength(double length);		
		void setBondLength(std::string name1, std::string name2, double length);
		
		void setAngleDegree(std::string name_a, std::string name_b, std::string name_c, double degree);
		void setDihedralDegree(std::string name, std::string name_b, std::string name_c, std::string name_d, double degree);	
		void setAngleDegree(unsigned int i,unsigned int j,unsigned int k, double degree);
		void setDihedralDegree(unsigned int i,unsigned int j,unsigned int k,unsigned int l, double degree);
		void setMass(double mass);
		void setMass(unsigned int i, double mass);		
		void setMass(std::string type, double mass);
		void setCharge(double charge);
		void setCharge(std::string type, double charge);
		void setCharge(unsigned int i, double charge);
		void setChargeFactor(double factor);		
		void setOrientation();
		void setOrientation(std::string type);
		void setOrientation(unsigned int i);
		void setInert(double inertx, double inerty, double inertz);
		void setInert(std::string type, double inertx, double inerty, double inertz);
		void setInert(unsigned int i, double inertx, double inerty, double inertz);		
		void setQuaternion();
		void setQuaternion(std::string type);
		void setQuaternion(unsigned int i);		
		void setDiameter(double diameter);
		void setDiameter(std::string type, double diameter);
		void setDiameter(unsigned int i, double diameter);
		void setCris(unsigned int cris);
		void setCris(std::string type, unsigned int cris);
		void setCris(unsigned int i, unsigned int cris);
		void setInit(unsigned int init);
		void setInit(std::string type, unsigned int init);
		void setInit(unsigned int i, unsigned int init);
		void setBody(unsigned int body);
		void setBody(std::string type, unsigned int body);
		void setBody(unsigned int i, unsigned int body);
		void setMolecule(unsigned int molecule);
		void setMolecule(std::string type, unsigned int molecule);
		void setMolecule(unsigned int i, unsigned int molecule);
		void setEllipsoidBondSpotVector(double spvx, double spvy, double spvz);		
		void setBox(double mol_Lx,double mol_Ly, double mol_Lz);
		void setBox(double mol_Lx_min, double mol_Lx_max, double mol_Ly_min, double mol_Ly_max, double mol_Lz_min, double mol_Lz_max);
		void setSphere(double ox,double oy, double oz, double r_min, double r_max);
		void setCylinder(double ox,double oy, double oz, double dx,double dy, double dz, double r_min, double r_max);
		void setBodyEvacuation();
		void setPutBox(double Lx,double Ly, double Lz);
		void generateAngle();
		void generateDihedral();		
		bool twoAnglesFixD(vec A, vec B, vec C, vec& D1, vec& D2, double lenthBD, double thetaABD,double thetaCBD);
		bool twoAnglesFixD(vec A, vec B, vec C, vec& D, double lenthBD, double thetaABD,double thetaCBD);
		bool twoAnglesFixE(vec A, vec B, vec C, vec D, vec& E1, vec& E2, double lenthBE, double lenthCE, double thetaABE, double thetaDCE);
		bool twoAnglesFixE(vec A, vec B, vec C, vec D, vec& E, double lenthBE, double lenthCE, double thetaABE, double thetaDCE);		
		bool threeAnglesFixE(vec A, vec B, vec C, vec D, vec& E, double lenthBE, double thetaABE,double thetaCBE,double thetaDBE);
		bool threeAnglesFixG(vec A, vec B, vec C, vec D, vec E, vec F, vec& G, 
								double lenthBG, double lenthCG,double lenthEG, 
								double thetaABG,double thetaDCG,double thetaFEG);
		bool arrayFixF(vec a, vec b, vec c, vec d, vec e, vec& F1, vec& F2);	
		void setDimention(unsigned int dimention);		
        std::vector<vec>& getPosition()
			{
			return m_xyz;
			}
        std::vector<std::string>& getType()
			{
			initData();
			return m_type;
			}		
        std::vector<Bond>& getBond()
			{
			initData();	
			return m_bond;
			}			
        std::vector<Angle>& getAngle()
			{
			return m_angle;
			}
        std::vector<Dihedral>& getDihedral()
			{
			return m_dihedral;
			}
        std::vector<Dihedral>& getVsite()
			{
			return m_vsite;
			}				
        std::vector<double>& getMass()
			{
			return m_mass;
			}
        std::vector<vec>& getInert()
			{
			return m_inert;
			}			
        std::vector<unsigned int >& getOrientation()
			{
			return m_orientation;
			}
        std::vector<unsigned int >& getQuaternion()
			{
			return m_quaternion;
			}				
        std::vector<double>& getCharge()
			{
			return m_charge;
			}
        std::vector<double>& getDiameter()
			{
			return m_diameter;
			}
        std::vector<unsigned int>& getBody()
			{
			return m_body;
			}
        std::vector<unsigned int>& getCris()
			{
			return m_cris;
			}
        std::vector<unsigned int>& getInit()
			{
			return m_init;
			}
        std::vector<unsigned int>& getMolecule()
			{
			return m_molecule;
			}
        std::vector<vec>& getOriVec()
			{
			return m_ori_vec;
			}
        std::vector<vec4>& getQuatVec()
			{
			return m_quat_vec;
			}
        std::vector<str_vec6>& getAsphere()
			{
			return m_asphere;
			}
        std::vector<str_vec6>& getPatch()
			{
			return m_patch;
			}
        std::vector<str_vec6>& getPatchNum()
			{
			return m_patch_num;
			}				
        unsigned int getBodyIdPlus()
			{
			return m_body_id_plus;
			}
        unsigned int getMolIdPlus()
			{
			return m_mol_id_plus;
			}
        void updatePos(std::vector<vec>& pos)
			{
			for(unsigned int i=0; i< pos.size(); i++)	
				m_pos_all.push_back(pos[i]);
			}
        void updateType(std::vector<std::string>& type)
			{
			for(unsigned int i=0; i< type.size(); i++)
				m_type_all.push_back(switchNametoType(type[i]));	
			}
        void updateBodyCom(std::vector<vec4>& body_com)
			{
			for(unsigned int i=0; i< body_com.size(); i++)	
				m_body_com_all.push_back(body_com[i]);
			}			
        unsigned int getNumParticle()
			{
			return m_NatomPerMole;
			}
        void setList(unsigned int* list, unsigned int* head)
			{
			m_list = list;
			m_head = head;
			}
        void setCell(vec_uint& dim,vec& width)
			{
			m_dim = dim;
			m_width = width;
			}	
        void setParam(std::vector<vec>& params, unsigned int NBtype, std::vector<string>& type_mapping_all, std::vector<double>& min_dis)
			{
			m_params = params;
			m_NBtype = NBtype;
			m_type_mapping_all = type_mapping_all;
			m_min_dis = min_dis;
			if(!m_set_testnum&&m_NatomPerMole>1)
				m_testnum=16;
			for (unsigned int i=0; i< m_min_dis.size(); i++)
				{
				if(m_min_dis[i]>0)
					m_check_distance = true;
				}
			}
        void setTestNum(unsigned int testnum)
			{
			m_testnum = testnum;
			m_set_testnum=true;
			}
        void genName();
		unsigned int cellid(int i, int j, int k);	
		void readData(const std::string& fname);
		bool existedDihedral(unsigned int a, unsigned int b, unsigned int c, unsigned int d);
		double existedAngleDigree(unsigned int a, unsigned int b, unsigned int c);	
		double isotactic(vec A, vec B, vec C, vec D);
		bool interMolCheck(unsigned int tag1, vec posa, double& boltz);
		bool intraMolCheck(unsigned int tag1, unsigned int tag2, vector<unsigned int>& tags, vec posa, double& boltz);
		bool checkdistance(unsigned int tag1, unsigned int tag2, vector<unsigned int>& tags, vector<vec>& testpos, vec& posa, unsigned int onoff);	
		unsigned int switchNametoType(const string& name);
		void initBond();
		void initType();
		void initData();
		void allocateData(unsigned int size);
		virtual void generate();

	protected:
        unsigned int m_NatomPerMole;
		unsigned int m_Ntypes;
		unsigned int m_dimention;
		vector<double> BondLength;
		vector<double> AngleRadian;
		vector<double> DihedralRadian;

		vector<double> RadianPerAngle;
		vector<double> RadianPerDihedral;		
		
		vector<std::string> m_type;
		vector<unsigned int> m_typeId;		
		vector<Bond> m_bond;					
		vector<Angle> m_angle;
		vector<Dihedral> m_dihedral;
		vector<Dihedral> m_vsite;		
		vector<double> m_mass;
		vector<unsigned int > m_orientation;
		vector<unsigned int > m_quaternion;		
		vector<double> m_charge;
		vector<double> m_diameter;
		vector<vec> m_inert;		
		vector<unsigned int> m_cris;
		vector<unsigned int> m_init;
		vector<unsigned int> m_body;
		vector<unsigned int> m_molecule;
		vector<vec> m_xyz;						// x , y , z
		vector<vec> m_xyz_read;						// x , y , z
		vector<vec> m_ori_vec;						// x , y , z
		vector<vec4> m_quat_vec;
		vector<vec> m_ori_vec_read;						// x , y , z
		vector<vec4> m_quat_vec_read;
		vector<str_vec6> m_asphere;
		vector<str_vec6> m_patch;
		vector<str_vec6> m_patch_num;			
		unsigned int m_body_id_plus;
		unsigned int m_mol_id_plus;

		unsigned int m_testnum;
		bool m_set_testnum;
        std::vector<vec> m_pos_all;	
        std::vector<unsigned int> m_type_all;
        std::vector< string > m_type_mapping_all;
        std::vector<vec4> m_body_com_all;		
		unsigned int m_NBtype;
        std::vector<vec> m_params;
        std::vector<double> m_min_dis;
		
		vector<unsigned int> m_nbond;
		vector<unsigned int> m_bondtable;
		unsigned int m_bondtableHight;
		vector<bool> m_be_generated;
		vector<bool> m_be_generated_read;
		vector<unsigned int> m_bond_init;
		Sphere m_sphere;
		Cylinder m_cylinder;
		bool m_set_sphere;
		bool m_set_cylinder;
		bool m_set_body_evacuation;
		double m_Lx;
		double m_Ly;
		double m_Lz;
		double m_mol_Lx;
		double m_mol_Ly;
		double m_mol_Lz;
		double m_shift_Lx;
		double m_shift_Ly;
		double m_shift_Lz;
		bool m_set_mol_box;
		unsigned int dimention;
		bool m_firststep;
		bool m_isotactic;
		bool m_initdata;
		bool m_check_distance;
		vec_uint m_dim;
		vec m_width;		
		unsigned int* m_list;
		unsigned int* m_head;
		
		unsigned int m_limit_ge;
		int str2num(string s);
		double R2S();
        std::vector<std::string> m_type_mapping;		
		unsigned int getTypeId(const std::string& name);
		std::string m_mol_name;
		unsigned int m_output_times;
		unsigned int m_Nread_particle;
		std::string m_type_str;
		std::string m_topology_str;
		vec m_eb_spv;
		unsigned int nwarning_delt;
		};
		

void export_Molecule(pybind11::module &m);
#endif
