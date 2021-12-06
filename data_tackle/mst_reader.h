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
#include <memory>
using namespace std;

#ifndef __MST_READER_H__
#define __MST_READER_H__

struct Bond
    {
    Bond(std::string bond_type, unsigned int tag_a, unsigned int tag_b) : type(bond_type), a(tag_a), b(tag_b), bc("b") { }
    Bond(std::string bond_type, unsigned int tag_a, unsigned int tag_b, std::string bond_bc) : type(bond_type), a(tag_a), b(tag_b), bc(bond_bc) { }	
    Bond(std::string bond_type, unsigned int tag_a, unsigned int tag_b, unsigned int bond_id) : type(bond_type), a(tag_a), b(tag_b), id(bond_id), bc("b") { }
    Bond(std::string bond_type, unsigned int tag_a, unsigned int tag_b, unsigned int bond_id, std::string bond_bc) : type(bond_type), a(tag_a), b(tag_b), id(bond_id), bc(bond_bc) { }	
    std::string type;  
    unsigned int a;    
    unsigned int b;
	unsigned int id;
	std::string bc;
    };
	
struct Angle
    {
    Angle(std::string angle_type, unsigned int tag_a, unsigned int tag_b, unsigned int tag_c) : type(angle_type), a(tag_a), b(tag_b), c(tag_c) { }
    Angle(std::string angle_type, unsigned int tag_a, unsigned int tag_b, unsigned int tag_c, unsigned int angle_id) : type(angle_type), a(tag_a), b(tag_b), c(tag_c), id(angle_id) { }	
    std::string type;  
    unsigned int a;    
    unsigned int b;
    unsigned int c; 
	unsigned int id;    
    };		
	
struct Dihedral
    {
    Dihedral(std::string dihedral_type, unsigned int tag_a, unsigned int tag_b, unsigned int tag_c, unsigned int tag_d) : type(dihedral_type), a(tag_a), b(tag_b), c(tag_c), d(tag_d) { }
    Dihedral(std::string dihedral_type, unsigned int tag_a, unsigned int tag_b, unsigned int tag_c, unsigned int tag_d, unsigned int dihedral_id) : type(dihedral_type), a(tag_a), b(tag_b), c(tag_c), d(tag_d),  id(dihedral_id) { }	
    std::string type;  
    unsigned int a;     
    unsigned int b;     
    unsigned int c;     
    unsigned int d;     
	unsigned int id;	
    };

struct BoxSize
    {
    BoxSize() : lx(0.0), ly(0.0), lz(0.0) { }		
    BoxSize(double len) : lx(len), ly(len), lz(len) { }
    BoxSize(double len_x, double len_y, double len_z) : lx(len_x), ly(len_y), lz(len_z) { } 		
		
    double lx;
    double ly;
    double lz;
    };
	
struct vec
	{
	vec() : x(0.0), y(0.0), z(0.0){ }
	vec(double xp, double yp, double zp) : x(xp), y(yp), z(zp){ }
	double x;  
	double y;  
	double z; 		
	};
	

struct vec4
	{
	vec4() : x(0.0), y(0.0), z(0.0), w(0.0){ }
	vec4(double xp, double yp, double zp, double wp) : x(xp), y(yp), z(zp), w(wp){ }
	double x;  
	double y;  
	double z;  
	double w;
	};

struct vec_int
	{
	vec_int() : x(0), y(0), z(0){ }
	vec_int(int xp, int yp, int zp) : x(xp), y(yp), z(zp){ }
	int x; 
	int y; 
	int z; 
	};            

struct vec_uint
	{
	vec_uint() : x(0), y(0), z(0){ }
	vec_uint(unsigned int xp, unsigned int yp, unsigned int zp) : x(xp), y(yp), z(zp){ }
	unsigned int x;   
	unsigned int y;   
	unsigned int z;   	
	};
	
struct str_vec6
	{
	str_vec6() : name(""), x(0.0), y(0.0), z(0.0), w(0.0), m(0.0), n(0.0){ }
	str_vec6(std::string namep, double xp, double yp, double zp, double wp, double mp, double np) : name(namep), x(xp), y(yp), z(zp), w(wp), m(mp), n(np){ }
	string name;
	double x;  
	double y;  
	double z;  
	double w;  
	double m;  
	double n;  	
	};
	
struct Str2
    {
	Str2() :  x(""), y(""){ }	
    Str2(std::string xp, std::string yp) :  x(xp), y(yp){ }	
    std::string x;  
    std::string y;  	
    };	
	

struct uint_2
	{
	//! Default construtor
	uint_2() : x(0), y(0)
		{
		}
	uint_2(unsigned int xp, unsigned int yp) : x(xp), y(yp)
		{
		}
	unsigned int x;  //!< x-component
	unsigned int y;  //!< y-component
	};
	
	
	
const unsigned int NO_BODY = 0xffffffff;
const unsigned int NO_INDEX = 0xffffffff;


class mst_reader 
    {
    public:
        mst_reader();
        unsigned int getNDimensions() const;
        unsigned int getNParticles() const;
        unsigned int getLastNParticles() const;		
        unsigned int getNParticleTypes() const;
        unsigned int getTimeStep() const;
        BoxSize getBox() const;
        unsigned int getNBondTypes() const;
        unsigned int getNAngleTypes() const;
        unsigned int getNDihedralTypes() const;
        unsigned int getNVsiteTypes() const;		
		
        const std::vector< Bond >& getBond() const { return m_bonds; }
        const std::vector< Angle >& getAngle() const { return m_angles; }
        const std::vector< Dihedral >& getDihedral() const { return m_dihedrals; }
        const std::vector< Dihedral >& getVsite() const { return m_vsites; }
        const std::vector< str_vec6 >& getAsphere() const { return m_asphere; } 
        const std::vector< str_vec6 >& getPatch() const { return m_patch; }
        const std::vector< str_vec6 >& getPatchNum() const { return m_patch_num; }			
        const std::vector< vec >& getPos() const { return m_pos; }
	
        const std::vector< vec >& getVel() const { return m_vel; }       
        const std::vector< vec_int >& getImage() const { return m_image; }
        const std::vector< double >& getMass() const { return m_mass; }			
        const std::vector< unsigned int >& getType() const { return m_type; }
        const std::vector< unsigned int >& getMolecule() const { return m_molecule; }
        const std::vector< double >& getCharge() const { return m_charge; }	
        const std::vector< unsigned int >& getBody() const { return m_body; }		
        const std::vector< double >& getDiameter() const { return m_diameter; }	
        const std::vector< std::string >& getTypeMap() const { return m_type_map; }	
        const std::vector< std::string >& getBondTypeMap() const { return m_bond_type_map; }
        const std::vector< std::string >& getAngleTypeMap() const { return m_angle_type_map; }
        const std::vector< std::string >& getDihedralTypeMap() const { return m_dihedral_type_map; }
        const std::vector< std::string >& getVsiteTypeMap() const { return m_vsite_type_map; }			
        const std::vector< vec >& getOrientation() const { return m_orientation; }
        const std::vector< vec4 >& getQuaternion() const { return m_quaternion; }
        const std::vector< vec >& getInert() const { return m_inert; }
        const std::vector< vec >& getRotation() const { return m_rotation; }		
        const std::vector< unsigned int >& getCris() const { return m_cris; }
        const std::vector< unsigned int >& getInit() const { return m_init; }
		const std::vector<vec>& getRotangle() const {return m_rotangle;}		
        virtual std::string getFilename() { return m_fname; }
        std::string getFirstFilename(){ return m_fname; }		
		virtual bool readDataFromMST(const string &fname);
		virtual bool ifchangedNp() {return m_if_changed_np; }		
        void outPutInfo();	
		const std::string& getObjectName() const { return m_object_name; }	
		bool iftrajectory() { return m_if_trajectory; }		
    protected:
        void reset_params();
		void clear_data();
        unsigned int getTypeId(const std::string& name);
        unsigned int getBondTypeId(const std::string& name);
        unsigned int getAngleTypeId(const std::string& name);
        unsigned int getDihedralTypeId(const std::string& name);
        unsigned int getVsiteTypeId(const std::string& name);		
        
        BoxSize m_box;
        std::string m_fname;
		
		unsigned int m_num_particles;
		unsigned int m_timestep;
		unsigned int m_dimension;			
          
        std::vector< vec > m_pos;           
        std::vector< vec_int > m_image;     
        std::vector< vec > m_vel;           
        std::vector< double > m_mass;        
        std::vector< unsigned int > m_type; 
        std::vector< unsigned int > m_body; 
        std::vector< double > m_charge;
        std::vector< double > m_diameter;		
        std::vector< unsigned int > m_init;       
        std::vector< unsigned int > m_cris; 
        std::vector< vec > m_force;		
        std::vector< double > m_virial;
              
        std::vector< Bond > m_bonds;              
        std::vector< Angle > m_angles;            
        std::vector< Dihedral > m_dihedrals;      
        std::vector< Dihedral > m_vsites;                         

        std::vector< vec > m_orientation;         
        std::vector< vec4 > m_quaternion;
        std::vector< vec > m_inert; 		
        std::vector<vec> m_rotangle;
        std::vector<vec> m_rotation; 		
        std::vector< str_vec6 > m_asphere;
        std::vector< str_vec6 > m_patch;
        std::vector< str_vec6 > m_patch_num;		
        std::vector<unsigned int> m_molecule;         
        std::vector<std::string> m_type_map;          
        std::vector<std::string> m_bond_type_map;     
        std::vector<std::string> m_angle_type_map;    
        std::vector<std::string> m_dihedral_type_map; 
        std::vector<std::string> m_vsite_type_map;
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
		std::string m_object_name;
		bool m_if_changed_np;
		bool m_if_trajectory;
		unsigned int m_last_np;
    };

#endif



