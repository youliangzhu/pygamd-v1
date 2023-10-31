#include "XMLParser.h"
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
#include <functional>
#include <cassert>
#include <iomanip> 
using namespace std;

#ifndef __XML_BULDER_H__
#define __XML_BULDER_H__

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
    unsigned int a;     //!< The tag of the first particle in the dihedral
    unsigned int b;     //!< The tag of the second particle in the dihedral
    unsigned int c;     //!< The tag of the third particle in the dihedral
    unsigned int d;     //!< The tag of the forth particle in the dihedral
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

		vec() : x(0.0), y(0.0), z(0.0)
		{
		}

		vec(double xp, double yp, double zp) : x(xp), y(yp), z(zp)
		{
		}
		double x;   //!< x-component
		double y;   //!< y-component
		double z;   //!< z-component		
	};
	

struct vec4
	{
	vec4() : x(0.0), y(0.0), z(0.0), w(0.0)
		{
		}

	vec4(double xp, double yp, double zp, double wp) : x(xp), y(yp), z(zp), w(wp)
		{
		}
		double x;   //!< x-component
		double y;   //!< y-component
		double z;   //!< z-component
		double w;
	};

struct vec_int
	{

		vec_int() : x(0), y(0), z(0)
		{
		}

		vec_int(int xp, int yp, int zp) : x(xp), y(yp), z(zp)
		{
		}
		int x;   //!< x-component
		int y;   //!< y-component
		int z;   //!< z-component		
	};            

struct vec_uint
	{

		vec_uint() : x(0), y(0), z(0)
		{
		}

		vec_uint(unsigned int xp, unsigned int yp, unsigned int zp) : x(xp), y(yp), z(zp)
		{
		}
		unsigned int x;   //!< x-component
		unsigned int y;   //!< y-component
		unsigned int z;   //!< z-component		
	};
	
struct str_vec6
	{

		str_vec6() : name(""), x(0.0), y(0.0), z(0.0), w(0.0), m(0.0), n(0.0)
		{
		}

		str_vec6(std::string namep, double xp, double yp, double zp, double wp, double mp, double np) : name(namep), x(xp), y(yp), z(zp), w(wp), m(mp), n(np)
		{
		}
		string name;
		double x;   //!< x-component
		double y;   //!< y-component
		double z;   //!< z-component
		double w;   //!< w-component
		double m;   //!< m-component
		double n;   //!< n-component		
	};
	
struct Str2
    {
	Str2() :  x(""), y(""){ }	
    Str2(std::string xp, std::string yp) :  x(xp), y(yp){ }	
    std::string x;  
    std::string y;  	
    };	
	
	
const unsigned int NO_BODY = 0xffffffff;
const unsigned int NO_INDEX = 0xffffffff;


class XMLBuilder 
    {
    public:
        //! Loads in the file and parses the data
        XMLBuilder(const std::string &fname);
		~XMLBuilder();
        //! Returns the number of particles to be initialized
        virtual unsigned int getNDimensions() const;
        
        //! Returns the number of particles to be initialized
        virtual unsigned int getNParticles() const;
        
        //! Returns the number of particle types to be initialized
        virtual unsigned int getNParticleTypes() const;
        
        //! Returns the timestep of the simulation
        virtual unsigned int getTimeStep() const;
        
        //! Sets the timestep of the simulation
        virtual void setTimeStep(unsigned int ts);
        
        //! Returns the box the particles will sit in
        virtual BoxSize getBox() const;

        //! Returns the number of bond types to be created
        virtual unsigned int getNBondTypes() const;
        
        //! Returns the number of angle types to be created
        virtual unsigned int getNAngleTypes() const;
 
        //! Returns the number of dihedral types to be created
        virtual unsigned int getNDihedralTypes() const;
		
        //! Returns the number of vsite types to be created
        virtual unsigned int getNVsiteTypes() const;		
		
        const std::vector< Bond >& getBond() const { return m_bonds; }
        const std::vector< Angle >& getAngle() const { return m_angles; }
        const std::vector< Dihedral >& getDihedral() const { return m_dihedrals; }
        const std::vector< Dihedral >& getVsite() const { return m_vsites; }
        const std::vector< str_vec6 >& getAsphere() const { return m_asphere; } 
        const std::vector< str_vec6 >& getPatch() const { return m_patch; }
        const std::vector< str_vec6 >& getPatchNum() const { return m_patch_num; }			
        //! Access the read particle positions
        const std::vector< vec >& getPos() const { return m_pos_array; }
	
         const std::vector< vec >& getVel() const { return m_vel_array; }       
        //! Access the read images
        const std::vector< vec_int >& getImage() const { return m_image_array; }
        const std::vector< double >& getMass() const { return m_mass_array; }			
        //! Access the read particle positions
        const std::vector< unsigned int >& getType() const { return m_type_array; }
        const std::vector< unsigned int >& getMolecule() const { return m_molecule_array; }
        const std::vector< double >& getCharge() const { return m_charge_array; }	
        const std::vector< unsigned int >& getBody() const { return m_body_array; }		
        const std::vector< double >& getDiameter() const { return m_diameter_array; }	
        const std::vector< std::string >& getTypeMap() const { return m_type_mapping; }	
        const std::vector< std::string >& getBondTypeMap() const { return m_bond_type_mapping; }
        const std::vector< std::string >& getAngleTypeMap() const { return m_angle_type_mapping; }
        const std::vector< std::string >& getDihedralTypeMap() const { return m_dihedral_type_mapping; }
        const std::vector< std::string >& getVsiteTypeMap() const { return m_vsite_type_mapping; }			
        //! Access the read particle positions
        const std::vector< vec >& getOrientation() const { return m_orientation; }
        const std::vector< vec4 >& getQuaternion() const { return m_quaternion; }
        const std::vector< vec >& getInert() const { return m_inert; }
        const std::vector< unsigned int >& getCris() const { return m_cris; }
        const std::vector< unsigned int >& getInit() const { return m_init; }		
        virtual std::string getFilename() { return m_fname; }
        std::string getFirstFilename(){ return m_fname; }		
		virtual bool updatePositionFromXML(const string &fname);
        void outPutInfo();	
		const std::string& getObjectName() const { return m_object_name; }
    protected:
        //! Helper function to read the input file
        void readFile(const std::string &fname);
        //! Helper function to parse the box node
        void parseBoxNode(const XMLNode& node);
        //! Helper function to parse the position node
        void parsePositionNode(const XMLNode& node);
        //! Helper function to parse the image node
        void parseImageNode(const XMLNode& node);
        //! Helper function to parse the velocity node
        void parseVelocityNode(const XMLNode& node);
        //! Helper function to parse the mass node
        void parseMassNode(const XMLNode& node);
        //! Helper function to parse diameter node
        void parseDiameterNode(const XMLNode& node);
        //! Helper function to parse the type node
        void parseTypeNode(const XMLNode& node);
        //! Helper function to parse the body node
        void parseBodyNode(const XMLNode& node);		
        //! Helper function to parse the bonds node
        void parseBondNode(const XMLNode& node);
        //! Helper function to parse the angle node
        void parseAngleNode(const XMLNode& node);
        //! Helper function to parse the dihedral node
        void parseDihedralNode(const XMLNode& node);
        //! Helper function to parse the improper node
        void parseMoleculeNode(const XMLNode& node);
        //! Parse charge node
        void parseChargeNode(const XMLNode& node);
        //! Parse orientation node
        void parseOrientationNode(const XMLNode& node);
        //! Parse Quaternion node
        void parseQuaternionNode(const XMLNode& node);
        void parseInertNode(const XMLNode& node);		
        void parseInitNode(const XMLNode &node);
        void parseCrisNode(const XMLNode &node);		
        void parseConstraintNode(const XMLNode& node);	
		
        void parseVsiteNode(const XMLNode& node);	
        void parseAsphereNode(const XMLNode& node);	
        void parsePatchNode(const XMLNode& node);			
        //! Parse position node for position update		
        void updatePositionNode(const XMLNode &node);
        //! Parse image node for position update		
        void updateImageNode(const XMLNode &node);	
        void updateVelocityNode(const XMLNode& node);
        void updateOrientationNode(const XMLNode& node);
        //! Parse Quaternion node
        void updateQuaternionNode(const XMLNode& node);		
        //! Helper function for identifying the particle type id
        unsigned int getTypeId(const std::string& name);
        //! Helper function for identifying the bond type id
        unsigned int getBondTypeId(const std::string& name);
        //! Helper function for identifying the angle type id
        unsigned int getAngleTypeId(const std::string& name);
        //! Helper function for identifying the dihedral type id
        unsigned int getDihedralTypeId(const std::string& name);
        //! Helper function for identifying the vsite type id
        unsigned int getVsiteTypeId(const std::string& name);		
        
        std::map< std::string, std::function< void (const XMLNode&) > > m_parser_map; //!< Map for dispatching parsers based on node type
        
        BoxSize m_box;   //!< Simulation box read from the file
        bool m_box_read;    //!< Stores the box we read in
        std::string m_fname;
        unsigned int m_num_dimensions;              //!< number of spatial dimensions
        std::vector< vec > m_pos_array;             //!< positions of all particles loaded
        std::vector< vec_int > m_image_array;       //!< images of all particles loaded
        std::vector< vec > m_vel_array;             //!< velocities of all particles loaded
        std::vector< double > m_mass_array;         //!< masses of all particles loaded
        std::vector< double > m_diameter_array;     //!< diameters of all particles loaded
        std::vector< unsigned int > m_type_array;   //!< type values for all particles loaded
        std::vector< unsigned int > m_body_array;   //!< body values for all particles loaded		
        std::vector< double > m_charge_array;       //!< charge of the particles loaded
        std::vector< unsigned int > m_init;       //!< init of the particles loaded
        std::vector< unsigned int > m_cris;       //!< cris of the particles loaded	
        std::vector< vec > m_inert;       //!< cris of the particles loaded			
        std::vector< Bond > m_bonds;                //!< Bonds read in from the file
        std::vector< Angle > m_angles;              //!< Angle read in from the file
        std::vector< Dihedral > m_dihedrals;        //!< Dihedral read in from the file
        std::vector< Dihedral > m_vsites;        //!< Vsite read in from the file		
        unsigned int m_timestep;                    //!< The time stamp

        std::vector< vec > m_orientation;             //!< orientation of all particles loaded	
        std::vector< vec4 > m_quaternion;             //!< quaternion of all particles loaded			
        std::vector< str_vec6 > m_asphere;
        std::vector< str_vec6 > m_patch;
        std::vector< str_vec6 > m_patch_num;		
        std::vector<unsigned int> m_molecule_array;         
        std::vector<std::string> m_type_mapping;          //!< The created mapping between particle types and ids
        std::vector<std::string> m_bond_type_mapping;     //!< The created mapping between bond types and ids
        std::vector<std::string> m_angle_type_mapping;    //!< The created mapping between angle types and ids
        std::vector<std::string> m_dihedral_type_mapping; //!< The created mapping between dihedral types and ids
        std::vector<std::string> m_vsite_type_mapping; //!< The created mapping between vsite types and ids		
		std::string m_node_name;
		std::string m_object_name;
        
    };



#endif



