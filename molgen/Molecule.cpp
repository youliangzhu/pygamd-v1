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
#include<time.h> 
#include<stdlib.h> 
using namespace std;

std::string etrim(std::string s)
	{
	unsigned int b=0;
	unsigned int e=0;
	for(unsigned int i=0;i<s.size();i++)
	{
		if(s[i]=='<')
		b = i;
		else if(s[i]=='>')
		e = i;
	}
	if(e>b)
		s=s.substr(b,e-b+1);
	return s;
	}

void Normalize(vec& v)
	{
	double r= sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
	v.x/=r;
	v.y/=r;
	v.z/=r;
	}

void RotateArbitraryLine(double m[4][4], vec& v1, vec& v2, double theta)
	{
    double a = v1.x;
    double b = v1.y;
    double c = v1.z;
	vec p;
	
    p.x = v2.x - v1.x;
	p.y = v2.y - v1.y;
    p.z = v2.z - v1.z;
	
    Normalize(p);
    double u = p.x;
    double v = p.y;
    double w = p.z;

    double uu = u * u;
    double uv = u * v;
    double uw = u * w;
    double vv = v * v;
    double vw = v * w;
    double ww = w * w;
    double au = a * u;
    double av = a * v;
    double aw = a * w;
    double bu = b * u;
    double bv = b * v;
    double bw = b * w;
    double cu = c * u;
    double cv = c * v;
    double cw = c * w;

    double costheta = cos(theta);
    double sintheta = sin(theta);

    m[0][0] = uu + (vv + ww) * costheta;
    m[0][1] = uv * (1 - costheta) + w * sintheta;
    m[0][2] = uw * (1 - costheta) - v * sintheta;
    m[0][3] = 0;

    m[1][0] = uv * (1 - costheta) - w * sintheta;
    m[1][1] = vv + (uu + ww) * costheta;
    m[1][2] = vw * (1 - costheta) + u * sintheta;
    m[1][3] = 0;

    m[2][0] = uw * (1 - costheta) + v * sintheta;
    m[2][1] = vw * (1 - costheta) - u * sintheta;
    m[2][2] = ww + (uu + vv) * costheta;
    m[2][3] = 0;

    m[3][0] = (a * (vv + ww) - u * (bv + cw)) * (1 - costheta) + (bw - cv) * sintheta;
    m[3][1] = (b * (uu + ww) - v * (au + cw)) * (1 - costheta) + (cu - aw) * sintheta;
    m[3][2] = (c * (uu + vv) - w * (au + bv)) * (1 - costheta) + (av - bu) * sintheta;
    m[3][3] = 1;
	}

void exyzFromQuaternion(vec4 &quat, vec4 &ex_space, vec4 &ey_space, vec4 &ez_space)
    {
    // ex_space
    ex_space.x = quat.x * quat.x + quat.y * quat.y - quat.z * quat.z - quat.w * quat.w;
    ex_space.y = double(2.0) * (quat.y * quat.z + quat.x * quat.w);
    ex_space.z = double(2.0) * (quat.y * quat.w - quat.x * quat.z);
    
    // ey_space
    ey_space.x = double(2.0) * (quat.y * quat.z - quat.x * quat.w);
    ey_space.y = quat.x * quat.x - quat.y * quat.y + quat.z * quat.z - quat.w * quat.w;
    ey_space.z = double(2.0) * (quat.z * quat.w + quat.x * quat.y);
    
    // ez_space
    ez_space.x = double(2.0) * (quat.y * quat.w + quat.x * quat.z);
    ez_space.y = double(2.0) * (quat.z * quat.w - quat.x * quat.y);
    ez_space.z = quat.x * quat.x - quat.y * quat.y - quat.z * quat.z + quat.w * quat.w;
    }
	
void quaternionFromEXYZ(vec4 &quat, vec4 &ex_space, vec4 &ey_space, vec4 &ez_space) 
	{
/* 	double qx = 0.5*sqrt(1.0 + ex_space.x + ey_space.y + ez_space.z);
	double qy = 0.5*sqrt(1.0 + ex_space.x - ey_space.y - ez_space.z);		
	double qz = 0.5*sqrt(1.0 - ex_space.x + ey_space.y - ez_space.z);		
	double qw = 0.5*sqrt(1.0 - ex_space.x - ey_space.y + ez_space.z); */
 	double trace = ex_space.x + ey_space.y + ez_space.z; // I removed + 1.0f; see discussion with Ethan
	if( trace > 0 ) 
		{// I changed M_EPSILON to 0
		double s = double(0.5) / sqrt(trace + double(1.0));
		quat.x = double(0.25) / s;
		quat.y = ( ey_space.z - ez_space.y ) * s;
		quat.z = ( ez_space.x - ex_space.z ) * s;
		quat.w = ( ex_space.y - ey_space.x ) * s;
		} 
	else{
		if ( ex_space.x > ey_space.y && ex_space.x > ez_space.z ) 
			{
			double s = double(2.0) * sqrt( double(1.0) + ex_space.x - ey_space.y - ez_space.z);
			quat.x = (ey_space.z - ez_space.y ) / s;
			quat.y = double(0.25) * s;
			quat.z = (ey_space.x + ex_space.y ) / s;
			quat.w = (ez_space.x + ex_space.z ) / s;
			} 
		else if (ey_space.y > ez_space.z) 
			{
			double s = double(2.0) * sqrt( double(1.0) + ey_space.y - ex_space.x - ez_space.z);
			quat.x = (ez_space.x - ex_space.z ) / s;
			quat.y = (ey_space.x + ex_space.y ) / s;
			quat.z = double(0.25) * s;
			quat.w = (ez_space.y + ey_space.z ) / s;
			} 
		else 
			{
			double s = double(2.0) * sqrt( double(1.0) + ez_space.z - ex_space.x - ey_space.y );
			quat.x = (ex_space.y - ey_space.x ) / s;
			quat.y = (ez_space.x + ex_space.z ) / s;
			quat.z = (ez_space.y + ey_space.z ) / s;
			quat.w = double(0.25) * s;
			}
		} 
		
//	cout<<quat.x*quat.x+quat.y*quat.y+quat.z*quat.z+quat.w*quat.w<<endl;		

	}		
	
Molecule::Molecule(unsigned int NatomPerMole): m_NatomPerMole(NatomPerMole)
    {
	allocateData(NatomPerMole);
    }

Molecule::Molecule(const std::string& fname, unsigned int NatomPerMole): m_NatomPerMole(NatomPerMole)
    {
	allocateData(NatomPerMole);
	readData(fname);
    }

void Molecule::allocateData(unsigned int size)
	{
	if(size==0) return;
	
	m_nbond.resize(size);
	m_xyz.resize(size);
	m_mass.resize(size);
	m_charge.resize(size);
	m_inert.resize(size);
	m_orientation.resize(size);
	m_quaternion.resize(size);
	m_diameter.resize(size);
	m_cris.resize(size);
	m_init.resize(size);
	m_molecule.resize(size);
	m_body.resize(size);
	m_be_generated.resize(size);
	m_xyz_read.resize(size);
	m_ori_vec.resize(size);
	m_quat_vec.resize(size);
	m_ori_vec_read.resize(size);
	m_quat_vec_read.resize(size);	
	m_be_generated_read.resize(size);		
	for(unsigned int i=0; i<size; i++)
		{
		m_mass[i] = 1.0;
		m_charge[i] = 0.0;
		m_inert[i] = vec(0.0, 0.0, 0.0);
		m_orientation[i] = 0;
		m_quaternion[i] = 0;
		m_diameter[i] = 0.0;
		m_body[i]=NO_INDEX;
		m_molecule[i]=NO_INDEX;
		m_ori_vec[i] = vec(0.0, 0.0, 0.0);
		m_quat_vec[i] = vec4(0.0, 0.0, 0.0, 0.0);
		m_ori_vec_read[i] = vec(0.0, 0.0, 0.0);
		m_quat_vec_read[i] = vec4(0.0, 0.0, 0.0, 0.0);		
		m_be_generated_read[i]=false;
		m_be_generated[i]=false;	
		}
	m_body_id_plus=0;
	m_mol_id_plus=0;
	m_bondtableHight = 0;
	m_firststep = true;	
	m_dimention = 3;
	m_limit_ge =0;
	m_output_times = 0;
	m_Nread_particle=0;
	m_isotactic = false;
	m_initdata = false;
	m_set_testnum = false;
	m_set_mol_box=false;
	m_check_distance = false;
	m_set_sphere = false;
	m_set_cylinder = false;
	m_set_body_evacuation = false;
	m_include_itself_in_angle = false;
	m_testnum =1;
	m_NBtype=50;
	m_shift_Lx=0.0;
	m_shift_Ly=0.0;
	m_shift_Lz=0.0;	
// initiate random generator	
	srand((int)time(0));
	m_eb_spv=vec(0.0, 0.0, 1.0);
	nwarning_delt = 0;
	}
	
void Molecule::readData(const std::string& fname)
	{
	mst_reader build;
	build.readDataFromMST(fname.c_str());
	unsigned int Np = build.getNParticles();
	std::vector< std::string > typeMapping = build.getTypeMap();
	std::vector<vec> pos = build.getPos();
	std::vector<vec_int> image = build.getImage();
	if(image.size()==0)
		{
		image.resize(pos.size());
		}
	BoxSize box = build.getBox();
	double Lx = double(box.lx);
	double Ly = double(box.ly);		 
	double Lz = double(box.lz);		
	std::vector< double > mass = build.getMass();
	std::vector<unsigned int> type = build.getType();
	std::vector<unsigned int> body = build.getBody();
	std::vector<unsigned int> molecule = build.getMolecule();
	std::vector<double> charge = build.getCharge();
	std::vector<double> diameter = build.getDiameter();	
	std::vector<vec> ori = build.getOrientation();
	std::vector<vec4> quat = build.getQuaternion();
	std::vector<vec> inert = build.getInert();
	std::vector<unsigned int> cris = build.getCris();
	std::vector<unsigned int> init = build.getInit();
	std::vector<str_vec6> asphere = build.getAsphere();	
	std::vector<str_vec6> patch = build.getPatch();
	std::vector<str_vec6> patch_num = build.getPatchNum();		
	m_Nread_particle = 0;
	if(m_NatomPerMole<Np)
		{
		cerr << endl << "***Error! The number of particles read from the file " << Np << " is greater than the initialized number " << m_NatomPerMole<<" !"<< endl << endl;
		throw runtime_error("Molecule::readData error");
		}
	if(type.size()!=Np)
		{
		cerr << endl << "***Error! The number of particle types " << type.size() << " is not equal to the number of particles " << Np<<" !"<< endl << endl;
		throw runtime_error("Molecule::readData error");
		}
	for(unsigned int i=0; i<Np; i++)
		{
		vec posi;
		posi.x = pos[i].x + double(image[i].x)*Lx;
		posi.y = pos[i].y + double(image[i].y)*Ly;
		posi.z = pos[i].z + double(image[i].z)*Lz;			
		m_xyz_read[m_Nread_particle]=posi;
		m_be_generated_read[m_Nread_particle]=true;				
		m_Nread_particle += 1;
		std::string name = typeMapping[type[i]];
		m_type.push_back(name);
		unsigned int id = getTypeId(name);
		m_typeId.push_back(id);
		}
	for(unsigned int i=0; i<mass.size(); i++)
		setMass(i,mass[i]);		
	for(unsigned int i=0; i<body.size(); i++)
		setBody(i,body[i]);
	for(unsigned int i=0; i<molecule.size(); i++)
		setMolecule(i,molecule[i]);	
	for(unsigned int i=0; i<charge.size(); i++)
		setCharge(i,charge[i]);
	for(unsigned int i=0; i<diameter.size(); i++)
		setDiameter(i,diameter[i]);
	for(unsigned int i=0; i<init.size(); i++)
		setInit(i,init[i]);	
	for(unsigned int i=0; i<cris.size(); i++)
		setCris(i,cris[i]);	
	for(unsigned int i=0; i<inert.size(); i++)
		setInert(i,inert[i].x,inert[i].y,inert[i].z);	
	for(unsigned int i=0; i<ori.size(); i++)
		{
		m_ori_vec[i] = ori[i];
		m_ori_vec_read[i] = ori[i];		
		m_orientation[i] = 2;
		}			

	for(unsigned int i=0; i<quat.size(); i++)
		{
		m_quat_vec[i] = quat[i];
		m_quat_vec_read[i] = quat[i];		
		m_quaternion[i] = 2;
		}
	for(unsigned int i=0; i<asphere.size(); i++)
		m_asphere.push_back(asphere[i]);

	for(unsigned int i=0; i<patch.size(); i++)
		m_patch.push_back(patch[i]);

	for(unsigned int i=0; i<patch_num.size(); i++)
		m_patch_num.push_back(patch_num[i]);	

	std::vector<Bond> bonds = build.getBond();
	std::vector<Angle> angles = build.getAngle();
	std::vector<Dihedral> dihedrals = build.getDihedral();
	std::vector<Dihedral> vsites = build.getVsite();	
	for(unsigned int i=0; i<bonds.size(); i++)
		{
		m_bond.push_back(Bond(bonds[i].type, bonds[i].a, bonds[i].b, bonds[i].bc));
		// m_nbond[bonds[i].a] += 1;
		// m_nbond[bonds[i].b] += 1;
		}
	for(unsigned int i=0; i<angles.size(); i++)
		{
		m_angle.push_back(Angle(angles[i].type,angles[i].a,angles[i].b,angles[i].c));
		RadianPerAngle.push_back(0.0);
		}
	for(unsigned int i=0; i<dihedrals.size(); i++)
		{
		m_dihedral.push_back(Dihedral(dihedrals[i].type,dihedrals[i].a,dihedrals[i].b,dihedrals[i].c,dihedrals[i].d));
		}
	for(unsigned int i=0; i<vsites.size(); i++)
		{
		m_vsite.push_back(Dihedral(vsites[i].type, vsites[i].a, vsites[i].b, vsites[i].c, vsites[i].d));
		}
	}		
void Molecule::setParticleTypes(std::string type_str)
	{
	m_type_str = type_str;
	}
void Molecule::initType()
	{
	std::string temp;
	std::string numTypes;
	unsigned column =0;
	bool multiply = false;
    for(unsigned int i =0; i< m_type_str.size(); i++)
		{    
	    if(m_type_str.at(i)!=','&&m_type_str.at(i)!=' '&&!multiply)
			temp.push_back(m_type_str.at(i));
		if(multiply)
			numTypes.push_back(m_type_str.at(i));

		if(m_type_str.at(i)==','||i ==(m_type_str.size()-1))
			{
			column += 1;		
			if(temp.size()==0)
				{
				cout<<"Warning! The void particle type input at column "<<column<<endl;
				if(multiply)
					multiply = false;
				}
			else if(multiply)
				{
				unsigned int num = str2num(numTypes);
				for(unsigned int n = 0; n<num ;n++)
					{
					m_type.push_back(temp);
					unsigned int id = getTypeId(temp);
					m_typeId.push_back(id);			
					}
				numTypes.clear();
				temp.clear();				
				multiply = false;
				}
			else
				{
				m_type.push_back(temp);
				unsigned int id = getTypeId(temp);
				m_typeId.push_back(id);
				temp.clear();
				}
			}
	    if(m_type_str.at(i)=='*')
			{
			multiply = true;
			unsigned int pos = temp.size()-1;
			temp.erase(pos,1);
			}	
		}

	if(m_type.size()!=m_NatomPerMole)
		{
		cerr << endl << "***Error! The number of particle types "<< m_type.size() << " is different from the initialized particle number " << m_NatomPerMole<<" !"<< endl << endl;
		throw runtime_error("Molecule::initType error");	
		}
		
	m_Ntypes = m_type_mapping.size();
	BondLength.resize(m_Ntypes*m_Ntypes);	
	AngleRadian.resize(m_Ntypes*m_Ntypes*m_Ntypes, -1.0);	
	DihedralRadian.resize(m_Ntypes*m_Ntypes*m_Ntypes*m_Ntypes, -1000.0);
	}

void Molecule::setTopology(std::string topology_str)
	{	
	m_topology_str = topology_str;
	}
	
void Molecule::initBond()
	{
	std::string temp_before;
	std::string temp_after;	
	unsigned int column =0;
	bool before = true;
	if( m_type.size()==0)
		{
		cerr << endl << "***Error! Please set particles types first! "  << endl << endl;
		throw runtime_error("Molecule::setTopology error");
		}
	if( m_type.size()!=m_NatomPerMole)
		{
		cerr << endl << "***Error! The number of types "<<m_type.size()<<" is not equal to target number "<< m_NatomPerMole << endl << endl;
		throw runtime_error("Molecule::setTopology error");		
		}			
    for(unsigned int i =0; i< m_topology_str.size(); i++)
		{   
		if(m_topology_str.at(i)=='-')
			before = false;
		if(m_topology_str.at(i)!=','&&m_topology_str.at(i)!=' '&&m_topology_str.at(i)!='-')
			{	
			if(before)
				temp_before.push_back(m_topology_str.at(i));
			else
				temp_after.push_back(m_topology_str.at(i));			
			}
			
		if(m_topology_str.at(i)==','|| i ==(m_topology_str.size()-1))
			{
			column += 1;
			if(temp_before.size()==0||temp_after.size()==0)
				{
				cout<<"Warning, the void topology input at column"<<column <<temp_before<<"-"<<temp_after<< endl;			
				}
			else
				{
				unsigned int bonda = str2num(temp_before);
				unsigned int bondb = str2num(temp_after);
				if(bonda>=m_NatomPerMole||bondb>=m_NatomPerMole||bonda==bondb)
					{
					cerr << endl << "***Error! The wrong particle number at topology input at column " << column <<", "<<temp_before<<"-"<<temp_after<<endl << endl;
					throw runtime_error("Molecule::setTopology error");			
					}
				
				unsigned int type_a = m_typeId[bonda];
				unsigned int type_b = m_typeId[bondb];
			
				string bondname;
				if(type_a<type_b)
					{
					bondname += m_type[bonda];
					bondname.push_back('-');
					bondname += m_type[bondb];
					}
				else
					{
					bondname += m_type[bondb];
					bondname.push_back('-');
					bondname += m_type[bonda];		  
					}
				m_bond.push_back(Bond(bondname,bonda,bondb));
				temp_before.clear();
				temp_after.clear();				
				before = true;
				}
			}
		}

	for(unsigned int i =0; i<m_bond.size();i++)
		{
		unsigned int bonda = m_bond[i].a;
		unsigned int bondb = m_bond[i].b;
		m_nbond[bonda] += 1;
		m_nbond[bondb] += 1;
		}


	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_nbond[i]>m_bondtableHight)
			m_bondtableHight = m_nbond[i];
		}

	vector<unsigned int> m_nbond_temp;
	m_nbond_temp.resize(m_NatomPerMole);
	m_bondtable.resize(m_bondtableHight*m_NatomPerMole);

	for(unsigned int i =0; i<m_bond.size();i++)
		{
		unsigned int bonda = m_bond[i].a;
		unsigned int bondb = m_bond[i].b;
		unsigned int h_a = m_nbond_temp[bonda];
		unsigned int h_b = m_nbond_temp[bondb];
		m_bondtable[bonda*m_bondtableHight+h_a] = bondb;
		m_bondtable[bondb*m_bondtableHight+h_b] = bonda;
		m_nbond_temp[bonda] = h_a +1;
		m_nbond_temp[bondb] = h_b +1;		
		}
	for (unsigned int i =0; i<m_NatomPerMole;i++)
		if(m_nbond_temp[i]!=m_nbond[i])
			cerr << endl << "***Error! bond number conflict!" << endl;
	}
		
void Molecule::initData()
	{
	if(m_initdata)
		return;
	initType();
	initBond();
	for(unsigned int i =0; i<m_bond.size();i++)
		{
		unsigned int bonda = m_bond[i].a;
		unsigned int bondb = m_bond[i].b;
		bool existeda=false;
		bool existedb=false;
		for(unsigned int j=0; j<m_bond_init.size(); j++)
			{
			if(m_bond_init[j]==bonda)
				existeda = true;
			if(m_bond_init[j]==bondb)
				existedb = true;
			}

		if(m_be_generated_read[bonda]&&!existeda)
			m_bond_init.push_back(bonda);
		if(m_be_generated_read[bondb]&&!existedb)
			m_bond_init.push_back(bondb);			
		}
	m_initdata=true;
	}

unsigned int Molecule::switchNametoType(const string& name)
	{
	for(unsigned int i=0; i<m_type_mapping_all.size(); i++)
		{
		if(m_type_mapping_all[i]==name)
			return i;
		}
    cerr << endl << "***Error! Type " << name << " do not exist!" << endl;
    throw runtime_error("Error Molecule switchNametoType");
    return 0;
	}

unsigned int Molecule::getTypeId(const std::string& name) 
    {
    // search for the type mapping
    for (unsigned int i = 0; i < m_type_mapping.size(); i++)
        {
        if (m_type_mapping[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_type_mapping.push_back(name);
    return m_type_mapping.size()-1;
    }

unsigned int Molecule::cellid(int i, int j, int k)
	{
	i = (i + (int)m_dim.x)%(int)m_dim.x;
	j = (j + (int)m_dim.y)%(int)m_dim.y;
	k = (k + (int)m_dim.z)%(int)m_dim.z;	
	return (unsigned int) (i + j*m_dim.x + k*m_dim.x*m_dim.y);
	}

bool Molecule::interMolCheck(unsigned int tag1, vec posa, double& boltz)
	{
	double pair_eng = 0.0;
	unsigned int typi = switchNametoType(m_type[tag1]);
 	double px = posa.x;
	double py = posa.y;
	double pz = posa.z;
		
	double shiftx = 0.0;
	double shifty = 0.0;
	double shiftz = 0.0;		
	if(m_Lx > 0.0)
		{
		shiftx = rint(px/m_Lx);	
		px -=  m_Lx *shiftx;
		}
	if(m_Ly > 0.0)
		{
		shifty = rint(py/m_Ly);
		py -=  m_Ly *shifty;
		}
	if(m_Lz > 0.0)
		{
		shiftz = rint(pz/m_Lz);				
		pz -=  m_Lz *shiftz;
		}
		
	int ix = int((px+0.5*m_Lx)/m_width.x);
	int iy = int((py+0.5*m_Ly)/m_width.y);
	int iz = int((pz+0.5*m_Lz)/m_width.z);			
//	unsigned int cid = cellid(ix, iy, iz);		
	for(int k =-1;  k <=1 ; k++)
		{
		for(int j =-1; j<=1; j++)
			{
			for(int i =-1; i<=1; i++)
				{
				unsigned int jcell = cellid(ix+i, iy+j, iz+k);
				unsigned int id = m_head[jcell];
				while(id!=NO_INDEX)
					{
					double dx = px - m_pos_all[id].x;
					double dy = py - m_pos_all[id].y;
					double dz = pz - m_pos_all[id].z;
					double shiftx = 0.0;
					double shifty = 0.0;
					double shiftz = 0.0;
					if(m_Lx > 0.0)
						{
						shiftx = rint(dx/m_Lx);
						dx -=  m_Lx *shiftx;
						}
					if(m_Ly > 0.0)
						{
						shifty = rint(dy/m_Ly);
						dy -=  m_Ly *shifty;
						}
					if(m_Lz > 0.0)
						{
						shiftz = rint(dz/m_Lz);				
						dz -=  m_Lz *shiftz;
						}						
					double dsq = dx*dx + dy*dy + dz*dz;
					unsigned int typ_pair = typi * m_NBtype + m_type_all[id];
					double min_dis=m_min_dis[typ_pair];					
					if(dsq<min_dis*min_dis)
						return false;

					double lj1 = m_params[typ_pair].x;
					double lj2 = m_params[typ_pair].y;
					double rcutsq = m_params[typ_pair].z;

					if (dsq < rcutsq)
						{
						double r2inv = 1.0/dsq;
						double r6inv = r2inv * r2inv * r2inv;
						pair_eng += r6inv * (lj1 * r6inv - lj2);				
						}
					id = m_list[id];		
					}	
				}
			}
		}
	boltz += exp(-pair_eng/3.741);
	return true;
	}
 
bool Molecule::intraMolCheck(unsigned int tag1, unsigned int tag2, vector<unsigned int>& tags, vec posa, double& boltz)
	{
	double pair_eng = 0.0;	
	unsigned int typi = switchNametoType(m_type[tag1]);
	for (unsigned int i = 0; i < m_NatomPerMole; i++)
		{
		bool exist =false;		
		for (unsigned int j=0; j< tags.size();j++)
			{
			if (i==tags[j])
				exist = true;
			}
		if(m_be_generated[i]&&i!=tag1&&i!=tag2&&!exist)
			{
			double dx = posa.x - m_xyz[i].x;
			double dy = posa.y - m_xyz[i].y;
			double dz = posa.z - m_xyz[i].z;
				
			double shiftx = 0.0;
			double shifty = 0.0;
			double shiftz = 0.0;
						
			if(m_Lx > 0.0)
				{
				shiftx = rint(dx/m_Lx);
				dx -=  m_Lx *shiftx;
				}
			if(m_Ly > 0.0)
				{
				shifty = rint(dy/m_Ly);
				dy -=  m_Ly *shifty;
				}
			if(m_Lz > 0.0)
				{
				shiftz = rint(dz/m_Lz);				
				dz -=  m_Lz *shiftz;
				}				
			double dsq = dx*dx + dy*dy + dz*dz;
			unsigned int typj = switchNametoType(m_type[i]);
			unsigned int typ_pair = typi * m_NBtype + typj;
			double min_dis=m_min_dis[typ_pair];
			if(dsq<min_dis*min_dis)
				return false;			
			double lj1 = m_params[typ_pair].x;
			double lj2 = m_params[typ_pair].y;
			double rcutsq = m_params[typ_pair].z;
			if (dsq < rcutsq)
				{
				double r2inv = 1.0/dsq;
				double r6inv = r2inv * r2inv * r2inv;
				pair_eng += r6inv * (lj1 * r6inv - lj2);				
				}
			}
		}
	boltz += exp(-pair_eng/3.741);		
	return true;
	}

bool Molecule::checkdistance(unsigned int tag1, unsigned int tag2, vector<unsigned int>& tags, vector<vec>& testpos, vec& posa, unsigned int onoff)
	{
	double min_Lx = m_shift_Lx-m_mol_Lx/2.0;
	double max_Lx = m_shift_Lx+m_mol_Lx/2.0;
	
	double min_Ly = m_shift_Ly-m_mol_Ly/2.0;
	double max_Ly = m_shift_Ly+m_mol_Ly/2.0;
	
	double min_Lz = m_shift_Lz-m_mol_Lz/2.0;
	double max_Lz = m_shift_Lz+m_mol_Lz/2.0;	
	
	if (!m_check_distance)
		{
		for(unsigned int i=0; i<testpos.size(); i++ )
			{
			vec posai = testpos[i];
			bool inbox = true;
			if(m_set_mol_box)
				{
				if(m_mol_Lx!=m_Lx&&(posai.x<min_Lx||posai.x>max_Lx))
					inbox=false;
				if(m_mol_Ly!=m_Ly&&(posai.y<min_Ly||posai.y>max_Ly))
					inbox=false;			
				if(m_mol_Lz!=m_Lz&&(posai.z<min_Lz||posai.z>max_Lz))
					inbox=false;
				}
			if(m_set_sphere)
				{
				double dx = posai.x - m_sphere.origin_x;
				double dy = posai.y - m_sphere.origin_y;
				double dz = posai.z - m_sphere.origin_z;

				if(m_Lx > 0.0)
					{
					double img = rint(dx/m_Lx);
					dx -= m_Lx * img;
					}
				if(m_Ly > 0.0)
					{
					double img = rint(dy/m_Ly);
					dy -= m_Ly * img;
					}
				if(m_Lz > 0.0)
					{
					double img = rint(dz/m_Lz);
					dz -= m_Lz * img;
					}

				double rsq = dx*dx + dy*dy + dz*dz;
				double r =sqrt(rsq);
				
				if(r>=m_sphere.radius_max||r<=m_sphere.radius_min)
					inbox=false;
				}
			if(m_set_cylinder)
				{
				double distFromOri = m_cylinder.direction_x * (posai.x - m_cylinder.origin_x) 
								   + m_cylinder.direction_y * (posai.y - m_cylinder.origin_y)
								   + m_cylinder.direction_z * (posai.z - m_cylinder.origin_z);
					
					// use the distance to create a vector pointing from the central line of cylinder to the particle
				double dx = posai.x - m_cylinder.origin_x - m_cylinder.direction_x * distFromOri;
				double dy = posai.y - m_cylinder.origin_y - m_cylinder.direction_y * distFromOri;
				double dz = posai.z - m_cylinder.origin_z - m_cylinder.direction_z * distFromOri;
				
				if(m_Lx > 0.0)
					{
					double img = rint(dx/m_Lx);
					dx -= m_Lx * img;
					}
				if(m_Ly > 0.0)
					{
					double img = rint(dy/m_Ly);
					dy -= m_Ly * img;
					}
				if(m_Lz > 0.0)
					{
					double img = rint(dz/m_Lz);
					dz -= m_Lz * img;
					}
				
				double rsq = dx*dx + dy*dy + dz*dz;
				double r =sqrt(rsq);

				if(r>=m_cylinder.radius_max||r<=m_cylinder.radius_min)
					inbox=false;
				}	

			if(m_set_body_evacuation)
				{
				for(unsigned int i=0; i< m_body_com_all.size(); i++)	
					{
					vec4 bc = m_body_com_all[i];
					double dx = posai.x - bc.x;
					double dy = posai.y - bc.y;
					double dz = posai.z - bc.z;
					
					if(m_Lx > 0.0)
						{
						double img = rint(dx/m_Lx);
						dx -= m_Lx * img;
						}
					if(m_Ly > 0.0)
						{
						double img = rint(dy/m_Ly);
						dy -= m_Ly * img;
						}
					if(m_Lz > 0.0)
						{
						double img = rint(dz/m_Lz);
						dz -= m_Lz * img;
						}
								
					double rsq = dx*dx + dy*dy + dz*dz;
					double r =sqrt(rsq);
					if(r<bc.w)
						{
						inbox=false;
						break;
						}
					}
				}					
				
			if(inbox)
				{
				posa.x = posai.x;
				posa.y = posai.y;
				posa.z = posai.z;
				return true;
				}
			}
		return false;
		}
	vector<vec> goodpos;
	vector<double> boltz;

	for(unsigned int i=0; i<testpos.size(); i++ )
		{
		double boltzi=0.0;
		vec posai = testpos[i];
		if(m_set_mol_box)
			{
			if(m_mol_Lx!=m_Lx&&(posai.x<min_Lx||posai.x>max_Lx))
				continue;
			if(m_mol_Ly!=m_Ly&&(posai.y<min_Ly||posai.y>max_Ly))
				continue;			
			if(m_mol_Lz!=m_Lz&&(posai.z<min_Lz||posai.z>max_Lz))
				continue;			
			}
		if(m_set_sphere)
			{
			double dx = posai.x - m_sphere.origin_x;
			double dy = posai.y - m_sphere.origin_y;
			double dz = posai.z - m_sphere.origin_z;
			
			if(m_Lx > 0.0)
				{
				double img = rint(dx/m_Lx);
				dx -= m_Lx * img;
				}
			if(m_Ly > 0.0)
				{
				double img = rint(dy/m_Ly);
				dy -= m_Ly * img;
				}
			if(m_Lz > 0.0)
				{
				double img = rint(dz/m_Lz);
				dz -= m_Lz * img;
				}
			
			double rsq = dx*dx + dy*dy + dz*dz;
			double r =sqrt(rsq);
			if(r>=m_sphere.radius_max||r<=m_sphere.radius_min)
				continue;
			}
		if(m_set_cylinder)
			{
			double distFromOri = m_cylinder.direction_x * (posai.x - m_cylinder.origin_x) 
							   + m_cylinder.direction_y * (posai.y - m_cylinder.origin_y)
							   + m_cylinder.direction_z * (posai.z - m_cylinder.origin_z);
					
					// use the distance to create a vector pointing from the central line of cylinder to the particle
			double dx = posai.x - m_cylinder.origin_x - m_cylinder.direction_x * distFromOri;
			double dy = posai.y - m_cylinder.origin_y - m_cylinder.direction_y * distFromOri;
			double dz = posai.z - m_cylinder.origin_z - m_cylinder.direction_z * distFromOri;
			
			if(m_Lx > 0.0)
				{
				double img = rint(dx/m_Lx);
				dx -= m_Lx * img;
				}
			if(m_Ly > 0.0)
				{
				double img = rint(dy/m_Ly);
				dy -= m_Ly * img;
				}
			if(m_Lz > 0.0)
				{
				double img = rint(dz/m_Lz);
				dz -= m_Lz * img;
				}
						
			
			double rsq = dx*dx + dy*dy + dz*dz;
			double r =sqrt(rsq);
				
			if(r>=m_cylinder.radius_max||r<=m_cylinder.radius_min)
				continue;
			}

		if(m_set_body_evacuation)
			{
			bool in_body=false;
			for(unsigned int i=0; i< m_body_com_all.size(); i++)	
				{
				vec4 bc = m_body_com_all[i];

				double dx = posai.x - bc.x;
				double dy = posai.y - bc.y;
				double dz = posai.z - bc.z;
				
				if(m_Lx > 0.0)
					{
					double img = rint(dx/m_Lx);
					dx -= m_Lx * img;
					}
				if(m_Ly > 0.0)
					{
					double img = rint(dy/m_Ly);
					dy -= m_Ly * img;
					}
				if(m_Lz > 0.0)
					{
					double img = rint(dz/m_Lz);
					dz -= m_Lz * img;
					}

				double rsq = dx*dx + dy*dy + dz*dz;
				double r =sqrt(rsq);
				if(r<bc.w)
					{
					in_body=true;
					break;
					}
				}
			if(in_body)
				continue;
			}				

		bool check0=true;
		bool check1=true;
		if(onoff==0)
			check0 = interMolCheck(tag1, posai, boltzi);
		else if(onoff==1)
			check1 = intraMolCheck(tag1, tag2, tags, posai, boltzi);
		else if(onoff==2)
			{
			check0 = interMolCheck(tag1, posai, boltzi);
			check1 = intraMolCheck(tag1, tag2, tags, posai, boltzi);
			}
		if(check0&&check1)
			{
			goodpos.push_back(posai);
			boltz.push_back(boltzi);
			}
		}
	if(goodpos.size()==0)
		return false;
	unsigned int max_id=0;
	double max_value =0.0;
	for(unsigned int i=0; i<boltz.size(); i++ )
		{
		double boltzi=boltz[i];
//		cout<<boltzi<<endl;
		if(boltzi>=max_value)
			{
			max_value = boltzi;
			max_id = i;
			}
		}	
	posa.x=goodpos[max_id].x;
	posa.y=goodpos[max_id].y;
	posa.z=goodpos[max_id].z;
		return true;
	}

void Molecule::setPosition(unsigned int i, double px,double py, double pz)
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set position for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setPosition error");		
		}	
	m_xyz_read[i]=vec(px, py, pz);
	m_be_generated_read[i]=true;
	}

void Molecule::setCharge(double charge)
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_charge[i] = charge;
	}	

void Molecule::setCharge(std::string type, double charge)
	{
	if(m_type_str.size()==0&&m_type.size()==0)
		{
		cerr << endl << "***Error! No type is given before! " << endl << endl;
		throw runtime_error("Molecule::setCharge error");			
		}			
	initData();
	unsigned int id = getTypeId(type);	
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_typeId[i] == id)
			m_charge[i] = charge;	
		}
	}

void Molecule::setCharge(unsigned int i, double charge)
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set charge for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setCharge error");		
		}	
	m_charge[i] = charge;
	}	

void Molecule::setChargeFactor(double factor)
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_charge[i] *= factor;
	}		
	
void Molecule::setMass(double mass)
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_mass[i] = mass;
	}

void Molecule::setMass(std::string type, double mass)
	{
	if(m_type_str.size()==0&&m_type.size()==0)
		{
		cerr << endl << "***Error! No type is given before! " << endl << endl;
		throw runtime_error("Molecule::setMass error");			
		}			
	initData();
	unsigned int id = getTypeId(type);	
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_typeId[i] == id)
			m_mass[i] = mass;	
		}
	}

void Molecule::setMass(unsigned int i, double mass)
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set mass for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setMass error");		
		}	
	m_mass[i] = mass;
	}

void Molecule::setInert(double inertx, double inerty, double inertz)
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_inert[i] = vec(inertx, inerty, inertz);
	}

void Molecule::setInert(std::string type, double inertx, double inerty, double inertz)
	{
	if(m_type_str.size()==0&&m_type.size()==0)
		{
		cerr << endl << "***Error! No type is given before! " << endl << endl;
		throw runtime_error("Molecule::setInert error");			
		}			
	initData();
	unsigned int id = getTypeId(type);	
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_typeId[i] == id)
			m_inert[i] = vec(inertx, inerty, inertz);	
		}
	}
	
void Molecule::setInert(unsigned int i, double inertx, double inerty, double inertz)		
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set orientation for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setOrientation error");		
		}	
	m_inert[i] = vec(inertx, inerty, inertz);	
	}	
	
void Molecule::setOrientation()
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_orientation[i] = 1;
	}

void Molecule::setOrientation(std::string type)
	{
	if(m_type_str.size()==0&&m_type.size()==0)
		{
		cerr << endl << "***Error! No type is given before! " << endl << endl;
		throw runtime_error("Molecule::setOrientation error");			
		}			
	initData();
	unsigned int id = getTypeId(type);	
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_typeId[i] == id)
			m_orientation[i] = 1;	
		}
	}

void Molecule::setOrientation(unsigned int i)
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set orientation for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setOrientation error");		
		}	
	m_orientation[i] = 1;
	}
	
void Molecule::setQuaternion()
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_quaternion[i] = 1;
	}

void Molecule::setQuaternion(std::string type)
	{
	if(m_type_str.size()==0&&m_type.size()==0)
		{
		cerr << endl << "***Error! No type is given before! " << endl << endl;
		throw runtime_error("Molecule::setQuaternion error");			
		}			
	initData();
	unsigned int id = getTypeId(type);	
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_typeId[i] == id)
			m_quaternion[i] = 1;
		}
	}

void Molecule::setQuaternion(unsigned int i)
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set orientation for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setQuaternion error");		
		}	
	m_quaternion[i] = 1;
	}	

void Molecule::setDiameter( double diameter)
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_diameter[i] = diameter;	
	}

void Molecule::setDiameter(std::string type, double diameter)
	{
	if(m_type_str.size()==0&&m_type.size()==0)
		{
		cerr << endl << "***Error! No type is given before! " << endl << endl;
		throw runtime_error("Molecule::setDiameter error");			
		}			
	initData();
	unsigned int id = getTypeId(type);	
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_typeId[i] == id)
			m_diameter[i] = diameter;	
		}
	}

void Molecule::setDiameter(unsigned int i, double diameter)
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set diameter for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setDiameter error");		
		}	
	m_diameter[i] = diameter;
	}

void Molecule::setCris( unsigned int cris)
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_cris[i] = cris;	
	}

void Molecule::setCris(std::string type, unsigned int cris)
	{
	if(m_type_str.size()==0&&m_type.size()==0)
		{
		cerr << endl << "***Error! No type is given before! " << endl << endl;
		throw runtime_error("Molecule::setCris error");			
		}			
	initData();
	unsigned int id = getTypeId(type);	
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_typeId[i] == id)
			m_cris[i] = cris;	
		}
	}

void Molecule::setCris(unsigned int i, unsigned int cris)
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set cris for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setCris error");		
		}	
	m_cris[i] = cris;
	}

void Molecule::setInit( unsigned int init)
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_init[i] = init;	
	}

void Molecule::setInit(std::string type, unsigned int init)
	{
	if(m_type_str.size()==0&&m_type.size()==0)
		{
		cerr << endl << "***Error! No type is given before! " << endl << endl;
		throw runtime_error("Molecule::setInit error");			
		}			
	initData();
	unsigned int id = getTypeId(type);	
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_typeId[i] == id)
			m_init[i] = init;	
		}
	}

void Molecule::setInit(unsigned int i, unsigned int init)
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set init for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setInit error");		
		}	
	m_init[i] = init;
	}

void Molecule::setBody(unsigned int body)
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_body[i] = body;
	if(body!=NO_BODY&&body+1>m_body_id_plus)
		m_body_id_plus=body+1;
	}
	
void Molecule::setBody(std::string type, unsigned int body)
	{
	if(m_type_str.size()==0&&m_type.size()==0)
		{
		cerr << endl << "***Error! No type is given before! " << endl << endl;
		throw runtime_error("Molecule::setBody error");			
		}			
	initData();
	unsigned int id = getTypeId(type);	
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_typeId[i] == id)
			m_body[i] = body;	
		}
	if(body!=NO_BODY&&body+1>m_body_id_plus)
		m_body_id_plus=body+1;
	}

void Molecule::setBody( unsigned int i, unsigned int body)
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set init for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setBody error");		
		}
	m_body[i] = body;
	if(body!=NO_BODY&&body+1>m_body_id_plus)
		m_body_id_plus=body+1;
	}

void Molecule::setMolecule(unsigned int molecule)
	{
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		m_molecule[i] = molecule;
	if(molecule!=NO_INDEX&&molecule+1>m_mol_id_plus)
		m_mol_id_plus=molecule+1;
	}
	
void Molecule::setMolecule(std::string type, unsigned int molecule)
	{
	if(m_type_str.size()==0&&m_type.size()==0)
		{
		cerr << endl << "***Error! No type is given before! " << endl << endl;
		throw runtime_error("Molecule::setMolecule error");			
		}		
	initData();
	unsigned int id = getTypeId(type);	
	for(unsigned int i =0; i<m_NatomPerMole; i++)
		{
		if(m_typeId[i] == id)
			m_molecule[i] = molecule;	
		}
	if(molecule!=NO_INDEX&&molecule+1>m_mol_id_plus)
		m_mol_id_plus=molecule+1;
	}

void Molecule::setMolecule( unsigned int i, unsigned int molecule)
	{
	if(i>=m_NatomPerMole)
		{
		cerr << endl << "***Error! set init for a non-existed particle " << i << endl << endl;
		throw runtime_error("Molecule::setMolecule error");		
		}
	m_molecule[i] = molecule;
	if(molecule!=NO_INDEX&&molecule+1>m_mol_id_plus)
		m_mol_id_plus=molecule+1;
	}
	
void Molecule::setEllipsoidBondSpotVector(double spvx, double spvy, double spvz)
	{
	m_eb_spv=vec(spvx, spvy, spvz);
	}
	
void Molecule::setIsotactic(bool iso)
	{
	m_isotactic = iso;
	}

void Molecule::setDimention(unsigned int dimention)
	{
	if(dimention!=1&&dimention!=2&&dimention!=3)
		{
		cerr << endl << "***Error! Dimention should be 1 or 2 or 3! " << dimention << endl << endl;
		throw runtime_error("Molecule::setDimention error");		
		}
	m_dimention = dimention;
	}

void Molecule::setPutBox(double Lx, double Ly, double Lz)
	{
	m_Lx = Lx;
	m_Ly = Ly;
	m_Lz = Lz;
	}

void Molecule::setBox(double mol_Lx,double mol_Ly, double mol_Lz)
	{
	m_mol_Lx = mol_Lx;
	m_mol_Ly = mol_Ly;
	m_mol_Lz = mol_Lz;
	m_set_mol_box=true;
	}

void Molecule::setBox(double mol_Lx_min, double mol_Lx_max, double mol_Ly_min, double mol_Ly_max, double mol_Lz_min, double mol_Lz_max)
	{
	m_mol_Lx = mol_Lx_max - mol_Lx_min;
	m_mol_Ly = mol_Ly_max - mol_Ly_min;
	m_mol_Lz = mol_Lz_max - mol_Lz_min;
	m_shift_Lx = (mol_Lx_max + mol_Lx_min)/2;
	m_shift_Ly = (mol_Ly_max + mol_Ly_min)/2;
	m_shift_Lz = (mol_Lz_max + mol_Lz_min)/2;
	if(m_mol_Lx < 0.0)
		{
		cerr << endl << "***Error! Trying to set mol_Lx_max = " <<mol_Lx_max <<" less than mol_Lx_min = "<<mol_Lx_min<< endl << endl;
		throw runtime_error("Molecule::setBox error");		
		}
	if(m_mol_Ly < 0.0)
		{
		cerr << endl << "***Error! Trying to set mol_Ly_max = " <<mol_Ly_max <<" less than mol_Ly_min = "<<mol_Ly_min<< endl << endl;
		throw runtime_error("Molecule::setBox error");		
		}
	if(m_mol_Lz < 0.0)
		{
		cerr << endl << "***Error! Trying to set mol_Lz_max = " <<mol_Lz_max <<" less than mol_Lz_min = "<<mol_Lz_min<< endl << endl;
		throw runtime_error("Molecule::setBox error");		
		}		
	m_set_mol_box=true;
	}

void Molecule::setSphere(double ox,double oy, double oz, double r_min, double r_max)
	{
	if(r_min > r_max)
		{
		cerr << endl << "***Error! Trying to set sphere with inner radius " <<r_min <<" great than outer radius "<<r_max<< endl << endl;
		throw runtime_error("Molecule::setSphere error");		
		}	
	m_sphere = Sphere(ox, oy, oz, r_min, r_max);
	m_set_sphere = true;
	}

void Molecule::setCylinder(double ox,double oy, double oz, double dx,double dy, double dz, double r_min, double r_max)
	{
	if(r_min > r_max)
		{
		cerr << endl << "***Error! Trying to set cylinder with inner radius " <<r_min <<" great than outer radius "<<r_max<< endl << endl;
		throw runtime_error("Molecule::setSphere error");		
		}	
	m_cylinder = Cylinder(ox, oy, oz, dx, dy, dz, r_min, r_max);
	m_set_cylinder = true;
	}
	
void Molecule::setBodyEvacuation()
	{
	m_set_body_evacuation = true;
	}	

void Molecule::setBondLength( double length)
	{
	if(m_topology_str.size()==0)
		{
		cerr << endl << "***Error! No topology is given before! " << endl << endl;
		throw runtime_error("Molecule::setBondLength error");			
		}
		
	initData();
	if(length <= 0.0)
		{
		cerr << endl << "***Error! Trying to set bond length less than or equal to zero! " << length << endl << endl;
		throw runtime_error("Molecule::setBondLength error");		
		}
	for(unsigned int i=0; i<m_Ntypes*m_Ntypes; i++)
		BondLength[i] = length;
	}

void Molecule::setBondLength(std::string name1, std::string name2, double length)
	{
	if(m_topology_str.size()==0)
		{
		cerr << endl << "***Error! No topology is given before!" << endl << endl;
		throw runtime_error("Molecule::setBondLength error");			
		}			
		
	initData();
	unsigned int typ1= getTypeId(name1);
	unsigned int typ2= getTypeId(name2);

	if (typ1 >= m_Ntypes || typ2 >= m_Ntypes)
		{
		cerr << endl << "***Error! Trying to set setBondLength for a non existant type! " << name1 << "," << name2 << endl << endl;
		throw runtime_error("Molecule::setBondLength error");
		}
	if(length <= 0.0)
		{
		cerr << endl << "***Error! Trying to set bond length less than or equal to zero! " << length << endl << endl;
		throw runtime_error("Molecule::setBondLength error");		
		}
	
	BondLength[typ1 + typ2*m_Ntypes] = length;
	BondLength[typ2 + typ1*m_Ntypes] = length;
	}

void Molecule::setAngleDegree(std::string name_a, std::string name_b, std::string name_c, double degree)
	{
	if(m_topology_str.size()==0)
		{
		cerr << endl << "***Error! No topology is given before! " << endl << endl;
		throw runtime_error("Molecule::setAngleDegree error");			
		}		
	initData();
	unsigned int typ_a= getTypeId(name_a);
	unsigned int typ_b= getTypeId(name_b);
	unsigned int typ_c= getTypeId(name_c);
	if (typ_a >= m_Ntypes || typ_b >= m_Ntypes|| typ_c >= m_Ntypes)
		{
		cerr << endl << "***Error! Trying to set setAngleDegree for a non existant type! " << name_a << "," << name_b << "," << name_c <<endl << endl;
		throw runtime_error("Molecule::setAngleDegree error");
		}
	if(degree < 0.0||degree>180.0)
		{
		cerr << endl << "***Error! Trying to set angle degree not in range (0 - 180]! " << degree << endl << endl;
		throw runtime_error("Molecule::setAngleDegree error");		
		}
	double radian =  M_PI*degree/180.0;		
	AngleRadian[typ_a + typ_b*m_Ntypes + typ_c*m_Ntypes*m_Ntypes] = radian;
	AngleRadian[typ_c + typ_b*m_Ntypes + typ_a*m_Ntypes*m_Ntypes] = radian;
	}

void Molecule::setAngleDegree(unsigned int i,unsigned int j,unsigned int k, double degree)
	{
	if(m_topology_str.size()==0)
		{
		cerr << endl << "***Error! No topology is given before! " << endl << endl;
		throw runtime_error("Molecule::setAngleDegree error");			
		}		
	initData();
	if (i >= m_NatomPerMole || j >= m_NatomPerMole|| k >= m_NatomPerMole)
		{
		cerr << endl << "***Error! Trying to set setAngleDegree for a non existant particle! " << i << " , " << j << " , " << k <<endl << endl;
		throw runtime_error("Molecule::setAngleDegree error");
		}
	if(i==j||i==k||i==k)
		{
		cerr << endl << "***Error! Trying to set setAngleDegree for a non existant angle! " << i << " , " << j << " , " << k <<endl << endl;
		throw runtime_error("Molecule::setAngleDegree error");
		}
	if(degree < 0.0||degree>180.0)
		{
		cerr << endl << "***Error! Trying to set angle degree not in range (0 - 180]! " << degree << endl << endl;
		throw runtime_error("Molecule::setAngleDegree error");		
		}	
	unsigned int typ_aa = m_typeId[i];
	unsigned int typ_cc = m_typeId[k];
	std::string anglename;
	if(typ_aa<typ_cc)
		{
		anglename += m_type[i];
		anglename.push_back('-');
		anglename += m_type[j];
		anglename.push_back('-');
		anglename += m_type[k];			
		}
	else
		{
		anglename += m_type[k];
		anglename.push_back('-');
		anglename += m_type[j];
		anglename.push_back('-');
		anglename += m_type[i];			
		}
	double radian =  M_PI*degree/180.0;	
	RadianPerAngle.push_back(radian);
	m_angle.push_back(Angle(anglename,i,j,k));		
	}

void Molecule::setDihedralDegree(std::string name_a, std::string name_b, std::string name_c, std::string name_d, double degree)
	{
	if(m_topology_str.size()==0)
		{
		cerr << endl << "***Error! No topology is given before! " << endl << endl;
		throw runtime_error("Molecule::setDihedralDegree error");			
		}		
	initData();
	unsigned int typ_a= getTypeId(name_a);
	unsigned int typ_b= getTypeId(name_b);
	unsigned int typ_c= getTypeId(name_c);
	unsigned int typ_d= getTypeId(name_d);	   
	if (typ_a >= m_Ntypes || typ_b >= m_Ntypes|| typ_c >= m_Ntypes|| typ_d >= m_Ntypes)
		{
		cerr << endl << "***Error! Trying to set setDihedralDegree for a non existant type! " << name_a << "," << name_b << "," << name_c<< "," << name_d <<endl << endl;
		throw runtime_error("Molecule::setDihedralDegree error");
		}
	if(degree <= -180||degree>180.0)
		{
		cerr << endl << "***Error! Trying to set dihedral degree not in range (-180 - 180] ! " << degree << endl << endl;
		throw runtime_error("Molecule::setDihedralDegree error");		
		}
	double radian =  M_PI*degree/180.0;		
	DihedralRadian[typ_a + typ_b*m_Ntypes + typ_c*m_Ntypes*m_Ntypes + typ_d*m_Ntypes*m_Ntypes*m_Ntypes] = radian;
	DihedralRadian[typ_d + typ_c*m_Ntypes + typ_b*m_Ntypes*m_Ntypes + typ_a*m_Ntypes*m_Ntypes*m_Ntypes] = radian;
	}	

void Molecule::setDihedralDegree(unsigned int i,unsigned int j,unsigned int k,unsigned int l, double degree)
	{
	if(m_topology_str.size()==0)
		{
		cerr << endl << "***Error! No topology is given before! " << endl << endl;
		throw runtime_error("Molecule::setDihedralDegree error");			
		}		
	initData();
	if (i >= m_NatomPerMole || j >= m_NatomPerMole|| k >= m_NatomPerMole|| l>=m_NatomPerMole)
		{
		cerr << endl << "***Error! Trying to set setDihedralDegree for a non existant particle! " << i << " , " << j << " , " << k<< " , " << l <<endl << endl;
		throw runtime_error("Molecule::setDihedralDegree error");
		}
	if(i==j||j==k||k==l||i==k||i==l||j==l)
		{
		cerr << endl << "***Error! Trying to set setDihedralDegree for a non existant angle! " << i << " , " << j << " , " << k<< " , " << l <<endl << endl;
		throw runtime_error("Molecule::setDihedralDegree error");
		}
	if(degree <= -180||degree>180.0)
		{
		cerr << endl << "***Error! Trying to set dihedral degree not in range (-180 - 180] ! " << degree << endl << endl;
		throw runtime_error("Molecule::setDihedralDegree error");		
		}
	unsigned int typ_aa = m_typeId[i];
	unsigned int typ_dd = m_typeId[l];
	std::string dihedralname;
	if(typ_aa<typ_dd)
		{
		dihedralname += m_type[i];
		dihedralname.push_back('-');
		dihedralname += m_type[j];
		dihedralname.push_back('-');
		dihedralname += m_type[k];
		dihedralname.push_back('-');
		dihedralname += m_type[l];										
		}
	else
		{
		dihedralname += m_type[l];
		dihedralname.push_back('-');
		dihedralname += m_type[k];
		dihedralname.push_back('-');
		dihedralname += m_type[j];
		dihedralname.push_back('-');
		dihedralname += m_type[i];		  
		}
	double radian =  M_PI*degree/180.0;			
	m_dihedral.push_back(Dihedral(dihedralname,i,j,k,l));
	RadianPerDihedral.push_back(radian);
	}

double Molecule::existedAngleDigree(unsigned int a, unsigned int b, unsigned int c)
	{
	for(unsigned int i =0; i< m_angle.size(); i++)
		{
		if(m_angle[i].a==a&&m_angle[i].b==b&&m_angle[i].c==c)
			return RadianPerAngle[i];
		else if(m_angle[i].a==c&&m_angle[i].b==b&&m_angle[i].c==a)
			return RadianPerAngle[i];
		}
	return 0;
	}	

int Molecule::str2num(string s)
	{   
	int num;
	stringstream ss(s);
	ss>>num;
	return num;
	}
	
double Molecule::R2S()
	{
	int ran = rand();
	double fran = (double)ran/(double)RAND_MAX;		
	return fran;
	}

bool Molecule::threeAnglesFixE(vec A, vec B, vec C, vec D, vec& E, double lenthBE, double thetaABE,double thetaCBE,double thetaDBE)	
	{
	vec E1,E2,E3,E4;
	if(twoAnglesFixD(A,B,C,E1,E2,lenthBE,thetaABE,thetaCBE)&&twoAnglesFixD(A,B,D,E3,E4,lenthBE,thetaABE,thetaDBE))
		{
		double R13 = (E1.x-E3.x)*(E1.x-E3.x) + (E1.y-E3.y)*(E1.y-E3.y) + (E1.z-E3.z)*(E1.z-E3.z);
		double R14 = (E1.x-E4.x)*(E1.x-E4.x) + (E1.y-E4.y)*(E1.y-E4.y) + (E1.z-E4.z)*(E1.z-E4.z);
		double R23 = (E2.x-E3.x)*(E2.x-E3.x) + (E2.y-E3.y)*(E2.y-E3.y) + (E2.z-E3.z)*(E2.z-E3.z);
		double R24 = (E2.x-E4.x)*(E2.x-E4.x) + (E2.y-E4.y)*(E2.y-E4.y) + (E2.z-E4.z)*(E2.z-E4.z);
		if(R13<0.001||R14<0.001)
			{
			E = E1;
			return true;
			}
		else if(R23<0.001||R24<0.001)
			{
			E = E2;
			return true;
			}
		else
			return false;
		}
	else
		return false;
	}

bool Molecule::threeAnglesFixG(vec A, vec B, vec C, vec D, vec E, vec F, vec& G, 
								double lenthBG, double lenthCG,double lenthEG, 
								double thetaABG,double thetaDCG,double thetaFEG)	
	{
	vec G1,G2,G3,G4;
	if(twoAnglesFixE(A,B,C,D,G1,G2,lenthBG,lenthCG,thetaABG,thetaDCG)&&twoAnglesFixE(A,B,E,F,G3,G4,lenthBG,lenthEG,thetaABG,thetaFEG))
		{
		double R13 = (G1.x-G3.x)*(G1.x-G3.x) + (G1.y-G3.y)*(G1.y-G3.y) + (G1.z-G3.z)*(G1.z-G3.z);
		double R14 = (G1.x-G4.x)*(G1.x-G4.x) + (G1.y-G4.y)*(G1.y-G4.y) + (G1.z-G4.z)*(G1.z-G4.z);
		double R23 = (G2.x-G3.x)*(G2.x-G3.x) + (G2.y-G3.y)*(G2.y-G3.y) + (G2.z-G3.z)*(G2.z-G3.z);
		double R24 = (G2.x-G4.x)*(G2.x-G4.x) + (G2.y-G4.y)*(G2.y-G4.y) + (G2.z-G4.z)*(G2.z-G4.z);
		if(R13<0.001||R14<0.001)
			{
			G = G1;
			return true;
			}
		else if(R23<0.001||R24<0.001)
			{
			G = G2;
			return true;
			}
		else
			return false;
		}
	else
		return false;
	}

bool Molecule::twoAnglesFixE(vec A, vec B, vec C, vec D, vec& E1, vec& E2, double lenthBE, double lenthCE, double thetaABE, double thetaDCE)	
	{
	vec BA, CD;
	BA.x = A.x-B.x;
	BA.y = A.y-B.y;
	BA.z = A.z-B.z;
	double lenthBA = sqrt(BA.x*BA.x+BA.y*BA.y+BA.z*BA.z);
	
	CD.x = D.x-C.x;
	CD.y = D.y-C.y;
	CD.z = D.z-C.z;	
	double lenthCD = sqrt(CD.x*CD.x+CD.y*CD.y+CD.z*CD.z);
	if (lenthBA==0 || lenthCD==0)
		{
		cout<<"lenth = 0"<<endl;
		return false;
		}
	vec a,b,c,d,e;
	a.x=BA.x; a.y=CD.x; a.z=1.0;
	b.x=BA.y; b.y=CD.y; b.z=1.0;	
	c.x=BA.z; c.y=CD.z; c.z=1.0;	
	d.x=cos(thetaABE)*lenthBA*lenthBE+B.x*BA.x+B.y*BA.y+B.z*BA.z;
	d.y=cos(thetaDCE)*lenthCD*lenthCE+C.x*CD.x+C.y*CD.y+C.z*CD.z;
	d.z=lenthBE*lenthBE;
	e.x=B.x; e.y=B.y; e.z=B.z;
	if(arrayFixF(a, b, c, d, e, E1, E2))
		{
		return true;
		}
	return false;
	}
	
bool Molecule::twoAnglesFixE(vec A, vec B, vec C, vec D, vec& E, double lenthBE, double lenthCE, double thetaABE, double thetaDCE)	
	{
	vec E1, E2;
	if(twoAnglesFixE(A, B, C, D, E1, E2, lenthBE,  lenthCE,  thetaABE, thetaDCE))
		{
		if(m_isotactic)
			{
			double value1=isotactic(A,B,D,E1);
			double value2=isotactic(A,B,D,E2);			
			if (value1>=0.0)
				E=E1;
			else if (value2>=0.0)
				E=E2;
			else
				throw runtime_error("Molecule::twoAnglesFixE: wrong isotactic compute!");
			return true;
			}
		
		if(R2S()<0.5)
			E=E1;
		else
			E=E2;
		return true;
		}
	else
		return false;
	}

bool Molecule::twoAnglesFixD(vec A, vec B, vec C, vec& D1, vec& D2, double lenthBD, double thetaABD,double thetaCBD)
	{	
	vec BA, BC;
	BA.x = A.x-B.x;
	BA.y = A.y-B.y;
	BA.z = A.z-B.z;
	double lenthBA = sqrt(BA.x*BA.x+BA.y*BA.y+BA.z*BA.z);
	
	BC.x = C.x-B.x;
	BC.y = C.y-B.y;
	BC.z = C.z-B.z;	
	double lenthBC = sqrt(BC.x*BC.x+BC.y*BC.y+BC.z*BC.z);
	if (lenthBA==0 || lenthBC==0)
		{
		cout<<"lenth = 0"<<endl;
		return false;
		}
	vec a,b,c,d,e;
	a.x=BC.x; a.y=BA.x; a.z=1.0;
	b.x=BC.y; b.y=BA.y; b.z=1.0;	
	c.x=BC.z; c.y=BA.z; c.z=1.0;	
	d.x=cos(thetaCBD)*lenthBC*lenthBD+B.x*BC.x+B.y*BC.y+B.z*BC.z;
	d.y=cos(thetaABD)*lenthBA*lenthBD+B.x*BA.x+B.y*BA.y+B.z*BA.z;
	d.z=lenthBD*lenthBD;
	e.x=B.x; e.y=B.y; e.z=B.z;
	if(arrayFixF(a, b, c, d, e, D1, D2))
		return true;
	return false;
	}

double Molecule::isotactic(vec A, vec B, vec C, vec D)
	{
	vec DA, DB, DC;
	DA.x = A.x -D.x;
	DA.y = A.y -D.y;	
	DA.z = A.z -D.z;	

	DB.x = B.x -D.x;
	DB.y = B.y -D.y;	
	DB.z = B.z -D.z;
	
	DC.x = C.x -D.x;
	DC.y = C.y -D.y;	
	DC.z = C.z -D.z;
	vec PV;
	PV.x = DA.y*DC.z-DA.z*DC.y;
	PV.y = DA.z*DC.x-DA.x*DC.z;
	PV.z = DA.x*DC.y-DA.y*DC.x;

	double value = PV.x*DB.x+PV.y*DB.y+PV.z*DB.z;
	return value;
	}

bool Molecule::twoAnglesFixD(vec A, vec B, vec C, vec& D, double lenthBD, double thetaABD,double thetaCBD)	
	{
	vec D1, D2;
	if(twoAnglesFixD(A,B,C,D1,D2,lenthBD,thetaABD,thetaCBD))
		{
		if(R2S()<0.5)
			D=D1;
		else
			D=D2;
		return true;
		}
	else
		return false;
	}

//a.x*x+b.x*y+c.x*z=d.x;  
//a.y*x+b.y*y+c.y*z=d.y;
//a.z*(x-e.x)^2+b.z*(y-e.y)^2+c.z*(z-e.z)^2=d.z;
bool Molecule::arrayFixF(vec a, vec b, vec c, vec d, vec e, vec& F1, vec& F2)       
	{
	double Ox=b.x*c.y-b.y*c.x;
	double Oy=a.x*c.y-a.y*c.x;
	double Oz=a.x*b.y-a.y*b.x;
	
	if (Ox==0&&Oy==0&&Oz==0)
		{
		cout<<"O == 0"<<endl;
		return false;
		}
	else if(Ox!=0)
		{
//		cout<<"x"<<endl;
		double M=a.y*c.x-a.x*c.y;
		double N=c.y*d.x-c.x*d.y;

		double M1=a.x*b.y-a.y*b.x;
		double N1=b.x*d.y-b.y*d.x;

		double q = M/Ox;
		double p = N/Ox;

		double q1 = M1/Ox;
		double p1 = N1/Ox;

		double A = a.z+b.z*q*q+c.z*q1*q1;
		double B = 2.0*(-a.z*e.x+b.z*q*(p-e.y)+c.z*q1*(p1-e.z));
		double C = a.z*e.x*e.x+b.z*(p-e.y)*(p-e.y)+c.z*(p1-e.z)*(p1-e.z)-d.z;
		
		double delt = B*B - 4.0*A*C;
		if(delt>-0.0001&&delt<0)
			delt=0;
		if(delt<0)
			{
			if(nwarning_delt<10)
				cout<<"delt = "<<delt<<" at Ox"<<endl;
			nwarning_delt += 1;
			return false;
			}
		else
			{
			F1.x = (-B+sqrt(delt))/(2.0*A);
			F2.x = (-B-sqrt(delt))/(2.0*A);
			F1.y = p+q*F1.x;
			F2.y = p+q*F2.x;	
			F1.z = p1+q1*F1.x;
			F2.z = p1+q1*F2.x;
				return true;
			}
		}
	else if(Oy!=0)
		{
		//cout<<"y"<<endl;		
		double M=b.y*c.x-b.x*c.y;
		double N=c.y*d.x-c.x*d.y;

		double M1=a.y*b.x-a.x*b.y;
		double N1=a.x*d.y-a.y*d.x;

		double q = M/Oy;
		double p = N/Oy;

		double q1 = M1/Oy;
		double p1 = N1/Oy;

		double A = a.z*q*q+b.z+c.z*q1*q1;
		double B = 2.0*(a.z*q*(p-e.x)-b.z*e.y+c.z*q1*(p1-e.z));
		double C = a.z*(p-e.x)*(p-e.x)+b.z*e.y*e.y+c.z*(p1-e.z)*(p1-e.z)-d.z;
		
		double delt = B*B - 4.0*A*C;
		if(delt>-0.0001&&delt<0)
			delt=0;		
		if(delt<0)
			{
			if(nwarning_delt<10)	
				cout<<"delt = "<<delt<<" at Oy"<<endl;
			nwarning_delt += 1;
			return false;
			}
		else
			{
			F1.y = (-B+sqrt(delt))/(2.0*A);
			F2.y = (-B-sqrt(delt))/(2.0*A);
			F1.x = p+q*F1.y;
			F2.x = p+q*F2.y;	
			F1.z = p1+q1*F1.y;
			F2.z = p1+q1*F2.y;
				return true;
			}
		}
	else if(Oz!=0)
		{
		//cout<<"z"<<endl;
		double M=b.x*c.y-b.y*c.x;
		double N=b.y*d.x-b.x*d.y;

		double M1=a.y*c.x-a.x*c.y;
		double N1=a.x*d.y-a.y*d.x;

		double q = M/Oz;
		double p = N/Oz;

		double q1 = M1/Oz;
		double p1 = N1/Oz;

		double A = a.z*q*q+b.z*q1*q1+c.z;
		double B = 2.0*(a.z*q*(p-e.x)+b.z*q1*(p1-e.y)-c.z*e.z);
		double C = a.z*(p-e.x)*(p-e.x)+b.z*(p1-e.y)*(p1-e.y)+c.z*e.z*e.z-d.z;
		
		double delt = B*B - 4.0*A*C;
		if(delt>-0.0001&&delt<0)
			delt=0;
		if(delt<0)
			{
			if(nwarning_delt<10)
				cout<<"delt = "<<delt<<" at Oz"<<endl;
			nwarning_delt += 1;
			return false;
			}
		else
			{
			F1.z = (-B+sqrt(delt))/(2.0*A);
			F2.z = (-B-sqrt(delt))/(2.0*A);
			F1.x = p+q*F1.z;
			F2.x = p+q*F2.z;	
			F1.y = p1+q1*F1.z;
			F2.y = p1+q1*F2.z;
				return true;
			}
		}
	return false;	
	}

void Molecule::genName()
	{
	unsigned int molmark0=0;
	unsigned int molmark1=0;
	vector<unsigned int> NEtype;	
	NEtype.resize(m_Ntypes);
	for(unsigned int i=0; i<m_NatomPerMole; i++)
		{
		unsigned int typi = m_typeId[i];
		unsigned int Nb = m_nbond[i];		
		molmark0 += Nb*typi;
		molmark1 += Nb*(typi+1);	
		NEtype[typi] += 1;
		}
	stringstream s0,s1,s2;
	s0<<m_NatomPerMole;
	s1<<molmark0;
	s2<<molmark1;
		
	string typemark;
	for(unsigned int j=0; j<m_Ntypes; j++)
		{
		unsigned int num = NEtype[j];
		stringstream st;				
		if(num>0)
			{
			st<<num;
			typemark += m_type_mapping[j]+"["+st.str()+"]";
			}
		}
	m_mol_name = s0.str()+"-"+typemark+"-"+s1.str()+"-"+s2.str();
	}
	
void Molecule::generate()
	{
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
		vector<std::string> bond_mapping;
		for(unsigned int i=0; i<m_bond.size(); i++)
			{
			std::string typi = m_bond[i].type;
			bool exist= false;
			for(unsigned int j=0; j<bond_mapping.size(); j++)
				{
				if(typi==bond_mapping[j])
					{
					exist=true;
					break;
					}
				}
			if(!exist)
				bond_mapping.push_back(typi);		
			}
	cout<<"The number of types of bonds: "<<bond_mapping.size()<<endl;
		for(unsigned int j=0; j<bond_mapping.size(); j++)
			cout<<bond_mapping[j]<<endl;		
		generateAngle();
		generateDihedral();	
		cout<<"generating ..."<<endl;		
		m_firststep = false;
		if(m_set_mol_box)
			{
			if(m_mol_Lx==m_Lx&&m_mol_Ly==m_Ly&&m_mol_Lz==m_Lz)
				m_set_mol_box = false;
			}
		else
			{
			m_mol_Lx=m_Lx;
			m_mol_Ly=m_Ly;
			m_mol_Lz=m_Lz;
			}
		}
	m_xyz.clear();
	m_xyz.resize(m_NatomPerMole);

	vector<unsigned int> last_list;
	vector<unsigned int> now_list;
	vector<unsigned int> next_list;	
	vector<unsigned int> last_pre_bond;
	vector<unsigned int> now_pre_bond;	
	vector<unsigned int> next_pre_bond;
	
	for(unsigned int i=0; i<m_NatomPerMole; i++)
		{
		if(m_be_generated_read[i])
			{
			m_be_generated[i]= true;
			m_xyz[i]=m_xyz_read[i];
			}
		else
			{
			m_be_generated[i]= false;
			}
		}
	for(unsigned int idx=0; idx<m_NatomPerMole; idx++)
		{
		if(!m_be_generated[idx])
			{
			vector<vector<unsigned int> > list;
			vector<vector<unsigned int> > pre_bond;
			now_list.clear();
			if(m_bond_init.size()>0)
				now_list=m_bond_init;
			else
				now_list.push_back(idx);
			bool go =true;
			unsigned int times=0;
			while(go)
				{
				next_list.clear();
				next_pre_bond.clear();
				for(unsigned int i =0; i<now_list.size();i++)
					{
					unsigned int taga = now_list[i];
					if (m_be_generated[taga]&&times>0)
						continue;				
					m_be_generated[taga] = true;						
					unsigned int Nneibor = m_nbond[taga];
					for(unsigned int j =0; j<Nneibor; j++)
						{
						unsigned int nextTag = m_bondtable[taga*m_bondtableHight + j];
						bool exist = false;
						if(m_be_generated[nextTag])
							exist =true;
						if(!exist)
							{
							next_list.push_back(nextTag);
							next_pre_bond.push_back(taga);							
							}
						}
					}						
				list.push_back(now_list);
				pre_bond.push_back(now_pre_bond);
			
				now_list = next_list;
				now_pre_bond = next_pre_bond;
				times +=1;
				if(now_list.size()==0)
					go = false;
				}
			unsigned int startfloor=0;
			if(m_bond_init.size()>0)	
				startfloor=1;
			for(unsigned int i=startfloor;i<list.size();i++)
				{
				for(unsigned int j=0;j<list[i].size();j++)
					{
					unsigned int taga = list[i][j];
					m_be_generated[taga] = false;
					}
				}
			unsigned int bound =0;
			for(int floor=0; floor<int(list.size()); floor++)
				{
				if(floor>=1)
					{
					last_list = list[floor-1];
					last_pre_bond = pre_bond[floor-1];
					}
				else
					{
					last_list.clear();
					last_pre_bond.clear();					
					}
				unsigned int Last_N = last_list.size();
				now_list = list[floor];
				now_pre_bond = pre_bond[floor];
				unsigned int N_now = now_list.size();
				for(unsigned int i =0; i<N_now;i++)
					{
					bool success = true;					
					unsigned int taga = now_list[i];
					if (m_be_generated[taga])
						continue;
					unsigned int bondb =  NO_INDEX;
					vector<unsigned int> anglec;
					vector<unsigned int> angleb;
					vector<double> bond_len;					
					vector<double> angle_radian;					
					unsigned int typ_a = m_typeId[taga];
					unsigned int typ_b = NO_INDEX;
//					unsigned int typ_c = NO_INDEX;
					double bond_length =0.0;
					
					if(now_pre_bond.size()!=0)
						{
						bondb = now_pre_bond[i];
						typ_b = m_typeId[bondb];
						bond_length = BondLength[typ_a + typ_b*m_Ntypes];
 						if(bond_length==0.0)
							{
							cout<<"Warning! the bond length have not been set between "<<m_type[taga]<<" and "<<m_type[bondb]<<endl;
							cout<<"The bond length is set 1.0 as default value!"<<endl;
							bond_length = 1.0;
							BondLength[typ_a + typ_b*m_Ntypes]=bond_length;
							BondLength[typ_b + typ_a*m_Ntypes]=bond_length;
							} 
						if(last_pre_bond.size()!=0)
							{						
							for(unsigned int k = 0; k <Last_N; k++)
								{
							   unsigned int lasttag = last_list[k];
							   if(bondb == lasttag)
									{
									unsigned int tagk = last_pre_bond[k];
//									typ_c = m_typeId[tagk];
									double a_r = existedAngleDigree(taga,bondb,tagk);
									if(a_r > 0.0)
										{
										anglec.push_back(tagk);
										angleb.push_back(bondb);
										bond_len.push_back(bond_length);
										angle_radian.push_back(a_r);
										}
									}
								}
							}
						for(unsigned int k =0; k<N_now;k++)
							{								
							unsigned int tagk = now_list[k];
							unsigned int bondk = now_pre_bond[k];
							if(tagk!=taga&&bondk==bondb&&m_be_generated[tagk])
								{								
//								typ_c = m_typeId[tagk];
								double a_r = existedAngleDigree(taga,bondb,tagk);								
								if(a_r > 0.0)
									{
									anglec.push_back(tagk);
									angleb.push_back(bondb);
									bond_len.push_back(bond_length);
									angle_radian.push_back(a_r);
									}
								}
							}
						for(unsigned int k =0; k<Last_N;k++)
							{
							unsigned int lasttag = last_list[k];
							unsigned int Nneibor = m_nbond[bondb];
							bool exist = false;				
							for(unsigned int j =0; j<Nneibor; j++)
								{
								unsigned int nextTag = m_bondtable[bondb*m_bondtableHight + j];
								if(nextTag==lasttag)
									exist=true;
								}
							if(exist)
								{								
//								typ_c = m_typeId[lasttag];
								double a_r = existedAngleDigree(taga,bondb,lasttag);								
								if(a_r > 0.0)
									{
									anglec.push_back(lasttag);
									angleb.push_back(bondb);
									bond_len.push_back(bond_length);
									angle_radian.push_back(a_r);
									}
								}
							}							
						if(last_pre_bond.size()!=0)
							{
							for(unsigned int j =0; j<N_now;j++)
								{
								unsigned int tagj = now_list[j];
								unsigned int bondj = now_pre_bond[j];
								typ_b = m_typeId[bondj];
								double bond_length1 = BondLength[typ_a + typ_b*m_Ntypes];
								if(tagj==taga&&i!=j)
									{
									for(unsigned int k = 0; k <Last_N; k++)
										{
									   unsigned int lasttag = last_list[k];
									   if(bondj == lasttag)
											{
											unsigned int tagk = last_pre_bond[k];
//											typ_c = m_typeId[tagk];
											double a_r = existedAngleDigree(taga,bondb,tagk);
											if(a_r > 0.0)
												{
												anglec.push_back(tagk);
												angleb.push_back(bondj);
												bond_len.push_back(bond_length1);
												angle_radian.push_back(a_r);
												}
											}
										}
									}
								}
							}
							
						bool fitangle = false;
						unsigned int turns =0;
						vec posa;						
						while(!fitangle)
							{
							vector<vec> testpos;
							if(anglec.size()==0)
								{
								vec posb=m_xyz[bondb];
								for(unsigned int i=0; i<m_testnum; i++)
									{
									double drX_tmp = ( R2S() - 0.5 )*m_mol_Lx;
									double drY_tmp = ( R2S() - 0.5 )*m_mol_Ly;
									double drZ_tmp = ( R2S() - 0.5 )*m_mol_Lz;
									double r_dist  = sqrt( drX_tmp*drX_tmp + drY_tmp*drY_tmp + drZ_tmp*drZ_tmp);
									double d_scale = bond_length/r_dist; 
									posa.x = posb.x + drX_tmp * d_scale; 
									posa.y = posb.y + drY_tmp * d_scale; 
									posa.z = posb.z + drZ_tmp * d_scale;
									testpos.push_back(posa);
									}
								vector<unsigned int> vd;
								fitangle = checkdistance(taga, bondb, vd, testpos, posa, 2);
								}
							else if(anglec.size()==1)
								{							
								vec posb = m_xyz[bondb];
								vec posc = m_xyz[anglec[0]];
								
								double drX_tmp = ( R2S() - 0.5 )*m_mol_Lx;
								double drY_tmp = ( R2S() - 0.5 )*m_mol_Ly;
								double drZ_tmp = ( R2S() - 0.5 )*m_mol_Lz;
								double r_dist  = sqrt( drX_tmp*drX_tmp + drY_tmp*drY_tmp + drZ_tmp*drZ_tmp);
								double d_scale = bond_length/r_dist;
								
								double dxab = drX_tmp * d_scale; 
								double dyab = drY_tmp * d_scale; 
								double dzab = drZ_tmp * d_scale;								

								posa.x = posb.x + dxab;
								posa.y = posb.y + dyab;
								posa.z = posb.z + dzab;
								
								double dxcb = posc.x - posb.x;
								double dycb = posc.y - posb.y;
								double dzcb = posc.z - posb.z;		
						
								double rsqab = dxab*dxab+dyab*dyab+dzab*dzab;
								double rab = sqrt(rsqab);
								double rsqcb = dxcb*dxcb+dycb*dycb+dzcb*dzcb;
								double rcb = sqrt(rsqcb);
								double c_abbc = dxab*dxcb+dyab*dycb+dzab*dzcb;
								c_abbc = c_abbc /(rab*rcb);
				
								if (c_abbc > 1.0)
									c_abbc = 1.0;
								if (c_abbc < -1.0)
									c_abbc = -1.0;
					
								double theta = acos(c_abbc);
					
								double fx = dyab*dzcb - dzab*dycb;
								double fy = dzab*dxcb - dxab*dzcb;
								double fz = dxab*dycb - dyab*dxcb;

								double dth = angle_radian[0] - theta;
								
								{
								vec vert = vec(posb.x + fx, posb.y + fy, posb.z + fz);
								vec temp;
								double m[4][4];
								RotateArbitraryLine(m, vert, posb, dth);
								temp.x = m[0][0]*posa.x + m[1][0]*posa.y+m[2][0]*posa.z+m[3][0];
								temp.y = m[0][1]*posa.x + m[1][1]*posa.y+m[2][1]*posa.z+m[3][1];
								temp.z = m[0][2]*posa.x + m[1][2]*posa.y+m[2][2]*posa.z+m[3][2];
								posa = temp;
								}

								if(m_mol_Lx>0.0&&m_mol_Ly>0.0&&m_mol_Lz>0.0)
									{
									for(unsigned int i=0; i<m_testnum; i++)
										{
										vec posai;
										double m[4][4];
										RotateArbitraryLine(m, posc, posb, 2.0*M_PI*double(i)/double(m_testnum));
										posai.x = m[0][0]*posa.x + m[1][0]*posa.y+m[2][0]*posa.z+m[3][0];
										posai.y = m[0][1]*posa.x + m[1][1]*posa.y+m[2][1]*posa.z+m[3][1];
										posai.z = m[0][2]*posa.x + m[1][2]*posa.y+m[2][2]*posa.z+m[3][2];
										testpos.push_back(posai);
										}
									}
								else
									{
									testpos.push_back(posa);									
									}
								fitangle = checkdistance(taga, bondb, anglec, testpos, posa, 2);
								}
							else if(anglec.size()==2)
								{							
								vec A = m_xyz[anglec[0]];
								vec D = m_xyz[anglec[1]];
								vec B = m_xyz[angleb[0]];
								vec C = m_xyz[angleb[1]];
//								cout<<" aa "<<A.x<<" "<<A.y<<" "<<A.z<<" "<<B.x<<" "<<B.y<<" "<<B.z<<" "<<C.x<<" "<<C.y<<" "<<C.z<<" "<<D.x<<" "<<D.y<<" "<<D.z<<endl;
								fitangle=twoAnglesFixE(A, B, C, D, posa, bond_len[0], bond_len[1], angle_radian[0],angle_radian[1]);
//									cout<<" F1 "<<taga<<" "<<m_xyz[taga].x<<" "<<m_xyz[taga].y<<" "<<m_xyz[taga].z<<endl;
								testpos.push_back(posa);
								if (fitangle)
									fitangle = checkdistance(taga, bondb, anglec, testpos, posa, 2);									
								}
							else if(anglec.size()>=3)
								{								
								vec A = m_xyz[anglec[0]];
								vec D = m_xyz[anglec[1]];
								vec F = m_xyz[anglec[2]];
//								cout<<" F "<<F.x<<" "<<F.y<<" "<<F.z<<" "<<anglec[2]<<" "<<m_xyz[anglec[2]].x<<" "<<m_xyz[anglec[2]].y<<" "<<m_xyz[anglec[2]].z<<endl;
								vec B = m_xyz[angleb[0]];
								vec C = m_xyz[angleb[1]];
								vec E = m_xyz[angleb[2]];
								fitangle=threeAnglesFixG(A, B, C, D, E, F, posa, bond_len[0], bond_len[1], bond_len[2],
								angle_radian[0], angle_radian[1], angle_radian[2]);
								testpos.push_back(posa);								
								if (fitangle)
									fitangle = checkdistance(taga, bondb, anglec, testpos, posa, 2);
								}								
							else
								{
								fitangle =true;
								}
							if(turns>100)
								{
								success = false;
								fitangle = true;	
								bound +=1;
								}
							turns +=1;								
							}
						m_xyz[taga].x = posa.x;
						m_xyz[taga].y = posa.y; 
						m_xyz[taga].z = posa.z;	
						
						if(m_mol_Lx==0.0)
							m_xyz[taga].x = 0.0;
						if(m_mol_Ly==0.0)
							m_xyz[taga].y = 0.0;
						if(m_mol_Lz==0.0)
							m_xyz[taga].z = 0.0;

						vec posb=m_xyz[bondb];
						if(m_orientation[taga]==1)
							{
							m_ori_vec[taga].x = posa.x - posb.x;
							m_ori_vec[taga].y = posa.y - posb.y;
							m_ori_vec[taga].z = posa.z - posb.z;
							if(m_mol_Lz==0.0)
								m_ori_vec[taga].z = 0.0;
							}

						if(m_orientation[bondb]==1)
							{
							m_ori_vec[bondb].x = posa.x - posb.x;
							m_ori_vec[bondb].y = posa.y - posb.y;
							m_ori_vec[bondb].z = posa.z - posb.z;
							if(m_mol_Lz==0.0)
								m_ori_vec[taga].z = 0.0;
							}

						if(m_quaternion[taga]==1)
							{	
							if(m_mol_Lz==0.0&&m_eb_spv.x==0.0&&m_eb_spv.y==0.0&&m_eb_spv.z==1.0)
								m_eb_spv=vec(1.0, 0.0, 0.0);							
								
							double vabx = posa.x - posb.x;
							double vaby = posa.y - posb.y;
							double vabz = posa.z - posb.z;	
							
							double fx = m_eb_spv.y*vabz - m_eb_spv.z*vaby;
							double fy = m_eb_spv.z*vabx - m_eb_spv.x*vabz;
							double fz = m_eb_spv.x*vaby - m_eb_spv.y*vabx;
							
							double f = sqrt(fx*fx + fy*fy + fz*fz);
							fx /= f;
							fy /= f;
							fz /= f;
							
							if(m_mol_Lz==0.0&&(fx!=0.0||fy!=0.0))
								{	
								cerr << endl << "***Wrong ellipsoid bond vector!" << endl << endl;
								throw runtime_error("Error Molecule generate quaternion!");
								}	

							double vabsq = vabx*vabx + vaby*vaby + vabz*vabz;
							double mul = (m_eb_spv.x*vabx + m_eb_spv.y*vaby + m_eb_spv.z*vabz)/sqrt(vabsq);
							double theta = acos(mul);

							double q0 = cos(theta/2.0);
							double q1 = sin(theta/2.0)*fx;
							double q2 = sin(theta/2.0)*fy;
							double q3 = sin(theta/2.0)*fz;					

							m_quat_vec[taga] = vec4(q0, q1, q2, q3);
							}

						if(m_quaternion[bondb]==1)
							{
							if(m_mol_Lz==0.0&&m_eb_spv.x==0.0&&m_eb_spv.y==0.0&&m_eb_spv.z==1.0)
								m_eb_spv=vec(1.0, 0.0, 0.0);
							
							double vabx = posa.x - posb.x;
							double vaby = posa.y - posb.y;
							double vabz = posa.z - posb.z;	
							
							double fx = m_eb_spv.y*vabz - m_eb_spv.z*vaby;
							double fy = m_eb_spv.z*vabx - m_eb_spv.x*vabz;
							double fz = m_eb_spv.x*vaby - m_eb_spv.y*vabx;
							
							double f = sqrt(fx*fx + fy*fy + fz*fz);
							fx /= f;
							fy /= f;
							fz /= f;
							
							if(m_mol_Lz==0.0&&(fx!=0.0||fy!=0.0))
								{	
								cerr << endl << "***Wrong ellipsoid bond vector!" << endl << endl;
								throw runtime_error("Error Molecule generate quaternion!");
								}	
						
							double vabsq = vabx*vabx + vaby*vaby + vabz*vabz;
							double mul = (m_eb_spv.x*vabx + m_eb_spv.y*vaby + m_eb_spv.z*vabz)/sqrt(vabsq);
							double theta = acos(mul);

							double q0 = cos(theta/2.0);
							double q1 = sin(theta/2.0)*fx;
							double q2 = sin(theta/2.0)*fy;
							double q3 = sin(theta/2.0)*fz;					

							m_quat_vec[bondb] = vec4(q0, q1, q2, q3);
							}								
						}
					else
						{
						bool fitpoint = false;
						unsigned int turns =0;
						vec posa;
						while(!fitpoint)
							{
							vector<vec> testpos;
							for(unsigned int i=0; i<m_testnum; i++)
								{
								posa.x = (R2S() - 0.5)*m_mol_Lx + m_shift_Lx;
								posa.y = (R2S() - 0.5)*m_mol_Ly + m_shift_Ly;
								posa.z = (R2S() - 0.5)*m_mol_Lz + m_shift_Lz;								
								testpos.push_back(posa);
								}
							vector<unsigned int> vd;							
							fitpoint = checkdistance(taga, taga, vd, testpos, posa, 2);
							if(turns>100)
								{
								success = false;
								fitpoint = true;	
								bound +=1;	
								}
							turns +=1;								
							}
							
						m_xyz[taga].x = posa.x;
						m_xyz[taga].y = posa.y; 
						m_xyz[taga].z = posa.z;

						if(m_orientation[taga]==1)
							{
							m_ori_vec[taga].x = (R2S() - 0.5);
							m_ori_vec[taga].y = (R2S() - 0.5);
							m_ori_vec[taga].z = (R2S() - 0.5);
							if(m_mol_Lz==0.0)
								m_ori_vec[taga].z = 0.0;
							}	
							
						if(m_quaternion[taga]==1)
							{
							double theta = R2S()*1.0*M_PI;
							double phi = R2S()*2.0*M_PI;
							double psi = R2S()*2.0*M_PI;
							if(m_mol_Lz==0.0)
								{
								theta = 0.0;
								psi  = 0.0;
								}
							double q0 = cos(0.5*theta)*cos(0.5*(phi+psi));
							double q1 = sin(0.5*theta)*cos(0.5*(phi-psi));
							double q2 = sin(0.5*theta)*sin(0.5*(phi-psi));
							double q3 = cos(0.5*theta)*sin(0.5*(phi+psi));

							if(m_mol_Lz==0.0)
								{
								q1 = 0.0;
								q2 = 0.0;
								}
								
							m_quat_vec[taga] = vec4(q0, q1, q2, q3);
							}
							
						}

					if(success)
						{
						m_be_generated[taga] = true;
						}
					else
						{
						if(floor>3&&bound>=10)
							{
							for(unsigned int i=0;i<N_now;i++)
								{
								unsigned int taga = now_list[i];
								m_be_generated[taga]=false;	
								}									
							for(unsigned int i=0;i<list[floor-1].size();i++)
								{
								unsigned int taga = list[floor-1][i];
								m_be_generated[taga]=false;	
								}
							for(unsigned int i=0;i<list[floor-2].size();i++)
								{
								unsigned int taga = list[floor-2][i];
								m_be_generated[taga]=false;	
								}
							for(unsigned int i=0;i<list[floor-3].size();i++)
								{
								unsigned int taga = list[floor-3][i];
								m_be_generated[taga]=false;	
								}									
							floor -= 4;	
							if(m_output_times<50)
								{
								cout<<"Go back three particles "<<bound<<" th!"<<endl;
								m_output_times +=1;
								}
							}						
						else if(floor>2&&bound>=5)
							{
							for(unsigned int i=0;i<N_now;i++)
								{
								unsigned int taga = now_list[i];
								m_be_generated[taga]=false;	
								}									
							for(unsigned int i=0;i<list[floor-1].size();i++)
								{
								unsigned int taga = list[floor-1][i];
								m_be_generated[taga]=false;	
								}
							for(unsigned int i=0;i<list[floor-2].size();i++)
								{
								unsigned int taga = list[floor-2][i];
								m_be_generated[taga]=false;	
								}								
							floor -= 3;
							if(m_output_times<50)
								{
								cout<<"Go back two particles "<<bound<<" th!"<<endl;
								m_output_times +=1;
								}								
							}
						else if(floor>1&&bound>=3)
							{
							for(unsigned int i=0;i<N_now;i++)
								{
								unsigned int taga = now_list[i];
								m_be_generated[taga]=false;	
								}									
							for(unsigned int i=0;i<list[floor-1].size();i++)
								{
								unsigned int taga = list[floor-1][i];
								m_be_generated[taga]=false;	
								}
							floor -= 2;
							if(m_output_times<50)
								{
								cout<<"Go back one particle "<<bound<<" th!"<<endl;
								m_output_times +=1;
								}
							}
						else
							{
							for(unsigned int i=0;i<N_now;i++)
								{
								unsigned int taga = now_list[i];
								m_be_generated[taga]=false;	
								}
							floor -=1;
							if(m_output_times<50)
								{							
								cout<<"Re generate the particles "<<bound<<" th!"<<endl;
								m_output_times +=1;
								}
							}
						break;
						}
					}
				if(bound>15)
					{
					m_limit_ge += 1;
					if(m_output_times<50)
						{						
						cout<<"Re generate the molecule "<<m_limit_ge<<" th!"<<endl;
						m_output_times +=1;
						}
					if(m_output_times==50)
						{
						cout<<"It is a hard work!"<<endl;
						cout<<"I will not output the detailed information anymore!"<<endl;
						cout<<"Save energy to generate your molecules!"<<endl;
						m_output_times +=1;
						}
					if (m_limit_ge > 1000)
						{
						cerr << endl << "***Sorry, can not generate the configuration!" << endl << endl;
						throw runtime_error("Error Molecule generate");
						}
					generate();							
					return;								
					}
				}
			}
		}
	m_limit_ge = 0;
	}

void Molecule::generateAngle()
	{
	m_firststep = false;
	if(m_angle.size()==0)
		{	
		for(unsigned int i =0; i<m_NatomPerMole;i++)
			{
			unsigned int Nn = m_nbond[i];
			vector<unsigned int > tempb; 
			unsigned int typ_aa = m_typeId[i];
			for (unsigned int j =0; j< Nn; j++)
				{
				unsigned int tag = m_bondtable[i*m_bondtableHight + j];
				tempb.push_back(tag); 
				}
			if (m_include_itself_in_angle)
				tempb.push_back(i);
			if(tempb.size()>1)
				{
				for(unsigned int ii =0; ii<tempb.size()-1;ii++)
					{
					for(unsigned int jj =ii+1; jj<tempb.size();jj++)
						{
						unsigned int typ_bb = m_typeId[tempb[ii]];
						unsigned int typ_cc = m_typeId[tempb[jj]];
						double a_r = AngleRadian[typ_bb + typ_aa*m_Ntypes + typ_cc*m_Ntypes*m_Ntypes];
	//					cout<<a_r<<endl;
						if (a_r>=0.0)
							{
							std::string anglename;
							if(typ_bb<typ_cc)
								{
								anglename += m_type[tempb[ii]];
								anglename.push_back('-');
								anglename += m_type[i];
								anglename.push_back('-');
								anglename += m_type[tempb[jj]];						  
								}
							else
								{
								anglename += m_type[tempb[jj]];
								anglename.push_back('-');
								anglename += m_type[i];
								anglename.push_back('-');
								anglename += m_type[tempb[ii]];			  
								}
							if(tempb[ii]<tempb[jj])	
								m_angle.push_back(Angle(anglename,tempb[ii],i,tempb[jj]));
							else
								m_angle.push_back(Angle(anglename,tempb[jj],i,tempb[ii]));
							RadianPerAngle.push_back(a_r);
							}
						}								
					}
				}		
			}
		}
	cout<<"The number of angles in a molecule:  "<<m_angle.size()<<endl;
	vector<std::string> angle_mapping;
	for(unsigned int i=0; i<m_angle.size(); i++)
		{
		std::string typi = m_angle[i].type;
		bool exist= false;
		for(unsigned int j=0; j<angle_mapping.size(); j++)
			{
			if(typi==angle_mapping[j])
				{
				exist=true;
				break;
				}
			}
		if(!exist)
			angle_mapping.push_back(typi);		
		}
	cout<<"The number of types of angles: "<<angle_mapping.size()<<endl;
	for(unsigned int j=0; j<angle_mapping.size(); j++)
		cout<<angle_mapping[j]<<endl;
	}

bool Molecule::existedDihedral(unsigned int a, unsigned int b, unsigned int c, unsigned int d)
	{
	for(unsigned int i =0; i< m_dihedral.size(); i++)
		{
		if(m_dihedral[i].a==a&&m_dihedral[i].b==b&&m_dihedral[i].c==c&&m_dihedral[i].d==d)
			return true;
		else if(m_dihedral[i].a==d&&m_dihedral[i].b==c&&m_dihedral[i].c==b&&m_dihedral[i].d==a)
			return true;			
		}
	return false;
	
	}

void Molecule::generateDihedral()
	{
	if(m_dihedral.size()==0)
		{
		for(unsigned int a =0; a<m_NatomPerMole;a++)
			{
			unsigned int Na = m_nbond[a];
			for(unsigned int j=0;j < Na; j++)
				{
				unsigned int b = m_bondtable[a*m_bondtableHight + j];
				unsigned int Nb = m_nbond[b];
				for(unsigned int k=0;k < Nb; k++)
					{
					unsigned int c = m_bondtable[b*m_bondtableHight + k];
					if(c!=a)
						{
						unsigned int Nc = m_nbond[c];
						for(unsigned int l=0;l < Nc; l++)
							{
							unsigned int d = m_bondtable[c*m_bondtableHight + l];
							if(d!=b)
								{
								unsigned int typ_a = m_typeId[a];							
								unsigned int typ_b = m_typeId[b];
								unsigned int typ_c = m_typeId[c];
								unsigned int typ_d = m_typeId[d];
								double a_r = DihedralRadian[typ_a + typ_b*m_Ntypes + typ_c*m_Ntypes*m_Ntypes + typ_d*m_Ntypes*m_Ntypes*m_Ntypes];
								if (a_r > -M_PI && a_r <= M_PI)
									{
									std::string dihedralname;
									if(typ_a<typ_d)
										{
										dihedralname += m_type[a];
										dihedralname.push_back('-');
										dihedralname += m_type[b];
										dihedralname.push_back('-');
										dihedralname += m_type[c];
										dihedralname.push_back('-');
										dihedralname += m_type[d];										
										}
									else
										{
										dihedralname += m_type[d];
										dihedralname.push_back('-');
										dihedralname += m_type[c];
										dihedralname.push_back('-');
										dihedralname += m_type[b];
										dihedralname.push_back('-');
										dihedralname += m_type[a];		  
										}
									if(!existedDihedral(a,b,c,d))
										m_dihedral.push_back(Dihedral(dihedralname,a,b,c,d));
									RadianPerDihedral.push_back(a_r);
									}
								}
							}
						}
					}
				}
			}
		}
	cout<<"The number of dihedrals in a molecule: "<<m_dihedral.size()<<endl;
	vector<std::string> dihedral_mapping;
	for(unsigned int i=0; i<m_dihedral.size(); i++)
		{
		std::string typi = m_dihedral[i].type;
		bool exist= false;
		for(unsigned int j=0; j<dihedral_mapping.size(); j++)
			{
			if(typi==dihedral_mapping[j])
				{
				exist=true;
				break;
				}
			}
		if(!exist)
			dihedral_mapping.push_back(typi);		
		}
	cout<<"The number of types of dihedrals: "<<dihedral_mapping.size()<<endl;
	for(unsigned int j=0; j<dihedral_mapping.size(); j++)
		cout<<dihedral_mapping[j]<<endl;	
	}

void export_Molecule(pybind11::module &m)
	{
	pybind11::class_<Molecule>(m, "Molecule")
        .def(pybind11::init<unsigned int >())
        .def(pybind11::init<const std::string&, unsigned int >())		
		.def("setParticleTypes", &Molecule::setParticleTypes)	
		.def("setTopology", &Molecule::setTopology)
		.def("setIsotactic", &Molecule::setIsotactic)
		.def("setPosition", &Molecule::setPosition)
		.def("setBox", static_cast< void (Molecule::*)(double, double, double) >(&Molecule::setBox)) 
		.def("setBox", static_cast< void (Molecule::*)(double, double, double, double, double, double)>(&Molecule::setBox))		
		.def("setBondLength", static_cast< void (Molecule::*)(double) >(&Molecule::setBondLength)) 
		.def("setBondLength", static_cast< void (Molecule::*)(std::string, std::string, double)>(&Molecule::setBondLength))
		.def("setMass", static_cast< void (Molecule::*)(double) >(&Molecule::setMass)) 
		.def("setMass", static_cast< void (Molecule::*)(std::string, double)>(&Molecule::setMass))
		.def("setMass", static_cast< void (Molecule::*)(unsigned int, double)>(&Molecule::setMass))
		.def("setAngleDegree", static_cast< void (Molecule::*)(std::string, std::string, std::string, double)>(&Molecule::setAngleDegree))
		.def("setAngleDegree", static_cast< void (Molecule::*)(unsigned int, unsigned int, unsigned int, double) >(&Molecule::setAngleDegree))
		.def("setDihedralDegree", static_cast< void (Molecule::*)(std::string, std::string, std::string, std::string, double)>(&Molecule::setDihedralDegree))
		.def("setDihedralDegree", static_cast< void (Molecule::*)(unsigned int, unsigned int, unsigned int, unsigned int, double) >(&Molecule::setDihedralDegree))
		.def("setCharge", static_cast< void (Molecule::*)(double)>(&Molecule::setCharge))
		.def("setCharge", static_cast< void (Molecule::*)(std::string, double)>(&Molecule::setCharge))
		.def("setCharge", static_cast< void (Molecule::*)(unsigned int, double) >(&Molecule::setCharge))
		.def("setInert", static_cast< void (Molecule::*)(double, double, double) >(&Molecule::setInert))
		.def("setInert", static_cast< void (Molecule::*)(std::string, double, double, double)>(&Molecule::setInert))
		.def("setInert", static_cast< void (Molecule::*)(unsigned int, double, double, double)>(&Molecule::setInert))		
		.def("setOrientation", static_cast< void (Molecule::*)() >(&Molecule::setOrientation))
		.def("setOrientation", static_cast< void (Molecule::*)(std::string)>(&Molecule::setOrientation))
		.def("setOrientation", static_cast< void (Molecule::*)(unsigned int)>(&Molecule::setOrientation))
		.def("setQuaternion", static_cast< void (Molecule::*)() >(&Molecule::setQuaternion))
		.def("setQuaternion", static_cast< void (Molecule::*)(std::string)>(&Molecule::setQuaternion))
		.def("setQuaternion", static_cast< void (Molecule::*)(unsigned int)>(&Molecule::setQuaternion))		
		.def("setDiameter", static_cast< void (Molecule::*)(double) >(&Molecule::setDiameter))		
		.def("setDiameter", static_cast< void (Molecule::*)(std::string, double)>(&Molecule::setDiameter))
		.def("setDiameter", static_cast< void (Molecule::*)(unsigned int, double)>(&Molecule::setDiameter))
		.def("setInit", static_cast< void (Molecule::*)(unsigned int) >(&Molecule::setInit))		
		.def("setInit", static_cast< void (Molecule::*)(std::string, unsigned int)>(&Molecule::setInit))
		.def("setInit", static_cast< void (Molecule::*)(unsigned int, unsigned int)>(&Molecule::setInit))
		.def("setCris", static_cast< void (Molecule::*)(unsigned int) >(&Molecule::setCris))		
		.def("setCris", static_cast< void (Molecule::*)(std::string, unsigned int)>(&Molecule::setCris))
		.def("setCris", static_cast< void (Molecule::*)(unsigned int, unsigned int)>(&Molecule::setCris))
		.def("setBody", static_cast< void (Molecule::*)(unsigned int) >(&Molecule::setBody))		
		.def("setBody", static_cast< void (Molecule::*)(std::string, unsigned int)>(&Molecule::setBody))
		.def("setBody", static_cast< void (Molecule::*)(unsigned int, unsigned int)>(&Molecule::setBody))
		.def("setMolecule", static_cast< void (Molecule::*)(unsigned int) >(&Molecule::setMolecule))		
		.def("setMolecule", static_cast< void (Molecule::*)(std::string, unsigned int)>(&Molecule::setMolecule))
		.def("setMolecule", static_cast< void (Molecule::*)(unsigned int, unsigned int)>(&Molecule::setMolecule))
		.def("setChargeFactor", &Molecule::setChargeFactor)	
		.def("setEllipsoidBondSpotVector", &Molecule::setEllipsoidBondSpotVector)		
		.def("setTestNum", &Molecule::setTestNum)
		.def("setSphere", &Molecule::setSphere)
		.def("setCylinder", &Molecule::setCylinder)
		.def("setBodyEvacuation", &Molecule::setBodyEvacuation)
		.def("setIncludeItselfInAngle", &Molecule::setIncludeItselfInAngle)		
		;
	} 
