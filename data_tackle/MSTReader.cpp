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

#include "MSTReader.h"

std::vector<std::string> split(std::string str, std::string pattern)
	{
	istringstream i_stream(str);
	std::vector<std::string> q;
	std::string s;
	while (i_stream>>s) {q.push_back(s);}
    return q;
	}

MSTReader::MSTReader() : XMLBuilder()
    {
	m_mst_read = false;
	m_invariant_data = false;
	m_variant_data = false;
	m_read_indicator = {{"bond", false}, {"angle", false}, {"dihedral", false}, {"vsite", false}, {"box", false}, {"position", false}, {"type", false},
						{"image", false}, {"mass", false}, {"velocity", false}, {"charge", false}, {"body", false}, {"diameter", false}, {"rotangle", false},
						{"force", false}, {"virial", false}, {"molecule", false}, {"init", false}, {"cris", false}, {"orientation", false}, {"quaternion", false}, 
						{"rotation", false}, {"inert", false}, {"asphere", false}, {"patch", false} };
	m_sp = ios::beg;
    m_fname = "XXXXXXXXX";
	m_object_name = "MSTReader";
    }

void MSTReader::reset_params()
	{
	m_num_particles_read = false;
	m_timestep_read = false;
	m_dimension_read = false;
	m_bond_read = false;
	m_angle_read = false;
	m_dihedral_read = false;
	m_vsite_read = false;		
	m_box_read = false;	
	m_position_read = false;
	m_type_read = false;
	m_image_read = false;
	m_mass_read = false;
	m_velocity_read = false;
	m_charge_read = false;
	m_body_read = false;
	m_diameter_read = false;	
	m_rotangle_read = false;
	m_force_read = false;	
	m_virial_read = false;
	m_molecule_read = false;
	m_init_read = false;
	m_cris_read = false;
	m_orientation_read = false;
	m_quaternion_read = false;	
	m_rotation_read = false;
	m_inert_read = false;
	m_asphere_read = false;
	m_patch_read = false;		
	}
	
void MSTReader::clear_data()
	{
	if (!m_read_indicator["bond"])
		m_bonds.clear();
	if (!m_read_indicator["angle"])	
		m_angles.clear();
	if (!m_read_indicator["dihedral"])
		m_dihedrals.clear();
	if (!m_read_indicator["vsite"])	
		m_vsites.clear();	
	if (!m_read_indicator["box"])				
		m_box=BoxSize(0.0, 0.0, 0.0);	
	if (!m_read_indicator["position"])				
		m_pos.clear();
	if (!m_read_indicator["type"])		
		m_type.clear();
	if (!m_read_indicator["image"])		
		m_image.clear();
	if (!m_read_indicator["mass"])			
		m_mass.clear();
	if (!m_read_indicator["velocity"])			
		m_vel.clear();
	if (!m_read_indicator["charge"])
		m_charge.clear();
	if (!m_read_indicator["body"])
		m_body.clear();
	if (!m_read_indicator["diameter"])
		m_diameter.clear();	
	if (!m_read_indicator["rotangle"])		
		m_rotangle.clear();
	if (!m_read_indicator["force"])			
		m_force.clear();
	if (!m_read_indicator["virial"])
		m_virial.clear();
	if (!m_read_indicator["molecule"])		
		m_molecule.clear();
	if (!m_read_indicator["init"])			
		m_init.clear();
	if (!m_read_indicator["cris"])			
		m_cris.clear();
	if (!m_read_indicator["orientation"])
		m_orientation.clear();
	if (!m_read_indicator["quaternion"])
		m_quaternion.clear();
	if (!m_read_indicator["rotation"])
		m_rotation.clear();	
	if (!m_read_indicator["inert"])		
		m_inert.clear();
	if (!m_read_indicator["asphere"])			
		m_asphere.clear();
	if (!m_read_indicator["patch"])
		m_patch.clear();		
	}


bool MSTReader::readDataFromMST(const string &fname)
    {
	m_if_changed_np = false;
	reset_params();
	bool read_frame = false;
	ifstream file;
	file.open(fname.c_str());	

	if (!file.good())
		{
		cerr << endl << "Unable to open file " << fname.c_str() << endl << endl;
		throw runtime_error("Error reading mst file");
		}

	if (m_fname == fname)
		{
		file.seekg(m_sp);
		}
	else
		{
		// cout<<"read '"<<fname.c_str()<<"'"<<endl;
		m_fname = fname;		
		m_mst_read = false;
		m_invariant_data = false;
		m_variant_data = false;	
		m_if_trajectory = false;
		for (std::map<std::string, bool>::iterator it = m_read_indicator.begin(); it != m_read_indicator.end();it++)
			it->second = false;
		
		file.seekg(0, ios::beg);
		clear_data();		
		}
		
	std::string line;

	while(getline(file, line))
		{
		if (line.find("mst_version") != line.npos && line.find("1.0") != line.npos )
			{
			// cout<<"read mst file with version 1.0"<<endl;
			m_mst_read = true;
			continue;
			}			
		// cout<<line<<endl;
		if (m_mst_read)
			{
			if (line.find("mst_end") != line.npos)
				{	
				read_frame = true;
				m_sp=file.tellg();				
				break;
				}	
				
			if (line.find("invariant_data") != line.npos)
				{
				m_invariant_data = true;
				m_variant_data = false;						
				continue;
				}
				
			if (line.find("variant_data") != line.npos)
				{
				m_invariant_data = false;
				m_variant_data = true;	
				continue;
				}
			
			if (line.find("frame_end") != line.npos)
				{
				// file.seekg(-line.size()-1, ios::cur);
				read_frame = true;
				m_sp=file.tellg();	
				break;
				}
			else if (line.find("frame") != line.npos)
				{
				// check node 
					{
					istringstream parser;
					parser.str(line);
					std::string name;
					int frame_id = -1;
					parser >> name;
					parser >> frame_id;
					
					if (name != "frame"||frame_id==-1)
						throw runtime_error("Error! mst file with wrong format for the indicator of 'frame'");
					
					m_if_trajectory = true;
					}
					
				if (m_variant_data)
					clear_data();
				else
					throw runtime_error("Error! mst files with multiple frames without the label of 'variant_data'");
				continue;
				}				

			if (line.find("num_particles") != line.npos)
				{
				reset_params();
				m_num_particles_read=true;
				continue;
				}
				
			if (line.find("timestep") != line.npos)
				{
				reset_params();
				m_timestep_read=true;
				continue;
				}

			if (line.find("dimension") != line.npos)
				{
				reset_params();
				m_dimension_read=true;
				continue;
				}	

			if (line.find("bond") != line.npos)
				{
				reset_params();
				m_bond_read=true;	
				if (m_invariant_data)
					m_read_indicator["bond"] = true;
				continue;
				}

			if (line.find("angle") != line.npos)
				{				
				reset_params();
				m_angle_read=true;
				if (m_invariant_data)
					m_read_indicator["angle"] = true;
				continue;
				}
				
			if (line.find("dihedral") != line.npos)
				{				
				reset_params();
				m_dihedral_read=true;
				if (m_invariant_data)
					m_read_indicator["dihedral"] = true;
				continue;
				}

			if (line.find("vsite") != line.npos)
				{				
				reset_params();
				m_vsite_read=true;
				if (m_invariant_data)
					m_read_indicator["vsite"] = true;
				continue;
				}				

			if (line.find("box") != line.npos)
				{
				reset_params();
				m_box_read=true;
				if (m_invariant_data)
					m_read_indicator["box"] = true;
				continue;
				}

			if (line.find("position") != line.npos)
				{		
				reset_params();			
				m_position_read=true;
				if (m_invariant_data)
					m_read_indicator["position"] = true;					
				continue;
				}

			if (line.find("type") != line.npos)
				{
				reset_params();			
				m_type_read=true;
				if (m_invariant_data)
					m_read_indicator["type"] = true;					
				continue;
				}

			if (line.find("image") != line.npos)
				{				
				reset_params();
				m_image_read=true;
				if (m_invariant_data)
					m_read_indicator["image"] = true;
				continue;
				}

			if (line.find("mass") != line.npos)
				{
				reset_params();
				m_mass_read=true;
				if (m_invariant_data)
					m_read_indicator["mass"] = true;						
				continue;
				}
				
			if (line.find("velocity") != line.npos)
				{
				reset_params();
				m_velocity_read=true;
				if (m_invariant_data)
					m_read_indicator["velocity"] = true;
				continue;
				}
				
			if (line.find("charge") != line.npos)
				{
				reset_params();
				m_charge_read=true;
				if (m_invariant_data)
					m_read_indicator["charge"] = true;
				continue;
				}

			if (line.find("body") != line.npos)
				{
				reset_params();
				m_body_read=true;
				if (m_invariant_data)
					m_read_indicator["body"] = true;
				continue;
				}

			if (line.find("diameter") != line.npos)
				{
				reset_params();
				m_diameter_read=true;
				if (m_invariant_data)
					m_read_indicator["diameter"] = true;
				continue;
				}				

			if (line.find("rotangle") != line.npos)
				{
				reset_params();
				m_rotangle_read=true;
				if (m_invariant_data)
					m_read_indicator["rotangle"] = true;
				continue;
				}
				
			if (line.find("force") != line.npos)
				{
				reset_params();
				m_force_read=true;
				if (m_invariant_data)
					m_read_indicator["force"] = true;
				continue;
				}

			if (line.find("virial") != line.npos)
				{
				reset_params();
				m_virial_read=true;
				if (m_invariant_data)
					m_read_indicator["virial"] = true;
				continue;
				}

			if (line.find("molecule") != line.npos)
				{
				reset_params();
				m_molecule_read=true;
				if (m_invariant_data)
					m_read_indicator["molecule"] = true;
				continue;
				}

			if (line.find("init") != line.npos)
				{
				reset_params();
				m_init_read=true;
				if (m_invariant_data)
					m_read_indicator["init"] = true;
				continue;
				}	

			if (line.find("cris") != line.npos)
				{
				reset_params();
				m_cris_read=true;
				if (m_invariant_data)
					m_read_indicator["cris"] = true;
				continue;
				}	

			if (line.find("orientation") != line.npos)
				{
				reset_params();
				m_orientation_read=true;
				if (m_invariant_data)
					m_read_indicator["orientation"] = true;
				continue;
				}	

			if (line.find("quaternion") != line.npos)
				{
				reset_params();
				m_quaternion_read=true;
				if (m_invariant_data)
					m_read_indicator["quaternion"] = true;
				continue;
				}	

			if (line.find("rotation") != line.npos)
				{
				reset_params();
				m_rotation_read=true;
				if (m_invariant_data)
					m_read_indicator["rotation"] = true;
				continue;
				}	

			if (line.find("inert") != line.npos)
				{
				reset_params();
				m_inert_read=true;
				if (m_invariant_data)
					m_read_indicator["inert"] = true;
				continue;
				}	

			if (line.find("asphere") != line.npos)
				{
				reset_params();
				m_asphere_read=true;
				if (m_invariant_data)
					m_read_indicator["asphere"] = true;
				continue;
				}	

			if (line.find("patch") != line.npos)
				{
				reset_params();
				m_patch_read=true;
				if (m_invariant_data)
					m_read_indicator["patch"] = true;
				continue;
				}					

			// read data				
			std::vector<std::string> line_array = split(line, "	");
			istringstream parser;
			parser.str(line);

			if (m_num_particles_read && line_array.size() >= 1)
				{
				unsigned int np;
				parser >> np;
				
				if(np==0)
					{	
					cerr << endl
					<< "***Error! number of particles is zero"
					<< endl << endl;
					throw runtime_error("Error extracting data from MST file");
					}
				
				if (np!=m_num_particles)
					{
					m_if_changed_np =  true;
					m_last_np = m_num_particles;
					}
				m_num_particles = np;

				}
				
			if (m_timestep_read && line_array.size() == 1)
				parser >> m_timestep;
				
			if (m_dimension_read && line_array.size() == 1)
				parser >> m_ndimension;		
				
			if (m_bond_read && line_array.size() == 3)
				{
				string name;
				unsigned int a;
				unsigned int b;
				parser>> name >> a >> b;
				m_bonds.push_back(Bond(name, a, b, getBondTypeId(name)));
				}
		
			if (m_angle_read && line_array.size() == 4)
				{
				string name;
				unsigned int a;
				unsigned int b;
				unsigned int c;
				parser>> name >> a >> b >> c;
				m_angles.push_back(Angle(name, a, b, c, getAngleTypeId(name)));
				}

			if (m_dihedral_read && line_array.size() == 5)
				{
				string name;
				unsigned int a;
				unsigned int b;
				unsigned int c;
				unsigned int d;				
				parser>> name >> a >> b >> c >> d;
				m_dihedrals.push_back(Dihedral(name, a, b, c, d, getDihedralTypeId(name)));
				}
		
			if (m_vsite_read && line_array.size() == 5)
				{
				string name;
				unsigned int a;
				unsigned int b;
				unsigned int c;
				unsigned int d;				
				parser>> name >> a >> b >> c >> d;
				m_vsites.push_back(Dihedral(name, a, b, c, d, getVsiteTypeId(name)));
				}					
				
			if (m_box_read && line_array.size() == 3)
				{
				double lx, ly, lz;
				parser>> lx >> ly >> lz;
				m_box = BoxSize(lx, ly, lz);
				}
				
			if (m_position_read && line_array.size() == 3)
				{
				double px, py, pz;				
				parser>> px >> py >> pz;
				// cout<<px<<" "<<py<<" "<<pz<<endl;				
				m_pos.push_back(vec(px, py, pz));
				}

			if (m_type_read && line_array.size() == 1)
				{
				// cout<<line_array[0]<<endl;
				string type;
				parser>> type;
				m_type.push_back(getTypeId(type));
				}
				
			if (m_image_read && line_array.size() == 3)
				{
				int ix, iy, iz;
				parser>> ix >> iy >> iz;
				m_image.push_back(vec_int(ix, iy, iz));
				}
				
			if (m_mass_read && line_array.size() == 1)
				{
				double mass;
				parser>> mass;
				m_mass.push_back(mass);
				}
				
			if (m_velocity_read && line_array.size() == 3)
				{
				double vx, vy, vz;
				parser>> vx >> vy >> vz;
				m_vel.push_back(vec(vx, vy, vz));
				}
				
			if (m_charge_read && line_array.size() == 1)
				{
				double c;
				parser>> c;
				m_charge.push_back(c);
				}

			if (m_body_read && line_array.size() == 1)
				{
				int b;
				parser>> b;
				if (b == -1)
					m_body.push_back(NO_BODY);
				else
					m_body.push_back(b);
				}

			if (m_diameter_read && line_array.size() == 1)
				{
				double d;
				parser>> d;
				m_diameter.push_back(d);
				}				
				
			if (m_rotangle_read && line_array.size() == 3)
				{
				double rx, ry, rz;
				parser>> rx >> ry >> rz;
				m_rotangle.push_back(vec(rx, ry, rz));
				}	
				
			if (m_force_read && line_array.size() == 3)
				{
				double fx, fy, fz;
				parser>> fx >> fy >> fz;
				m_force.push_back(vec(fx, fy, fz));
				}	

			if (m_virial_read && line_array.size() == 1)
				{
				double v;
				parser>> v;
				m_virial.push_back(v);
				}

			if (m_molecule_read && line_array.size() == 1)
				{
				int mol;
				parser>> mol;
				if (mol == -1)
					m_molecule.push_back(NO_INDEX);
				else
					m_molecule.push_back(mol);
				}

			if (m_init_read && line_array.size() == 1)
				{
				unsigned int init;
				parser>> init;
				m_init.push_back(init);
				}

			if (m_cris_read && line_array.size() == 1)
				{
				unsigned int cris;
				parser>> cris;
				m_cris.push_back(cris);
				}

			if (m_orientation_read && line_array.size() == 3)
				{
				double x, y, z;
				parser>> x >> y >> z;
				m_orientation.push_back(vec(x, y, z));
				}	

			if (m_quaternion_read && line_array.size() == 4)
				{
				double x, y, z, w;
				parser>> x >> y >> z >> w;
				m_quaternion.push_back(vec4(x, y, z, w));
				}
				
			if (m_rotation_read && line_array.size() == 3)
				{
				double x, y, z;
				parser>> x >> y >> z;
				m_rotation.push_back(vec(x, y, z));
				}

			if (m_inert_read && line_array.size() == 3)
				{
				double x, y, z;
				parser>> x >> y >> z;
				m_inert.push_back(vec(x, y, z));
				}

			if (m_asphere_read && line_array.size() == 7)
				{
				string na;
				double x, y, z, w, m, n;
				parser>> na >> x >> y >> z >> w >> m >> n;
				m_asphere.push_back(str_vec6(na, x, y, z, w, m, n));
				}

			if (m_patch_read && line_array.size() == 5)
				{
				string na;
				double x, y, z, w;
				parser>> na >> x >> y >> z >> w;
				m_patch.push_back(str_vec6(na, x, y, z, w, 0.0, 0.0));
				}				
			}
		}

	if(!m_mst_read)
		{
        cerr << endl
             << "***Error! This is not a MST file"
             << endl << endl;
        throw runtime_error("Error extracting data from MST file");			
		}		
		
    if (m_box.lx==0.0&&m_box.ly==0.0&&m_box.lz==0.0)
        {
        cerr << endl
             << "***Error! A box is required to define simulation box"
             << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }
		

	if(m_ndimension==2&&m_box.lz>0.0)
		{
        cerr << "***Error! two dimensions of the simulation box should be with Lz = 0.0, the Lz = " << m_box.lz <<" in MST files"<< endl << endl;
        throw runtime_error("Error extracting data from MST file");
		}
	
	if(m_ndimension==3&&m_box.lz<1.0e-6)
		{
        cerr << "***Error! Lz = 0.0 should be with two dimensions of the simulation, three dimensions defind in MST files "<< endl << endl;
        throw runtime_error("Error extracting data from MST file");
		}
	
    if (m_pos.size() == 0)
        {
        cerr << endl << "***Error! No particles defined in position node" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }
		
    if (m_type.size() == 0)
        {
        cerr << endl << "***Error! No particles defined in type node" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }

    if (m_pos.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_pos.size() << " positions != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }
		
    if (m_type.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_type.size() << " types != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }	

    if (m_image.size() != 0 && m_image.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_image.size() << " images != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }	

    if (m_mass.size() != 0 && m_mass.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_mass.size() << " masses != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }		
		
    if (m_vel.size() != 0 && m_vel.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_vel.size() << " velocities != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }		
		
    if (m_molecule.size() != 0 && m_molecule.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_molecule.size() << " molecule != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }

    if (m_charge.size() != 0 && m_charge.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_charge.size() << " charge values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }

    if (m_body.size() != 0 && m_body.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_body.size() << " body values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }			

    if (m_diameter.size() != 0 && m_diameter.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_diameter.size() << " diameters != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }		

    if (m_rotangle.size() != 0 && m_rotangle.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_rotangle.size() << " rotangle values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }

    if (m_force.size() != 0 && m_force.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_force.size() << " force values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }	

    if (m_virial.size() != 0 && m_virial.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_virial.size() << " virial values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }

    if (m_init.size() != 0 && m_init.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_init.size() << " init values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }

    if (m_cris.size() != 0 && m_cris.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_cris.size() << " cris values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }

    if (m_orientation.size() != 0 && m_orientation.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_orientation.size() << " orientation values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }

    if (m_quaternion.size() != 0 && m_quaternion.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_orientation.size() << " quaternion values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }

    if (m_rotation.size() != 0 && m_rotation.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_rotation.size() << " rotation values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }

    if (m_inert.size() != 0 && m_inert.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_inert.size() << " inert values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from MST file");
        }		
		
	//-check bonds, angles and dihedrals
	for (unsigned int i=0; i<m_bonds.size();i++)
		{
		Bond bi = m_bonds[i];
		if(bi.a>=m_num_particles||bi.b>=m_num_particles)
			{
			cerr << endl << "***Error! bond '" << bi.type <<" "<<bi.a<<" "<<bi.b<< "' with the particle index not in the range [0, N-1], with N = "<< m_num_particles<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from MST file");				
			}
		}
		
	for (unsigned int i=0; i<m_angles.size();i++)
		{
		Angle ai = m_angles[i];
		if(ai.a>=m_num_particles||ai.b>=m_num_particles||ai.c>=m_num_particles)
			{
			cerr << endl << "***Error! angle '" << ai.type <<" "<<ai.a<<" "<<ai.b<<" "<<ai.c<< "' with the particle index not in the range [0, N-1], with N = "<< m_num_particles<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from MST file");				
			}
		}	
		
	for (unsigned int i=0; i<m_dihedrals.size();i++)
		{
		Dihedral di = m_dihedrals[i];
		if(di.a>=m_num_particles||di.b>=m_num_particles||di.c>=m_num_particles||di.d>=m_num_particles)
			{
			cerr << endl << "***Error! dihedral '" << di.type <<" "<<di.a<<" "<<di.b<<" "<<di.c<<" "<<di.d<< "' with the particle index not in the range [0, N-1], with N = "<< m_num_particles<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from MST file");				
			}
		}

	for (unsigned int i=0; i<m_vsites.size();i++)
		{
		Dihedral di = m_vsites[i];
		if(di.a>=m_num_particles||di.b>=m_num_particles||di.c>=m_num_particles||di.d>=m_num_particles)
			{
			cerr << endl << "***Error! vsite '" << di.type <<" "<<di.a<<" "<<di.b<<" "<<di.c<<" "<<di.d<< "' with the particle index not in the range [0, N-1], with N = "<< m_num_particles<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from MST file");				
			}
		}
		
	if (file.eof())
		return false;
	else
		return read_frame;
    }
	
void MSTReader::outPutInfo()
	{
    // notify the user of what we have accomplished
    cout <<"--- MST file read summary" << endl;
    cout <<" "<< getNParticles() << " particles at timestep " << m_timestep << endl;
    cout <<" "<< getNParticleTypes() <<  " particle types" << endl;	

    if (m_image.size() > 0)
        cout <<" "<< m_image.size() << " images" << endl;
    if (m_vel.size() > 0)
        cout <<" "<< m_vel.size() << " velocities" << endl;
    if (m_mass.size() > 0)
        cout <<" "<< m_mass.size() << " masses" << endl;
    if (m_charge.size() > 0)
        cout <<" "<< m_charge.size() << " charges" << endl;
    if (m_body.size() > 0)
        cout <<" "<< m_body.size() << " particle body values" << endl; 	
    if (m_diameter.size() > 0)
        cout <<" "<< m_diameter.size() << " diameters" << endl;	
    if (m_rotangle.size() > 0)
        cout <<" "<< m_rotangle.size() << " rotangles" << endl;		
    if (m_force.size() > 0)
        cout <<" "<< m_force.size() << " forces" << endl;	
    if (m_virial.size() > 0)
        cout <<" "<< m_virial.size() << " virials" << endl;	
    if (m_molecule.size() > 0)
        cout <<" "<< m_molecule.size() << " molecules" << endl;
    if (m_init.size() > 0)
        cout <<" "<< m_init.size() << " inits" << endl;
    if (m_cris.size() > 0)
        cout <<" "<< m_cris.size() << " crises" << endl;
    if (m_quaternion.size() > 0)
        cout <<" "<< m_quaternion.size() << " quaternions" << endl;
    if (m_orientation.size() > 0)
        cout <<" "<< m_orientation.size() << " orientations" << endl;	
    if (m_rotation.size() > 0)
        cout <<" "<< m_rotation.size() << " rotations" << endl;	
    if (m_inert.size() > 0)
        cout <<" "<< m_inert.size() << " inerts" << endl;	
    if (m_asphere.size() > 0)
        cout <<" "<< m_asphere.size() << " aspheres" << endl;	
    if (m_patch.size() > 0)
        cout <<" "<< m_patch.size() << " patches" << endl;		
    if (m_bonds.size() > 0)
        cout <<" "<< m_bonds.size() << " bonds" << endl;
    if (m_angles.size() > 0)
        cout <<" "<< m_angles.size() << " angles" << endl;
    if (m_dihedrals.size() > 0)
        cout <<" "<< m_dihedrals.size() << " dihedrals" << endl;
    if (m_vsites.size() > 0)
        cout <<" "<< m_vsites.size() << " vsites" << endl;		
	}


