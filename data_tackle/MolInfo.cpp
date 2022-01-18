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

#include "MolInfo.h"

MolInfo::MolInfo(mst_reader& build): m_build(build)
	{
	}
	
void MolInfo::initialize()
	{
	m_mol_id_particle_device = true;
	m_bond_ex = false;
	m_angle_ex = false;
	m_dihedral_ex = false;
	m_build_ex_host_table = false;
	m_build_ex_device_table = true;
	m_N = m_build.getNParticles();
	m_dimension = m_build.getNDimensions();
    std::vector<Bond> bonds = m_build.getBond();
    std::vector<Angle> angles = m_build.getAngle();
    std::vector<Dihedral> dihedrals = m_build.getDihedral();
	std::vector<unsigned int> n_bond;	
	n_bond.resize(m_N);
	std::vector<unsigned int> n_angle;	
	n_angle.resize(m_N);
	std::vector<unsigned int> n_dihedral;	
	n_dihedral.resize(m_N);
	for (unsigned int i=0; i< bonds.size(); i++)
		{
		Bond bond = bonds[i];
		n_bond[bond.a] +=1;
		n_bond[bond.b] +=1;
		}
		
	for (unsigned int i=0; i< angles.size(); i++)
		{
		Angle angle = angles[i];
		n_angle[angle.a] +=1;
		n_angle[angle.b] +=1;
		n_angle[angle.c] +=1;
		}		
		
	for (unsigned int i=0; i< dihedrals.size(); i++)
		{
		Dihedral dihedral = dihedrals[i];
		n_dihedral[dihedral.a] +=1;
		n_dihedral[dihedral.b] +=1;
		n_dihedral[dihedral.c] +=1;
		n_dihedral[dihedral.d] +=1;
		}
	m_n_max_bond = 0;
	m_n_max_angle = 0;
	m_n_max_dihedral = 0;
	m_n_max_ex = 0;
	for (unsigned int i=0; i< m_N; i++)
		{
		if(n_bond[i]>m_n_max_bond)
			m_n_max_bond = n_bond[i];
		}
	for (unsigned int i=0; i< m_N; i++)
		{
		if(n_angle[i]>m_n_max_angle)
			m_n_max_angle = n_angle[i];
		}
	for (unsigned int i=0; i< m_N; i++)
		{
		if(n_dihedral[i]>m_n_max_dihedral)
			m_n_max_dihedral = n_dihedral[i];
		}	
	
	m_bonds.resize(m_N*m_n_max_bond);
	m_n_bond.resize(m_N);
	m_n_exclusion.resize(m_N);
	
	for(unsigned int i =0; i< m_N; i++)
		{
		m_n_bond[i] = 0;
		m_n_exclusion[i] = 0;
		}

	for (unsigned int i=0; i< bonds.size(); i++)
		{
		Bond bond = bonds[i];
		unsigned int Na = m_n_bond[bond.a];
		unsigned int Nb = m_n_bond[bond.b];		
		m_bonds[Na*m_N+ bond.a] = uint_2(bond.b, bond.id);
		m_bonds[Nb*m_N+ bond.b] = uint_2(bond.a, bond.id);		
		m_n_bond[bond.a] +=1;
		m_n_bond[bond.b] +=1;
		}

	buildmap();
	buildMol();	
    }

void MolInfo::bondExclude()
	{
	m_bond_ex = true;
	m_build_ex_host_table = true;
	}
	
void MolInfo::angleExclude()
	{
	m_angle_ex = true;
	m_build_ex_host_table = true;
	}
	
void MolInfo::dihedralExclude()
	{
	m_dihedral_ex = true;
	m_build_ex_host_table = true;
	}	
	
void MolInfo::buildExclusionList()
	{
	if(m_build_ex_host_table)
		{
		if(m_bond_ex)
			m_n_max_ex += m_n_max_bond;
		if(m_angle_ex)
			m_n_max_ex += m_n_max_angle;
		if(m_dihedral_ex)
			m_n_max_ex += m_n_max_dihedral;
		m_exclusion_list.resize(m_N*m_n_max_ex);
		if(m_bond_ex)
			{
			std::vector<Bond> bonds = m_build.getBond();
			for (unsigned int i=0; i< bonds.size(); i++)
				{
				Bond bond = bonds[i];
				unsigned int ena = m_n_exclusion[bond.a]; 
				unsigned int enb = m_n_exclusion[bond.b]; 		
				m_exclusion_list[bond.a+ena*m_N]=bond.b;
				m_exclusion_list[bond.b+enb*m_N]=bond.a;
				m_n_exclusion[bond.a] += 1;
				m_n_exclusion[bond.b] += 1;
				}
			}
		if(m_angle_ex)
			{
			std::vector<Angle> angles = m_build.getAngle();
			for (unsigned int i=0; i< angles.size(); i++)
				{
				Angle angle = angles[i];
				unsigned int ena = m_n_exclusion[angle.a];
				unsigned int enc = m_n_exclusion[angle.c]; 				
				m_exclusion_list[angle.a+ena*m_N]=angle.c;
				m_exclusion_list[angle.c+enc*m_N]=angle.a;			
				m_n_exclusion[angle.a] += 1;
				m_n_exclusion[angle.c] += 1;
				}
			}
		if(m_dihedral_ex)
			{
			std::vector<Dihedral> dihedrals = m_build.getDihedral();
			for (unsigned int i=0; i< dihedrals.size(); i++)
				{
				Dihedral dihedral = dihedrals[i];
				unsigned int ena = m_n_exclusion[dihedral.a];
				unsigned int end = m_n_exclusion[dihedral.d]; 				
				m_exclusion_list[dihedral.a+ena*m_N]=dihedral.d;
				m_exclusion_list[dihedral.d+end*m_N]=dihedral.a;			
				m_n_exclusion[dihedral.a] += 1;
				m_n_exclusion[dihedral.d] += 1;
				}
			}
		m_build_ex_host_table = false;
		if(m_n_max_ex==0)
			{
			cerr << endl << "***Error, No exclusion information are given in read file!" << endl << endl;	
			throw runtime_error("Error MolInfo::buildExclusionList");
			}
		}
	}

void MolInfo::buildExclusionGPUList()
	{
	buildExclusionList();
	if(m_build_ex_device_table)
		{
		cudaMalloc(&d_n_exclusion, sizeof(unsigned int)*m_N);
		cudaMalloc(&d_exclusion_list, sizeof(unsigned int)*m_N*m_n_max_ex);
		cudaMemcpy(d_n_exclusion, &m_n_exclusion[0], sizeof(unsigned int)*m_N, cudaMemcpyHostToDevice);
		cudaMemcpy(d_exclusion_list, &m_exclusion_list[0], sizeof(unsigned int)*m_N*m_n_max_ex, cudaMemcpyHostToDevice);
		m_build_ex_device_table = false;
		}
	}
	
void MolInfo::buildmap()
	{
	m_map.clear();
	m_map.push_back(vec_int(1,0,0));
	m_map.push_back(vec_int(1,1,0));
	m_map.push_back(vec_int(0,1,0));
	m_map.push_back(vec_int(-1,1,0));
	m_map.push_back(vec_int(1,0,-1));
	m_map.push_back(vec_int(1,1,-1));
	m_map.push_back(vec_int(0,1,-1));
	m_map.push_back(vec_int(-1,1,-1));
	m_map.push_back(vec_int(1,0,1));
	m_map.push_back(vec_int(1,1,1));
	m_map.push_back(vec_int(0,1,1));
	m_map.push_back(vec_int(-1,1,1));
	m_map.push_back(vec_int(0,0,1));	
	}
	
void MolInfo::buildMol()
	{
	std::vector<unsigned int> h_type = m_build.getType();	
	std::vector< std::string> type_map  = m_build.getTypeMap();
	unsigned int Ntype = type_map.size();
	m_box = m_build.getBox();
	m_pos = m_build.getPos();
	std::vector<unsigned int> h_body = m_build.getBody();		
	m_pos0.resize(m_N);
	
	float Lx = m_box.lx;
	float Ly = m_box.ly;
	float Lz = m_box.lz;

	float Lxinv = 0.0;      
	float Lyinv = 0.0; 
	float Lzinv = 0.0;

	if(Lx > 0.0)
		Lxinv = 1.0/Lx;
	
	if(Ly > 0.0)
		Lyinv = 1.0/Ly;
	
	if(Lz > 0.0)
		Lzinv = 1.0/Lz;
	
	m_mol_id_particle.resize(m_N);
	m_mol_stat_end.clear();
	std::vector<unsigned int> molecule = m_build.getMolecule();

	if(molecule.size()==m_N)
		{
		std::vector<std::vector<unsigned int> > mollist;
		mollist.resize(m_N);
		int mol_max=-1;
		for(unsigned int idx=0; idx<m_N; idx++)
			{
			unsigned int mid = molecule[idx];
			if(mid!=NO_INDEX)
				{
				if(mid>=m_N)
					{
					cerr << endl << "***Error, the mol id "<<mid<< " greater than the particle number " << m_N<< endl << endl;	
					throw runtime_error("Error MolInfo::buildMol!");						
					}
				if((int)mid>mol_max)
					mol_max=(int)mid;
				mollist[mid].push_back(idx);
				}
			m_mol_id_particle[idx] = mid;
			}
		if(mol_max>=0)
			m_n_mol = (unsigned int)mol_max+1;
		else
			m_n_mol = 0;
		for(unsigned int mid=0; mid<m_n_mol; mid++)
			{
			unsigned int npmol = mollist[mid].size();
			if(npmol==0)
				{
				cerr << endl << "***Error, the mol id "<<mid<< " has no particles !"<< endl << endl;	
				throw runtime_error("Error MolInfo::buildMol!");				
				}
			
			unsigned int tag0 = mollist[mid][0];
			vec posi = m_pos[tag0];
			m_pos0[tag0] = posi;			
			for(unsigned int idx=1; idx<npmol; idx++)
				{
				unsigned int nextTag = mollist[mid][idx];
				float dx = m_pos[nextTag].x - posi.x;
				float dy = m_pos[nextTag].y - posi.y;
				float dz = m_pos[nextTag].z - posi.z;
				dx -= Lx*rint(dx*Lxinv);
				dy -= Ly*rint(dy*Lyinv);
				dz -= Lz*rint(dz*Lzinv);
				m_pos0[nextTag].x = posi.x + dx;
				m_pos0[nextTag].y = posi.y + dy;
				m_pos0[nextTag].z = posi.z + dz;
				posi = m_pos0[nextTag];
				}

			unsigned int p_start = mollist[mid][0];
			unsigned int p_end = mollist[mid][npmol-1];		
			m_mol_stat_end.push_back(uint_2(p_start, p_end));				
			}				
		
		}
	else
		{
		vector<unsigned int> templist;
		vector<unsigned int> next_templist;	
		m_n_mol = 0;
		unsigned int p_start;
		unsigned int p_end;
		for(unsigned int i=0; i<m_N; i++)
			m_mol_id_particle[i]=NO_INDEX;			
			
		for(unsigned int idx=0; idx<m_N; idx++)
			{
			unsigned int N_b = m_n_bond[idx];
			if(m_mol_id_particle[idx]==NO_INDEX&&N_b>0)
				{
				m_mol_id_particle[idx] = m_n_mol;
				p_start = idx;
				p_end = idx;
				templist.clear();
				templist.push_back(idx);
				m_pos0[idx] = m_pos[idx];
				unsigned int N_temp = templist.size();
	
	
				bool go =true;
				unsigned int bound = 0;
	
				while(go)
					{
					bound += 1;
					next_templist.clear();
					for(unsigned int i =0; i<N_temp;i++)
						{
						unsigned int tag = templist[i];
						vec posi = m_pos0[tag];
						unsigned int Nneibor = m_n_bond[tag];
						for(unsigned int j =0; j<Nneibor; j++)
							{
							unsigned int nextTag =(unsigned int) m_bonds[j*m_N + tag].x;
							bool exist = false;
							if(m_mol_id_particle[nextTag]!=NO_INDEX)
								exist =true;
							if(!exist)
								{
								next_templist.push_back(nextTag);
								m_mol_id_particle[nextTag] = m_n_mol;
								float dx = m_pos[nextTag].x - posi.x;
								float dy = m_pos[nextTag].y - posi.y;
								float dz = m_pos[nextTag].z - posi.z;
								dx -= Lx*rint(dx*Lxinv);
								dy -= Ly*rint(dy*Lyinv);
								dz -= Lz*rint(dz*Lzinv);
								m_pos0[nextTag].x = posi.x + dx;
								m_pos0[nextTag].y = posi.y + dy;
								m_pos0[nextTag].z = posi.z + dz;
								}
							}
						}
					if(next_templist.size()==0)
						p_end = templist[templist.size()-1];	
					templist = next_templist;
					N_temp = templist.size();
					if(N_temp ==0)
						{
						go = false;
						}
					if(bound >= 100000)
						{
						go = false;
						cout<<"Warning! inital molecule error!"<<endl;
						}	
					}
				m_n_mol += 1;
				m_mol_stat_end.push_back(uint_2(p_start,p_end));
				}
			}
	
		if(h_body.size()==m_N)
			{
			int maxbody = -1;
			for(unsigned int idx=0; idx<m_N; idx++)
				{
				unsigned int body = h_body[idx];
				if(body!=NO_INDEX)
					{
					m_mol_id_particle[idx] = m_n_mol+body;
					if ((int) body>maxbody)
						maxbody=(int) body;
					}
				}
			
			if(maxbody>=0)
				{
				maxbody += 1;	
				std::vector<vec> temp_pos;
				std::vector<uint_2> temp_mol_stat_end;
				std::vector<unsigned int> body_size;
				
				temp_pos.resize(maxbody);
				temp_mol_stat_end.resize(maxbody);
				body_size.resize(maxbody);
			
				for(unsigned int idx=0; idx<m_N; idx++)
					{
					unsigned int body = h_body[idx];
					if(body!=NO_INDEX)
						{
						if (body_size[body]==0)
							{
							temp_mol_stat_end[body].x=idx;
							temp_pos[body] = m_pos[idx];
							m_pos0[idx] = m_pos[idx];
							}
						else
							{
							temp_mol_stat_end[body].y=idx;
							vec posi = temp_pos[body];
							float dx = m_pos[idx].x - posi.x;
							float dy = m_pos[idx].y - posi.y;
							float dz = m_pos[idx].z - posi.z;
							dx -= Lx*rint(dx*Lxinv);
							dy -= Ly*rint(dy*Lyinv);
							dz -= Lz*rint(dz*Lzinv);
							m_pos0[idx].x = posi.x + dx;
							m_pos0[idx].y = posi.y + dy;
							m_pos0[idx].z = posi.z + dz;					
							}
						body_size[body] += 1;
						}
					}
				m_n_mol += maxbody;
				for (unsigned int i=0; i<temp_mol_stat_end.size(); i++)
					m_mol_stat_end.push_back(temp_mol_stat_end[i]);
	
		//update the posion0 of grafed polymer
				vector<unsigned int> mol_id_particle;
				mol_id_particle.resize(m_N);			
				unsigned int n_mol = 0;
				for(unsigned int i=0; i<m_N; i++)
					mol_id_particle[i]=NO_INDEX;	
	
				for(unsigned int idx=0; idx<m_N; idx++)
					{
					unsigned int N_b = m_n_bond[idx];
					unsigned int body = h_body[idx];				
					if(mol_id_particle[idx]==NO_INDEX&&body!=NO_INDEX&&N_b>0)
						{
						mol_id_particle[idx] = n_mol;						
						templist.clear();
						templist.push_back(idx);
						unsigned int N_temp = templist.size();
						
						bool go =true;
						unsigned int bound = 0;
	
						while(go)
							{
							bound += 1;
							next_templist.clear();
							for(unsigned int i =0; i<N_temp;i++)
								{
								unsigned int tag = templist[i];
								vec posi = m_pos0[tag];
								unsigned int Nneibor = m_n_bond[tag];
								for(unsigned int j =0; j<Nneibor; j++)
									{
									unsigned int nextTag =(unsigned int) m_bonds[j*m_N + tag].x;
									bool exist = false;
									if(mol_id_particle[nextTag]!=NO_INDEX||h_body[nextTag]!=NO_INDEX)
										exist =true;
									if(!exist)
										{
										next_templist.push_back(nextTag);
										mol_id_particle[nextTag] = n_mol;									
										float dx = m_pos[nextTag].x - posi.x;
										float dy = m_pos[nextTag].y - posi.y;
										float dz = m_pos[nextTag].z - posi.z;
										dx -= Lx*rint(dx*Lxinv);
										dy -= Ly*rint(dy*Lyinv);
										dz -= Lz*rint(dz*Lzinv);
										m_pos0[nextTag].x = posi.x + dx;
										m_pos0[nextTag].y = posi.y + dy;
										m_pos0[nextTag].z = posi.z + dz;
										}
									}
								}	
							templist = next_templist;
							N_temp = templist.size();
							if(N_temp ==0)
								{
								go = false;
								}
							if(bound >= 100000)
								{
								go = false;
								cout<<"Warning! inital molecule error!"<<endl;
								}	
							}
						n_mol += 1;					
						}
					}
				}
			}		
		}
	
	for(unsigned int idx=0; idx<m_N; idx++)
		{
		unsigned int N_b = m_n_bond[idx];
		if(m_mol_id_particle[idx]==NO_INDEX)
			{
			if(N_b!=0&&molecule.size()!=m_N)
				{
				cerr << endl << "***Error, the particle "<<idx<< " without molecule index, but with " << N_b<<" bonds!" << endl << endl;	
				throw runtime_error("Error MolInfo::buildMol!");
				}
			m_pos0[idx] = m_pos[idx];
			addFreeParticleType(type_map[h_type[idx]]);
			}
		}		

	vector<unsigned int> molmark0;
	vector<unsigned int> molmark1;
	vector< vector<unsigned int> > mol_Ntype;	

	m_mol_size.resize(m_n_mol);
	for(unsigned int i =0; i<m_n_mol; i++)
		m_mol_size[i] = 0;
	molmark0.resize(m_n_mol);
	molmark1.resize(m_n_mol);
	mol_Ntype.resize(m_n_mol);
	m_mol_type.resize(m_n_mol);
	m_mol_type_id.resize(m_n_mol);
	
	for(unsigned int i=0; i<m_n_mol; i++)
		{
		mol_Ntype[i].resize(Ntype);
		}
	for(unsigned int i=0; i<m_N; i++)
		{
	    unsigned int molid = m_mol_id_particle[i];	
		if(molid!=NO_INDEX)
			{
			m_mol_size[molid] += 1;
			unsigned int typi = h_type[i];
			unsigned int Nb = m_n_bond[i];		
			molmark0[molid] += Nb*typi;
			molmark1[molid] += Nb*(typi+1);	
			mol_Ntype[molid][typi] += 1;
			}
		}
		
	for(unsigned int i=0; i<m_n_mol; i++)
		{
		stringstream s0,s1,s2;
		s0<<m_mol_size[i];
		s1<<molmark0[i];
		s2<<molmark1[i];
		
		string typemark;
		for(unsigned int j=0; j<Ntype; j++)
			{
			unsigned int num = mol_Ntype[i][j];
			stringstream st;				
			if(num>0)
				{
				st<<num;
				typemark += type_map[j]+"["+st.str()+"]";
				}
			}
		string mark;
		mark = s0.str()+"-"+typemark+"-"+s1.str()+"-"+s2.str();
		m_mol_type[i]=mark;
		}		
	for(unsigned int i=0; i<m_n_mol; i++)
		{
		m_mol_type_id[i] = switchNameToIndex(m_mol_type[i]);
		}

	m_n_mol_per_kind.resize(m_mol_type_exchmap.size());
	for(unsigned int i=0; i<m_mol_type_exchmap.size(); i++)
		m_n_mol_per_kind[i]=0;
	
	for(unsigned int i=0; i<m_n_mol; i++)
		{
	    unsigned int moltid = m_mol_type_id[i];
		m_n_mol_per_kind[moltid] += 1;
		}
	}
	
void MolInfo::outPutInfo()
	{
	if(m_mol_type_exchmap.size()>0)	
		cout<<"--- Molecules statistics"<<endl;		
	for(unsigned int i=0; i<m_mol_type_exchmap.size(); i++)
		{
		cout<<" "<<m_n_mol_per_kind[i]<<" Mol"<<i<<" with assigned name "<<m_mol_type_exchmap[i]<<endl;
		}
	}	
	

void MolInfo::updatePosition0()	
	{
	m_box = m_build.getBox();
	m_pos = m_build.getPos();

	float Lx = m_box.lx;
	float Ly = m_box.ly;
	float Lz = m_box.lz;

	double Lxinv = 0.0;
	double Lyinv = 0.0;
	double Lzinv = 0.0;
	
	if(Lx!=0.0)
		Lxinv = 1.0/Lx;
	if(Ly!=0.0)
		Lyinv = 1.0/Ly;
	if(Lz!=0.0)
		Lzinv = 1.0/Lz;
	std::vector<unsigned int> molecule = m_build.getMolecule();
	if(molecule.size()==m_N)
		{
		std::vector<std::vector<unsigned int> > mollist;
		mollist.resize(m_N);
		unsigned int mol_max=0;
		for(unsigned int idx=0; idx<m_N; idx++)
			{
			unsigned int mid = m_mol_id_particle[idx];
			if(mid!=NO_INDEX)
				{
				if(mid>=m_N)
					{
					cerr << endl << "***Error, the mol id "<<mid<< " greater than the particle number " << m_N<< endl << endl;	
					throw runtime_error("Error MolInfo::updatePosition0!");						
					}
				if(mid>mol_max)
					mol_max=mid;
				mollist[mid].push_back(idx);
				}
			}	
		for(unsigned int mid=0; mid<mol_max+1; mid++)
			{
			unsigned int npmol = mollist[mid].size();
			if(npmol==0)
				{
				cerr << endl << "***Error, the mol id "<<mid<< " has no particles !"<< endl << endl;	
				throw runtime_error("Error MolInfo::updatePosition0!");				
				}
			
			unsigned int tag0 = mollist[mid][0];
			vec posi = m_pos[tag0];
			m_pos0[tag0] = posi;			
			for(unsigned int idx=1; idx<npmol; idx++)
				{
				unsigned int nextTag = mollist[mid][idx];
				float dx = m_pos[nextTag].x - posi.x;
				float dy = m_pos[nextTag].y - posi.y;
				float dz = m_pos[nextTag].z - posi.z;
				dx -= Lx*rint(dx*Lxinv);
				dy -= Ly*rint(dy*Lyinv);
				dz -= Lz*rint(dz*Lzinv);
				m_pos0[nextTag].x = posi.x + dx;
				m_pos0[nextTag].y = posi.y + dy;
				m_pos0[nextTag].z = posi.z + dz;
				posi = m_pos0[nextTag];
				}			
			}				
		
		}
	else
		{
		std::vector<unsigned int> h_body = m_build.getBody();			
		std::vector<unsigned int > mol_id_particle;	
		mol_id_particle.resize(m_N);
		vector<unsigned int> templist;
		vector<unsigned int> next_templist;	
		unsigned int n_mol = 0;
	
		for(unsigned int i=0; i<m_N; i++)
			mol_id_particle[i]=NO_INDEX;
			
		for(unsigned int idx=0; idx<m_N; idx++)
			{
			unsigned int N_b = m_n_bond[idx];
			if(mol_id_particle[idx]==NO_INDEX&&N_b>0)
				{
				mol_id_particle[idx] = n_mol;
	
				templist.clear();
				templist.push_back(idx);
				m_pos0[idx] = m_pos[idx];
				unsigned int N_temp = templist.size();
	
	
				bool go =true;
				unsigned int bound = 0;
	
				while(go)
					{
					bound += 1;
					next_templist.clear();
					for(unsigned int i =0; i<N_temp;i++)
						{
						unsigned int tag = templist[i];
						vec posi = m_pos0[tag];
						unsigned int Nneibor = m_n_bond[tag];
						for(unsigned int j =0; j<Nneibor; j++)
							{
							unsigned int nextTag =(unsigned int) m_bonds[j*m_N + tag].x;
							bool exist = false;
							if(mol_id_particle[nextTag]!=NO_INDEX)
								exist =true;
							if(!exist)
								{
								next_templist.push_back(nextTag);
								mol_id_particle[nextTag] = n_mol;
								float dx = m_pos[nextTag].x - posi.x;
								float dy = m_pos[nextTag].y - posi.y;
								float dz = m_pos[nextTag].z - posi.z;
								dx -= Lx*rint(dx*Lxinv);
								dy -= Ly*rint(dy*Lyinv);
								dz -= Lz*rint(dz*Lzinv);
								m_pos0[nextTag].x = posi.x + dx;
								m_pos0[nextTag].y = posi.y + dy;
								m_pos0[nextTag].z = posi.z + dz;
								}
							}
						}
		
					templist = next_templist;
					N_temp = templist.size();
					if(N_temp ==0)
						{
						go = false;
						}
					if(bound >= 100000)
						{
						go = false;
						cout<<"Warning! inital molecule error!"<<endl;
						}	
					}
				n_mol += 1;
				}
			}
	
		if(h_body.size()==m_N)
			{
			int maxbody = -1;
			for(unsigned int idx=0; idx<m_N; idx++)
				{
				unsigned int body = h_body[idx];
				if(body!=NO_INDEX)
					{
					mol_id_particle[idx] = n_mol+body;
					if ((int) body>maxbody)
						maxbody=(int) body;
					}
				}
			
			if(maxbody>=0)
				{
				maxbody += 1;	
				std::vector<vec> temp_pos;
				std::vector<unsigned int> body_size;
				
				temp_pos.resize(maxbody);
				body_size.resize(maxbody);
			
				for(unsigned int idx=0; idx<m_N; idx++)
					{
					unsigned int body = h_body[idx];
					if(body!=NO_INDEX)
						{
						if (body_size[body]==0)
							{
							temp_pos[body] = m_pos[idx];
							m_pos0[idx] = m_pos[idx];
							}
						else
							{
							vec posi = temp_pos[body];
							float dx = m_pos[idx].x - posi.x;
							float dy = m_pos[idx].y - posi.y;
							float dz = m_pos[idx].z - posi.z;
							dx -= Lx*rint(dx*Lxinv);
							dy -= Ly*rint(dy*Lyinv);
							dz -= Lz*rint(dz*Lzinv);
							m_pos0[idx].x = posi.x + dx;
							m_pos0[idx].y = posi.y + dy;
							m_pos0[idx].z = posi.z + dz;					
							}
						body_size[body] += 1;
						}
					}	
				n_mol += maxbody;
				
		//update the posion0 of grafed polymer
		
				n_mol = 0;
				for(unsigned int i=0; i<m_N; i++)
					mol_id_particle[i]=NO_INDEX;	
				for(unsigned int idx=0; idx<m_N; idx++)
					{
					unsigned int N_b = m_n_bond[idx];
					unsigned int body = h_body[idx];				
					if(mol_id_particle[idx]==NO_INDEX&&body!=NO_INDEX&&N_b>0)
						{
						mol_id_particle[idx] = n_mol;
						templist.clear();
						templist.push_back(idx);
						unsigned int N_temp = templist.size();
						
						bool go =true;
						unsigned int bound = 0;
	
						while(go)
							{
							bound += 1;
							next_templist.clear();
							for(unsigned int i =0; i<N_temp;i++)
								{
								unsigned int tag = templist[i];
								vec posi = m_pos0[tag];
								unsigned int Nneibor = m_n_bond[tag];
								for(unsigned int j =0; j<Nneibor; j++)
									{
									unsigned int nextTag =(unsigned int) m_bonds[j*m_N + tag].x;
									bool exist = false;
									if(mol_id_particle[nextTag]!=NO_INDEX||h_body[nextTag]!=NO_INDEX)
										exist =true;
									if(!exist)
										{
										next_templist.push_back(nextTag);
										mol_id_particle[nextTag] = n_mol;
										float dx = m_pos[nextTag].x - posi.x;
										float dy = m_pos[nextTag].y - posi.y;
										float dz = m_pos[nextTag].z - posi.z;
										dx -= Lx*rint(dx*Lxinv);
										dy -= Ly*rint(dy*Lyinv);
										dz -= Lz*rint(dz*Lzinv);
										m_pos0[nextTag].x = posi.x + dx;
										m_pos0[nextTag].y = posi.y + dy;
										m_pos0[nextTag].z = posi.z + dz;
										}
									}
								}	
							templist = next_templist;
							N_temp = templist.size();
							if(N_temp ==0)
								{
								go = false;
								}
							if(bound >= 100000)
								{
								go = false;
								cout<<"Warning! inital molecule error!"<<endl;
								}	
							}
						n_mol += 1;	
						}
					}
				}
			}
		}
		
	for(unsigned int idx=0; idx<m_N; idx++)
		{
		unsigned int N_b = m_n_bond[idx];
		if(m_mol_id_particle[idx]==NO_INDEX)
			{
			if(N_b!=0&&molecule.size()!=m_N)
				{
				cerr << endl << "***Error, the particle "<<idx<< " without molecule index, but with " << N_b<<" bonds!" << endl << endl;	
				throw runtime_error("Error MolInfo::updatePosition0!");
				}
			m_pos0[idx] = m_pos[idx];
			}
		}		
	}
unsigned int MolInfo::switchNameToIndex(const std::string &name)
    {
    for (unsigned int i = 0; i < m_mol_type_exchmap.size(); i++)
        {
        if (m_mol_type_exchmap[i] == name)
            return i;
        }
	m_mol_type_exchmap.push_back(name);       
 
    return m_mol_type_exchmap.size()-1;
    }

unsigned int MolInfo::cellid(int i, int j, int k)
	{
	i = (i + (int)m_dim.x)%(int)m_dim.x;
	j = (j + (int)m_dim.y)%(int)m_dim.y;
	k = (k + (int)m_dim.z)%(int)m_dim.z;	
	return (unsigned int) (i + j*m_dim.x + k*m_dim.x*m_dim.y);
	}

void MolInfo::computeList(double r_cut)
	{
	double Lx = m_box.lx;
	double Ly = m_box.ly;
	double Lz = m_box.lz;

	double Lxinv = 0.0;
	double Lyinv = 0.0;
	double Lzinv = 0.0;
	
	if(Lx!=0.0)
		Lxinv = 1.0/Lx;
	if(Ly!=0.0)
		Lyinv = 1.0/Ly;
	if(Lz!=0.0)
		Lzinv = 1.0/Lz;
	
	if(r_cut<0)
		throw runtime_error("Error MolInfo computeList, negative rcut!");	
	if(r_cut>Lx/2.0||r_cut>Ly/2.0||(r_cut>Lz/2.0&&m_dimension==3))
		throw runtime_error("Error MolInfo computeList, rcut larger than half box!");

	m_dim.x = (unsigned int)(Lx/double(r_cut));
	m_dim.y = (unsigned int)(Ly/double(r_cut));
	m_dim.z = (unsigned int)(Lz/double(r_cut));
	
	
	if(m_dim.x==0)
		m_dim.x = 1;
	if(m_dim.y==0)
		m_dim.y = 1;
	if(m_dim.z==0)
		m_dim.z = 1;	
	
	m_width.x = Lx/double(m_dim.x);
	m_width.y = Ly/double(m_dim.y);
	m_width.z = Lz/double(m_dim.z);

	if(m_width.x==0)
		m_width.x = 1;
	if(m_width.y==0)
		m_width.y = 1;
	if(m_width.z==0)
		m_width.z = 1;	
	
	unsigned int m_N = m_pos.size(); 
	unsigned int ncell = m_dim.x*m_dim.y*m_dim.z;	
	m_head.resize(ncell);
	m_list.resize(m_N);
	for(unsigned int i=0; i< ncell; i++)
		m_head[i] = NO_INDEX;


	for(unsigned int i =0; i< m_N; i++)
		{
		vec posi = m_pos[i];
		double shiftx = 0.0;
		double shifty = 0.0;
		double shiftz = 0.0;		
		if(Lx != 0.0&&Ly !=0.0 &&Lz!=0.0)
			{
			shiftx = rint(posi.x*Lxinv);
			shifty = rint(posi.y*Lyinv);
			shiftz = rint(posi.z*Lzinv);
							
			posi.x -=  Lx *shiftx;
			posi.y -=  Ly *shifty;				
			posi.z -=  Lz *shiftz;
			}	
		int ix = int((posi.x+0.5*Lx)/m_width.x);
		int iy = int((posi.y+0.5*Ly)/m_width.y);
		int iz = int((posi.z+0.5*Lz)/m_width.z);		
		unsigned int cid = cellid(ix, iy, iz);
		if (cid >= ncell)
			{
			cerr << endl << "***Error, cell id unnormal!" << endl << endl;
			throw runtime_error("Error Generator updatePos");		
			}
		m_list[i]  = m_head[cid];
		m_head[cid] = i;
		}	
	}	
bool MolInfo::ifexclude(unsigned int a, unsigned int b)
    {
	buildExclusionList();
	unsigned int nex = m_n_exclusion[a];

    for (unsigned int i = 0; i < nex; i++)
        {
		if(i>=m_n_max_ex)
			{
			cerr << endl << "***Error, the number of excluded particles "<<nex<<" of particle "<<a<<" greater than the uplimited "<<m_n_max_ex<< endl << endl;
			throw runtime_error("Error ifexclude");	
			}
		unsigned int eid = m_exclusion_list[a+m_N*i];
        if (eid==b)
            return true;
        }
    return false;
    }


	