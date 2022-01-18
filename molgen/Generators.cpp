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


#include "Generators.h"
#include<time.h> 
#include<stdlib.h> 

using namespace std;

Generators::Generators(double Lx, double Ly, double Lz) : m_Lx(Lx), m_Ly(Ly), m_Lz(Lz)
    {
	m_edge = 0.0;
	m_N = 0;
	m_Num_bonds = 0;
	m_Num_mol = 0;
	srand((int)time(0)); 
	m_dimension = 3;
	m_nprecision = 10;
	m_nhead = 7;
	m_rcut_max=1.0;
	m_NBtype=50;
	m_generated = false;
	m_params.resize(m_NBtype*m_NBtype);
	m_min_dis.resize(m_NBtype*m_NBtype);
    }

Generators::~Generators()
    {
	free(m_list);
	free(m_head);
    }	
	
void Generators::addMolecule(Molecule* Mol, unsigned int Nm)	
	{
	m_molecules.push_back(Mol);
	m_Nmol.push_back(Nm);
	Mol->setPutBox(m_Lx, m_Ly, m_Lz);
	unsigned int numb = Mol->getBond().size()*Nm;
	unsigned int num = Mol->getNumParticle()*Nm;
	m_N += num;
	m_Num_bonds += numb;
	m_Num_mol += Nm;
	}

unsigned int Generators::switchNametoType(const string& name)
	{
	for(unsigned int i=0; i<m_type_mapping_all.size(); i++)
		{
		if(m_type_mapping_all[i]==name)
			return i;
		}
	m_type_mapping_all.push_back(name);
	if(m_type_mapping_all.size()>m_NBtype)
		{
		throw runtime_error("Error Generators switchNametoType, m_NBtype is too small!");		
		}	
	return m_type_mapping_all.size() -1;
	}

void Generators::setParam(const string& name1, const string& name2, double epsilon, double sigma, double rcut)
	{
	unsigned int typ1 = switchNametoType(name1);
	unsigned int typ2 = switchNametoType(name2);
	if(typ1>=m_NBtype||typ2>=m_NBtype)
		{
		throw runtime_error("Error Generators setParam, m_NBtype is too small!");		
		}
	double lj1 = 4.0 * epsilon * pow(sigma, int(12));
	double lj2 = 4.0 * epsilon * pow(sigma, int(6));	
	m_params[typ1*m_NBtype+typ2] = vec(lj1, lj2, rcut*rcut);
	m_params[typ2*m_NBtype+typ1] = vec(lj1, lj2, rcut*rcut);
	
	if(rcut>m_rcut_max)
		m_rcut_max=rcut;
	}

double Generators::R2S()	
	{
	int ran = rand();
	double fran = (double)ran/(double)RAND_MAX;		
	return fran;
	}

unsigned int Generators::cellid(int i, int j, int k)
	{
	i = (i + (int)m_dim.x)%(int)m_dim.x;
	j = (j + (int)m_dim.y)%(int)m_dim.y;
	k = (k + (int)m_dim.z)%(int)m_dim.z;	
	return (unsigned int) (i + j*m_dim.x + k*m_dim.x*m_dim.y);
	}

void Generators::initiateList()
	{
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		std::vector<std::string> type = m_molecules[j]->getType();
		for (unsigned int k = 0; k < type.size(); k++)
			{
			unsigned int typi = switchNametoType(type[k]);
			if(typi>=m_NBtype)
				{
				throw runtime_error("Error Generators updateType, m_NBtype is too small!");				
				}			
			}
		}
	m_dim.x=1;
	m_dim.y=1;
	m_dim.z=1;
	if(m_Lx>0)
		m_dim.x = (unsigned int)(m_Lx/double(m_rcut_max));
	if(m_Ly>0)		
		m_dim.y = (unsigned int)(m_Ly/double(m_rcut_max));
	if(m_Lz>0)	
		m_dim.z = (unsigned int)(m_Lz/double(m_rcut_max));
	
	m_width.x = m_Lx/double(m_dim.x);
	m_width.y = m_Ly/double(m_dim.y);
	m_width.z = m_Lz/double(m_dim.z);
	
	unsigned int ncell = m_dim.x*m_dim.y*m_dim.z;	
	m_head = (unsigned int*) malloc(sizeof(unsigned int) * ncell);
	m_list = (unsigned int*) malloc(sizeof(unsigned int) * m_N);
	for(unsigned int i=0; i< ncell; i++)
		m_head[i] = NO_INDEX;

	for(unsigned int i=0; i<m_molecules.size(); i++)
		{
		m_molecules[i]->setCell(m_dim, m_width);
		m_molecules[i]->setList(m_list, m_head);
		m_molecules[i]->setParam(m_params, m_NBtype, m_type_mapping_all, m_min_dis);
		}
	}

void Generators::setMinimumDistance(double mini_dis)
	{
	for(unsigned int i=0;i<m_NBtype*m_NBtype;i++)
		m_min_dis[i]=mini_dis;
		
	if(mini_dis>m_rcut_max)
		m_rcut_max=mini_dis;
	}

void Generators::setMinimumDistance(const std::string& name1, const std::string& name2, double mini_dis)
	{
	unsigned int typ1 = switchNametoType(name1);
	unsigned int typ2 = switchNametoType(name2);
	if(typ1>=m_NBtype||typ2>=m_NBtype)
		{
		throw runtime_error("Error Generators setMinimumDistance, m_NBtype is too small!");		
		}

	m_min_dis[typ1*m_NBtype+typ2] = mini_dis;
	m_min_dis[typ2*m_NBtype+typ1] = mini_dis;

	if(mini_dis>m_rcut_max)
		m_rcut_max=mini_dis;	
	}	

void Generators::updatePos(std::vector<vec>& pos, std::vector<std::string>& type, std::vector<vec>& ori, std::vector<vec4>& quat,
							std::vector<unsigned int>& body, std::vector<unsigned int>& molecule)
	{
	for(unsigned int i=0; i<m_molecules.size(); i++)
		{
		m_molecules[i]->updatePos(pos);			
		m_molecules[i]->updateType(type);		
		}
	unsigned int offset = m_pos_all.size();
	
	for(unsigned int i=0; i< pos.size(); i++)
		m_pos_all.push_back(pos[i]);
		
	for(unsigned int i=0; i< ori.size(); i++)
		m_ori_all.push_back(ori[i]);
	
	for(unsigned int i=0; i< quat.size(); i++)
		m_quat_all.push_back(quat[i]);	
		
	for(unsigned int i=0; i< body.size(); i++)
		m_body_all.push_back(body[i]);
		
	for(unsigned int i=0; i< molecule.size(); i++)
		m_molecule_all.push_back(molecule[i]);
		
	unsigned int ncell = m_dim.x*m_dim.y*m_dim.z;
	unsigned int np = pos.size(); 
	for(unsigned int i =0; i< np; i++)
		{
		vec posi = pos[i];
		double shiftx = 0.0;
		double shifty = 0.0;
		double shiftz = 0.0;		
		if(m_Lx > 0.0)
			{
			shiftx = rint(posi.x/m_Lx);
			posi.x -=  m_Lx *shiftx;
			}	
		if(m_Ly > 0.0)
			{
			shifty = rint(posi.y/m_Ly);
			posi.y -=  m_Ly *shifty;
			}	
		if(m_Lz > 0.0)
			{
			shiftz = rint(posi.z/m_Lz);			
			posi.z -=  m_Lz *shiftz;
			}
			
		int ix = int((posi.x+0.5*m_Lx)/m_width.x);
		int iy = int((posi.y+0.5*m_Ly)/m_width.y);
		int iz = int((posi.z+0.5*m_Lz)/m_width.z);		
		unsigned int cid = cellid(ix, iy, iz);
		if (cid >= ncell)
			{
			cerr << endl << "***Error, cell id unnormal!" << endl << endl;
			throw runtime_error("Error Generators updatePos");		
			}
		m_list[i+offset]  = m_head[cid];
		m_head[cid] = i+offset;
		}
		
// for body evacuation		
	int body_min=1e6;
	int body_max=-1;
	for(unsigned int i =0; i< body.size(); i++)
		{
		unsigned int body_id = body[i];
		if (body_id == NO_INDEX) continue;
		if(int(body_id)<body_min)
			body_min = int(body_id);
		if(int(body_id)>body_max)
			body_max = int(body_id);
		}
		
	if(body_max== -1)
		return;

	unsigned int nbody = body_max - body_min +1;
	std::vector<vec4> body_com;
	std::vector<double> body_mass;
	body_com.resize(nbody);
	body_mass.resize(nbody);
	
	for(unsigned int i =0; i< body.size(); i++)
		{
		unsigned int body_id = body[i];
		if (body_id == NO_INDEX) continue;
		unsigned int body_shift = body_id - body_min;
		body_mass[body_shift] += 1.0;

		body_com[body_shift].x += pos[i].x;
		body_com[body_shift].y += pos[i].y;
		body_com[body_shift].z += pos[i].z;
		}

	for (unsigned int idx = 0; idx < nbody; idx++)
		{
		double mass = body_mass[idx];
		body_com[idx].x /= mass;
		body_com[idx].y /= mass;
		body_com[idx].z /= mass;
		}	

	for(unsigned int i =0; i< body.size(); i++)
		{
		unsigned int body_id = body[i];
		if (body_id == NO_INDEX) continue;
		unsigned int body_shift = body_id - body_min;
	
		double dx = pos[i].x - body_com[body_shift].x;
		double dy = pos[i].y - body_com[body_shift].y;
		double dz = pos[i].z - body_com[body_shift].z;
		double dis = sqrt(dx*dx + dy*dy + dz*dz);
		
		if(dis>body_com[body_shift].w)
			body_com[body_shift].w = dis;
		}	
	for(unsigned int i=0; i<m_molecules.size(); i++)
		{
		m_molecules[i]->updateBodyCom(body_com);		
		}			
	}
	
void Generators::generate()
	{
	if(m_dimension==2&&m_Lz>0.0)
		{
        cerr << "***Error! two dimensions of the simulation box should be with Lz = 0.0, the Lz = " << m_Lz <<" in xml files"<< endl << endl;
        throw runtime_error("Error Generators::generate");
		}
	
	if(m_dimension==3&&m_Lz<1.0e-6)
		{
        cerr << "***Error! Lz = 0.0 should be with two dimensions of the simulation, three dimensions defind in xml files!"<< endl << endl;
        cerr << "Please use setDimension() to change the dimension!"<< endl << endl;		
        throw runtime_error("Error Generators::generate");
		}

	if(m_dimension>3||m_dimension<1)
		{
		cerr << endl << "***Error! wrong dimension "<< m_dimension << endl << endl;
		throw runtime_error("Error Generators::generate");					
		}		
		
	if(m_generated)
		return;
	initiateList();
	unsigned int body_all_plus=0;
	unsigned int mol_all_plus=0;
	unsigned int indicate=0;
	unsigned int segm = m_Num_mol/100+1;
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		unsigned int Nm = m_Nmol[j];
//		unsigned int Np = m_molecules[j]->getNumParticle();
		for(unsigned int i =0; i< Nm; i++)
			{
			m_molecules[j]->generate();
			std::vector<vec> pos = m_molecules[j]->getPosition();
			std::vector<std::string> type = m_molecules[j]->getType();
			std::vector<vec> ori = m_molecules[j]->getOriVec();
			std::vector<vec4> quat = m_molecules[j]->getQuatVec();
			std::vector<unsigned int> body = m_molecules[j]->getBody();
			for(unsigned int k = 0; k< body.size(); k++)
				{
				if(body[k]!=NO_BODY)
					body[k]+=body_all_plus;
				}
			std::vector<unsigned int> molecule = m_molecules[j]->getMolecule();
			for(unsigned int k = 0; k< molecule.size(); k++)
				{
				if(molecule[k]!=NO_BODY)
					molecule[k]+=mol_all_plus;
				}
			unsigned int bip = m_molecules[j]->getBodyIdPlus();
			body_all_plus += bip;
			unsigned int mip = m_molecules[j]->getMolIdPlus();
			mol_all_plus += mip;			
			updatePos(pos, type, ori, quat, body, molecule);
			indicate += 1;
			if(indicate%segm==0)
				cout<<"complete "<<indicate*100/m_Num_mol<< "%"<<endl;
			}
		}
	m_generated=true;
	}

void Generators::outPutXml(std::string fname)
	{
	generate();
    std::string m_fname = fname + ".xml";
    ofstream f(m_fname.c_str());
	cout<<"Output xml ..."<<endl;
    if (!f.good())
        {
        cerr << endl << "***Error! Unable to open dump file for writing: " << fname << endl << endl;
        throw runtime_error("Error writing galamost_xml dump file");
        }
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
    f << "<galamost_xml version=\"1.3\">" << "\n";
    f << "<configuration time_step=\"" << "0" << "\" "
      << "dimensions=\"" << m_dimension << "\" "
      << "natoms=\"" << m_N << "\" "
      << ">" << "\n";
    f << "<box " << "lx=\""<< m_Lx << "\" ly=\""<< m_Ly << "\" lz=\""<< m_Lz << "\"/>" << "\n";	
	std::vector<vec_int> image;		
    f << "<position num=\"" << m_N << "\">" << "\n";
	for (unsigned int k = 0; k < m_pos_all.size(); k++)
		{
		double px = m_pos_all[k].x;
		double py = m_pos_all[k].y;
		double pz = m_pos_all[k].z;
					
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
		image.push_back(vec_int(int(shiftx),int(shifty),int(shiftz)));
		f <<setiosflags(ios::fixed)<<setprecision(m_nprecision)<<setw(m_nprecision+m_nhead)<< px <<setw(m_nprecision+m_nhead)<< py <<setw(m_nprecision+m_nhead)<< pz << "\n";
		if (!f.good())
			{
			cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
			throw runtime_error("Error writing galamost dump file");
			}
		}
	f <<"</position>" << "\n";
	
	f << "<image num=\"" << m_N << "\">" << "\n";
	for (unsigned int j = 0; j < m_N; j++)
		{
		f << image[j].x << " " << image[j].y << " "<< image[j].z << "\n";
		if (!f.good())
			{
			cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
			throw runtime_error("Error writing galamost dump file");
			}
		}			
    f <<"</image>" << "\n";

    f <<"<mass num=\"" << m_N << "\">" << "\n";
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		unsigned int Nm = m_Nmol[j];
		unsigned int Np = m_molecules[j]->getNumParticle();
		std::vector<double> mass = m_molecules[j]->getMass();
			
		for(unsigned int i =0; i< Nm; i++)
			{
			for(unsigned int k =0; k<Np; k++)
				f << mass[k] << "\n";
			if (!f.good())
				{
				cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
				throw runtime_error("Error writing galamost dump file");
				}			
			}
		}
    f <<"</mass>" << "\n";

    f <<"<type num=\"" << m_N << "\">" << "\n";
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		unsigned int Nm = m_Nmol[j];
		unsigned int Np = m_molecules[j]->getNumParticle();
		std::vector<std::string> type = m_molecules[j]->getType();
			
		for(unsigned int i =0; i< Nm; i++)
			{
			for(unsigned int k =0; k<Np; k++)
				f << type[k] << "\n";
            if (!f.good())
                {
                cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
                throw runtime_error("Error writing galamost dump file");
                }			
			}
		}
    f <<"</type>" << "\n";
	
	double tcharge = 0.0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<double> charge = m_molecules[j]->getCharge();
		for(unsigned int k=0;k<charge.size();k++)
			tcharge += fabs(charge[k]);
		}
	if(tcharge>0.0)
		{
		f <<"<charge num=\"" << m_N << "\">" << "\n";
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<double> charge = m_molecules[j]->getCharge();
				
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					f << charge[k] << "\n";
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}			
				}
			}
		f <<"</charge>" << "\n";
		}	
		
	bool outinert = false;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<vec> inert = m_molecules[j]->getInert();
		for(unsigned int k=0; k<inert.size();k++)
			{
			if(inert[k].x!=0.0||inert[k].y!=0.0||inert[k].z!=0.0)
				outinert=true;
			}
		}	

	if(outinert)
		{
		f <<"<inert num=\"" << m_N << "\">" << "\n";
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<vec> inert = m_molecules[j]->getInert();	
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					{
					vec dia = inert[k];
					f << dia.x << "  " <<dia.y << "  "<< dia.z << "\n";

					if (!f.good())
						{
						cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
						throw runtime_error("Error writing galamost dump file");
						}
					}
				}
			}
		f <<"</inert>" << "\n";
		}		
		

	unsigned int taniso = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<unsigned int> orientation = m_molecules[j]->getOrientation();
		for(unsigned int k=0;k<orientation.size();k++)
			taniso += orientation[k];
		}
	if(taniso>0)
		{
		f <<"<orientation num=\"" << m_N << "\">" << "\n";
		unsigned int count=0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<unsigned int> orientation = m_molecules[j]->getOrientation();	
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					{
					unsigned int aniso = orientation[k];
					vec ori = m_ori_all[count];
					double val = sqrt(ori.x*ori.x + ori.y*ori.y + ori.z*ori.z);
					if(aniso!=0&&val>0)
						{
						double orix = ori.x/val;
						double oriy = ori.y/val;
						double oriz = ori.z/val;
						f << orix<<"  "<<oriy<<"  "<<oriz<< "\n";
						}
					else
						{
						double orix = 0.0;
						double oriy = 0.0;
						double oriz = 0.0;							
						f << orix<<"  "<<oriy<<"  "<<oriz<< "\n";
						}
					if (!f.good())
						{
						cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
							throw runtime_error("Error writing galamost dump file");
						}
					count += 1;
					}
				}
			}
		f <<"</orientation>" << "\n";
		}			

	unsigned int tquat = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<unsigned int> quaternion = m_molecules[j]->getQuaternion();
		for(unsigned int k=0;k<quaternion.size();k++)
			tquat += quaternion[k];
		}
	if(tquat>0)
		{
		f <<"<quaternion num=\"" << m_N << "\">" << "\n";
		unsigned int count=0;		
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<unsigned int> quaternion = m_molecules[j]->getQuaternion();	
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					{
					unsigned int aniso = quaternion[k];
					vec4 quat = m_quat_all[count];
					double val = sqrt(quat.x*quat.x + quat.y*quat.y + quat.z*quat.z + quat.w*quat.w);					
					if(aniso!=0&&val>0)
						{
						double q0 = quat.x/val;
						double q1 = quat.y/val;
						double q2 = quat.z/val;
						double q3 = quat.w/val;						
						f << q0<<"  "<<q1<<"  "<<q2<<"  "<<q3<< "\n";
						}
					else
						{
						double q0 = 0.0;
						double q1 = 0.0;
						double q2 = 0.0;
						double q3 = 0.0;						
						f << q0<<"  "<<q1<<"  "<<q2<<"  "<<q3<< "\n";
						}

					if (!f.good())
						{
						cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
							throw runtime_error("Error writing galamost dump file");
						}
					count += 1;
					}
				}
			}
		f <<"</quaternion>" << "\n";
		}			

	bool outdiameter = false;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<double> diameter = m_molecules[j]->getDiameter();
		for(unsigned int k=0; k<diameter.size();k++)
			{
			if(diameter[k]!=0.0)
				outdiameter=true;
			}
		}	

	if(outdiameter)
		{
		f <<"<diameter num=\"" << m_N << "\">" << "\n";
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<double> diameter = m_molecules[j]->getDiameter();	
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					{
					double dia = diameter[k];
					f << dia << "\n";

					if (!f.good())
						{
						cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
						throw runtime_error("Error writing galamost dump file");
						}
					}
				}
			}
		f <<"</diameter>" << "\n";
		}

	unsigned int tinit = 0.0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<unsigned int> init = m_molecules[j]->getInit();
		for(unsigned int k=0;k<init.size();k++)
			tinit += fabs(init[k]);
		}

	if(tinit>0.0)
		{
		f <<"<h_init num=\"" << m_N << "\">" << "\n";
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<unsigned int> init = m_molecules[j]->getInit();
				
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					f << init[k] << "\n";
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}			
				}
			}
		f <<"</h_init>" << "\n";
		
		f <<"<h_cris num=\"" << m_N << "\">" << "\n";
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<unsigned int> cris = m_molecules[j]->getCris();
				
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					f << cris[k] << "\n";
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}			
				}
			}
		f <<"</h_cris>" << "\n";
		}
		
	bool tbody = false;
    for (unsigned int j = 0; j < m_body_all.size(); j++)
		{
		if(m_body_all[j]!=NO_INDEX)
			tbody=true;
		}
	if(tbody)
		{
		f <<"<body num=\"" << m_N << "\">" << "\n";
		for (unsigned int j = 0; j < m_body_all.size(); j++)
			{
			if(m_body_all[j]!=NO_INDEX)
				f << m_body_all[j]<< "\n";
			else 
				f << "-1"<< "\n";
			if (!f.good())
				{
				cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
				throw runtime_error("Error writing galamost dump file");
				}
			}
		f <<"</body>" << "\n";
		}

	bool outmolecule = false;
	for(unsigned int k=0; k<m_molecule_all.size();k++)
		{
		if(m_molecule_all[k]!=NO_INDEX)
			outmolecule=true;
		}
	if(outmolecule)
		{
		f <<"<molecule num=\"" << m_N << "\">" << "\n";
		for (unsigned int j = 0; j < m_molecule_all.size(); j++)
			{
			unsigned int mol = m_molecule_all[j];
			if(mol!=NO_INDEX)
				f << mol<< "\n";
			else 
				f << "-1"<< "\n";
			if (!f.good())
				{
				cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
				throw runtime_error("Error writing galamost dump file");
				}
			}
		f <<"</molecule>" << "\n";
		}

	unsigned int num_bond = 0;
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		unsigned int nb = 0;
		std::vector<Bond> bonds = m_molecules[j]->getBond();
		for(unsigned int k=0; k<bonds.size(); k++)
			if(bonds[k].bc=="b")
				nb +=1;
		num_bond += nb*m_Nmol[j];
		}	
	if(num_bond!=0)
		{
		f <<"<bond num=\"" << num_bond << "\">" << "\n";
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();			
			std::vector<Bond>& bonds = m_molecules[j]->getBond();
			unsigned int size = bonds.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					if(bonds[k].bc=="b")
						f << bonds[k].type << " " << bonds[k].a+count << " " << bonds[k].b+count<< "\n";	
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}
				count += Np;
				}
			}
		f << "</bond>" << "\n";
		}
	unsigned int num_constaint = 0;
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		unsigned int nc = 0;
		std::vector<Bond> bonds = m_molecules[j]->getBond();
		for(unsigned int k=0; k<bonds.size(); k++)
			if(bonds[k].bc=="c")
				nc +=1;
		num_constaint += nc*m_Nmol[j];
		}	
	if(num_constaint!=0)
		{
		f <<"<constraint num=\"" << num_constaint << "\">" << "\n";
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();			
			std::vector<Bond>& bonds = m_molecules[j]->getBond();
			unsigned int size = bonds.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					if(bonds[k].bc=="c")
						f << bonds[k].type << " " << bonds[k].a+count << " " << bonds[k].b+count<< "\n";	
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}
				count += Np;
				}
			}
		f << "</constraint>" << "\n";
		}		
	unsigned int num_angle = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<Angle> angles = m_molecules[j]->getAngle();
		num_angle += angles.size()*m_Nmol[j];
		}
	if(num_angle!=0)
		{
		f <<"<angle num=\"" << num_angle << "\">" << "\n";
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();				
			std::vector<Angle>& angles = m_molecules[j]->getAngle();
			unsigned int size = angles.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					f << angles[k].type << " " << angles[k].a+count << " " << angles[k].b+count<<" " << angles[k].c+count<< "\n";	
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}
				count += Np;
				}
			}
		f << "</angle>" << "\n";
		}

	unsigned int num_dihedral = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<Dihedral>& dihedrals = m_molecules[j]->getDihedral();
		num_dihedral += dihedrals.size()*m_Nmol[j];
		}	
	if(num_dihedral!=0)
		{
		f <<"<dihedral num=\"" << num_dihedral << "\">" << "\n";		
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();				
			std::vector<Dihedral>& dihedrals = m_molecules[j]->getDihedral();
			unsigned int size = dihedrals.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					f << dihedrals[k].type << " " << dihedrals[k].a+count << " " << dihedrals[k].b+count<<" " << dihedrals[k].c+count<<" " << dihedrals[k].d+count<< "\n";	
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}
				count += Np;
				}
			}
		f << "</dihedral>" << "\n";
		}
		
	unsigned int num_vsite = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<Dihedral>& vsites = m_molecules[j]->getVsite();
		num_vsite += vsites.size()*m_Nmol[j];
		}	
	if(num_vsite!=0)
		{
		f <<"<vsite num=\"" << num_vsite << "\">" << "\n";		
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();				
			std::vector<Dihedral>& vsites = m_molecules[j]->getVsite();
			unsigned int size = vsites.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					f << vsites[k].type << " " << vsites[k].a+count << " " << vsites[k].b+count<<" " << vsites[k].c+count<<" " << vsites[k].d+count<< "\n";	
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}
				count += Np;
				}
			}
		f << "</vsite>" << "\n";
		}		
		
	unsigned int num_asphere = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<str_vec6>& asphere = m_molecules[j]->getAsphere();
		num_asphere += asphere.size();
		}	
	if(num_asphere!=0)
		{
		f <<"<Aspheres>" << "\n";		
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{		
			std::vector<str_vec6>& asphere = m_molecules[j]->getAsphere();
			for(unsigned int k =0; k<asphere.size(); k++)
				f << asphere[k].name << " " << asphere[k].x << " " << asphere[k].y<<" " << asphere[k].z<<" " << asphere[k].w<<" " << asphere[k].m<<" " << asphere[k].n<< "\n";	
			}
		f << "</Aspheres>" << "\n";
		}

	unsigned int num_patch = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {		
		std::vector<str_vec6>& patch = m_molecules[j]->getPatch();
		num_patch += patch.size();
		}	
	if(num_patch!=0)
		{
		f <<"<Patchs>" << "\n";		
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{		
			std::vector<str_vec6>& patch_num = m_molecules[j]->getPatchNum();
			std::vector<str_vec6>& patch = m_molecules[j]->getPatch();			
			for(unsigned int k =0; k<patch_num.size(); k++)
				{
				unsigned int num = (unsigned int)patch_num[k].x;
				unsigned int start = (unsigned int)patch_num[k].y;
				unsigned int end = (unsigned int)patch_num[k].z;
				f << patch_num[k].name << " " << num<< "\n";			
				for(unsigned int l = start; l<end; l++)
					f << patch[l].name << " " << patch[l].x << " " << patch[l].y<<" " << patch[l].z<<" " << patch[l].w<< "\n";				
				}
	
			}
		f << "</Patchs>" << "\n";
		}		
    f << "</configuration>" << "\n";
    f << "</galamost_xml>" << "\n";
    
    if (!f.good())
        {
        cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
        throw runtime_error("Error writing galamost dump file");
        }
    f.close();
	cout<<"Success"<<endl;
	}
	
void Generators::outPutMol2(std::string fname)
	{
	generate();
    std::string m_fname = fname + ".mol2";
    ofstream f(m_fname.c_str());
	cout<<"Output mol2 ..."<<endl;
    if (!f.good())
        {
        cerr << endl << "***Error! Unable to open dump file for writing: " << fname << endl << endl;
        throw runtime_error("Error writing galamost_mol2 dump file");
        }
    f << "@<TRIPOS>MOLECULE" << "\n";
    f << "Generated by galamost" << "\n";		
    f << m_N << " " << m_Num_bonds << "\n";
    f << "NO_CHARGES" << "\n";
    f << "@<TRIPOS>ATOM" << "\n";
	unsigned int indicate=0;
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		unsigned int Nm = m_Nmol[j];
		unsigned int Np = m_molecules[j]->getNumParticle();
		std::vector<std::string> type = m_molecules[j]->getType();
		for(unsigned int i =0; i< Nm; i++)
			{
			for(unsigned int k =0; k<Np; k++)
				{
				double px = m_pos_all[indicate].x;
				double py = m_pos_all[indicate].y;
				double pz = m_pos_all[indicate].z;
							
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
				f << indicate+1 << " " << type[k] << " " << px << " " << py << " "<< pz << " " << type[k] << "\n";
				indicate += 1;
				}
			}
		}

	indicate=0;
    f << "@<TRIPOS>BOND" << "\n";
	if(m_Num_bonds>0)
		{
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
            {
            unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();			
			std::vector<Bond> bonds = m_molecules[j]->getBond();
			unsigned int size = bonds.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					{
					f << indicate+1 << " " << bonds[k].a+count+1 << " " << bonds[k].b+count+1 << " 1" << "\n";
					indicate +=1;
					}
				count += Np;
				}
            }
		}
    else
        {
        // write a dummy bond since VMD doesn't like loading mol2 files without bonds
        f << "1 1 2 1" << "\n";
        }
        
    if (!f.good())
        {
        cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
        throw runtime_error("Error writing mol2 dump file");
        }
        
    f.close();
	cout<<"Success"<<endl;
	}		

void Generators::outPutMST(std::string fname)
	{
	generate();
    std::string m_fname = fname + ".mst";
    ofstream f(m_fname.c_str());
	cout<<"Output mst ..."<<endl;
    if (!f.good())
        {
        cerr << endl << "***Error! Unable to open dump file for writing: " << fname << endl << endl;
        throw runtime_error("Error writing mst dump file");
        }

    f << "mst_version 1.0 #Copyright You-Liang Zhu" << "\n";
    f << "\tnum_particles" << "\n";	
    f << "\t\t" <<  m_N << "\n";	
    f << "\ttimestep" << "\n";		
    f << "\t\t0" << "\n";		
    f << "\tdimension" << "\n";		
    f << "\t\t" <<  m_dimension << "\n";		
    f << "\tbox " << "\n";	
    f << "\t\t" <<  m_Lx << "\t" <<  m_Ly << "\t" <<  m_Lz << "\n";	
	std::vector<vec_int> image;		
    f << "\tposition" << "\n";
	for (unsigned int k = 0; k < m_pos_all.size(); k++)
		{
		double px = m_pos_all[k].x;
		double py = m_pos_all[k].y;
		double pz = m_pos_all[k].z;
					
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
		image.push_back(vec_int(int(shiftx),int(shifty),int(shiftz)));
		f <<"\t\t"<<setiosflags(ios::fixed)<<setprecision(m_nprecision)<<setw(m_nprecision+m_nhead)<< px <<setw(m_nprecision+m_nhead)<< py <<setw(m_nprecision+m_nhead)<< pz << "\n";
		if (!f.good())
			{
			cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
			throw runtime_error("Error writing galamost dump file");
			}
		}

	f << "\timage" << "\n";
	for (unsigned int j = 0; j < m_N; j++)
		{
		f <<"\t\t"<<image[j].x << "\t" << image[j].y << "\t"<< image[j].z << "\n";
		if (!f.good())
			{
			cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
			throw runtime_error("Error writing galamost dump file");
			}
		}			

    f <<"\tmass" << "\n";
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		unsigned int Nm = m_Nmol[j];
		unsigned int Np = m_molecules[j]->getNumParticle();
		std::vector<double> mass = m_molecules[j]->getMass();
			
		for(unsigned int i =0; i< Nm; i++)
			{
			for(unsigned int k =0; k<Np; k++)
				f <<"\t\t"<< mass[k] << "\n";
			if (!f.good())
				{
				cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
				throw runtime_error("Error writing galamost dump file");
				}			
			}
		}

    f <<"\ttype" << "\n";
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		unsigned int Nm = m_Nmol[j];
		unsigned int Np = m_molecules[j]->getNumParticle();
		std::vector<std::string> type = m_molecules[j]->getType();
			
		for(unsigned int i =0; i< Nm; i++)
			{
			for(unsigned int k =0; k<Np; k++)
				f <<"\t\t"<<type[k] << "\n";
            if (!f.good())
                {
                cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
                throw runtime_error("Error writing galamost dump file");
                }			
			}
		}

	double tcharge = 0.0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<double> charge = m_molecules[j]->getCharge();
		for(unsigned int k=0;k<charge.size();k++)
			tcharge += fabs(charge[k]);
		}
	if(tcharge>0.0)
		{
		f <<"\tcharge" << "\n";
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<double> charge = m_molecules[j]->getCharge();
				
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					f <<"\t\t"<<charge[k] << "\n";
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}			
				}
			}
		}	
		
	bool outinert = false;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<vec> inert = m_molecules[j]->getInert();
		for(unsigned int k=0; k<inert.size();k++)
			{
			if(inert[k].x!=0.0||inert[k].y!=0.0||inert[k].z!=0.0)
				outinert=true;
			}
		}	

	if(outinert)
		{
		f <<"\tinert"<< "\n";
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<vec> inert = m_molecules[j]->getInert();	
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					{
					vec dia = inert[k];
					f <<"\t\t"<< dia.x << "\t" <<dia.y << "\t"<< dia.z << "\n";

					if (!f.good())
						{
						cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
						throw runtime_error("Error writing galamost dump file");
						}
					}
				}
			}
		}		

	unsigned int taniso = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<unsigned int> orientation = m_molecules[j]->getOrientation();
		for(unsigned int k=0;k<orientation.size();k++)
			taniso += orientation[k];
		}
	if(taniso>0)
		{
		f <<"\torientation"<< "\n";
		unsigned int count=0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<unsigned int> orientation = m_molecules[j]->getOrientation();	
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					{
					unsigned int aniso = orientation[k];
					vec ori = m_ori_all[count];
					double val = sqrt(ori.x*ori.x + ori.y*ori.y + ori.z*ori.z);
					if(aniso!=0&&val>0)
						{
						double orix = ori.x/val;
						double oriy = ori.y/val;
						double oriz = ori.z/val;
						f <<"\t\t"<< orix<<"\t"<<oriy<<"\t"<<oriz<< "\n";
						}
					else
						{
						double orix = 0.0;
						double oriy = 0.0;
						double oriz = 0.0;							
						f <<"\t\t"<< orix<<"\t"<<oriy<<"\t"<<oriz<< "\n";
						}
					if (!f.good())
						{
						cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
							throw runtime_error("Error writing galamost dump file");
						}
					count += 1;
					}
				}
			}
		}			

	unsigned int tquat = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<unsigned int> quaternion = m_molecules[j]->getQuaternion();
		for(unsigned int k=0;k<quaternion.size();k++)
			tquat += quaternion[k];
		}
	if(tquat>0)
		{
		f <<"\tquaternion"<< "\n";
		unsigned int count=0;		
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<unsigned int> quaternion = m_molecules[j]->getQuaternion();	
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					{
					unsigned int aniso = quaternion[k];
					vec4 quat = m_quat_all[count];
					double val = sqrt(quat.x*quat.x + quat.y*quat.y + quat.z*quat.z + quat.w*quat.w);					
					if(aniso!=0&&val>0)
						{
						double q0 = quat.x/val;
						double q1 = quat.y/val;
						double q2 = quat.z/val;
						double q3 = quat.w/val;						
						f <<"\t\t"<<q0<<"\t"<<q1<<"\t"<<q2<<"\t"<<q3<< "\n";
						}
					else
						{
						double q0 = 0.0;
						double q1 = 0.0;
						double q2 = 0.0;
						double q3 = 0.0;						
						f <<"\t\t" << q0<<"\t"<<q1<<"\t"<<q2<<"\t"<<q3<< "\n";
						}

					if (!f.good())
						{
						cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
							throw runtime_error("Error writing galamost dump file");
						}
					count += 1;
					}
				}
			}
		}			

	bool outdiameter = false;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<double> diameter = m_molecules[j]->getDiameter();
		for(unsigned int k=0; k<diameter.size();k++)
			{
			if(diameter[k]!=0.0)
				outdiameter=true;
			}
		}	

	if(outdiameter)
		{
		f <<"\tdiameter"<< "\n";
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<double> diameter = m_molecules[j]->getDiameter();	
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					{
					double dia = diameter[k];
					f <<"\t\t"<< dia << "\n";

					if (!f.good())
						{
						cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
						throw runtime_error("Error writing galamost dump file");
						}
					}
				}
			}
		}

	unsigned int tinit = 0.0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<unsigned int> init = m_molecules[j]->getInit();
		for(unsigned int k=0;k<init.size();k++)
			tinit += fabs(init[k]);
		}

	if(tinit>0.0)
		{
		f <<"\tinit" << "\n";
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<unsigned int> init = m_molecules[j]->getInit();
				
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					f <<"\t\t"<< init[k] << "\n";
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}			
				}
			}
		
		f <<"<\tcris" << "\n";
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();
			std::vector<unsigned int> cris = m_molecules[j]->getCris();
				
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<Np; k++)
					f <<"\t\t"<<cris[k] << "\n";
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}			
				}
			}
		}
		
	bool tbody = false;
    for (unsigned int j = 0; j < m_body_all.size(); j++)
		{
		if(m_body_all[j]!=NO_INDEX)
			tbody=true;
		}
	if(tbody)
		{
		f <<"\tbody" << "\n";
		for (unsigned int j = 0; j < m_body_all.size(); j++)
			{
			if(m_body_all[j]!=NO_INDEX)
				f <<"\t\t"<<m_body_all[j]<< "\n";
			else 
				f <<"\t\t"<<"-1"<< "\n";
			if (!f.good())
				{
				cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
				throw runtime_error("Error writing galamost dump file");
				}
			}
		}

	bool outmolecule = false;
	for(unsigned int k=0; k<m_molecule_all.size();k++)
		{
		if(m_molecule_all[k]!=NO_INDEX)
			outmolecule=true;
		}
	if(outmolecule)
		{
		f <<"\tmolecule"<< "\n";
		for (unsigned int j = 0; j < m_molecule_all.size(); j++)
			{
			unsigned int mol = m_molecule_all[j];
			if(mol!=NO_INDEX)
				f <<"\t\t"<<mol<< "\n";
			else 
				f <<"\t\t"<< "-1"<< "\n";
			if (!f.good())
				{
				cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
				throw runtime_error("Error writing galamost dump file");
				}
			}
		}

	unsigned int num_bond = 0;
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		unsigned int nb = 0;
		std::vector<Bond> bonds = m_molecules[j]->getBond();
		for(unsigned int k=0; k<bonds.size(); k++)
			if(bonds[k].bc=="b")
				nb +=1;
		num_bond += nb*m_Nmol[j];
		}	
	if(num_bond!=0)
		{
		f <<"\tbond"<< "\n";
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();			
			std::vector<Bond>& bonds = m_molecules[j]->getBond();
			unsigned int size = bonds.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					if(bonds[k].bc=="b")
						f <<"\t\t"<< bonds[k].type << "\t" << bonds[k].a+count << "\t" << bonds[k].b+count<< "\n";	
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}
				count += Np;
				}
			}
		}
	unsigned int num_constaint = 0;
	for (unsigned int j = 0; j < m_molecules.size(); j++)
		{
		unsigned int nc = 0;
		std::vector<Bond> bonds = m_molecules[j]->getBond();
		for(unsigned int k=0; k<bonds.size(); k++)
			if(bonds[k].bc=="c")
				nc +=1;
		num_constaint += nc*m_Nmol[j];
		}	
	if(num_constaint!=0)
		{
		f <<"\tconstraint"<< "\n";
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();			
			std::vector<Bond>& bonds = m_molecules[j]->getBond();
			unsigned int size = bonds.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					if(bonds[k].bc=="c")
						f <<"\t\t"<< bonds[k].type << "\t" << bonds[k].a+count << "\t" << bonds[k].b+count<< "\n";	
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}
				count += Np;
				}
			}
		}		
	unsigned int num_angle = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<Angle> angles = m_molecules[j]->getAngle();
		num_angle += angles.size()*m_Nmol[j];
		}
	if(num_angle!=0)
		{
		f <<"\tangle" << "\n";
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();				
			std::vector<Angle>& angles = m_molecules[j]->getAngle();
			unsigned int size = angles.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					f <<"\t\t"<< angles[k].type << "\t" << angles[k].a+count << "\t" << angles[k].b+count<<"\t" << angles[k].c+count<< "\n";	
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}
				count += Np;
				}
			}
		}

	unsigned int num_dihedral = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<Dihedral>& dihedrals = m_molecules[j]->getDihedral();
		num_dihedral += dihedrals.size()*m_Nmol[j];
		}	
	if(num_dihedral!=0)
		{
		f <<"\tdihedral"<< "\n";		
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();				
			std::vector<Dihedral>& dihedrals = m_molecules[j]->getDihedral();
			unsigned int size = dihedrals.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					f <<"\t\t"<< dihedrals[k].type << "\t" << dihedrals[k].a+count << "\t" << dihedrals[k].b+count<<"\t" << dihedrals[k].c+count<<"\t" << dihedrals[k].d+count<< "\n";	
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}
				count += Np;
				}
			}
		}
		
	unsigned int num_vsite = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<Dihedral>& vsites = m_molecules[j]->getVsite();
		num_vsite += vsites.size()*m_Nmol[j];
		}	
	if(num_vsite!=0)
		{
		f <<"\tvsite" << "\n";		
		unsigned int count =0;
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{
			unsigned int Nm = m_Nmol[j];
			unsigned int Np = m_molecules[j]->getNumParticle();				
			std::vector<Dihedral>& vsites = m_molecules[j]->getVsite();
			unsigned int size = vsites.size();
			for(unsigned int i =0; i< Nm; i++)
				{
				for(unsigned int k =0; k<size; k++)
					f <<"\t\t"<<vsites[k].type << "\t" << vsites[k].a+count << "\t" << vsites[k].b+count<<"\t" << vsites[k].c+count<<"\t" << vsites[k].d+count<< "\n";	
				if (!f.good())
					{
					cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
					throw runtime_error("Error writing galamost dump file");
					}
				count += Np;
				}
			}
		}		
		
	unsigned int num_asphere = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {
		std::vector<str_vec6>& asphere = m_molecules[j]->getAsphere();
		num_asphere += asphere.size();
		}	
	if(num_asphere!=0)
		{
		f <<"\tAspheres" << "\n";		
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{		
			std::vector<str_vec6>& asphere = m_molecules[j]->getAsphere();
			for(unsigned int k =0; k<asphere.size(); k++)
				f <<"\t\t"<< asphere[k].name << "\t" << asphere[k].x << "\t" << asphere[k].y<<"\t" << asphere[k].z<<"\t" << asphere[k].w<<"\t" << asphere[k].m<<"\t" << asphere[k].n<< "\n";	
			}
		}

	unsigned int num_patch = 0;
    for (unsigned int j = 0; j < m_molecules.size(); j++)
        {		
		std::vector<str_vec6>& patch = m_molecules[j]->getPatch();
		num_patch += patch.size();
		}	
	if(num_patch!=0)
		{
		f <<"\tPatchs" << "\n";		
		for (unsigned int j = 0; j < m_molecules.size(); j++)
			{		
			std::vector<str_vec6>& patch_num = m_molecules[j]->getPatchNum();
			std::vector<str_vec6>& patch = m_molecules[j]->getPatch();			
			for(unsigned int k =0; k<patch_num.size(); k++)
				{
				unsigned int num = (unsigned int)patch_num[k].x;
				unsigned int start = (unsigned int)patch_num[k].y;
				unsigned int end = (unsigned int)patch_num[k].z;
				f <<"\t\t"<<patch_num[k].name << "\t" << num<< "\n";			
				for(unsigned int l = start; l<end; l++)
					f <<"\t\t"<<patch[l].name << "\t" << patch[l].x << "\t" << patch[l].y<<"\t" << patch[l].z<<"\t" << patch[l].w<< "\n";				
				}
	
			}
		}		
    f <<"mst_end" << "\n";	
    if (!f.good())
        {
        cerr << endl << "***Error! Unexpected error writing galamost dump file" << endl << endl;
        throw runtime_error("Error writing galamost dump file");
        }
    f.close();
	cout<<"Success"<<endl;
	}
	
void export_Generators(pybind11::module &m)
	{
	pybind11::class_<Generators>(m, "Generators")	
        .def(pybind11::init<double, double, double >())	
		.def("addMolecule", &Generators::addMolecule)
		.def("setMinimumDistance", static_cast< void (Generators::*)(double) >(&Generators::setMinimumDistance))
		.def("setMinimumDistance", static_cast< void (Generators::*)(const std::string&, const std::string&, double)>(&Generators::setMinimumDistance))
		.def("outPutXml", &Generators::outPutXml)
		.def("outPutMol2", &Generators::outPutMol2)
		.def("outPutMST", &Generators::outPutMST)		
		.def("setParam", &Generators::setParam)
		.def("setDimension", &Generators::setDimension)		
		;
	}

