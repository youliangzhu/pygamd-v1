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


#include "DNAchain.h"

DNAchain::DNAchain(unsigned int Nbp, DNAchain::Strand s): Molecule((Nbp*3-1)*(int(s)+1)), m_s(s) 
    {
	allocate_DNAdata();
	m_xyz_temp.resize(m_NatomPerMole);
	m_mol_id_plus =0;
	}

DNAchain::DNAchain(const std::string& fname, unsigned int NatomPerMole, DNAchain::Strand s): Molecule(fname, NatomPerMole), m_s(s) 
    {
	allocate_DNAdata();
	m_xyz_temp.resize(m_NatomPerMole);
	m_mol_id_plus =0;
    }

DNAchain::DNAchain(unsigned int Nbp, DNAchain::Strand s, std::string form): Molecule((Nbp*3)*(int(s)+1)), m_s(s), m_form(form)
    {
	allocate_DNAdata();
	m_xyz_temp.resize(m_NatomPerMole);
	if(form=="ring")
		m_circle = true;
	}

DNAchain::DNAchain(const std::string& fname, unsigned int NatomPerMole, DNAchain::Strand s, std::string form): Molecule(fname, NatomPerMole), m_s(s), m_form(form)
    {
	allocate_DNAdata();
	m_xyz_temp.resize(m_NatomPerMole);
	if(form=="ring")
		m_circle=true;
    }
	
	
void DNAchain::allocate_DNAdata()
	{
	m_site_data1.resize(6);
	m_site_data2.resize(6);
	m_site_name.resize(6);	
	
	m_site_data1[sites::P] = vec4(2.186,8.918,94.038,94.97);			//the z ,r ,degree, mass	
	m_site_data1[sites::S] = vec4(1.280,6.981,70.197,83.11);			//the z ,r ,degree, mass
	m_site_data1[sites::Ab] = vec4(0.051,0.773,41.905,134.1);			//the z ,r ,degree, mass
	m_site_data1[sites::Tb] = vec4(0.191,2.349,86.119,125.1);			//the z ,r ,degree, mass	
	m_site_data1[sites::Cb] = vec4(0.187,2.296,85.027,110.1);			//the z ,r ,degree, mass	
	m_site_data1[sites::Gb] = vec4(0.053,0.828,40.691,150.1);			//the z ,r ,degree, mass

	m_site_data2[sites::P] = vec4(-2.186,8.918,265.962,94.97);			//the z ,r ,degree, mass
	m_site_data2[sites::S] = vec4(-1.280,6.981,289.803,83.11);			//the z ,r ,degree, mass
	m_site_data2[sites::Ab] = vec4(-0.051,0.773,318.095,134.1);			//the z ,r ,degree, mass
	m_site_data2[sites::Tb] = vec4(-0.191,2.349,273.881,125.1);			//the z ,r ,degree, mass	
	m_site_data2[sites::Cb] = vec4(-0.187,2.296,274.973,110.1);			//the z ,r ,degree, mass	
	m_site_data2[sites::Gb] = vec4(-0.053,0.828,319.309,150.1);			//the z ,r ,degree, mass

	m_site_name[sites::P] = "P";			
	m_site_name[sites::S] = "S";		
	m_site_name[sites::Ab] = "Ab";		
	m_site_name[sites::Tb] = "Tb";		
	m_site_name[sites::Cb] = "Cb";		
	m_site_name[sites::Gb] = "Gb";

	m_scale = 10.0;
	m_delt_degree = 36.0;
	m_delt_z = 3.38;

	m_circle = false;
	m_radius = 100;
	m_set_start_point=false;
	m_set_direction = false;
	}

void DNAchain::setSequences(std::string types, std::string fname)
	{
	readSequences(fname);
	if(m_circle)
		{
		double lenth = double(m_N_bps)*m_delt_z;
		m_radius = 0.5*lenth/M_PI;
		unsigned int res = m_N_bps%10;
		m_delt_degree = double(m_N_bps - res)/double(m_N_bps)*36.0;		
		}
	generateType();
	if(types.size()!=0&&types!=" "&&types!="")
		types +=",";
	types+=m_site_name[m_site_id[0]];
	for(unsigned int i=1; i< m_site_id.size(); i++)
		types+=","+m_site_name[m_site_id[i]];
	setParticleTypes(types);
	}

void DNAchain::setSequences(std::string fname)
	{
	readSequences(fname);
	if(m_circle)
		{
		double lenth = double(m_N_bps)*m_delt_z;
		m_radius = 0.5*lenth/M_PI;
		unsigned int res = m_N_bps%10;
		m_delt_degree = double(m_N_bps - res)/double(m_N_bps)*36.0;		
		}
	generateType();
	std::string types;
	types+=m_site_name[m_site_id[0]];
	for(unsigned int i=1; i< m_site_id.size(); i++)
		types+=","+m_site_name[m_site_id[i]];
	setParticleTypes(types);
	}

void DNAchain::readSequences(std::string fname)
	{
	ifstream file;
	file.open(fname.c_str());
	file.seekg(0,ios::beg);		
	if (!file.good())
        {
        cerr << endl << "Unable to open file " << fname << endl << endl;
        throw runtime_error("Error reading DNAchain::readSequences imput file");
        }

	std::string line;		
	std::string aim1 = "<sequence>";
	std::string aim2 = "</sequence>";
	unsigned int count = 0;
		
	while(getline(file, line) && line != aim1)
		{
		}
		
	if (!file.eof())
		{	
		cout << "INFO : read: " << line << "\n";
		while(getline(file, line) && line != aim2)
			{
			for(unsigned int i= 0; i< line.size();i++)
				{
				string temp;
				temp.push_back(line.at(i));
				m_sequence.push_back(temp);
				temp.clear();
				count +=1;
				}
			}
		}		
	m_N_bps = m_sequence.size();
	if(m_s==ss)
		{
		if(m_circle)
			m_N_sites = m_N_bps*3;
		else
			m_N_sites = m_N_bps*3-1;
		}
	else if(m_s==ds)
		{
		if(m_circle)
			m_N_sites = m_N_bps*3*2;
		else
			m_N_sites = m_N_bps*3*2-2;
		}
	else
        {
        cerr << endl << "Unable to recognise the strand type " << m_s << endl << endl;
        throw runtime_error("Error DNAchain::DNAchain");
        }
	cout << "INFO : Sequences statistics " << count<<" bp"<< endl;
	}

void DNAchain::generateType()
	{
	unsigned int count1 = 0;
	unsigned int count2 = m_N_sites-3;
	double z = 0.0;
	double r = 0.0;
	double d = 0.0;					//degree
	int t = 0;					//type
	vector<vec> rzd;
	rzd.resize(m_N_sites);
	m_site_id.resize(m_N_sites);
	m_site_xyz.resize(m_N_sites);
	for(unsigned int i = 0; i< m_N_bps ; i++)
		{
		std::string nucleo = m_sequence[i];

// the base of first strand
		if(m_circle||i!=m_N_bps-1)
			{
			z = m_site_data1[sites::P].x + double(i)*m_delt_z;
			r = m_site_data1[sites::P].y;
			d = m_site_data1[sites::P].z + double(i)*m_delt_degree;
			d = d - int(d/360.0)*360.0;
			t = int(sites::P);
			rzd[count1] = vec(z,r,d);
			m_site_id[count1] = t;
			}
	
		z = m_site_data1[sites::S].x + double(i)*m_delt_z;
		r = m_site_data1[sites::S].y;
		d = m_site_data1[sites::S].z + double(i)*m_delt_degree;
		d = d - int(d/360.0)*360.0;
		t = int(sites::S);
		rzd[count1+1] = vec(z,r,d);
		m_site_id[count1+1] = t;

		if(nucleo=="A")
			{
			z = m_site_data1[sites::Ab].x + double(i)*m_delt_z;
			r = m_site_data1[sites::Ab].y;
			d = m_site_data1[sites::Ab].z + double(i)*m_delt_degree;
			d = d - int(d/360.0)*360.0;
			t = int(sites::Ab);
			rzd[count1+2] = vec(z,r,d);
			m_site_id[count1+2] = t;
			}
		else if(nucleo=="G")
			{
			z = m_site_data1[sites::Gb].x + double(i)*m_delt_z;
			r = m_site_data1[sites::Gb].y;
			d = m_site_data1[sites::Gb].z + double(i)*m_delt_degree;
			d = d - int(d/360.0)*360.0;
			t =int(sites::Gb);
			rzd[count1+2] = vec(z,r,d);
			m_site_id[count1+2] = t;		
			}		
		else if(nucleo=="C")
			{
			z = m_site_data1[sites::Cb].x + double(i)*m_delt_z;
			r = m_site_data1[sites::Cb].y;
			d = m_site_data1[sites::Cb].z + double(i)*m_delt_degree;
			d = d - int(d/360.0)*360.0;
			t =int(sites::Cb);
			rzd[count1+2] = vec(z,r,d);
			m_site_id[count1+2] = t;		
			}		
		else if(nucleo=="T")
			{
			z = m_site_data1[sites::Tb].x + double(i)*m_delt_z;
			r = m_site_data1[sites::Tb].y;
			d = m_site_data1[sites::Tb].z + double(i)*m_delt_degree;
			d = d - int(d/360.0)*360.0;
			t =int(sites::Tb);
			rzd[count1+2] = vec(z,r,d);
			m_site_id[count1+2] = t;		
			}		
		else
			{
			cerr << endl << "Unable to recognise the inputed name " << nucleo << endl << endl;
			throw runtime_error("Error DNAchain::generateType");
			}
		if(m_circle)
			count1 +=3;
		else if(i==m_N_bps-2)
			count1 += 2;
		else
			count1 +=3;
// the base of second strand	
		if(m_s==ds)
			{
			if(m_circle||i!=0)
				{
				z = m_site_data2[sites::P].x + double(i)*m_delt_z;
				r = m_site_data2[sites::P].y;
				d = m_site_data2[sites::P].z + double(i)*m_delt_degree;
				d = d - int(d/360.0)*360.0;
				t =int(sites::P);
				rzd[count2] = vec(z,r,d);
				m_site_id[count2] = t;
				}

			z = m_site_data2[sites::S].x + double(i)*m_delt_z;
			r = m_site_data2[sites::S].y;
			d = m_site_data2[sites::S].z + double(i)*m_delt_degree;
			d = d - int(d/360.0)*360.0;
			t =int(sites::S);
			rzd[count2+1] = vec(z,r,d);
			m_site_id[count2+1] = t;	
			
			if(nucleo=="A")
				{
				z = m_site_data2[sites::Tb].x + double(i)*m_delt_z;
				r = m_site_data2[sites::Tb].y;
				d = m_site_data2[sites::Tb].z + double(i)*m_delt_degree;
				d = d - int(d/360.0)*360.0;
				t =int(sites::Tb);
				rzd[count2+2] = vec(z,r,d);
				m_site_id[count2+2] = t;				
				}
			else if(nucleo=="G")
				{
				z = m_site_data2[sites::Cb].x + double(i)*m_delt_z;
				r = m_site_data2[sites::Cb].y;
				d = m_site_data2[sites::Cb].z + double(i)*m_delt_degree;
				d = d - int(d/360.0)*360.0;
				t =int(sites::Cb);
				rzd[count2+2] = vec(z,r,d);
				m_site_id[count2+2] = t;			
				}		
			else if(nucleo=="C")
				{
				z = m_site_data2[sites::Gb].x + double(i)*m_delt_z;
				r = m_site_data2[sites::Gb].y;
				d = m_site_data2[sites::Gb].z + double(i)*m_delt_degree;
				d = d - int(d/360.0)*360.0;
				t =int(sites::Gb);
				rzd[count2+2] = vec(z,r,d);
				m_site_id[count2+2] = t;			
				}		
			else if(nucleo=="T")
				{
				z = m_site_data2[sites::Ab].x + double(i)*m_delt_z;
				r = m_site_data2[sites::Ab].y;
				d = m_site_data2[sites::Ab].z + double(i)*m_delt_degree;
				d = d - int(d/360.0)*360.0;
				t =int(sites::Ab);
				rzd[count2+2] = vec(z,r,d);
				m_site_id[count2+2] = t;		
				}		
			else
				{
				cerr << endl << "Unable to recognise the inputed name " << nucleo << endl << endl;
				throw runtime_error("Error DNAchain::generateType");
				}
			if(m_circle)
				count2 -=3;
			else if(i==0)
				count2 -= 2;
			else
				count2 -=3;
			}
		}
	double x = 0.0;
	double y = 0.0;
	
	for(unsigned int i = 0; i< m_N_sites ; i++)
		{
		z = rzd[i].x;
		r = rzd[i].y;
		d = rzd[i].z;
		
		x = cos(M_PI*d/180.0)*r;
		y = sin(M_PI*d/180.0)*r;

		m_site_xyz[i] = vec(x,y,z);
		}
	}

void DNAchain::generateSites()
	{
	if(m_N_sites>m_NatomPerMole)
        {
        cerr << endl << "The particle number "<< m_NatomPerMole << " is less than the generated sites number "<<m_N_sites<<" for "<<m_N_bps<<" base pairs!"<< endl << endl;
        throw runtime_error("Error DNAchain::DNAchain");
        }
// circle the linear DNA double helix		
	if(m_circle)
		{
	  for(unsigned int i = 0; i< m_N_sites; i++)
			{
			double x = m_site_xyz[i].x;
			double y = m_site_xyz[i].y;
			double z = m_site_xyz[i].z;
			
			double theta = z/m_radius;
			double cir_z = (m_radius - y)*sin(theta);
			double cir_y = m_radius - (m_radius - y)*cos(theta);

			m_site_xyz[i] = vec(x,cir_y,cir_z);
			}
		}
//transpose the molecular in box
	unsigned int count=0;
	vec pos0 = m_site_xyz[0];
	for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
		{
		if(!m_be_generated[taga])
			{
			vec posa = m_site_xyz[count];
			int t = m_site_id[count];
			m_mass[taga] = m_site_data1[t].w;

			m_xyz_temp[taga].x = (posa.x - pos0.x)/m_scale;
			m_xyz_temp[taga].y = (posa.y - pos0.y)/m_scale;
			m_xyz_temp[taga].z = (posa.z - pos0.z)/m_scale;
			count += 1;
			}
		}
	}

void DNAchain::generateTopology()
	{
//generate bonds
	vector<Bond> bonds;
	std::string bond_name;
	if(m_circle)
		{
		unsigned int offset = m_N_bps*3;
		for(unsigned int j = 0; j < m_N_bps; j++)
			{
			bond_name ="S-" + m_site_name[m_site_id[j*3+2]];
			bonds.push_back(Bond("S5-P", j*3, j*3+1));
			bonds.push_back(Bond(bond_name, j*3+1, j*3+2));
			if(j==m_N_bps-1)
				bonds.push_back(Bond("S3-P", j*3, 1));
			else
				bonds.push_back(Bond("S3-P", j*3, j*3+4));
			}
		if(m_s==ds)
			{
			for(unsigned int j = 0; j < m_N_bps; j++)
				{
				bond_name ="S-" + m_site_name[m_site_id[j*3+2+offset]];
				bonds.push_back(Bond("S5-P", j*3+offset, j*3+1+offset));
				bonds.push_back(Bond(bond_name, j*3+1+offset, j*3+2+offset));
				if(j==m_N_bps-1)
					bonds.push_back(Bond("S3-P", j*3+offset, 1+offset));
				else
					bonds.push_back(Bond("S3-P", j*3+offset, j*3+4+offset));
				}
			}
		}
	else
		{
		unsigned int offset = m_N_bps*3-1;
		for(unsigned int j = 0; j < m_N_bps-1; j++)
			{
			bond_name ="S-" + m_site_name[m_site_id[j*3+2]];
			bonds.push_back(Bond("S5-P", j*3, j*3+1));
			bonds.push_back(Bond(bond_name, j*3+1, j*3+2));
			if(j==m_N_bps-2)
				bonds.push_back(Bond("S3-P", j*3, j*3+3));
			else
				bonds.push_back(Bond("S3-P", j*3, j*3+4));			
			}
		bond_name ="S-" + m_site_name[m_site_id[(m_N_bps-1)*3+1]];
		bonds.push_back(Bond(bond_name, (m_N_bps-1)*3, (m_N_bps-1)*3+1));

		if(m_s==ds)
			{
			for(unsigned int j = 0; j < m_N_bps-1; j++)
				{
				bond_name ="S-" + m_site_name[m_site_id[j*3+2+offset]];
				bonds.push_back(Bond("S5-P", j*3+offset, j*3+1+offset));
				bonds.push_back(Bond(bond_name, j*3+1+offset, j*3+2+offset));
				if(j==m_N_bps-2)
					bonds.push_back(Bond("S3-P", j*3+offset, j*3+3+offset));
				else
					bonds.push_back(Bond("S3-P", j*3+offset, j*3+4+offset));			
				}
			bond_name ="S-" + m_site_name[m_site_id[(m_N_bps-1)*3+1+offset]];
			bonds.push_back(Bond(bond_name, (m_N_bps-1)*3+offset, (m_N_bps-1)*3+1+offset));
			}
		}
//generate angles
	vector<Angle> angles;
	std::string angle_name;
	if(m_circle)
		{
		unsigned int offset = m_N_bps*3;
		for(unsigned int j = 0; j < m_N_bps-1; j++)
			{
			angle_name ="P-5S-" + m_site_name[m_site_id[j*3+2]];
			angles.push_back(Angle(angle_name, j*3, j*3+1, j*3+2));			
			angles.push_back(Angle("S5-P-3S", j*3+1, j*3, j*3+4));
			angles.push_back(Angle("P-5S3-P", j*3, j*3+4, j*3+3));
			angle_name ="P-3S-" + m_site_name[m_site_id[j*3+5]];
			angles.push_back(Angle(angle_name, j*3, j*3+4, j*3+5));			
			}   
		angle_name ="P-5S-" + m_site_name[m_site_id[(m_N_bps-1)*3+2]];
		angles.push_back(Angle(angle_name, (m_N_bps-1)*3, (m_N_bps-1)*3+1, (m_N_bps-1)*3+2));
		angles.push_back(Angle("S5-P-3S", (m_N_bps-1)*3+1, (m_N_bps-1)*3, 1));
		angles.push_back(Angle("P-5S3-P", (m_N_bps-1)*3, 1, 0));
		angle_name ="P-3S-" + m_site_name[m_site_id[2]];
		angles.push_back(Angle(angle_name, (m_N_bps-1)*3, 1, 2));
		if(m_s==ds)
			{
			for(unsigned int j = 0; j < m_N_bps-1; j++)
				{
				angle_name ="P-5S-" + m_site_name[m_site_id[j*3+2+offset]];
				angles.push_back(Angle(angle_name, j*3+offset, j*3+1+offset, j*3+2+offset));			
				angles.push_back(Angle("S5-P-3S", j*3+1+offset, j*3+offset, j*3+4+offset));
				angles.push_back(Angle("P-5S3-P", j*3+offset, j*3+4+offset, j*3+3+offset));
				angle_name ="P-3S-" + m_site_name[m_site_id[j*3+5+offset]];
				angles.push_back(Angle(angle_name, j*3+offset, j*3+4+offset, j*3+5+offset));			
				}
					   
			angle_name ="P-5S-" + m_site_name[m_site_id[(m_N_bps-1)*3+2+offset]];
			angles.push_back(Angle(angle_name, (m_N_bps-1)*3+offset, (m_N_bps-1)*3+1+offset, (m_N_bps-1)*3+2+offset));
			angles.push_back(Angle("S5-P-3S", (m_N_bps-1)*3+1+offset, (m_N_bps-1)*3+offset, 1+offset));
			angles.push_back(Angle("P-5S3-P", (m_N_bps-1)*3+offset, 1+offset, 0+offset));
			angle_name ="P-3S-" + m_site_name[m_site_id[2+offset]];
			angles.push_back(Angle(angle_name, (m_N_bps-1)*3+offset, 1+offset, 2+offset));
			}	
		}
	else
		{
		unsigned int offset = m_N_bps*3-1;
		for(unsigned int j = 0; j < m_N_bps-1; j++)
			{
			angle_name ="P-5S-" + m_site_name[m_site_id[j*3+2]];
			angles.push_back(Angle(angle_name, j*3, j*3+1, j*3+2));
			if(j==m_N_bps-2)
				{
				angles.push_back(Angle("S5-P-3S", j*3+1, j*3, j*3+3));
				angle_name ="P-3S-" + m_site_name[m_site_id[j*3+4]];
				angles.push_back(Angle(angle_name, j*3, j*3+3, j*3+4));
				}
			else
				{
				angles.push_back(Angle("S5-P-3S", j*3+1, j*3, j*3+4));	
				angles.push_back(Angle("P-5S3-P", j*3, j*3+4, j*3+3));
				angle_name ="P-3S-" + m_site_name[m_site_id[j*3+5]];
				angles.push_back(Angle(angle_name, j*3, j*3+4, j*3+5));
				}			
			}
		if(m_s==ds)
			{			
			for(unsigned int j = 0; j < m_N_bps-1; j++)
				{
				angle_name ="P-5S-" + m_site_name[m_site_id[j*3+2+offset]];
				angles.push_back(Angle(angle_name, j*3+offset, j*3+1+offset, j*3+2+offset));
				if(j==m_N_bps-2)
					{
					angles.push_back(Angle("S5-P-3S", j*3+1+offset, j*3+offset, j*3+3+offset));
					angle_name ="P-3S-" + m_site_name[m_site_id[j*3+4+offset]];
					angles.push_back(Angle(angle_name, j*3+offset, j*3+3+offset, j*3+4+offset));
					}
				else
					{
					angles.push_back(Angle("S5-P-3S", j*3+1+offset, j*3+offset, j*3+4+offset));	
					angles.push_back(Angle("P-5S3-P", j*3+offset, j*3+4+offset, j*3+3+offset));
					angle_name ="P-3S-" + m_site_name[m_site_id[j*3+5+offset]];
					angles.push_back(Angle(angle_name, j*3+offset, j*3+4+offset, j*3+5+offset));
					}			
				}
			}				
		}
//generate dihedrals
	vector<Dihedral> dihedrals;
	std::string dihedral_name;
	if(m_circle)
		{
		unsigned int offset = m_N_bps*3;
		for(unsigned int j = 0; j < m_N_bps-2; j++)
			{
			dihedral_name ="S3-P-5S-" + m_site_name[m_site_id[j*3+2]];
			dihedrals.push_back(Dihedral(dihedral_name, j*3+4, j*3, j*3+1, j*3+2));				
			dihedrals.push_back(Dihedral("P-5S3-P-5S", j*3+3, j*3+4, j*3, j*3+1));
			dihedral_name =m_site_name[m_site_id[j*3+5]]+"-S3-P-5S";		
			dihedrals.push_back(Dihedral(dihedral_name, j*3+5, j*3+4, j*3, j*3+1));	
			dihedrals.push_back(Dihedral("S3-P-5S3-P", j*3+7, j*3+3, j*3+4, j*3));			
			}   
		dihedral_name ="S3-P-5S-" + m_site_name[m_site_id[(m_N_bps-2)*3+2]];
		dihedrals.push_back(Dihedral(dihedral_name, (m_N_bps-2)*3+4, (m_N_bps-2)*3, (m_N_bps-2)*3+1, (m_N_bps-2)*3+2));
		dihedrals.push_back(Dihedral("P-5S3-P-5S", (m_N_bps-2)*3+3, (m_N_bps-2)*3+4, (m_N_bps-2)*3, (m_N_bps-2)*3+1));
		dihedral_name =m_site_name[m_site_id[(m_N_bps-2)*3+5]]+"-S3-P-5S"; 			
		dihedrals.push_back(Dihedral(dihedral_name, (m_N_bps-2)*3+5, (m_N_bps-2)*3+4, (m_N_bps-2)*3, (m_N_bps-2)*3+1));			

		dihedrals.push_back(Dihedral("S3-P-5S3-P", 1, (m_N_bps-2)*3+3, (m_N_bps-2)*3+4, (m_N_bps-2)*3));
		dihedral_name ="S3-P-5S-" + m_site_name[m_site_id[(m_N_bps-1)*3+2]];
		dihedrals.push_back(Dihedral(dihedral_name, 1, (m_N_bps-1)*3, (m_N_bps-1)*3+1, (m_N_bps-1)*3+2));
		dihedrals.push_back(Dihedral("P-5S3-P-5S", 0, 1, (m_N_bps-1)*3, (m_N_bps-1)*3+1));
		dihedral_name =m_site_name[m_site_id[2]]+"-S3-P-5S";		
		dihedrals.push_back(Dihedral(dihedral_name, 2, 1, (m_N_bps-1)*3, (m_N_bps-1)*3+1));
		dihedrals.push_back(Dihedral("S3-P-5S3-P", 4, 0, 1, (m_N_bps-1)*3));
		if(m_s==ds)
			{
			for(unsigned int j = 0; j < m_N_bps-2; j++)
				{
				dihedral_name ="S3-P-5S-" + m_site_name[m_site_id[j*3+2+offset]];
				dihedrals.push_back(Dihedral(dihedral_name, j*3+4+offset, j*3+offset, j*3+1+offset, j*3+2+offset));				
				dihedrals.push_back(Dihedral("P-5S3-P-5S", j*3+3+offset, j*3+4+offset, j*3+offset, j*3+1+offset));
				dihedral_name =m_site_name[m_site_id[j*3+5+offset]]+"-S3-P-5S";		
				dihedrals.push_back(Dihedral(dihedral_name, j*3+5+offset, j*3+4+offset, j*3+offset, j*3+1+offset));	
				dihedrals.push_back(Dihedral("S3-P-5S3-P", j*3+7+offset, j*3+3+offset, j*3+4+offset, j*3+offset));				
				}   
			dihedral_name ="S3-P-5S-" + m_site_name[m_site_id[(m_N_bps-2)*3+2+offset]];
			dihedrals.push_back(Dihedral(dihedral_name, (m_N_bps-2)*3+4+offset, (m_N_bps-2)*3+offset, (m_N_bps-2)*3+1+offset, (m_N_bps-2)*3+2+offset));
			dihedrals.push_back(Dihedral("P-5S3-P-5S", (m_N_bps-2)*3+3+offset, (m_N_bps-2)*3+4+offset, (m_N_bps-2)*3+offset, (m_N_bps-2)*3+1+offset));
			dihedral_name =m_site_name[m_site_id[(m_N_bps-2)*3+5+offset]]+"-S3-P-5S"; 			
			dihedrals.push_back(Dihedral(dihedral_name, (m_N_bps-2)*3+5+offset, (m_N_bps-2)*3+4+offset, (m_N_bps-2)*3+offset, (m_N_bps-2)*3+1+offset));			

			dihedrals.push_back(Dihedral("S3-P-5S3-P", 1+offset, (m_N_bps-2)*3+3+offset, (m_N_bps-2)*3+4+offset, (m_N_bps-2)*3+offset));
			dihedral_name ="S3-P-5S-" + m_site_name[m_site_id[(m_N_bps-1)*3+2+offset]];
			dihedrals.push_back(Dihedral(dihedral_name, 1+offset, (m_N_bps-1)*3+offset, (m_N_bps-1)*3+1+offset, (m_N_bps-1)*3+2+offset));
			dihedrals.push_back(Dihedral("P-5S3-P-5S", 0+offset, 1+offset, (m_N_bps-1)*3+offset, (m_N_bps-1)*3+1+offset));
			dihedral_name =m_site_name[m_site_id[2+offset]]+"-S3-P-5S";		
			dihedrals.push_back(Dihedral(dihedral_name, 2+offset, 1+offset, (m_N_bps-1)*3+offset, (m_N_bps-1)*3+1+offset));
			dihedrals.push_back(Dihedral("S3-P-5S3-P", 4+offset, 0+offset, 1+offset, (m_N_bps-1)*3+offset));
			}
		}
	else
		{
		unsigned int offset = m_N_bps*3-1;
		for(unsigned int j = 0; j < m_N_bps-2; j++)
			{
			dihedral_name ="S3-P-5S-" + m_site_name[m_site_id[j*3+2]];
			dihedrals.push_back(Dihedral(dihedral_name, j*3+4, j*3, j*3+1, j*3+2));					
			dihedrals.push_back(Dihedral("P-5S3-P-5S", j*3+3, j*3+4, j*3, j*3+1));
			dihedral_name =m_site_name[m_site_id[j*3+5]]+"-S3-P-5S";		
			dihedrals.push_back(Dihedral(dihedral_name, j*3+5, j*3+4, j*3, j*3+1));
			if(j==m_N_bps-3)
				dihedrals.push_back(Dihedral("S3-P-5S3-P", j*3+6, j*3+3, j*3+4, j*3));
			else
				dihedrals.push_back(Dihedral("S3-P-5S3-P", j*3+7, j*3+3, j*3+4, j*3));
			}
				   
		dihedral_name ="S3-P-5S-" + m_site_name[m_site_id[(m_N_bps-2)*3+2]];
		dihedrals.push_back(Dihedral(dihedral_name, (m_N_bps-2)*3+3, (m_N_bps-2)*3, (m_N_bps-2)*3+1, (m_N_bps-2)*3+2));
		dihedral_name =m_site_name[m_site_id[(m_N_bps-2)*3+4]]+"-S3-P-5S"; 			
		dihedrals.push_back(Dihedral(dihedral_name, (m_N_bps-2)*3+4, (m_N_bps-2)*3+3, (m_N_bps-2)*3, (m_N_bps-2)*3+1));	
		if(m_s==ds)
			{
			for(unsigned int j = 0; j < m_N_bps-2; j++)
				{
				dihedral_name ="S3-P-5S-" + m_site_name[m_site_id[j*3+2+offset]];
				dihedrals.push_back(Dihedral(dihedral_name, j*3+4+offset, j*3+offset, j*3+1+offset, j*3+2+offset));
				dihedrals.push_back(Dihedral("P-5S3-P-5S", j*3+3+offset, j*3+4+offset, j*3+offset, j*3+1+offset));
				dihedral_name =m_site_name[m_site_id[j*3+5+offset]]+"-S3-P-5S";		
				dihedrals.push_back(Dihedral(dihedral_name, j*3+5+offset, j*3+4+offset, j*3+offset, j*3+1+offset));
				if(j==m_N_bps-3)
					dihedrals.push_back(Dihedral("S3-P-5S3-P", j*3+6+offset, j*3+3+offset, j*3+4+offset, j*3+offset));
				else
					dihedrals.push_back(Dihedral("S3-P-5S3-P", j*3+7+offset, j*3+3+offset, j*3+4+offset, j*3+offset));		
				}
					   
			dihedral_name ="S3-P-5S-" + m_site_name[m_site_id[(m_N_bps-2)*3+2+offset]];
			dihedrals.push_back(Dihedral(dihedral_name, (m_N_bps-2)*3+3+offset, (m_N_bps-2)*3+offset, (m_N_bps-2)*3+1+offset, (m_N_bps-2)*3+2+offset));
			dihedral_name =m_site_name[m_site_id[(m_N_bps-2)*3+4+offset]]+"-S3-P-5S"; 			
			dihedrals.push_back(Dihedral(dihedral_name, (m_N_bps-2)*3+4+offset, (m_N_bps-2)*3+3+offset, (m_N_bps-2)*3+offset, (m_N_bps-2)*3+1+offset));	
			}
		}
	for(unsigned int i=0; i<bonds.size(); i++)
		{
		bonds[i].a += m_Nread_particle;
		bonds[i].b += m_Nread_particle;
		m_bond.push_back(bonds[i]);
		}
	for(unsigned int i=0; i<angles.size(); i++)
		{
		angles[i].a +=m_Nread_particle;
		angles[i].b +=m_Nread_particle;
		angles[i].c +=m_Nread_particle;
		m_angle.push_back(angles[i]);
		}
	for(unsigned int i=0; i<dihedrals.size(); i++)
		{
		dihedrals[i].a +=m_Nread_particle;
		dihedrals[i].b +=m_Nread_particle;
		dihedrals[i].c +=m_Nread_particle;
		dihedrals[i].d +=m_Nread_particle;
		m_dihedral.push_back(dihedrals[i]);
		}
	}

void DNAchain::placeSites()
	{
	bool fitpoint = false;
	bool success = true;
	unsigned int turns =0;
	vec posa, centermass, angle;
	vec4 q, ar0, ar1, ar2;
	if(m_set_direction)
		{
		while(!fitpoint)
			{
			if(m_set_start_point)
				{
				centermass.x = m_start_point.x;
				centermass.y = m_start_point.y;
				centermass.z = m_start_point.z;
				}
			else
				{
				centermass.x = ( R2S() - 0.5 )*m_mol_Lx + m_shift_Lx;
				centermass.y = ( R2S() - 0.5 )*m_mol_Ly + m_shift_Ly;
				centermass.z = ( R2S() - 0.5 )*m_mol_Lz + m_shift_Lz;
				}
			double rsq = m_direction.x*m_direction.x + m_direction.y*m_direction.y + m_direction.z*m_direction.z;
			double r = sqrt(rsq);
			m_direction.x /= r;
			m_direction.y /= r;
			m_direction.z /= r;
			ar2.x = m_direction.x;
			ar2.y = m_direction.y;
			ar2.z = m_direction.z;
			vec random_point;
			random_point.x = ( R2S() - 0.5 );
			random_point.y = ( R2S() - 0.5 );
			random_point.z = ( R2S() - 0.5 );
			ar0.x = random_point.y*m_direction.z - random_point.z*m_direction.y;
			ar0.y = random_point.z*m_direction.x - random_point.x*m_direction.z;
			ar0.z = random_point.x*m_direction.y - random_point.y*m_direction.x;
			rsq = ar0.x*ar0.x + ar0.y*ar0.y + ar0.z*ar0.z;
			r = sqrt(rsq);		
			ar0.x /= r;
			ar0.y /= r;
			ar0.z /= r;
			ar1.x = ar2.y*ar0.z - ar2.z*ar0.y;
			ar1.y = ar2.z*ar0.x - ar2.x*ar0.z;
			ar1.z = ar2.x*ar0.y - ar2.y*ar0.x;
			for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
				{
				if(!m_be_generated[taga])
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
				}
			if(turns>10000)
				{
				success = false;
				fitpoint = true;
				}
			turns +=1;
			}		
		
		}
	else
		{
		while(!fitpoint)
			{
			if(m_set_start_point)
				{
				centermass.x = m_start_point.x;
				centermass.y = m_start_point.y;
				centermass.z = m_start_point.z;
				}
			else
				{
				centermass.x = ( R2S() - 0.5 )*m_mol_Lx + m_shift_Lx;
				centermass.y = ( R2S() - 0.5 )*m_mol_Ly + m_shift_Ly;
				centermass.z = ( R2S() - 0.5 )*m_mol_Lz + m_shift_Lz;
				}
			angle.x = R2S()*M_PI*2.0;
			angle.y = R2S()*M_PI*2.0;
			angle.z = R2S()*M_PI*2.0;
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
				if(!m_be_generated[taga])
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
				}
			if(turns>10000)
				{
				success = false;
				fitpoint = true;
				}
			turns +=1;
			}
		}
	if(!success)
		{
		cerr << endl << "***Error! Can not generate this DNAchain!"<<endl << endl;
		throw runtime_error("DNAchain::generate error");	
		}
	for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
		{
		if(!m_be_generated[taga])
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
			}
		}
	}

void DNAchain::generate()
	{
	m_xyz.clear();
	m_xyz.resize(m_NatomPerMole);

	m_ori_vec.clear();
	m_ori_vec.resize(m_NatomPerMole);
	
	if (m_firststep)
		{
		generateTopology();
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
		generateSites();
		
	for(unsigned int taga=0; taga<m_NatomPerMole; taga++)
		{
		if(m_s==ds)
			{
			if(taga>=m_Nread_particle&&taga-m_Nread_particle<m_N_sites/2)
				setMolecule(taga, 0);
			else
				setMolecule(taga, 1);
			}
		else if(m_s==ss)
			{
			if(taga>=m_Nread_particle)
				setMolecule(taga, 0);		
			}
		}		
		m_firststep = false;
		}
	placeSites();
	}	
	
void export_DNAchain(pybind11::module &m)
	{
	pybind11::class_<DNAchain>(m, "DNAchain")
		.def(pybind11::init<unsigned int, DNAchain::Strand >())
		.def(pybind11::init<const std::string&, unsigned int, DNAchain::Strand >())
		.def(pybind11::init<unsigned int, DNAchain::Strand, std::string>())
		.def(pybind11::init<const std::string&, unsigned int, DNAchain::Strand, std::string >())
		.def("setScale", &DNAchain::setScale)
		.def("setStartPoint", &DNAchain::setStartPoint)
		.def("setDirection", &DNAchain::setDirection)
		.def("setSequences", static_cast< void (DNAchain::*)(std::string) >(&DNAchain::setSequences))		
		.def("setSequences", static_cast< void (DNAchain::*)(std::string, std::string)>(&DNAchain::setSequences))		
		;
    // enum_<DNAchain::Strand>("Strand")
    // .value("ss",DNAchain::ss)
	// .value("ds",DNAchain::ds)
	// ;
	}

