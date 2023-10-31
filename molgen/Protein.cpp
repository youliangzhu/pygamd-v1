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


#include "Protein.h"

Protein::Protein(const std::string& fname): Molecule(0)
    {
	m_cg_nmax=5;
	m_n_amino_acid_types=20;
	m_atom_nmax=10;	
	allocate_amino_acid_data();
	
	m_NatomPerMole = readSequences(fname);
	readPos(fname);
	allocateData(m_NatomPerMole);
	generateType();
	generateTopology();
	cout<<"Martini CG model"<<endl;
	}
	
	
Protein::Protein(const std::string& fname, Protein::Model m): Molecule(0)
    {
	if(m==martini)
		{
		m_cg_nmax=5;
		m_n_amino_acid_types=20;
		m_atom_nmax=10;	
		allocate_amino_acid_data();
		
		m_NatomPerMole = readSequences(fname);
		readPos(fname);
		allocateData(m_NatomPerMole);
		generateType();
		generateTopology();	
		cout<<"Martini CG model"<<endl;		
		}
	else if (m==one_bead_one_amino_acid)
		{
		m_cg_nmax=1;
		m_n_amino_acid_types=20;
		m_atom_nmax=10;	
		allocate_amino_acid_data2();

		m_NatomPerMole = readSequences(fname);
		readPos(fname);
		allocateData(m_NatomPerMole);
		generateType();
		generateTopology();		
		cout<<"One bead one amino acid CG model"<<endl;		
		}

	}	

void Protein::allocate_amino_acid_data()
	{
	m_name.resize(m_n_amino_acid_types);
	m_cg_np.resize(m_n_amino_acid_types);
	m_cg_ptype.resize(m_n_amino_acid_types*m_cg_nmax);
	m_natom.resize(m_n_amino_acid_types*m_cg_nmax);
	
	m_name[0]  = "ALA";
	m_name[1]  = "ARG";
	m_name[2]  = "ASN";
	m_name[3]  = "ASP";
	m_name[4]  = "CYS";
	m_name[5]  = "GLN";
	m_name[6]  = "GLU";
	m_name[7]  = "GLY";
	m_name[8]  = "HIS";
	m_name[9]  = "ILE";
	m_name[10] = "LEU";
	m_name[11] = "LYS";
	m_name[12] = "MET";
	m_name[13] = "PHE";
	m_name[14] = "PRO";
	m_name[15] = "SER";
	m_name[16] = "THR";
	m_name[17] = "TRP";
	m_name[18] = "TYR";
	m_name[19] = "VAL";
	

	m_cg_np[0]  = 1;	
	m_cg_np[1]  = 3;	
	m_cg_np[2]  = 2;	
	m_cg_np[3]  = 2;	
	m_cg_np[4]  = 2;	
	m_cg_np[5]  = 2;	
	m_cg_np[6]  = 2;	
	m_cg_np[7]  = 1;	
	m_cg_np[8]  = 4;	
	m_cg_np[9]  = 2;	
	m_cg_np[10] = 2;	
	m_cg_np[11] = 3;	
	m_cg_np[12] = 2;	
	m_cg_np[13] = 4;	
	m_cg_np[14] = 2;	
	m_cg_np[15] = 2;	
	m_cg_np[16] = 2;	
	m_cg_np[17] = 5;	
	m_cg_np[18] = 4;	
	m_cg_np[19] = 2;
	
	std::string positive = "K";
	std::string negtive = "S";

	m_cg_ptype[0*m_cg_nmax+0] = "G";
	m_cg_ptype[1*m_cg_nmax+0] = "G";
	m_cg_ptype[1*m_cg_nmax+1] = "C";	
	m_cg_ptype[1*m_cg_nmax+2] = positive;
	m_cg_ptype[2*m_cg_nmax+0] = "G";
	m_cg_ptype[2*m_cg_nmax+1] = "G";
	m_cg_ptype[3*m_cg_nmax+0] = "G";
	m_cg_ptype[3*m_cg_nmax+1] = negtive;	
	m_cg_ptype[4*m_cg_nmax+0] = "G";
	m_cg_ptype[4*m_cg_nmax+1] = "C";
	m_cg_ptype[5*m_cg_nmax+0] = "G";
	m_cg_ptype[5*m_cg_nmax+1] = "G";
	m_cg_ptype[6*m_cg_nmax+0] = "G";
	m_cg_ptype[6*m_cg_nmax+1] = negtive;
	m_cg_ptype[7*m_cg_nmax+0] = "G";
	m_cg_ptype[8*m_cg_nmax+0] = "G";
	m_cg_ptype[8*m_cg_nmax+1] = "C";
	m_cg_ptype[8*m_cg_nmax+2] = "G";
	m_cg_ptype[8*m_cg_nmax+3] = "G";
	m_cg_ptype[9*m_cg_nmax+0] = "G";
	m_cg_ptype[9*m_cg_nmax+1] = "C";
	m_cg_ptype[10*m_cg_nmax+0] = "G";
	m_cg_ptype[10*m_cg_nmax+1] = "C";
	m_cg_ptype[11*m_cg_nmax+0] = "G";
	m_cg_ptype[11*m_cg_nmax+1] = "C";
	m_cg_ptype[11*m_cg_nmax+2] = positive;
	m_cg_ptype[12*m_cg_nmax+0] = "G";
	m_cg_ptype[12*m_cg_nmax+1] = "C";
	m_cg_ptype[13*m_cg_nmax+0] = "G";
	m_cg_ptype[13*m_cg_nmax+1] = "C";
	m_cg_ptype[13*m_cg_nmax+2] = "C";
	m_cg_ptype[13*m_cg_nmax+3] = "C";
	m_cg_ptype[14*m_cg_nmax+0] = "G";
	m_cg_ptype[14*m_cg_nmax+1] = "C";	
	m_cg_ptype[15*m_cg_nmax+0] = "G";
	m_cg_ptype[15*m_cg_nmax+1] = "G";
	m_cg_ptype[16*m_cg_nmax+0] = "G";
	m_cg_ptype[16*m_cg_nmax+1] = "G";
	m_cg_ptype[17*m_cg_nmax+0] = "G";
	m_cg_ptype[17*m_cg_nmax+1] = "C";
	m_cg_ptype[17*m_cg_nmax+2] = "G";
	m_cg_ptype[17*m_cg_nmax+3] = "C";	
	m_cg_ptype[17*m_cg_nmax+4] = "C";
	m_cg_ptype[18*m_cg_nmax+0] = "G";
	m_cg_ptype[18*m_cg_nmax+1] = "C";
	m_cg_ptype[18*m_cg_nmax+2] = "C";
	m_cg_ptype[18*m_cg_nmax+3] = "G";
	m_cg_ptype[19*m_cg_nmax+0] = "G";
	m_cg_ptype[19*m_cg_nmax+1] = "C";	

	m_natom[0*m_cg_nmax+0] = 5;
	m_natom[1*m_cg_nmax+0] = 4;
	m_natom[1*m_cg_nmax+1] = 3;	
	m_natom[1*m_cg_nmax+2] = 4;
	m_natom[2*m_cg_nmax+0] = 4;
	m_natom[2*m_cg_nmax+1] = 4;
	m_natom[3*m_cg_nmax+0] = 4;
	m_natom[3*m_cg_nmax+1] = 4;	
	m_natom[4*m_cg_nmax+0] = 4;
	m_natom[4*m_cg_nmax+1] = 2;
	m_natom[5*m_cg_nmax+0] = 4;
	m_natom[5*m_cg_nmax+1] = 5;
	m_natom[6*m_cg_nmax+0] = 4;
	m_natom[6*m_cg_nmax+1] = 5;
	m_natom[7*m_cg_nmax+0] = 4;
	m_natom[8*m_cg_nmax+0] = 4;
	m_natom[8*m_cg_nmax+1] = 2;
	m_natom[8*m_cg_nmax+2] = 2;
	m_natom[8*m_cg_nmax+3] = 2;
	m_natom[9*m_cg_nmax+0] = 4;
	m_natom[9*m_cg_nmax+1] = 4;
	m_natom[10*m_cg_nmax+0] = 4;
	m_natom[10*m_cg_nmax+1] = 4;
	m_natom[11*m_cg_nmax+0] = 4;
	m_natom[11*m_cg_nmax+1] = 3;
	m_natom[11*m_cg_nmax+2] = 2;
	m_natom[12*m_cg_nmax+0] = 4;
	m_natom[12*m_cg_nmax+1] = 4;
	m_natom[13*m_cg_nmax+0] = 4;
	m_natom[13*m_cg_nmax+1] = 3;
	m_natom[13*m_cg_nmax+2] = 2;
	m_natom[13*m_cg_nmax+3] = 2;
	m_natom[14*m_cg_nmax+0] = 4;
	m_natom[14*m_cg_nmax+1] = 3;	
	m_natom[15*m_cg_nmax+0] = 4;
	m_natom[15*m_cg_nmax+1] = 2;
	m_natom[16*m_cg_nmax+0] = 4;
	m_natom[16*m_cg_nmax+1] = 3;
	m_natom[17*m_cg_nmax+0] = 4;
	m_natom[17*m_cg_nmax+1] = 3;
	m_natom[17*m_cg_nmax+2] = 3;
	m_natom[17*m_cg_nmax+3] = 2;	
	m_natom[17*m_cg_nmax+4] = 2;
	m_natom[18*m_cg_nmax+0] = 4;
	m_natom[18*m_cg_nmax+1] = 3;
	m_natom[18*m_cg_nmax+2] = 2;
	m_natom[18*m_cg_nmax+3] = 3;
	m_natom[19*m_cg_nmax+0] = 4;
	m_natom[19*m_cg_nmax+1] = 3;
	}
	
	
void Protein::allocate_amino_acid_data2()
	{
	m_name.resize(m_n_amino_acid_types);
	m_cg_np.resize(m_n_amino_acid_types);
	m_cg_ptype.resize(m_n_amino_acid_types*m_cg_nmax);
	m_natom.resize(m_n_amino_acid_types*m_cg_nmax);
	
	m_name[0]  = "ALA"; 
	m_name[1]  = "ARG"; 
	m_name[2]  = "ASN"; 
	m_name[3]  = "ASP"; 
	m_name[4]  = "CYS"; 
	m_name[5]  = "GLN";  
	m_name[6]  = "GLU";  
	m_name[7]  = "GLY";  
	m_name[8]  = "HIS"; 
	m_name[9]  = "ILE"; 
	m_name[10] = "LEU";  
	m_name[11] = "LYS"; 
	m_name[12] = "MET";  
	m_name[13] = "PHE"; 
	m_name[14] = "PRO";  
	m_name[15] = "SER"; 
	m_name[16] = "THR"; 
	m_name[17] = "TRP"; 
	m_name[18] = "TYR"; 
	m_name[19] = "VAL"; 
	

	m_cg_np[0]  = 1;	
	m_cg_np[1]  = 1;	
	m_cg_np[2]  = 1;	
	m_cg_np[3]  = 1;	
	m_cg_np[4]  = 1;	
	m_cg_np[5]  = 1;	
	m_cg_np[6]  = 1;	
	m_cg_np[7]  = 1;	
	m_cg_np[8]  = 1;	
	m_cg_np[9]  = 1;	
	m_cg_np[10] = 1;	
	m_cg_np[11] = 1;	
	m_cg_np[12] = 1;	
	m_cg_np[13] = 1;	
	m_cg_np[14] = 1;	
	m_cg_np[15] = 1;	
	m_cg_np[16] = 1;	
	m_cg_np[17] = 1;	
	m_cg_np[18] = 1;	
	m_cg_np[19] = 1;


	m_cg_ptype[0] =  "A";
	m_cg_ptype[1] =  "R";
	m_cg_ptype[2] =  "N";	
	m_cg_ptype[3] =  "D";
	m_cg_ptype[4] =  "C";
	m_cg_ptype[5] =  "Q";
	m_cg_ptype[6] =  "E";
	m_cg_ptype[7] =  "G";	
	m_cg_ptype[8] =  "H";
	m_cg_ptype[9] =  "I";
	m_cg_ptype[10] = "L";
	m_cg_ptype[11] = "K";
	m_cg_ptype[12] = "M";	
	m_cg_ptype[13] = "F";
	m_cg_ptype[14] = "P";
	m_cg_ptype[15] = "S";
	m_cg_ptype[16] = "T";
	m_cg_ptype[17] = "W";	
	m_cg_ptype[18] = "Y";
	m_cg_ptype[19] = "V";

	m_natom[0]  =  5;
	m_natom[1]  = 11;
	m_natom[2]  =  8;	
	m_natom[3]  =  8;
	m_natom[4]  =  6;
	m_natom[5]  =  9;
	m_natom[6]  =  9;
	m_natom[7]  =  4;	
	m_natom[8]  = 10;
	m_natom[9]  =  8;
	m_natom[10] =  8;
	m_natom[11] =  9;
	m_natom[12] =  8;	
	m_natom[13] = 11;
	m_natom[14] =  7;
	m_natom[15] =  6;
	m_natom[16] =  7;
	m_natom[17] = 14;	
	m_natom[18] = 12;
	m_natom[19] =  7;
	}
	

unsigned int Protein::getIndex(std::string name)
	{
	for (unsigned int i = 0; i < m_n_amino_acid_types; i++)
        {
        if (m_name[i] == name)
            return i;
        }
        
    cerr << endl << "***Error! Amino Acid Type " << name << " do not exist!" << endl;
    throw runtime_error("Error Protein::getIndex");
    return 0;
	}

unsigned int Protein::readSequences(std::string fname)
	{
	ifstream file;
	file.open(fname.c_str());
	file.seekg(0,ios::beg);		
	if (!file.good())
        {
        cerr << endl << "Unable to open file " << fname << endl << endl;
        throw runtime_error("Error reading Protein::readSequences imput file");
        }

	std::string line;		
	std::string aim1 = "<sequence>";
	std::string aim2 = "</sequence>";
		
	while(getline(file, line) && line != aim1)
		{
		}
		
	if (!file.eof())
		{	
		cout << "INFO : read: " << line << "\n";
		while(getline(file, line) && line != aim2)
			{
			istringstream parser(line);
			if(parser.good())
				{
				std::string aa_name;
				while(parser>>aa_name)
					{
					m_sequence.push_back(aa_name);
					}
				}
			else
				{
				cerr << endl << "Unable to parse line, parser.good() failed" << endl << endl;
				throw runtime_error("Error parser(line)");			
				}
			}
		}
	else
		{
		cout << endl << "Warning!!! Can not find sequence node!" << endl << endl;	
		}			
	m_N_amino_acid = m_sequence.size();
	unsigned int ncg=0;
	unsigned int natom=0;
	for(unsigned int i=0; i< m_N_amino_acid; i++)
		{
		string aa = m_sequence[i];
		unsigned int idx = getIndex(aa);
//		cout<<i+1<<" "<<aa<<endl;
		ncg += m_cg_np[idx];
		for(unsigned int j=0; j< m_cg_np[idx]; j++)
			{
			natom += m_natom[idx*m_cg_nmax+j];
			}
		}
	file.close();
	cout << "INFO : Sequences statistics " << m_N_amino_acid<<" bases, "<< ncg<<" coarse-grained particles, "<<natom<<" atoms "<<endl;
	return ncg;
	}


void Protein::readPos(std::string fname)
	{
    ifstream file;
    file.open(fname.c_str());	

	if (!file.good())
		{
		cerr << endl << "Unable to open file " << fname.c_str() << endl << endl;
		throw runtime_error("Error reading Protein::readPos file");
		}
	else
		{
		cout<<"INFO : Read the file "<<fname.c_str()<<endl;
		}

	file.seekg(0,ios::beg);			
	std::string line;		
	std::string aim1 = "<position>";
	std::string aim2 = "</position>";
	while(getline(file, line) && line != aim1)
		{
		}
	if (!file.eof())
		{	
		cout << "read " << line << endl;
//		cout<<"	"<<"type"<<", "<<"px"<<", "<<"py"<<", "<<"pz"<<", "<<"mol_id"<<", "<<"mol_name"<<endl;						
		
		while(getline(file, line) && line != aim2)
			{
			istringstream parser(line);
			if(parser.good())
				{
				std::string type, mol_name;
				unsigned int mol_id;
				double px, py, pz;
				while(parser>>type>>px>>py>>pz>>mol_id>>mol_name)
					{
					m_read_atom_type.push_back(type);
					m_read_atom_molid.push_back(mol_id);
					m_read_atom_molname.push_back(mol_name);
					m_read_atom_pos.push_back(vec(px, py, pz));
					}
				}
			else
				{
				cerr << endl << "Unable to parse line, parser.good() failed" << endl << endl;
				throw runtime_error("Error parser(line)");			
				}
			}
		}
	else
		{
		cout << endl << "Warning!!! Can not find position node!" << endl << endl;	
		}	
	file.close();
	cout << "INFO : Position statistics " << m_read_atom_pos.size()<<" atoms  "<<endl;	
	}

void Protein::generateType()
	{
	string type="";
	for(unsigned int j = 0; j < m_N_amino_acid; j++)
		{
		string base = m_sequence[j];
		unsigned int idx = getIndex(base);
		unsigned int np = m_cg_np[idx];
		for(unsigned int t = 0; t < np; t++)
			{
			type+=m_cg_ptype[idx*m_cg_nmax + t]+",";
			}
		}	
	unsigned int st=type.size();
	type=type.substr(0,st-1);
	setParticleTypes(type);		
	}

void Protein::generateSites()
	{
	unsigned int count_cg = 0;
	unsigned int count_atom = 0;

	for(unsigned int j = 0; j < m_N_amino_acid; j++)
		{
		string aa = m_sequence[j];
		unsigned int idx = getIndex(aa);
		unsigned int np = m_cg_np[idx];

		for(unsigned int k = 0; k < np; k++)
			{
			double px = 0.0;
			double py = 0.0;
			double pz = 0.0;
			unsigned int natom = m_natom[idx*m_cg_nmax+k];
			for(unsigned int n = 0; n < natom; n++)
				{
				unsigned int atom_idx = count_atom + n;
				string mol_name = m_read_atom_molname[atom_idx];
				string atom_type = m_read_atom_type[atom_idx];
				unsigned int molid = m_read_atom_molid[atom_idx];
				if((mol_name!=aa)||(molid!=(j+1)))
					{
					cerr << endl <<"Atom_idx "<< atom_idx<<" "<<m_read_atom_pos[atom_idx].x<<" "<<m_read_atom_pos[atom_idx].y<<" "<<m_read_atom_pos[atom_idx].z<<", atom_type "<< atom_type <<", mol_id " << molid << ", unmatched amino acid "<< mol_name<< " from sequence " << j<<" "<<aa<< endl << endl;
					throw runtime_error("Error parser(line)");	
					}

				px += m_read_atom_pos[atom_idx].x;
				py += m_read_atom_pos[atom_idx].y;
				pz += m_read_atom_pos[atom_idx].z;
				}
				
			unsigned int extra = 0;
			if(k==np-1)
				extra = 10;
			unsigned int extra_natom = 0;	
			for(unsigned int n = natom; n < extra+natom; n++)
					{
					unsigned int atom_idx = count_atom + n;
					if (atom_idx<m_read_atom_molname.size())
					{
					string mol_name = m_read_atom_molname[atom_idx];
					string atom_type = m_read_atom_type[atom_idx];
					unsigned int molid = m_read_atom_molid[atom_idx];
					if((mol_name==aa)&&(molid==(j+1)))
						{
						extra_natom += 1;
						}					
					}
				}	

			if(extra_natom>0)
				cout<<extra_natom<<" extra atoms for mol "<<m_read_atom_molid[count_atom]<<endl;
			
			px /= double(natom+extra_natom);
			py /= double(natom+extra_natom);
			pz /= double(natom+extra_natom);
			m_xyz[count_cg]=vec(px/10.0, py/10.0, pz/10.0);
			count_cg += 1;
			count_atom += natom + extra_natom;
			}
		}
	}

void Protein::generateTopology()
	{
//generate topology information
	vector<Bond> bonds;
	vector<Angle> angles;
	vector<Dihedral> dihedrals;
	int d0 = -1;
	int d1 = -1;
	int d2 = -1;
	int offset = 0;
	
	for(unsigned int j = 0; j < m_N_amino_acid; j++)
		{
		string aa = m_sequence[j];
		unsigned int np = m_cg_np[getIndex(aa)];
		
		if(d2>=0)
			bonds.push_back(Bond("B-B", d2, offset));
		if(d2>=0&&d1>=0)
			angles.push_back(Angle("B-B-B", d1, d2, offset));
		if(d2>=0&&d1>=0&&d0>=0)
			dihedrals.push_back(Dihedral("B-B-B", d0, d1, d2, offset));		

		if (np==2)
			{
			bonds.push_back(Bond("S-S", offset, offset+1));				
			}			
		else if (np==3)
			{
			bonds.push_back(Bond("S-S", offset, offset+1));
			bonds.push_back(Bond("S-S", offset+1, offset+2));
			}
		else if (np==4)
			{
			bonds.push_back(Bond("S-S", offset, offset+1));
			bonds.push_back(Bond("S-S", offset+1, offset+2));
			bonds.push_back(Bond("S-S", offset+1, offset+3));
			bonds.push_back(Bond("S-S", offset+2, offset+3));
			}
		else if (np==5)
			{
			bonds.push_back(Bond("B-B", d2, offset));
			bonds.push_back(Bond("S-S", offset, offset+1));
			bonds.push_back(Bond("S-S", offset+1, offset+2));
			bonds.push_back(Bond("S-S", offset+1, offset+3));
			bonds.push_back(Bond("S-S", offset+2, offset+4));
			bonds.push_back(Bond("S-S", offset+3, offset+4));
			}

		d0 = d1;
		d1 = d2;
		d2 = offset;
		offset += np;
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
	
	
void Protein::assignTypes()
	{
	for(unsigned int i=0; i<m_bond.size(); i++)
		{
		if(m_bond[i].type=="S-S")
			{
			double dxab = m_xyz[m_bond[i].a].x - m_xyz[m_bond[i].b].x;
			double dyab = m_xyz[m_bond[i].a].y - m_xyz[m_bond[i].b].y;		
			double dzab = m_xyz[m_bond[i].a].z - m_xyz[m_bond[i].b].z;		
	
	
			double rsqab = dxab*dxab+dyab*dyab+dzab*dzab;
			double rab = sqrt(rsqab);
			if(rab<0.2)
				m_bond[i].type="S-S1";
			else
				m_bond[i].type="S-S2";
			}
		}

	for(unsigned int i=0; i<m_angle.size();i++)
		{
		double dxab = m_xyz[m_angle[i].a].x - m_xyz[m_angle[i].b].x;
		double dyab = m_xyz[m_angle[i].a].y - m_xyz[m_angle[i].b].y;		
		double dzab = m_xyz[m_angle[i].a].z - m_xyz[m_angle[i].b].z;

		double dxcb = m_xyz[m_angle[i].c].x - m_xyz[m_angle[i].b].x;
		double dycb = m_xyz[m_angle[i].c].y - m_xyz[m_angle[i].b].y;		
		double dzcb = m_xyz[m_angle[i].c].z - m_xyz[m_angle[i].b].z;		


        double rsqab = dxab*dxab+dyab*dyab+dzab*dzab;
        double rab = sqrt(rsqab);
        double rsqcb = dxcb*dxcb+dycb*dycb+dzcb*dzcb;
        double rcb = sqrt(rsqcb);		
		
        double c_abbc = dxab*dxcb+dyab*dycb+dzab*dzcb;
        c_abbc /= rab*rcb;
        if (c_abbc > 1.0) c_abbc = 1.0;
        if (c_abbc < -1.0) c_abbc = -1.0;		
		double th = acos(c_abbc);
		if(th<2.0)
			m_angle[i].type="B-B-B-helix";
		else if(th<2.3)
			m_angle[i].type="B-B-B-coil";
		else
			m_angle[i].type="B-B-B-extended";
		}		
		
		
	for(unsigned int i=0; i<m_dihedral.size();i++)
		{
		double dxab = m_xyz[m_dihedral[i].a].x - m_xyz[m_dihedral[i].b].x;
		double dyab = m_xyz[m_dihedral[i].a].y - m_xyz[m_dihedral[i].b].y;		
		double dzab = m_xyz[m_dihedral[i].a].z - m_xyz[m_dihedral[i].b].z;

		double dxcb = m_xyz[m_dihedral[i].c].x - m_xyz[m_dihedral[i].b].x;
		double dycb = m_xyz[m_dihedral[i].c].y - m_xyz[m_dihedral[i].b].y;		
		double dzcb = m_xyz[m_dihedral[i].c].z - m_xyz[m_dihedral[i].b].z;	

		double dxdc = m_xyz[m_dihedral[i].d].x - m_xyz[m_dihedral[i].c].x;
		double dydc = m_xyz[m_dihedral[i].d].y - m_xyz[m_dihedral[i].c].y;		
		double dzdc = m_xyz[m_dihedral[i].d].z - m_xyz[m_dihedral[i].c].z;			

        double dxcbm = -dxcb;
        double dycbm = -dycb;
        double dzcbm = -dzcb;

        double aax = dyab*dzcbm - dzab*dycbm;
        double aay = dzab*dxcbm - dxab*dzcbm;
        double aaz = dxab*dycbm - dyab*dxcbm;
        
        double bbx = dydc*dzcbm - dzdc*dycbm;
        double bby = dzdc*dxcbm - dxdc*dzcbm;
        double bbz = dxdc*dycbm - dydc*dxcbm;
        
        double raasq = aax*aax + aay*aay + aaz*aaz;
        double rbbsq = bbx*bbx + bby*bby + bbz*bbz;
        double rgsq = dxcbm*dxcbm + dycbm*dycbm + dzcbm*dzcbm;
        double rg = sqrt(rgsq);
        
        double  raa2inv, rbb2inv;
        raa2inv = rbb2inv = 0.0f;
        if (raasq > 0.0f) raa2inv = 1.0f/raasq;
        if (rbbsq > 0.0f) rbb2inv = 1.0f/rbbsq;
        double rabinv = sqrt(raa2inv*rbb2inv);
        
        double c_abcd = (aax*bbx + aay*bby + aaz*bbz)*rabinv;
        double s_abcd = rg*rabinv*(aax*dxdc + aay*dydc + aaz*dzdc);
        
        if (c_abcd > 1.0f) c_abcd = 1.0f;
        if (c_abcd < -1.0f) c_abcd = -1.0f;

		double th = atan2f(s_abcd,c_abcd);		
		if (th<0)
			th += 2.0*M_PI;
		
		if(th<2.0)
			m_dihedral[i].type="B-B-B-B-helix";
		else
			m_dihedral[i].type="B-B-B-B-extended";
		}		
	
	}
	
void Protein::generate()
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
				}
			else
				{
				m_be_generated[i]= false;
				}
			}
		m_firststep = false;
		}
	
	generateSites();
	assignTypes();

	}

void export_Protein(pybind11::module &m)
	{
	pybind11::class_<Protein, Molecule, std::shared_ptr<Protein> >(m, "Protein")
		.def(pybind11::init<const std::string& >())
		.def(pybind11::init<const std::string&, Protein::Model >())
		;
    pybind11::enum_<Protein::Model>(m, "Model")
    .value("one_bead_one_amino_acid",Protein::Model::one_bead_one_amino_acid)
	.value("martini",Protein::Model::martini)
	.export_values()
	;
	}
