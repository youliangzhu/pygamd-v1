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

#include <cuda_runtime.h>
#include "MolInfo.h"
#include <unordered_map>

#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

class Function
	{
    public:
        Function() {};
        virtual ~Function() {};	
		void add(XMLBuilder* build, MolInfo* mol)
			{
			m_build = build;
			m_mol = mol;
			};
		void add_datalog(std::string datalog)
			{
            m_datalog = std::move(datalog);
			}
			
		virtual void compute(){};
		XMLBuilder* m_build;		
		MolInfo* m_mol;
		std::string m_datalog;
	};

//--- case 1
class Rg2 : public Function
    {
    public:
        Rg2(std::string filename): Function()
			{
			m_Nf =0;
			m_file.open(filename.c_str(), ios::app);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Rg2 dump");
				}
			}
		virtual ~Rg2() 
			{
			for(unsigned int i=0; i<m_av_rg2.size();i++)
				cout<<"Mol"<<i<<" rg2= "<<m_av_rg2[i]/double(m_Nf)<<endl;
			};
		virtual void compute();		
    private:
		std::ofstream m_file;
		std::vector<double> m_av_rg2;
		unsigned int m_Nf;
	};
//--- case 2	
class Ed2 : public Function
    {
    public:
        Ed2(std::string filename): Function()
			{
			m_Nf =0;			
			m_file.open(filename.c_str(), ios::app);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Ed2 dump");
				}
			}
		virtual ~Ed2()
			{
			for(unsigned int i=0; i<m_av_ed2.size();i++)
				cout<<"Mol"<<i<<" ed2 = "<<m_av_ed2[i]/double(m_Nf)<<endl;
			};		

		virtual void compute();			
    private:
		std::ofstream m_file;
		std::vector<double> m_av_ed2;
		unsigned int m_Nf;		
	};
//--- case 3
class RDF : public Function
    {
    public:
        RDF(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error RDF dump");
				}
			m_maxbin =100;
			m_block_size = 256;
			m_gpu_id = 0;
			m_rdf.resize(m_maxbin);
			m_r.resize(m_maxbin);
			m_Nf=0;
			m_rmax=0.0;
			m_exclusion_mol = false;
			m_exclusion_list = false;
			m_exclusion_angle = false;
			m_exclusion_bond = false;			
			}
		virtual ~RDF() 
			{
			for(unsigned int i=0; i<m_rdf.size(); i++)
				m_file<<m_r[i]<<"  "<<m_rdf[i]/double(m_Nf)<<"\n";
			m_file.close();
			}
		void setMaxbin(unsigned int maxbin)
			{
			m_maxbin=maxbin;
			m_rdf.clear();
			m_rdf.resize(m_maxbin);
			m_r.clear();
			m_r.resize(m_maxbin);
			}
		void setRmax(double rmax)
			{
			m_rmax = rmax;
			}
		void setGPU(unsigned int gpu_id)
			{
			m_gpu_id = gpu_id;
			cudaSetDevice(m_gpu_id);
			}
		void setBondEx(bool bond_ex)
			{
			m_exclusion_bond = bond_ex;
			}			
		void setAngleEx(bool angle_ex)
			{
			m_exclusion_angle = angle_ex;
			}	
		void setMolEx(bool mol_ex)
			{
			m_exclusion_mol = mol_ex;
			}				
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_maxbin;
		unsigned int m_block_size;
		int m_gpu_id;
		std::vector<double> m_rdf;
		std::vector<double> m_r;		
		unsigned int m_Nf;
		double m_rmax;
		bool m_exclusion_mol;
		bool m_exclusion_list;
		bool m_exclusion_angle;
		bool m_exclusion_bond;
		unsigned int* d_n_exclusion;
		unsigned int* d_exclusion_list;
		unsigned int* d_mol_id_per_particle;
	};
//--- case 4
class Bond_distr : public Function
    {
    public:
        Bond_distr(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Bond_distr dump");
				}
			m_Npot = 2001;
			m_Nf=0;			
			}
		virtual ~Bond_distr() 
			{
			for(unsigned int i=0; i<m_Nb;i++)
				{
				m_file<<m_bondMap[i]<<endl;
				for(unsigned int j=0;j<m_Npot;j++)
					{
					if(m_distr[i*m_Npot+j]>0)
						{
						double value = m_distr[i*m_Npot+j]/double(m_Nf);
						m_file<<double(j)*m_delt<<"  "<<value/m_delt<<"\n";
						}
					}
				cout<<"The averaged length of bond "<<m_bondMap[i]<<" is "<<m_bond_lenth[i]/double(m_Nf)<<endl;					
				}
			m_file.close();				
			}
		void setNpot(unsigned int npot)
			{
			m_Npot=npot;
			}
		virtual void compute();			
    private:
		std::ofstream m_file;
		double m_rcut;
		unsigned int m_Npot;
		unsigned int m_Nf;
		unsigned int m_Nb;
		double m_delt;		
		std::vector<double> m_distr;
		std::vector<double> m_bond_lenth;		
		std::vector<std::string> m_bondMap;
	};
//--- case 5	
class Angle_distr : public Function
    {
    public:
        Angle_distr(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Angle_distr dump");
				}
			m_Npot = 2001;	
			m_Nf=0;				
			}
		virtual ~Angle_distr() 
			{
			for(unsigned int i=0; i<m_Nb;i++)
				{
				m_file<<m_angleMap[i]<<endl;
				for(unsigned int j=0;j<m_Npot;j++)
					{
					if(m_distr[i*m_Npot+j]>0)
						{
						double value = m_distr[i*m_Npot+j]/double(m_Nf);
						m_file<<double(j)*m_delt<<"  "<<value/m_delt<<"\n";
						}
					}
				cout<<"The averaged radian of angle "<<m_angleMap[i]<<" is "<<m_angle_radian[i]/double(m_Nf)<<endl;							
				}
			m_file.close();						
			}
		void setNpot(unsigned int npot)
			{
			m_Npot=npot;
			}
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_Npot;
		unsigned int m_Nf;
		unsigned int m_Nb;
		double m_delt;
		std::vector<double > m_distr;
		std::vector<double> m_angle_radian;		
		std::vector<std::string> m_angleMap;		
	};
//--- case 6
class Dihedral_distr : public Function
    {
    public:
        Dihedral_distr(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Dihedral_distr dump");
				}
			m_Npot = 2001;	
			m_Nf=0;					
			}
		virtual ~Dihedral_distr() 
			{
			for(unsigned int i=0; i<m_Nb;i++)
				{
				m_file<<m_dihedralMap[i]<<endl;
				for(unsigned int j=0;j<m_Npot;j++)
					{
					if(m_distr[i*m_Npot+j]>0)
						{
						double value = m_distr[i*m_Npot+j]/double(m_Nf);
						m_file<<double(j)*m_delt<<"  "<<value/m_delt<<"\n";
						}
					}
				cout<<"The averaged radian of dihedral "<<m_dihedralMap[i]<<" is "<<m_dihedral_radian[i]/double(m_Nf)<<endl;					
				}
			m_file.close();						
			}
		void setNpot(unsigned int npot)
			{
			m_Npot=npot;
			}
		virtual void compute();			
    private:		
		std::ofstream m_file;
		unsigned int m_Npot;
		unsigned int m_Nf;
		unsigned int m_Nb;
		double m_delt;	
		std::vector<double > m_distr;
		std::vector<double> m_dihedral_radian;			
		std::vector<std::string> m_dihedralMap;
	};
//--- case 7
class StressTensor : public Function
    {
    public:
        StressTensor(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios::app);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error StressTensor dump");
				}
			m_Nf=0;
			m_bondex=true;
			m_bodyex=true;
			m_diameter = true;
			m_rcut=3.0;
			m_delt=0.0;
			}
		virtual ~StressTensor() {};
		void setParam();	
		virtual void compute();
		void setBondEx(bool bondex)
			{
			m_bondex = bondex;
			}
		void setBodyEx(bool bodyex)
			{
			m_bodyex = bodyex;
			}
		void setDiameterConsider(bool diameter)
			{
			m_diameter = diameter;
			}
    private:
		std::ofstream m_file;
		unsigned int m_Nf;
		std::vector<vec> m_pparams;
		std::vector<vec> m_bparams;
		double m_rcut;
		double m_delt;
		bool m_bondex;
		bool m_bodyex;
		bool m_diameter;
	};
//--- case 8
class Density : public Function
    {
    public:
        Density(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios::app);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Density dump");
				}
			m_Nf=0;
			}
		virtual ~Density() {};

		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_Nf;
	};	
//--- case 9
class Reimage : public Function
    {
    public:
        Reimage(): Function()
			{
			m_solute_center = false;
			m_label_free_particle = false;
			m_unwrap_molecule = true;
			m_keep_molecule_center_in_box = false;
			m_convert_constraints_to_bonds = false;
			m_remove_bond_cross_box = false;
			m_remove_image = false;
			m_body_keep = false;
			m_target_type = "";
			m_shiftx = 0.0;
			m_shifty = 0.0;
			m_shiftz = 0.0;
			m_nprecision = 10;
			m_nhead = 7;
			m_add_image_to_pos = true;
			m_image_integrate = false;
			m_Nf=0;
			m_sp = ios::beg;
			m_file="xxxxxxxx";
			}
		virtual ~Reimage() {};
		void setShiftX(double shiftx)
			{
			m_shiftx = shiftx;
			}
		void setShiftY(double shifty)
			{
			m_shifty = shifty;
			}
		void setShiftZ(double shiftz)
			{
			m_shiftz = shiftz;
			}
		void setLabelFreeParticle(std::string target_type)
			{
			m_target_type = target_type;
			m_label_free_particle = true;
			}
		void setUnwrapMolecule(bool unwrap_molecule)
			{
			m_unwrap_molecule = unwrap_molecule;
			}
		void setMoleculeCenterInBox(bool keep_molecule_center_in_box)
			{
			m_keep_molecule_center_in_box = keep_molecule_center_in_box;
			}	
		void setRemoveImage(bool remove_image)
			{
			m_remove_image = remove_image;
			m_add_image_to_pos=false;
			}
		void setConvertConstraintsToBonds(bool convert_constraints_to_bonds)
			{
			m_convert_constraints_to_bonds = convert_constraints_to_bonds;
			}
		void setRemoveBondCrossBox(bool remove_bond_cross_box)
			{
			m_remove_bond_cross_box = remove_bond_cross_box;
			}
		void addImageToPos(bool add_image_to_pos)
			{
			m_add_image_to_pos = add_image_to_pos;
			}
		void setBodyKeep(bool body_keep)
			{
			m_body_keep = body_keep;
			}
		void setImageIntegrate(bool image_integrate)
			{
			m_image_integrate = image_integrate;
			}				
		virtual void compute();			
    private:
		bool m_solute_center;
		bool m_label_free_particle;
		bool m_unwrap_molecule;
		bool m_keep_molecule_center_in_box;
		bool m_remove_image;
		bool m_convert_constraints_to_bonds;
		bool m_remove_bond_cross_box;
		bool m_add_image_to_pos;
		bool m_body_keep;
		bool m_image_integrate;
		double m_shiftx;
		double m_shifty;
		double m_shiftz;
		unsigned int m_nprecision;
		unsigned int m_nhead;
		std::string m_target_type;
		unsigned int m_Nf;
		ifstream::pos_type m_sp;
		std::string m_file;
	};	
//--- case 10
class MSD : public Function
    {
    public:
        MSD(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error MSD dump");
				}
			m_Nf=0;
			}
		virtual ~MSD() {};
		void setParam();
		void setDirection(std::string direction)
			{
			m_direction=direction;
			}
		virtual void compute();			
    private:
		unsigned int m_Nf;
		std::string m_direction;
		std::ofstream m_file;
		std::vector<vec> m_pos_offset;
		std::vector<vec> m_pos_cm;		
	};
//--- case 11
class RDFCM : public Function
    {
    public:
        RDFCM(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error RDFCM dump");
				}
			m_maxbin =100;
			m_block_size = 256;
			m_gpu_id = 0;
			m_rdf.resize(m_maxbin);
			m_r.resize(m_maxbin);			
			m_Nf=0;
			m_rmax = 0.0;
			m_exclusion_mol = false;
			m_exclusion_list = false;
			m_exclusion_angle = false;
			m_exclusion_bond = false;		
			}
		virtual ~RDFCM() 
			{
			for(unsigned int i=0; i<m_rdf.size(); i++)
				m_file<<m_r[i]<<"  "<<m_rdf[i]/double(m_Nf)<<"\n";
			m_file.close();
			}
		void setMaxbin(unsigned int maxbin)
			{
			m_maxbin=maxbin;
			m_rdf.clear();
			m_rdf.resize(m_maxbin);
			m_r.clear();
			m_r.resize(m_maxbin);
			}
		void setRmax(double rmax)
			{
			m_rmax = rmax;
			}			
		void setGPU(unsigned int gpu_id)
			{
			m_gpu_id = gpu_id;
			cudaSetDevice(m_gpu_id);
			}
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_maxbin;
		unsigned int m_block_size;
		int m_gpu_id;
		std::vector<double> m_rdf;
		std::vector<double> m_r;		
		unsigned int m_Nf;
		double m_rmax;
		bool m_exclusion_mol;
		bool m_exclusion_list;
		bool m_exclusion_angle;
		bool m_exclusion_bond;	
		unsigned int* d_n_exclusion;
		unsigned int* d_exclusion_list;
		unsigned int* d_mol_id_per_particle;		
	};	
//--- case 12
class MSDCM : public Function
    {
    public:
        MSDCM(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error MSDCM dump");
				}
			m_Nf=0;
			}
		virtual ~MSDCM();
		virtual void compute();	
		void setDirection(std::string direction)
			{
			m_direction=direction;
			}		
    private:
		unsigned int m_Nf;
		std::string m_direction;		
		std::ofstream m_file;
		std::vector<std::vector<vec> > m_pos_all;
		std::vector<unsigned int > m_mol_type_id;
		unsigned int m_n_kind_mol;	
	};

//--- case 13
class Entanglement : public Function
    {
    public:
        Entanglement(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Entanglement dump");
				}
			m_Nf=0;
			m_rcut_max = 0.0;
			m_totalnum = 0;
			m_fdelt = 1.0;
			m_fmax = 10000.0;
			m_fdistr.resize((unsigned int)(m_fmax/m_fdelt));
			}
		virtual ~Entanglement() 
			{
			for(unsigned int i=0; i<m_fdistr.size();i++)
				{
				if(m_fdistr[i]>0)
					{
					double distr = double(m_fdistr[i])/(double(m_totalnum)*m_fdelt);
					m_file<<(double(i)+0.5)*m_fdelt<<"  "<<distr<<"  "<<m_fdistr[i]<<endl;
					}
				}
			m_file.close();
			};
		void setParam();			
		virtual void compute();			
    private:
		unsigned int m_Nf;
		std::ofstream m_file;
		std::vector<double> m_rcutsq;
		std::vector<unsigned int> m_fdistr;
		unsigned int m_totalnum;
		double m_fdelt;
		double m_fmax;
		std::vector<vec> m_params;		
		double m_rcut_max;
	};	
//--- case 14
class STRFAC : public Function
    {
    public:
        STRFAC(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error STRFAC dump");
				}
			m_Kmax=80;
			m_LKmax=0;
			m_MKmax=0;
			m_NKmax=0;			
			m_Nf=0;
			m_Ntypes=0;
			m_block_size = 256;
			m_gpu_id = 0;
			m_L = 0.0;
			m_Qmax = 0.0;
			m_direction = "XYZ";
			m_2D=false;
			}
		virtual ~STRFAC() 
			{
			unsigned int Ksqmax = m_Kmax*m_Kmax;
			if(m_2D)
				{
				m_file<<"qx"<<"  "<<"qy"<<"  ";
				for(unsigned int typi =0; typi <m_Ntypes; typi++)
					{
					for(unsigned int typj =typi; typj <m_Ntypes; typj++)
						{
						m_file<<m_type_map[typi]+"-"+m_type_map[typj]+"_r"<<"  "<<m_type_map[typi]+"-"+m_type_map[typj]+"_i"<<"  ";
						}
					}
				m_file<<"\n";							

				for(unsigned int ii =1; ii<(unsigned int)(m_LKmax+1); ii++)
					{
					unsigned int i = m_LKmax + 1 - ii;
					double qi =double(i)*M_PI*2.0/m_L;						
					for(unsigned int jj =0; jj<(unsigned int)(2*m_MKmax+1); jj++)
						{
						unsigned int j = 2*m_MKmax -jj;
						double qj =double(int(j)-int(m_MKmax))*M_PI*2.0/m_L;
						m_file<<-qi<<"   "<<-qj<<"   ";
						unsigned int icn =0;
						for(unsigned int typi =0; typi < m_Ntypes; typi++)
							{	
							for(unsigned int typj =typi; typj < m_Ntypes; typj++)
								{	
								unsigned int idi = i+icn*(m_LKmax+1);
								unsigned int idj = j+icn*(2*m_MKmax+1);
								m_file<<m_sqreal2D[idi][idj]/double(m_Nf)<<"   "<<m_sqimage2D[idi][idj]/double(m_Nf)<<"  ";
								icn +=1;
								}
							}
						m_file<<"\n";	
						}
					m_file<<"\n";						
					}

				for(unsigned int i =0; i<(unsigned int)(m_LKmax+1); i++)
					{
					double qi =double(i)*M_PI*2.0/m_L;						
					for(unsigned int j =0; j<(unsigned int)(2*m_MKmax+1); j++)
						{
						double qj =double(int(j)-int(m_MKmax))*M_PI*2.0/m_L;
						m_file<<qi<<"   "<<qj<<"   ";
						unsigned int icn =0;
						for(unsigned int typi =0; typi < m_Ntypes; typi++)
							{	
							for(unsigned int typj =typi; typj < m_Ntypes; typj++)
								{	
								unsigned int idi = i+icn*(m_LKmax+1);
								unsigned int idj = j+icn*(2*m_MKmax+1);
								if(i==0&&int(j)<m_MKmax)
									idj = 2*m_MKmax-j+icn*(2*m_MKmax+1);
								m_file<<m_sqreal2D[idi][idj]/double(m_Nf)<<"   "<<m_sqimage2D[idi][idj]/double(m_Nf)<<"  ";
								icn +=1;
								}
							}
						m_file<<"\n";	
						}
					m_file<<"\n";	
					}						
				}
			else
				{
				m_file<<"q"<<"  ";
				for(unsigned int typi =0; typi <m_Ntypes; typi++)
					{
					for(unsigned int typj =typi; typj <m_Ntypes; typj++)
						{
						m_file<<m_type_map[typi]+"-"+m_type_map[typj]+"_r"<<"  "<<m_type_map[typi]+"-"+m_type_map[typj]+"_i"<<"  ";
						}
					}
				m_file<<"\n";					
					
				for(unsigned int j =0; j<Ksqmax; j++)
					{
					unsigned int icn = 0;
					float sumx=0.0;
					float sumy=0.0;
					for(unsigned int typi =0; typi < m_Ntypes; typi++)
						{
						for(unsigned int typj =typi; typj < m_Ntypes; typj++)
							{								
							unsigned int id = icn*Ksqmax + j;
							sumx += m_sqreal[id]/double(m_Nf);
							sumy += m_sqimage[id]/double(m_Nf);
							icn +=1;
							}
						}
					if(sumx!=0.0||sumy!=0.0)
						{						
						unsigned int icn = 0;
						double q =sqrt(double(j+1))*M_PI*2.0/m_L;
						m_file<<q<<"  ";
						for(unsigned int typi =0; typi < m_Ntypes; typi++)
							{
							for(unsigned int typj =typi; typj < m_Ntypes; typj++)
								{								
								unsigned int id = icn*Ksqmax + j;
								m_file<<m_sqreal[id]/double(m_Nf)<<"  "<<m_sqimage[id]/double(m_Nf)<<"  ";
								icn +=1;
								}
							}
						m_file<<"\n";
						}					
					}						
				}
			m_file.close();
			}
		void setQmax(float qmax)
			{
			m_Qmax=qmax;
			}
		void setGPU(unsigned int gpu_id)
			{
			m_gpu_id = gpu_id;
			cudaSetDevice(m_gpu_id);
			}
		void setDeltaq(float deltaq)
			{
			m_L = M_PI*2.0/deltaq;
			}
		void setDirection(string direction)
			{
			m_direction=direction;
			}	
		void set2D(bool d2)
			{
			m_2D=d2;
			}				
		virtual void compute();			
    private:
		std::ofstream m_file;
		int m_Kmax;
		int m_LKmax;
		int m_MKmax;
		int m_NKmax;		
		float m_Qmax;
		std::vector<double> m_sqreal;
		std::vector<double> m_sqimage;
		std::vector<std::vector<double> > m_sqreal2D;
		std::vector<std::vector<double> > m_sqimage2D;

		std::vector<string> m_type_map;
		unsigned int m_Ntypes;
		unsigned int m_Nf;
		unsigned int m_gpu_id;
		unsigned int m_block_size;
		float m_L;
		bool m_2D;
		string m_direction;		
	};
//--- case 15
class DOMAINSIZE : public Function
    {
    public:
        DOMAINSIZE(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error DOMAINSIZE dump");
				}
			m_Kmax=20;
			m_Nf=0;
			m_block_size = 256;
			m_gpu_id = 0;
			m_qc = 0.4690;     //scf=0.3708;0.4690 //dpd=0.6188;0.5809
			}
		virtual ~DOMAINSIZE() 
			{
			m_file.close();
			}
		void setKmax(unsigned int kmax)
			{
			m_Kmax=kmax;
			}
		void setQc(float qc)
			{
			m_qc=qc;
			}
		
		void setGPU(unsigned int gpu_id)
			{
			m_gpu_id = gpu_id;
			cudaSetDevice(m_gpu_id);
			}
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_Kmax;
		unsigned int m_Nf;
		unsigned int m_gpu_id;
		unsigned int m_block_size;
		float m_qc;
	};
//--- case 16
class DSTRFAC : public Function
    {
    public:
        DSTRFAC(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error DSTRFAC dump");
				}
			m_Kmax=vec_int(0, 0, 0);
			m_q=7.0;
			m_Nf=0;
			}
		virtual ~DSTRFAC() {};
		void setParam();
		virtual void compute();	
		void setKmax(unsigned int kmax)
			{
			m_Kmax=vec_int(kmax, kmax, kmax);
			}
		void setQ(double q)
			{
			m_q=q;
			}
    private:
		unsigned int m_Nf;
		vec_int m_Kmax;
		double m_q;
		std::vector<vec_int> m_qvec;	
		std::ofstream m_file;
		std::vector<vec> m_pos_offset;
		std::vector<vec> m_pos_cm;		
	};
//--- case 17
class ConfigCheck : public Function
    {
    public:
        ConfigCheck(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error ConfigCheck dump");
				}
			m_Nf=0;
			m_bondex=true;
			m_angleex=true;
			m_dihedralex=true;
			m_bodyex=true;
			m_rcut=2.0;
			}
		virtual ~ConfigCheck() {};
		virtual void compute();
		void setBondEx(bool bondex)
			{
			m_bondex = bondex;
			}
		void setAngleEx(bool angleex)
			{
			m_angleex = angleex;
			}
		void setDihedralEx(bool dihedralex)
			{
			m_dihedralex = dihedralex;
			}
		void setBodyEx(bool bodyex)
			{
			m_bodyex = bodyex;
			}
		void setRcut(double rcut)
			{
			m_rcut = rcut;
			}
    private:
		unsigned int m_Nf;	
		std::ofstream m_file;
		double m_rcut;
		bool m_bondex;
		bool m_bodyex;
		bool m_angleex;
		bool m_dihedralex;
	};
	
//--- case 18
class RDFBetweenTypes : public Function
    {
    public:
        RDFBetweenTypes(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error RDFBetweenTypes dump");
				}
			m_maxbin =100;
			m_block_size = 256;
			m_gpu_id = 0;
			m_rdf.resize(m_maxbin);
			m_r.resize(m_maxbin);
			m_Nf=0;
			m_rmax = 0.0;
			m_exclusion_mol = false;
			m_exclusion_list = false;
			m_exclusion_angle = false;
			m_exclusion_bond = false;	
			}
		virtual ~RDFBetweenTypes() 
			{
			m_file<<"r"<<"  ";
			for(unsigned int typi =0; typi <m_Ntype; typi++)
				{
				for(unsigned int typj =typi; typj <m_Ntype; typj++)
					{
					m_file<<m_type_map[typi]+"-"+m_type_map[typj]<<"  ";
					}
				}
			m_file<<"\n";
			for (unsigned int bin = 0; bin < m_maxbin; bin++ )
				{ 
				double r = m_r[bin];
				m_file<<r<<"  ";
				for(unsigned int typi =0; typi <m_Ntype; typi++)
					{
					for(unsigned int typj =typi; typj <m_Ntype; typj++)
						{
						double g = m_rdf[(typi*m_Ntype+typj)*m_maxbin+bin]/double(m_Nf);
						m_file<<g<<"  ";
						}  
					}
				m_file<<"\n";
				}
			m_file.close();
			}
		void setMaxbin(unsigned int maxbin)
			{
			m_maxbin=maxbin;
			m_rdf.clear();
			m_rdf.resize(m_maxbin);
			m_r.clear();
			m_r.resize(m_maxbin);
			}
		void setRmax(double rmax)
			{
			m_rmax = rmax;
			}				
		void setGPU(unsigned int gpu_id)
			{
			m_gpu_id = gpu_id;
			cudaSetDevice(m_gpu_id);
			}
		void setBondEx(bool bond_ex)
			{
			m_exclusion_bond = bond_ex;
			}			
		void setAngleEx(bool angle_ex)
			{
			m_exclusion_angle = angle_ex;
			}	
		void setMolEx(bool mol_ex)
			{
			m_exclusion_mol = mol_ex;
			}			
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_maxbin;
		unsigned int m_block_size;
		int m_gpu_id;
		std::vector<double> m_rdf;
		std::vector<double> m_r;		
		unsigned int m_Nf;
		std::vector<std::string> m_type_map;
		unsigned int m_Ntype;
		double m_rmax;
		bool m_exclusion_mol;
		bool m_exclusion_list;
		bool m_exclusion_angle;
		bool m_exclusion_bond;
		unsigned int* d_n_exclusion;
		unsigned int* d_exclusion_list;
		unsigned int* d_mol_id_per_particle;		
	};
//--- case 19
class FileConversion : public Function
    {
    public:
        FileConversion(): Function()
			{
			m_lammps=false;
			m_gromacs = false;
			m_xml=false;
			m_Nf = 0;
			m_nprecision=10;
			m_nhead = 7;
			}
		virtual ~FileConversion() 
			{
			}
		void setLammps(bool lammps)
			{
			m_lammps = lammps;
			}
		void setGromacs(bool gromacs)
			{
			m_gromacs = gromacs;
			}
		void setXML(bool xml)
			{
			m_xml = xml;
			}		
		virtual void compute();			
    private:
		std::ofstream m_file;
		bool m_lammps;
		bool m_gromacs;
		bool m_xml;
		unsigned int m_Nf;
		unsigned int m_nprecision;
		unsigned int m_nhead;		
	};	
//--- case 20	
class PatchToParticle : public Function
    {
    public:
        PatchToParticle(): Function()
			{
			m_nprecision=10;
			m_nhead = 7;
			m_separatedis = 0.1;
			m_scale = 0.97;
			m_filter_sphere = false;
			m_Nf = 0;
			}
		virtual ~PatchToParticle() 
			{
			}
		void setSeparateDis(double sd)
			{
			m_separatedis = sd;
			}
		void setFilterSphere(bool fs)
			{
			m_filter_sphere = fs;
			}
		void setPatchParticleScale(double sc)
			{
			m_scale = sc;
			}			
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_nprecision;
		unsigned int m_nhead;
		double m_separatedis;
		double m_scale;
		bool m_filter_sphere;
		unsigned int m_Nf;
	};		
//--- case 21 
class SSF : public Function
    {
    public:
        SSF(std::string filename): Function()
			{
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error SSF dump");
				}
            m_Nf=0;
			m_qnummax=40;
			}
        virtual ~SSF()
			{
			std::vector<double> q_all_sum(q_all[0].size()*q_all[0].size(), 0.0);
		    std::vector <double> ssf_all_sum(q_all[0].size()*q_all[0].size(), 0.0);
			for(unsigned int i=0; i<q_all[0].size(); i++)
				{	
				for(unsigned int j=0; j<m_Nf; j++)
					{
					q_all_sum[i] += q_all[j][i];
					ssf_all_sum[i] += ssf_all[j][i];
					}
				m_file << q_all_sum[i]/m_Nf << "  " << ssf_all_sum[i]/m_Nf << endl;					 
				}					
			cout << "21. Good Luck! Outputting results of the static structure factor (SSF) to 'ssf.log'." << endl;		
			};
		void setqnummax(unsigned int qnummax)
			{
			m_qnummax=qnummax;
			}
        virtual void compute();
    private:
        unsigned int m_Nf, pos_size, m_qnummax;
		double deltat;
        std::ofstream m_file;
		std::vector<std::vector<double> > q_all;
		std::vector<std::vector<double> > ssf_all;
    };

//--- case 22                                                                                                                                   
class ADF : public Function
    {
    public:
        ADF(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if(!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error ADF dump");
				}
			m_Nf=0;
			maxbin=200;
			m_rcut=1.0;
			m_beta=90.0;
			}
		virtual ~ADF() 
			{
			double sum_g_angular=0.0;
			std::vector<double> r_angular(maxbin, 0.0);
		    std::vector<double> g_angular(maxbin, 0.0);
			for(unsigned int i=0; i<maxbin; i++)
				{	
				for(unsigned int j=0; j<m_Nf ; j++)
					{
					r_angular[i] += r_all[j][i];
					g_angular[i] += g_all[j][i];				 
					}
				sum_g_angular += g_angular[i];				 			 
				}
            for(unsigned int i=0; i<maxbin; i++)
				{			
				m_file << r_angular[i]/m_Nf << "  " << g_angular[i]/(sum_g_angular) << endl;	
				}	
			cout << "22. Good Luck! Outputting results of the angular distribution function (ADF) to 'adf.log'." << endl;			
			m_file.close();	
            r_all.clear();
			g_all.clear();
			};
		void setBeta(double beta)
			{
			m_beta=beta;
			}
		void setRcut(double rcut)
			{
			m_rcut=rcut;
			}		
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_Nf;		
		unsigned int maxbin;
		double m_rcut;
		double m_beta;
		std::vector<std::vector<double> > r_all;
		std::vector<std::vector<double> > g_all;
	};

//--- case 23                                                                                                                                    //1
class CND : public Function
    {
    public:
        CND(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if(!m_file.good())
				{
				cerr<<endl<<"***Error! Error opening dump file "<<filename<<endl<<endl;
				throw runtime_error("Error CND dump");
				}
			m_Nf=0;
			maxbin=200;
			m_rcut=1.0;
			m_beta=90.0;
			}
		virtual ~CND() 
			{
			std::vector<unsigned int> score_s(maxbin,0);
		    std::vector <double> N_s(maxbin,0.0);
			for(unsigned int i=1; i<maxbin; i++)
				{	
				for(unsigned int j=0; j<m_Nf ; j++)
					{
					score_s[i] += score_s_all[j][i];
					N_s[i] += N_s_all[j][i];
					}
				m_file << score_s[i]/m_Nf << "  " << N_s[i]/m_Nf << endl;					 
				}				
			cout << "23. Good Luck! Outputting results of the contact number distribution (CND) to 'cnd.log'." << endl;
			m_file.close();	
            score_s_all.clear();
			N_s_all.clear();
			};
		void setBeta(double beta)
			{
			m_beta=beta;
			}
		void setRcut(double rcut)
			{
			m_rcut=rcut;
			}			
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int maxbin;
		unsigned int m_Nf;	
		double m_rcut;
	    double m_beta;
		std::vector<std::vector<unsigned int> > score_s_all;
		std::vector<std::vector <double> > N_s_all;
	};

//--- case 24	
class MSAD : public Function
    {
    public:
        MSAD(std::string filename): Function()
			{
            m_file.open(filename.c_str(), ios_base::out);
            if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file" << filename << endl << endl;
				throw runtime_error("Error MSAD dump");
				}
            m_Nf=0;
			m_dt=0.005;
			}
        virtual ~MSAD()
			{				
            std::vector<double> msad; msad.reserve(m_Nf);
			if(m_Nf<=1000)
				{
				n_ensemble=0.1*m_Nf;	
				}				
            else
				{
				n_ensemble=1000;	
				}
			pos_size=m_Rotangle_all[0].size();	
			deltat=(delta_t[1]-delta_t[0])*m_dt;			
            for(unsigned int i=1; i<m_Nf-n_ensemble; i++)                                                                                   
				{
				msad[i]=0.0;
				unsigned int count=0;
                for(unsigned int j=i; j<i+n_ensemble; j++)                                                
					{
                    for(unsigned int k=0; k<pos_size; k++)
						{
//						if(m_type_all[j-i][k]==0)
							{
							double dxtheta_t=m_Rotangle_all[j][k].x-m_Rotangle_all[j-i][k].x;
							double dytheta_t=m_Rotangle_all[j][k].y-m_Rotangle_all[j-i][k].y;
							double dztheta_t=m_Rotangle_all[j][k].z-m_Rotangle_all[j-i][k].z;
							msad[i] += dxtheta_t*dxtheta_t + dytheta_t*dytheta_t + dztheta_t*dztheta_t;	
							count += 1;
							}
						}
					}				
                msad[i] /= double(count);
				m_file << i*deltat << "  " << msad[i] << endl;
				}	
            cout << "24. Good Luck! Outputting results of the mean square angular displacement (MSAD) to 'msad.log'." << endl;			
			m_file.close();
			m_type_all.clear();
			m_Rotangle_all.clear();	
            delta_t.clear();
			msad.clear();
			};	
		void setDt(double dt)
			{
			m_dt=dt;
			}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat;
		double m_dt;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;
		std::vector< std::vector<unsigned int> > m_type_all;
        std::vector<std::vector<vec> > m_Rotangle_all;
    } ;	
	
//--- case 25	
class RMSAD : public Function
    {
    public:
        RMSAD(std::string filename): Function()
			{
            m_file.open(filename.c_str(), ios_base::out);
            if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file" << filename << endl << endl;
				throw runtime_error("Error RMSAD dump");
				}
            m_Nf=0;
			m_dt=0.005;
			}
        virtual ~RMSAD()
			{				
            std::vector<double> rotationmsd0; rotationmsd0.reserve(m_Nf);
			std::vector<double> rotationmsd1; rotationmsd1.reserve(m_Nf);
			std::vector<double> rotationmsd; rotationmsd.reserve(m_Nf);
			pos_size=m_ori_all[0].size();	
			deltat=(delta_t[1]-delta_t[0])*m_dt;
			if(m_Nf<=1000)
				{
				n_ensemble=0.1*m_Nf;	
				}				
            else
				{
				n_ensemble=1000;	
				}
			std::vector<vec> deltaRotation0_t; deltaRotation0_t.resize(pos_size);
			std::vector<vec> deltaRotation1_t; deltaRotation1_t.resize(pos_size);				
            for(unsigned int j=0; j<1; j++) 
				{   
				for(unsigned int i=0; i<m_Nf; i++) 	
					{
					for(unsigned int k=0; k<pos_size; k++)
						{					
						double dotproduct=m_Rotangle_all[i][k].x*m_ori_all[j][k].x + m_Rotangle_all[i][k].y*m_ori_all[j][k].y + m_Rotangle_all[i][k].z*m_ori_all[j][k].z;					 
						deltaRotation0_t[k].x=dotproduct*m_ori_all[j][k].x;	
						deltaRotation0_t[k].y=dotproduct*m_ori_all[j][k].y;
						deltaRotation0_t[k].z=dotproduct*m_ori_all[j][k].z;
						
						deltaRotation1_t[k].x=m_Rotangle_all[i][k].x - dotproduct*m_ori_all[j][k].x;	
						deltaRotation1_t[k].y=m_Rotangle_all[i][k].y - dotproduct*m_ori_all[j][k].y;
						deltaRotation1_t[k].z=m_Rotangle_all[i][k].z - dotproduct*m_ori_all[j][k].z;					
						}
					m_Rotation0_all.push_back(deltaRotation0_t);
					m_Rotation1_all.push_back(deltaRotation1_t);					
					}
				}														
            for(unsigned int i=1; i<m_Nf-n_ensemble; i++)                                                                                   
				{
				rotationmsd0[i]=0.0;
				rotationmsd1[i]=0.0;
				rotationmsd[i]=0.0;
				unsigned int count=0;
                for(unsigned int j=i; j<i+1; j++)                                                
					{
                    for(unsigned int k=0; k<pos_size; k++)
						{
//						if(m_type_all[i][k]==0)
							{
							
							double dxtheta0_t=m_Rotation0_all[j][k].x-m_Rotation0_all[j-i][k].x;
							double dytheta0_t=m_Rotation0_all[j][k].y-m_Rotation0_all[j-i][k].y;
							double dztheta0_t=m_Rotation0_all[j][k].z-m_Rotation0_all[j-i][k].z;
							rotationmsd0[i] += dxtheta0_t*dxtheta0_t+dytheta0_t*dytheta0_t+dztheta0_t*dztheta0_t;	
	
							double dxtheta1_t=m_Rotation1_all[j][k].x-m_Rotation1_all[j-i][k].x;
							double dytheta1_t=m_Rotation1_all[j][k].y-m_Rotation1_all[j-i][k].y;
							double dztheta1_t=m_Rotation1_all[j][k].z-m_Rotation1_all[j-i][k].z;
							rotationmsd1[i] += dxtheta1_t*dxtheta1_t+dytheta1_t*dytheta1_t+dztheta1_t*dztheta1_t;
							
							double dxtheta_t=m_Rotangle_all[j][k].x-m_Rotangle_all[j-i][k].x;
							double dytheta_t=m_Rotangle_all[j][k].y-m_Rotangle_all[j-i][k].y;
							double dztheta_t=m_Rotangle_all[j][k].z-m_Rotangle_all[j-i][k].z;
							rotationmsd[i] += dxtheta_t*dxtheta_t+dytheta_t*dytheta_t+dztheta_t*dztheta_t;
							count += 1;
							}												
						}
					}				
                rotationmsd0[i] /= double(count);
				rotationmsd1[i] /= double(count);
				rotationmsd[i] /= double(count);
				m_file << i*deltat << "  " << rotationmsd0[i] <<"   "<< rotationmsd1[i] << "  " << rotationmsd[i] << endl;
				}			
            cout << "25. Good Luck! Outputting results of the reorientational mean square angular displacement (RMSAD) to 'rmsad.log'."<<endl;			
			m_file.close();
			m_ori_all.clear();
			m_Rotation0_all.clear();
			m_Rotation1_all.clear();
			m_Rotangle_all.clear();	
            delta_t.clear();
			rotationmsd0.clear();
			rotationmsd1.clear();
			rotationmsd.clear();
			m_type_all.clear();
			};	
		void setDt(double dt)
			{
			m_dt=dt;
			}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat;
		double m_dt;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;
		std::vector<std::vector<vec> > m_ori_all;
		std::vector<std::vector<vec> > m_Rotation0_all;
		std::vector<std::vector<vec> > m_Rotation1_all;
		std::vector<std::vector<vec> > m_Rotangle_all;
		std::vector< std::vector<unsigned int> > m_type_all;
    } ;
	
//--- case 26
class ISF : public Function
    {
    public:
        ISF(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
            {
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error ISF dump");
            }
            m_Nf=0;
			m_q=6.02;
			m_dt=0.005;
        }
        virtual ~ISF()
		{						
            double c=2.0*3.1415926/Lx, cq=m_q/c;
            int Qnum=int(cq*cq+0.5);
            for(unsigned int i=0; i<=Qnum; i++)
            {
                if(qvec.size()>24) break;
                for(unsigned int j=0; j<=Qnum; j++)
                {
                    if(qvec.size()>24) break;
                    for(unsigned int k=0; k<=Qnum; k++)
                    {
                        unsigned int mm=i*i+j*j+k*k;
                        if(mm==Qnum)
                        {
                         double ic=i*c, jc=j*c, kc=k*c;
                         qvec.push_back(vec(ic,jc,kc));
						 qvec.push_back(vec(ic,-jc,kc));
						 qvec.push_back(vec(ic,jc,-kc));
						 qvec.push_back(vec(ic,-jc,-kc));							
                        }
                        Qcount=qvec.size();
                        if(qvec.size()>24) break;
                    }
                }
            }				
            std::vector<double> isf; isf.reserve(m_Nf);
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}			
			pos_size=m_pos_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;		
			for(unsigned int i=1; i<m_Nf-n_ensemble; i++) 
			{	
				unsigned int count=0;
				isf[i]=0.0;				
				if(Qcount<1) 
				{
				cout<<"***Wrong! q number is 0."<<endl; break;
				}
				for(unsigned int j=i; j<i+n_ensemble; j++)
				{
					for(unsigned int k=0; k<pos_size; k++)
					{
//						if(m_type_all[j-i][k]==0)
						{						
							double dx=m_pos_all[j][k].x-m_pos_all[j-i][k].x;
							double dy=m_pos_all[j][k].y-m_pos_all[j-i][k].y;
							double dz=m_pos_all[j][k].z-m_pos_all[j-i][k].z;
							for(unsigned int l=0; l<Qcount; l++)
							{
							 double theta=qvec[l].x*dx+qvec[l].y*dy+qvec[l].z*dz;
							 isf[i] += cosf(float(theta));
                             count += 1;							
							}						
						}						
					}
				}
				isf[i] /= double(count);
				m_file << i*deltat << " " << isf[i] << endl;
			}
			cout << "26. Good Luck! Outputting results of the self-part intermediate scattering function (ISF) to 'isf.log'." << endl;
			m_file.close();	
			m_type_all.clear();
	    	m_pos_all.clear();
			qvec.clear();
			delta_t.clear();
			isf.clear();
        }; 
		void setDt(double dt)
		{
		 m_dt=dt;
		}
		void setQ(double q)
		{
		 m_q=q;
		}	
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size, Qcount;
		double deltat, m_q, Lx, Ly, Lz;
		double m_dt;
        std::ofstream m_file;
		std::vector<vec> qvec;
		std::vector<unsigned int> delta_t;
		std::vector< std::vector<unsigned int> > m_type_all;
        std::vector<std::vector<vec> > m_pos_all;
    };	
	
//--- case 27
class OACF : public Function
    {
    public:
        OACF(std::string filename): Function()
		{
            m_file.open(filename.c_str(), ios_base::out);
            if (!m_file.good())
			{
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error OACF dump");
			}
            m_Nf=0;
			m_dt=0.005;			
		}
        virtual ~OACF()
		{            
            std::vector<double> oacf1; oacf1.reserve(m_Nf);
			std::vector<double> oacf2; oacf2.reserve(m_Nf);
			std::vector<double> oacf3; oacf3.reserve(m_Nf);
			std::vector<double> oacf4; oacf4.reserve(m_Nf);
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}			
			pos_size=m_ori_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;				
            for(unsigned int i=1; i<m_Nf-n_ensemble; i++)                                                                                
            {
				oacf1[i]=0.0;
				oacf2[i]=0.0;
				oacf3[i]=0.0;
				oacf4[i]=0.0;				
				unsigned int count=0;
				double oacfx=0.0;
                for(unsigned int j=i; j<i+n_ensemble; j++)                                                
				{
                    for(unsigned int k=0; k<pos_size; k++)
					{
//						if(m_type_all[j-i][k]==0)
						{
						 double dotprox=m_ori_all[j][k].x*m_ori_all[j-i][k].x;
						 double dotproy=m_ori_all[j][k].y*m_ori_all[j-i][k].y;
						 double dotproz=m_ori_all[j][k].z*m_ori_all[j-i][k].z;
						 oacfx=dotprox+dotproy+dotproz;
						 oacf1[i] += oacfx;
                         oacf2[i] += (3*oacfx*oacfx - 1)/2;
						 oacf3[i] += (5*oacfx*oacfx*oacfx - 3*oacfx)/2;
						 oacf4[i] += (35*oacfx*oacfx*oacfx*oacfx - 30*oacfx*oacfx + 3)/8;
						 count += 1;
						}
					}
                }
                oacf1[i] /= double(count);
				oacf2[i] /= double(count);
                oacf3[i] /= double(count);
				oacf4[i] /= double(count);				
				m_file << i*deltat << "  " << oacf1[i] << "  " << oacf2[i] << "  " << oacf3[i] << "  " << oacf4[i] << endl;
			}
            cout << "27. Good Luck! Outputting results of the orientational autocorrelation function (OACF) to 'oacf.log'." << endl;			
			m_file.close();
			m_type_all.clear();
			m_ori_all.clear();
			delta_t.clear();
			oacf1.clear();
			oacf2.clear();
			oacf3.clear();
			oacf4.clear();			
		};
		void setDt(double dt)
		{
		 m_dt=dt;
		}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat;
		double m_dt;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;
		std::vector< std::vector<unsigned int> > m_type_all;
        std::vector<std::vector<vec> > m_ori_all;
    } ;

//---case 28
class Q4Q6 : public Function
	{
    public:
        Q4Q6(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
            {
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error Q4Q6 dump");
            }			
            m_Nf=0;
			maxbin=100;
			q4_all=0.0;
			q4max=-1.0;
			q4min=1.0;
			q6_all=0.0;
			q6max=-1.0;
			q6min=1.0;
			m_rcut=1.647;
			m_Voronoi=false;			
        }
        virtual ~Q4Q6()
        {
			double dq4=(q4max-q4min)/double(maxbin);
			double dq6=(q6max-q6min)/double(maxbin);
			std::vector<double> layer4(maxbin+1, 0.0);
			std::vector<double> layer6(maxbin+1, 0.0);
			for(unsigned int i=0; i<m_Nf; i++)
			{	
				for(unsigned int k=0; k<pos_size; k++)
				{
				 unsigned int mth;
				 mth=int((q4_local_all[i][k]-q4min)/dq4); 
				 layer4[mth] += 1.0;
				 unsigned int lth;
				 lth=int((q6_local_all[i][k]-q6min)/dq6); 
				 layer6[lth] += 1.0;
				}
			}
			for(unsigned int bin=0; bin<maxbin+1; bin++)
			{
				if(bin==0) 
				 m_file << q4min+bin*dq4 << "  " << layer4[bin]/double(pos_size*m_Nf) << "  " << q6min+bin*dq6 << "  " << layer6[bin]/double(pos_size*m_Nf) << "  " << q4_all/m_Nf << "  " << q6_all/m_Nf << endl;
				else 
				 m_file << q4min+bin*dq4 << "  " << layer4[bin]/double(pos_size*m_Nf) << "  " << q6min+bin*dq6 << "  " << layer6[bin]/double(pos_size*m_Nf) << endl;
			}
           cout << "28. Good Luck! Outputting results of bond order parameters (Q4Q6) and outputting results to 'q4q6.log'." << endl;			
			m_file.close();	
			q4_local_all.clear();
			q6_local_all.clear();
     	}
		void setVoronoi(bool Voronoi)
		{
		 m_Voronoi = Voronoi;
		}
		void setRcut(double rcut)
		{
		 m_rcut=rcut;
		}		
        virtual void compute();
    private:
        std::ofstream m_file;
		bool m_Voronoi;
        unsigned int m_Nf, pos_size, maxbin;
		double m_rcut;
		double q4_all, q4max, q4min, q6_all, q6max, q6min;
		std::vector<std::vector<double> > q4_local_all;
		std::vector<std::vector<double> > q6_local_all;
	};

//---case 29
class VORONOI : public Function
	{
    public:
        VORONOI(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
            {
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error VORONOI dump");
            }			
            m_Nf=0;
			m_rcut=1.647;
			maxbin=100;
			voronoimax=0.0;
			voronoimin=9.0;
			voronoi_all=0.0;
        }
        virtual ~VORONOI()
        {			
			double dv=(voronoimax-voronoimin)/double(maxbin);
			std::vector<double> layer(maxbin+1, 0.0);
			for(unsigned int i=0; i<m_Nf; i++)
			{	
				for(unsigned int k=0; k<pos_size; k++)
				{
				 unsigned int lth;
				 lth=int((voronoi_local_all[i][k]-voronoimin)/dv); 
				 layer[lth] += 1.0;
				}
			}
			for(unsigned int bin=0; bin<maxbin+1; bin++)
			{
				if(bin==0)
				 m_file << "  " << voronoimin+bin*dv << "  " << layer[bin]/double(pos_size*m_Nf) << "  " << voronoi_all/m_Nf << "  " << Lx*Ly*Lz<< endl;
				else
				 m_file << "  " << voronoimin+bin*dv << "  " << layer[bin]/double(pos_size*m_Nf) << endl;
			}
           cout << "29. Good Luck! Outputting results of the volume of the Voronoi cells (VORONOI) and outputting results to 'voronoi.log'." << endl;				
			m_file.close();	
			voronoi_local_all.clear();		
     	}
		void setRcut(double rcut)
		{
		 m_rcut=rcut;
		}		
        virtual void compute();
    private:
        std::ofstream m_file;
		bool Voronoi;
        unsigned int m_Nf, pos_size, maxbin;
		double m_rcut;
		double Lx, Ly, Lz, LxINV, LyINV, LzINV;
		double voronoi_all, voronoimax, voronoimin;
		std::vector<std::vector<double> > voronoi_local_all;
	};

//--- case 30
class nonGauPar: public Function
    {
    public:
        nonGauPar(std::string filename): Function()
		{
            m_file.open(filename.c_str(), ios_base::out);
            if (!m_file.good())
			{
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error nonGauPar dump");
			}
            m_Nf=0;
			m_dt=0.005;
		}
        virtual ~nonGauPar()
		{
			std::ofstream msd_file;
			msd_file.open("msd.log", ios_base::app);
			if (!msd_file.good())
			{
			 cerr << endl << "***Error! Error opening dump file." << endl << endl;
			 throw runtime_error("Error rave dump");
			}
			std::vector<double> msd, msd2; msd.resize(m_Nf); msd2.resize(m_Nf);
			std::vector<double> alpha2t; alpha2t.resize(m_Nf);
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}			
			pos_size=m_pos_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;									
            for(unsigned int i=1; i<m_Nf-n_ensemble; i++)
			{
				unsigned int count=0;
				double r2modul=0.0;
                for(unsigned int j=i; j<i+n_ensemble; j++)
				{
                    for(unsigned int k=0; k<pos_size; k++)
					{
//						if(m_type_all[j-i][k]==0)
						{
						 double dx=m_pos_all[j][k].x-m_pos_all[j-i][k].x;
						 double dy=m_pos_all[j][k].y-m_pos_all[j-i][k].y;
						 double dz=m_pos_all[j][k].z-m_pos_all[j-i][k].z;
						 r2modul=dx*dx+dy*dy+dz*dz;
						 msd[i] += r2modul;
						 msd2[i]+= r2modul*r2modul;	
						 count += 1;
						}
					}
				}
                msd[i] /= double(count); 
				msd2[i] /= double(count);
				alpha2t[i]=0.6*msd2[i]/(msd[i]*msd[i])-1;
                m_file << i*deltat << "  " << alpha2t[i] << endl;
			}
			cout << "30. Good Luck! Outputting results of the non-Gaussian parameter (NGP) to 'nongaupar.log'." << endl;
            m_file.close();
			msd_file.close();
			m_type_all.clear();
			m_pos_all.clear();
			delta_t.clear();
			msd.clear();
			msd2.clear();
			alpha2t.clear();
		}
		void setDt(double dt)
		{
		 m_dt=dt;
		}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat;
		double m_dt;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;
		std::vector< std::vector<unsigned int> > m_type_all;
        std::vector<std::vector<vec> > m_pos_all;
    };

//--- case 31
class RnonGauPar: public Function
    {
    public:
        RnonGauPar(std::string filename): Function()
		{
            m_file.open(filename.c_str(), ios_base::out);
            if (!m_file.good())
			{
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error RnonGauPar dump");
			}
            m_Nf=0;
			m_dt=0.005;
		}
        virtual ~RnonGauPar()
		{
			std:: ofstream msad_file;
			msad_file.open("msad.log", ios_base::app);
			if (!msad_file.good())
			{
			 cerr << endl << "***Error! Error opening dump file." << endl << endl;
			 throw runtime_error("Error rave dump");
			}
			std::vector<double> msad, msad2; msad.resize(m_Nf); msad2.resize(m_Nf);
			std::vector<double> Ralpha2t; Ralpha2t.resize(m_Nf);				
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}			
            pos_size=m_Rotangle_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;								
            for(unsigned int i=1; i<m_Nf-n_ensemble; i++)
			{
				unsigned int count=0;
				double theta2modul=0.0;
				msad2[i]=0.0;
                for(unsigned int j=i; j<i+n_ensemble; j++)
				{
                    for(unsigned int k=0; k<pos_size; k++)
					{
//						if(m_type_all[j-i][k]==0)
						{									
						 double dxtheta_t=m_Rotangle_all[j][k].x-m_Rotangle_all[j-i][k].x;
						 double dytheta_t=m_Rotangle_all[j][k].y-m_Rotangle_all[j-i][k].y;
						 double dztheta_t=m_Rotangle_all[j][k].z-m_Rotangle_all[j-i][k].z;
						 theta2modul=dxtheta_t*dxtheta_t + dytheta_t*dytheta_t + dztheta_t*dztheta_t;
						 msad[i] += theta2modul;
						 msad2[i] += theta2modul*theta2modul;
						 count += 1;
						}
					}
				}
                msad[i] /= double(count);
				msad2[i] /= double(count);
				Ralpha2t[i]=0.60*msad2[i]/(msad[i]*msad[i])-1;
                m_file << i*deltat << "  "<< Ralpha2t[i] << endl;
			}
			cout << "31. Good Luck! Outputting results of the rotational non-Gaussian parameter RNGP to 'rnongaupar.log'." << endl;
            m_file.close();
			m_type_all.clear();
			msad_file.close();
			m_Rotangle_all.clear();
			delta_t.clear();
			msad.clear();
			msad2.clear();
			Ralpha2t.clear();
		}
		void setDt(double dt)
		{
		 m_dt=dt;
		}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat;
		double m_dt;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;
		std::vector< std::vector<unsigned int> > m_type_all;
        std::vector<std::vector<vec> > m_Rotangle_all;
    };	
	
//--- case 32
class SVH : public Function
    {
    public:	
       	SVH(std::string filename): Function()
		{
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
			{
             cerr << endl << "***Error! Error opening dump file" << filename << endl << endl;
             throw runtime_error("Error SVH dump");
			}			
            m_Nf=0;
			m_dt=0.005;
		}
        virtual ~SVH()
		{
			std::vector< std::vector<double> > rt;
			pos_size=m_pos_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}	
			unsigned int num_t=0;			
			std::vector<unsigned int> delta_time;
			for(unsigned int i=0; i!=10; i++)
			{
				for(unsigned int j=1; j!=10; j++)
				{
					unsigned int deltatime=j*pow(10,i);
					if((n_ensemble+deltatime)<m_Nf)
					{
					 delta_time.push_back(deltatime);
					 num_t += 1;
					}					
				}			
			}
			double maxdelr=0.0, mindelr=100.0;
			std::vector<unsigned int> count; count.resize(num_t);
			std::vector<double> m; m.resize(num_t);	
            for(unsigned int i=0; i!=num_t; i++)
			{
                unsigned int delt=delta_time[i];
				if((n_ensemble+delt)<m_Nf)
				{
					double r2temp=0.0;
					count[i]=0, m[i]=0.0;
					std::vector<double> r; r.resize(n_ensemble*pos_size);
					for(unsigned int j=0; j!=n_ensemble; j++)
					{
						for(unsigned int k=0; k!=pos_size; k++)
						{
//							if(m_type_all[i][k]==0)
							{
								double dx=m_pos_all[j+delt][k].x-m_pos_all[j][k].x; 
								double dy=m_pos_all[j+delt][k].y-m_pos_all[j][k].y;
								double dz=m_pos_all[j+delt][k].z-m_pos_all[j][k].z;
								double delr=sqrt(dx*dx + dy*dy + dz*dz);
								r[count[i]]=delr;
								r2temp += dx*dx + dy*dy + dz*dz;
								if(delr>=maxdelr)
								{
								 maxdelr=delr;	
								}					
								if(delr<=mindelr)
								{
								 mindelr=delr;	
								}						
								count[i] += 1;																
							}					
						}
					}
					rt.push_back(r);
					r.clear();
					m[i]=1.5*count[i]/r2temp;										
				} 
				else
				{	
				 cerr << endl << "***Error! Error delta_time" << endl;
				 throw runtime_error("Error delta_time dump");				 
				}				
			}
			double delr=0.01;	
			cout << delr << " " << (maxdelr-mindelr)/1000.0 << maxdelr << " " << mindelr << endl;	
			
			std::vector< std::vector<double> > nlayer;
            for(unsigned int i=0; i!=num_t; i++)
			{	
				std::vector<double> layer(n_ensemble*pos_size, 0.0);
                for(unsigned int j=0; j!=count[i]; j++)
				{   
				 unsigned int lth=int(rt[i][j]/delr); 
				 layer[lth] += 1.0;
				}				
				nlayer.push_back(layer);
				layer.clear();	
                for(unsigned int j=0; j!=count[i]; j++)
				{
				 nlayer[i][j] /= double(count[i]);					
                }				
			}
			
			m_file << "r";
			for(unsigned int i=0; i!=num_t; i++)
			{	
			 m_file << "  t=" << deltat*delta_time[i];
			}
			m_file << endl;	
            for(unsigned int j=1; j!=n_ensemble*pos_size; j++)
			{				
				double r1=delr*j;			
				if(r1<=100.0)
				{						
					m_file << r1;
					for(unsigned int i=0; i!=num_t; i++)
					{	
					 m_file << "  " << nlayer[i][j]/delr;
					}
					m_file << endl;				
				}
			}										
			cout << "32. Good Luck! Outputting results of the self van Hove fucntion (VHF) to 'selfvhf.log'." << endl;			
			m_file.close();	
			m_pos_all.clear();
		}
		void setDt(double dt)
		{
		 m_dt=dt;
		}			
		virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat;
		double m_dt;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;				
        std::vector<std::vector<vec> > m_pos_all;
		std::vector<std::vector<unsigned int> > m_type_all;
    };
	
//--- case 33
class RSVH : public Function
    {
    public:	
       	RSVH(std::string filename): Function()
		{
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
			{
             cerr << endl << "***Error! Error opening dump file" << filename << endl << endl;
             throw runtime_error("Error RSVH dump");
			}			
            m_Nf=0;
			m_dt=0.005;
		}
        virtual ~RSVH()
		{
			std::vector< std::vector<double> > thetat;
			pos_size=m_ori_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}			
			unsigned int num_t=0;			
			std::vector<unsigned int> delta_time;
			for(unsigned int i=0; i!=10; i++)
			{
				for(unsigned int j=1; j!=10; j++)
				{
					unsigned int deltatime=j*pow(10,i);
					if((n_ensemble+deltatime)<m_Nf)
					{
					 delta_time.push_back(deltatime);
					 num_t += 1;
					}					
				}			
			}
			double maxdeltheta=0.0, mindeltheta=100.0;
			std::vector<unsigned int> count; count.resize(num_t);
			std::vector<double> m; m.resize(num_t);									
            for(unsigned int i=0; i!=num_t; i++)
			{
                unsigned int delt=delta_time[i];
				if((n_ensemble+delt)<m_Nf)
				{
					double r2temp=0.0;
					count[i]=0, m[i]=0.0;
					std::vector<double> theta; theta.resize(n_ensemble*pos_size);
					for(unsigned int j=0; j!=n_ensemble; j++)
					{
						for(unsigned int k=0; k!=pos_size; k++)
						{
//							if(m_type_all[i][k]==0)
							{
								double dotprox=m_ori_all[j+delt][k].x*m_ori_all[j][k].x; 
								double dotproy=m_ori_all[j+delt][k].y*m_ori_all[j][k].y;
								double dotproz=m_ori_all[j+delt][k].z*m_ori_all[j][k].z;
								double dotpro=dotprox+dotproy+dotproz;
								double deltheta;
								if(fabs(dotpro)<1.0)
								{
								 deltheta=acos(dotpro);	
								}							
								else if(dotpro>=1.0)
								{
								 deltheta=0.0;	
								}						
								else
								{
								 deltheta=3.1415926;	
								}							 						 
								theta[count[i]]=deltheta;
								double dxtheta_t=m_Rotangle_all[j+delt][k].x-m_Rotangle_all[j][k].x;
								double dytheta_t=m_Rotangle_all[j+delt][k].y-m_Rotangle_all[j][k].y;
								double dztheta_t=m_Rotangle_all[j+delt][k].z-m_Rotangle_all[j][k].z;
								r2temp += dxtheta_t*dxtheta_t + dytheta_t*dytheta_t + dztheta_t*dztheta_t;					
								if(deltheta>=maxdeltheta)
								{
								 maxdeltheta=deltheta;	
								}					
								if(deltheta<=mindeltheta)
								{
								 mindeltheta=deltheta;	
								}						
								count[i] += 1;																
							}																				
						}
					}
					thetat.push_back(theta);
					theta.clear();
					m[i]=1.5*count[i]/r2temp;										
				} 
				else
				{	
				 cerr << endl << "***Error! Error delta_time" << endl;
				 throw runtime_error("Error delta_time dump");				 
				}				
			}
			double dtheta=0.005;
			cout << dtheta << " " << (maxdeltheta-mindeltheta)/1000.0 << "  " << maxdeltheta << " " << mindeltheta << endl;	

			std::vector< std::vector<double> > nlayer;
            for(unsigned int i=0; i!=num_t; i++)
			{	
				std::vector<double> layer(n_ensemble*pos_size, 0.0);
                for(unsigned int j=0; j!=count[i]; j++)
				{   
				 unsigned int lth=int(thetat[i][j]/dtheta); 
				 layer[lth] += 1.0;
				}				
				nlayer.push_back(layer);
				layer.clear();			
                for(unsigned int j=0; j!=count[i]; j++)
				{
				 nlayer[i][j] /= double(count[i]);					
                }				
			}
			m_file << "theta";
			for(unsigned int i=0; i!=num_t; i++)
			{	
			 m_file << "  t=" << m_dt*delta_time[i];
			}
			m_file << endl;				
            for(unsigned int j=1; j!=n_ensemble*pos_size; j++)
			{
				double theta1=dtheta*j;				
				if(theta1<=3.1415926)
				{						
					m_file << theta1;
					for(unsigned int i=0; i!=num_t; i++)
					{	
					 m_file << "  " << nlayer[i][j]/dtheta;
					}
					m_file << endl;				
				}
			}										
			cout << "33. Good Luck! Outputting results of the rotational self van Hove fucntion to 'rselfvhf.log'." << endl;			
			m_file.close();	
			m_ori_all.clear();
			m_Rotangle_all.clear();
			nlayer.clear();
		}
		void setDt(double dt)
		{
		 m_dt=dt;
		}			
		virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat;
		double m_dt;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;
        std::vector<std::vector<vec> > m_ori_all;
		std::vector<std::vector<vec> > m_Rotangle_all;
		std::vector<std::vector<unsigned int> > m_type_all;
    };	
	
//-- case 34
class fpSus : public Function
    {
    public:
       	fpSus(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if (!m_file.good())
			{
            cerr<<endl<<"***Error! Error opening dump file."<<filename<<endl<<endl;
            throw runtime_error("Error fpSus dump");
			}
            m_Nf=0;
			m_q=6.02;
			m_dt=0.005;
		}
        virtual ~fpSus()
		{
		    std::vector<double> isf; isf.resize(m_Nf);
			std::vector<double> SQISF; SQISF.resize(m_Nf);
			std::vector<double> chi4t; chi4t.resize(m_Nf);
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}
			pos_size=m_pos_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;
            double c=2.0*3.1415926/Lx, cq=m_q/c;
            int Qnum=int(cq*cq+0.5);
            for(unsigned int i=0; i<=Qnum; i++)
            {
                if(qvec.size()>24) break;
                for(unsigned int j=0; j<=Qnum; j++)
                {
                    if(qvec.size()>24) break;
                    for(unsigned int k=0; k<=Qnum; k++)
                    {
                        unsigned int mm=i*i+j*j+k*k;
                        if(mm==Qnum)
                        {
                         double ic=i*c, jc=j*c, kc=k*c;
                         qvec.push_back(vec(ic,jc,kc));
						 qvec.push_back(vec(ic,-jc,kc));
						 qvec.push_back(vec(ic,jc,-kc));
						 qvec.push_back(vec(ic,-jc,-kc));							
                        }
                        Qcount=qvec.size();
                        if(qvec.size()>24) break;
                    }
                }
            }
			
			for(unsigned int i=1; i<m_Nf-n_ensemble; i++) 
			{	
				unsigned int count=0;
				isf[i]=0.0;
				SQISF[i]=0.0;
				if(Qcount<1) 
				{
				 cout<<"***Wrong! q number is 0."<<endl; break;
				}
				for(unsigned int j=i; j<i+n_ensemble; j++)
				{    	 
					
					for(unsigned int l=0; l<Qcount; l++)
					{
						double cos_qdr_all=0.0;
						for(unsigned int k=0; k<pos_size; k++)
						{
						 double dx=m_pos_all[j][k].x-m_pos_all[j-i][k].x;
						 double dy=m_pos_all[j][k].y-m_pos_all[j-i][k].y;
						 double dz=m_pos_all[j][k].z-m_pos_all[j-i][k].z;
						 double theta=qvec[l].x*dx+qvec[l].y*dy+qvec[l].z*dz;	
						 cos_qdr_all += cosf(float(theta));     
						}
						cos_qdr_all /= double(pos_size);
						isf[i] += cos_qdr_all;
						SQISF[i] += cos_qdr_all*cos_qdr_all;	
						count += 1;					
					}
				}
				isf[i] /= double(count);			
				SQISF[i] /= double(count);
				chi4t[i]=pos_size*(SQISF[i]-isf[i]*isf[i]);
				m_file << deltat*i << "  " << chi4t[i] << "  "<< isf[i] <<endl;				
			}								
			cout<<"34. Good Luck! Outputting results of the four-point susceptibility to 'fpsus.log'."<<endl;
			m_file.close();
			m_pos_all.clear();
			isf.clear();
			SQISF.clear();
			chi4t.clear();			
		};
		void setDt(double dt)
		{
		 m_dt=dt;
		}
		void setQ(double q)
		{
		 m_q=q;
		}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size, Qcount;
		double deltat, Lx, Ly, Lz;
		double m_q;
		double m_dt;
		std::ofstream m_file;
		std::vector<unsigned int> delta_t;
		std::vector<vec> qvec;
        std::vector<std::vector<vec> > m_pos_all;
		std::vector<std::vector<unsigned int> > m_type_all;
    };	
	
//-- case 35
class RfpSus : public Function
    {
    public:
       	RfpSus(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if (!m_file.good())
			{
            cerr<<endl<<"***Error! Error opening dump file."<<filename<<endl<<endl;
            throw runtime_error("Error RfpSus dump");
			}
            m_Nf=0;
			m_dt=0.005;
		}
        virtual ~RfpSus()
		{
            std::vector<double> oacf1; oacf1.reserve(m_Nf);
			std::vector<double> oacf2; oacf2.reserve(m_Nf);
            std::vector<double> sqoacf1; sqoacf1.reserve(m_Nf);
			std::vector<double> sqoacf2; sqoacf2.reserve(m_Nf);
			std::vector<double> chi4t_1; chi4t_1.resize(m_Nf);
			std::vector<double> chi4t_2; chi4t_2.resize(m_Nf);
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}
			pos_size=m_ori_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;
			
            for(unsigned int i=1; i<m_Nf-n_ensemble; i++)                                                                                
            {
				oacf1[i]=0.0;
				oacf2[i]=0.0;
				sqoacf1[i]=0.0;
				sqoacf2[i]=0.0;				
				unsigned int count=0;
                for(unsigned int j=i; j<i+n_ensemble; j++)                                                
				{
					double ori_ijt1_all=0.0;
					double ori_ijt2_all=0.0;
                    for(unsigned int k=0; k<pos_size; k++)
					{
					 double dotprox=m_ori_all[j][k].x*m_ori_all[j-i][k].x;
					 double dotproy=m_ori_all[j][k].y*m_ori_all[j-i][k].y;
					 double dotproz=m_ori_all[j][k].z*m_ori_all[j-i][k].z;
					 double oacfx=dotprox+dotproy+dotproz;
					 ori_ijt1_all +=  oacfx;
					 ori_ijt2_all +=  (3*oacfx*oacfx - 1)/2;
					}
					ori_ijt1_all /= double(pos_size);
					ori_ijt2_all /= double(pos_size);
					oacf1[i] += ori_ijt1_all;
					sqoacf1[i] += ori_ijt1_all*ori_ijt1_all;
					oacf2[i] += ori_ijt2_all;
					sqoacf2[i] += ori_ijt2_all*ori_ijt2_all;					
					count += 1;					
                }
                oacf1[i] /= double(count);
				oacf2[i] /= double(count);				
                sqoacf1[i] /= double(count);
				sqoacf2[i] /= double(count);
				chi4t_1[i]=pos_size*(sqoacf1[i]-oacf1[i]*oacf1[i]);
				chi4t_2[i]=pos_size*(sqoacf2[i]-oacf2[i]*oacf2[i]);
				m_file << deltat*i << "  " << chi4t_1[i] << "  " << chi4t_2[i] << "  "<< oacf1[i] << "  "<< oacf2[i] <<endl;						
			}
			cout<<"35. Good Luck! Outputting results of the rotaional four-point susceptibility to 'rfpsus.log'."<<endl;
			m_file.close();
			m_ori_all.clear();
			oacf1.clear();
			oacf2.clear();
			sqoacf1.clear();
			sqoacf1.clear();
			chi4t_1.clear();
			chi4t_2.clear();
		};
		void setDt(double dt)
		{
		 m_dt=dt;
		}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat;
		double m_dt;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;
		std::vector< std::vector<unsigned int> > m_type_all;
        std::vector<std::vector<vec> > m_ori_all;
    };
	
//--- case 36
class OVLAF : public Function
    {
    public:
        OVLAF(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
            {
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error OVLAF dump");
            }
            m_Nf=0;
			m_dt=0.005;
			m_a=0.30;
        }
        virtual ~OVLAF()
		{										
            std::vector<double> ovlaf; ovlaf.reserve(m_Nf);
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}			
			pos_size=m_pos_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;			
			for(unsigned int i=1; i<m_Nf-n_ensemble; i++) 
			{	
				unsigned int count=0;
				ovlaf[i]=0.0;				
				for(unsigned int j=i; j<i+n_ensemble; j++)
				{
					for(unsigned int k=0; k<pos_size; k++)
					{
//						if(m_type_all[i][k]==0)
						{				
							double dx=m_pos_all[j][k].x-m_pos_all[j-i][k].x;
							double dy=m_pos_all[j][k].y-m_pos_all[j-i][k].y;
							double dz=m_pos_all[j][k].z-m_pos_all[j-i][k].z;
							double r2modul=dx*dx+dy*dy+dz*dz;
							if(r2modul<=m_a*m_a)
							{
							 ovlaf[i] += 1.0;	 
							}
							else
							{
							 ovlaf[i] += 0.0;
							}
                            count += 1;												
						}						
					}
				}
				ovlaf[i] /= double(count);
				m_file << i*deltat << " " << ovlaf[i] << endl;
			}
			cout << "36. Good Luck! Outputting results of the self-part overlap function (OVLAF) to 'ovlaf.log'." << endl;
			m_file.close();	
			m_type_all.clear();
	    	m_pos_all.clear();
			delta_t.clear();
			ovlaf.clear();
        }; 
		void setDt(double dt)
		{
		 m_dt=dt;
		}
		void setA(double a)
		{
		 m_a=a;
		}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat, Lx, Ly, Lz;
		double m_dt;
		double m_a;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;
		std::vector< std::vector<unsigned int> > m_type_all;
        std::vector<std::vector<vec> > m_pos_all;
    };

//--- case 37
class CISF : public Function
    {
    public:
        CISF(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
            {
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error CISF dump");
            }
            m_Nf=0;
			m_dt=0.005;
			m_q=6.02;
        }
        virtual ~CISF()
		{						
            double c=2.0*3.1415926/Lx, cq=m_q/c;
            int Qnum=int(cq*cq+0.5);
            for(unsigned int i=0; i<=Qnum; i++)
            {
                if(qvec.size()>48) break;
                for(unsigned int j=0; j<=Qnum; j++)
                {
                    if(qvec.size()>48) break;
                    for(unsigned int k=0; k<=Qnum; k++)
                    {
                        unsigned int mm=i*i+j*j+k*k;
                        if(mm==Qnum)
                        {
                         double ic=i*c, jc=j*c, kc=k*c;
                         qvec.push_back(vec(ic,jc,kc));
						 qvec.push_back(vec(ic,-jc,kc));
						 qvec.push_back(vec(ic,jc,-kc));
						 qvec.push_back(vec(ic,-jc,-kc));							
                        }
                        Qcount=qvec.size();
                        if(qvec.size()>48) break;
                    }
                }
            }				
            std::vector<double> cisf; cisf.reserve(m_Nf);
			std::vector<double> cssf; cssf.reserve(m_Nf);
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;			 
			}				
            else
			{
			 n_ensemble=1000;	
			}			
			pos_size=m_pos_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;
			std::vector<unsigned int> pow_time; pow_time.push_back(1);
			for(unsigned int i=0; i<1000; i++)
			{
				unsigned int pow_num = pow(10, i*0.01)+0.50;
				unsigned int num_t=pow_time.size();
				if(pow_num<m_Nf-n_ensemble && pow_time[num_t-1]!=pow_num)
				{
				 pow_time.push_back(pow_num);					
				}				
			}													
			for(unsigned int ii=0; ii<pow_time.size(); ii++) 
			{
				unsigned int i=pow_time[ii];
				unsigned int count=0;
				cisf[i]=0.0;
				cssf[i]=0.0;
				if(Qcount<1) 
				{
				 cout << "***Wrong! q number is 0." << endl; break;
				}
				for(unsigned int j=i; j<i+n_ensemble; j++)
				{
					for(unsigned int l=0; l<Qcount; l++)
					{						
						double cos_qrt_all=0.0;
						double cos_qr0_all=0.0;
						double sin_qrt_all=0.0;
						double sin_qr0_all=0.0;						
						for(unsigned int k=0; k<pos_size; k++)
						{												
//							if(m_type_all[i][k]==0)
							{												
							 cos_qrt_all += cosf(qvec[l].x*(m_pos_all[j][k].x)+qvec[l].y*(m_pos_all[j][k].y)+qvec[l].z*(m_pos_all[j][k].z));
							 cos_qr0_all += cosf(qvec[l].x*(m_pos_all[j-i][k].x)+qvec[l].y*(m_pos_all[j-i][k].y)+qvec[l].z*(m_pos_all[j-i][k].z));
							 sin_qrt_all += sinf(qvec[l].x*(m_pos_all[j][k].x)+qvec[l].y*(m_pos_all[j][k].y)+qvec[l].z*(m_pos_all[j][k].z));
							 sin_qr0_all += sinf(qvec[l].x*(m_pos_all[j-i][k].x)+qvec[l].y*(m_pos_all[j-i][k].y)+qvec[l].z*(m_pos_all[j-i][k].z));							 
							}						
						}
						cisf[i] += cos_qrt_all*cos_qr0_all + sin_qrt_all*sin_qr0_all;						
						cssf[i] += cos_qr0_all*cos_qr0_all + sin_qr0_all*sin_qr0_all;
						count += 1;
					}
				}
				cisf[i] /= double(count);
				cssf[i] /= double(count);
				m_file << i*deltat << " " << cisf[i]/cssf[i] << endl;
			}
					
			cout << "37. Good Luck! Outputting results of the coherent intermediate scattering function (CISF) to 'cisf.log'." << endl;
			m_file.close();	
			m_type_all.clear();
	    	m_pos_all.clear();
			qvec.clear();
			delta_t.clear();
			cisf.clear();			
        }; 
		void setDt(double dt)
		{
		 m_dt=dt;
		}
		void setQ(double q)
		{
		 m_q=q;
		}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size, Qcount;
		double deltat, m_q, Lx, Ly, Lz;
		double m_dt;		
        std::ofstream m_file;
		std::vector<vec> qvec;
		std::vector<unsigned int> delta_t;
		std::vector<std::vector<unsigned int> > m_type_all;
        std::vector<std::vector<vec> > m_pos_all;
    };

//--- case 38
class CAGEISF : public Function
    {
    public:
        CAGEISF(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
            {
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error CAGEIS dump");
            }
            m_Nf=0;
			m_dt=0.005;			
			m_q=6.02;
			m_rcut=1.65;
			m_Voronoi=false;
			
        }
        virtual ~CAGEISF()
		{						
            double c=2.0*3.1415926/Lx, cq=m_q/c;
            int Qnum=int(cq*cq+0.5);
            for(unsigned int i=0; i<=Qnum; i++)
            {
                if(qvec.size()>24) break;
                for(unsigned int j=0; j<=Qnum; j++)
                {
                    if(qvec.size()>24) break;
                    for(unsigned int k=0; k<=Qnum; k++)
                    {
                        unsigned int mm=i*i+j*j+k*k;
                        if(mm==Qnum)
                        {
                         double ic=i*c, jc=j*c, kc=k*c;
                         qvec.push_back(vec(ic,jc,kc));
						 qvec.push_back(vec(ic,-jc,kc));
						 qvec.push_back(vec(ic,jc,-kc));
						 qvec.push_back(vec(ic,-jc,-kc));							
                        }
                        Qcount=qvec.size();
                        if(qvec.size()>24) break;
                    }
                }
            }				
            std::vector<double> cageisf; cageisf.reserve(m_Nf);
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}			
			pos_size=m_pos_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;	
			for(unsigned int i=1; i<m_Nf-n_ensemble; i++) 
			{	
				unsigned int count=0;
				cageisf[i]=0.0;				
				if(Qcount<1) 
				{
				cout<<"***Wrong! q number is 0."<<endl; break;
				}
				for(unsigned int j=i; j<i+n_ensemble; j++)
				{
					for(unsigned int k=0; k<pos_size; k++)
					{
//						if(m_type_all[i][k]==0)
						{
							double nbdx=0.0, nbdy=0.0, nbdz=0.0;
							for(unsigned int n=0; n<m_num_all[j-i][k]; n++)
							{
							 unsigned int nb=m_nb_all[j-i][k][n];
							 nbdx += m_pos_all[j][nb].x - m_pos_all[j-i][nb].x;
							 nbdy += m_pos_all[j][nb].y - m_pos_all[j-i][nb].y;
							 nbdz += m_pos_all[j][nb].z - m_pos_all[j-i][nb].z;
							}			
							double dx=m_pos_all[j][k].x-m_pos_all[j-i][k].x-nbdx/m_num_all[j-i][k];
							double dy=m_pos_all[j][k].y-m_pos_all[j-i][k].y-nbdy/m_num_all[j-i][k];
							double dz=m_pos_all[j][k].z-m_pos_all[j-i][k].z-nbdz/m_num_all[j-i][k];								
							for(unsigned int l=0; l<Qcount; l++)
							{
							 double theta=qvec[l].x*dx+qvec[l].y*dy+qvec[l].z*dz;
							 cageisf[i] += cosf(float(theta));
                             count += 1;							
							}						
						}						
					}
				}
				cageisf[i] /= double(count);
				m_file << i*deltat << " " << cageisf[i] << endl;
			}
			cout << "38. Good Luck! Outputting results of the cage-relative self-part intermediate scattering function (CAGEISF) to 'cageisf.log'." << endl;
			m_file.close();	
			m_type_all.clear();
	    	m_pos_all.clear();
			qvec.clear();
			delta_t.clear();
			cageisf.clear();
        };
		void setDt(double dt)
		{
		 m_dt=dt;
		}
		void setQ(double q)
		{
		 m_q=q;
		}
		void setVoronoi(bool Voronoi)
		{
		 m_Voronoi = Voronoi;
		}
		void setRcut(double rcut)
		{
		 m_rcut=rcut;
		}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size, Qcount;
		double deltat, m_rcut, m_q, Lx, Ly, Lz, LxINV, LyINV, LzINV;
		double m_dt;
		bool m_Voronoi;
        std::ofstream m_file;
		std::vector<vec> qvec;
		std::vector<unsigned int> delta_t;
		std::vector< std::vector<unsigned int> > m_type_all;
		std::vector< std::vector<unsigned int> > m_num_all;
		std::vector< std::vector<std::vector<unsigned int> > > m_nb_all;	
        std::vector<std::vector<vec> > m_pos_all;
    };

//--- case 39
class CAGEMSD : public Function
    {
    public:
        CAGEMSD(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
            {
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error CAGEMSD dump");
            }
            m_Nf=0;
			m_dt=0.005;
			m_rcut=1.65;
			m_Voronoi=false;
			
        }
        virtual ~CAGEMSD()
		{						
            std::vector<double> cagemsd; cagemsd.reserve(m_Nf);
			if(m_Nf<=1000) 				
			{
			 n_ensemble=0.1*m_Nf;
            }					
            else 
			{
			 n_ensemble=1000;		
			}					
            pos_size=m_pos_all[0].size();
			deltat=(delta_t[1]-delta_t[0])*m_dt;			
            for(unsigned int i=1; i<m_Nf-n_ensemble; i++)                                                                                   
            {
				cagemsd[i]=0.0;
				unsigned int count=0;
                for(unsigned int j=i; j<i+n_ensemble; j++)                                                
				{
                    for(unsigned int k=0; k<pos_size; k++)
					{
//						if(m_type_all[i][k]==0)
						{
							double nbdx=0.0, nbdy=0.0, nbdz=0.0;
							for(unsigned int n=0; n<m_num_all[j-i][k]; n++)
							{
							 unsigned int nb=m_nb_all[j-i][k][n];
//							 if(i==1 && j==1) cout <<"  i="<< i <<"  j="<< j <<"  k="<< k <<"  m_num_all[j-i][k]="<< m_num_all[j-i][k] <<"  nb="<< nb << endl;
							 nbdx += m_pos_all[j][nb].x - m_pos_all[j-i][nb].x;
							 nbdy += m_pos_all[j][nb].y - m_pos_all[j-i][nb].y;
							 nbdz += m_pos_all[j][nb].z - m_pos_all[j-i][nb].z;
							}										
							double dx=m_pos_all[j][k].x-m_pos_all[j-i][k].x-nbdx/m_num_all[j-i][k];
							double dy=m_pos_all[j][k].y-m_pos_all[j-i][k].y-nbdy/m_num_all[j-i][k];
							double dz=m_pos_all[j][k].z-m_pos_all[j-i][k].z-nbdz/m_num_all[j-i][k];							
							cagemsd[i] += dx*dx + dy*dy + dz*dz;	
							count += 1;						
						}						
					}	
                }
				cagemsd[i] /= double(count);
				m_file << i*deltat << "  " << cagemsd[i] << endl;
			}

            cout << "39. Good Luck! Outputting results of the cage-relative mean square displacement CAGEMSD to 'cagemsd.log'." << endl;			
			m_file.close();
			delta_t.clear();
			m_type_all.clear();
			m_pos_all.clear();
			cagemsd.clear();
        };
		void setDt(double dt)
		{
		 m_dt=dt;
		}
		void setVoronoi(bool Voronoi)
		{
		 m_Voronoi = Voronoi;
		}
		void setRcut(double rcut)
		{
		 m_rcut=rcut;
		}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat, m_rcut, Lx, Ly, Lz, LxINV, LyINV, LzINV;
		double m_dt;
		bool m_Voronoi;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;
		std::vector< std::vector<unsigned int> > m_type_all;
		std::vector< std::vector<unsigned int> > m_num_all;
		std::vector< std::vector<std::vector<unsigned int> > > m_nb_all;	
        std::vector<std::vector<vec> > m_pos_all;
    };	

//--- case 40	
class RMSD : public Function
    {
    public:
        RMSD(std::string filename): Function()
		{
            m_file.open(filename.c_str(), ios_base::out);
            if (!m_file.good())
			{
             cerr << endl << "***Error! Error opening dump file" << filename << endl << endl;
             throw runtime_error("Error RMSD dump");
			}
            m_Nf=0;
			m_dt=0.005;
		}
        virtual ~RMSD()
		{				
            std::vector<double> rotationmsd0; rotationmsd0.reserve(m_Nf);
			std::vector<double> rotationmsd1; rotationmsd1.reserve(m_Nf);
			std::vector<double> rotationmsd; rotationmsd.reserve(m_Nf);
			pos_size=m_ori_all[0].size();	
			deltat=(delta_t[1]-delta_t[0])*m_dt;
			if(m_Nf<=1000)
			{
			 n_ensemble=0.1*m_Nf;	
			}				
            else
			{
			 n_ensemble=1000;	
			}
			std::vector<vec> deltaRotation0_t; deltaRotation0_t.resize(pos_size);
			std::vector<vec> deltaRotation1_t; deltaRotation1_t.resize(pos_size);				
            for(unsigned int j=0; j<1; j++) 
			{   
				for(unsigned int i=0; i<m_Nf; i++) 	
				{
					for(unsigned int k=0; k<pos_size; k++)
					{					
					 double dotproduct=m_Rotangle_all[i][k].x*m_ori_all[j][k].x + m_Rotangle_all[i][k].y*m_ori_all[j][k].y + m_Rotangle_all[i][k].z*m_ori_all[j][k].z;					 
					 deltaRotation0_t[k].x=dotproduct*m_ori_all[j][k].x;	
					 deltaRotation0_t[k].y=dotproduct*m_ori_all[j][k].y;
					 deltaRotation0_t[k].z=dotproduct*m_ori_all[j][k].z;
					 
					 deltaRotation1_t[k].x=m_Rotangle_all[i][k].x - dotproduct*m_ori_all[j][k].x;	
					 deltaRotation1_t[k].y=m_Rotangle_all[i][k].y - dotproduct*m_ori_all[j][k].y;
					 deltaRotation1_t[k].z=m_Rotangle_all[i][k].z - dotproduct*m_ori_all[j][k].z;					
					}
					m_Rotation0_all.push_back(deltaRotation0_t);
					m_Rotation1_all.push_back(deltaRotation1_t);					
				}
			}														
            for(unsigned int i=1; i<m_Nf-n_ensemble; i++)                                                                                   
            {
				rotationmsd0[i]=0.0;
				rotationmsd1[i]=0.0;
				rotationmsd[i]=0.0;
				unsigned int count=0;
                for(unsigned int j=i; j<i+1; j++)                                                
				{
                    for(unsigned int k=0; k<pos_size; k++)
					{
//						if(m_type_all[i][k]==0)
						{
						
						 double dxtheta0_t=m_Rotation0_all[j][k].x-m_Rotation0_all[j-i][k].x;
						 double dytheta0_t=m_Rotation0_all[j][k].y-m_Rotation0_all[j-i][k].y;
						 double dztheta0_t=m_Rotation0_all[j][k].z-m_Rotation0_all[j-i][k].z;
						 rotationmsd0[i] += dxtheta0_t*dxtheta0_t+dytheta0_t*dytheta0_t+dztheta0_t*dztheta0_t;	

						 double dxtheta1_t=m_Rotation1_all[j][k].x-m_Rotation1_all[j-i][k].x;
						 double dytheta1_t=m_Rotation1_all[j][k].y-m_Rotation1_all[j-i][k].y;
						 double dztheta1_t=m_Rotation1_all[j][k].z-m_Rotation1_all[j-i][k].z;
						 rotationmsd1[i] += dxtheta1_t*dxtheta1_t+dytheta1_t*dytheta1_t+dztheta1_t*dztheta1_t;
					 	
						 double dxtheta_t=m_Rotangle_all[j][k].x-m_Rotangle_all[j-i][k].x;
						 double dytheta_t=m_Rotangle_all[j][k].y-m_Rotangle_all[j-i][k].y;
						 double dztheta_t=m_Rotangle_all[j][k].z-m_Rotangle_all[j-i][k].z;
						 rotationmsd[i] += dxtheta_t*dxtheta_t+dytheta_t*dytheta_t+dztheta_t*dztheta_t;
						 count += 1;
						}												
					}
                }				
                rotationmsd0[i] /= double(count);
				rotationmsd1[i] /= double(count);
				rotationmsd[i] /= double(count);
				m_file << i*deltat << "  " << rotationmsd0[i] <<"   "<< rotationmsd1[i] << "  " << rotationmsd[i] << endl;
			}			
            cout << "40. Good Luck! Outputting results of the rotaional mean square displacement (RMSD) to 'rmsd.log'."<<endl;			
			m_file.close();
			m_ori_all.clear();
			m_Rotation0_all.clear();
			m_Rotation1_all.clear();
			m_Rotangle_all.clear();	
            delta_t.clear();
			rotationmsd0.clear();
			rotationmsd1.clear();
			rotationmsd.clear();
			m_type_all.clear();
		};
		void setDt(double dt)
		{
		 m_dt=dt;
		}		
        virtual void compute();
    private:
        unsigned int m_Nf, n_ensemble, pos_size;
		double deltat;
		double m_dt;
        std::ofstream m_file;
		std::vector<unsigned int> delta_t;
		std::vector<std::vector<vec> > m_ori_all;
		std::vector<std::vector<vec> > m_Rotation0_all;
		std::vector<std::vector<vec> > m_Rotation1_all;
		std::vector<std::vector<vec> > m_Rotangle_all;
		std::vector< std::vector<unsigned int> > m_type_all;
    } ;
	
//---case 41
class P2P4 : public Function
	{
    public:
        P2P4(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
            {
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error P2P4 dump");
            }			
            m_Nf=0;
			m_rcut=1.0;
			m_Voronoi=false;
			maxbin=100;
			p2_all=0.0;
			p2max=1.0;
			p2min=-1.0;
			p4_all=0.0;
			p4max=1.0;
			p4min=-1.0;
        }
        virtual ~P2P4()
        {
			double dp2=(p2max-p2min)/double(maxbin);
			double dp4=(p4max-p4min)/double(maxbin);
			std::vector<double> layer2(maxbin+1, 0.0);
			std::vector<double> layer4(maxbin+1, 0.0);
			for(unsigned int i=0; i<m_Nf; i++)
			{	
				for(unsigned int k=0; k<pos_size; k++)
				{
//				 if(k==27) cout << p2_local_all[i][k] << "  " << p2min << "  " << dp2 << "    " <<  p4_local_all[i][k] << "  " << p4min << "  " << dp4 << endl;
				 unsigned int mth;
				 mth=int((p2_local_all[i][k]-p2min)/dp2); 
				 layer2[mth] += 1.0;
				 unsigned int lth;
				 lth=int((p4_local_all[i][k]-p4min)/dp4); 
				 layer4[lth] += 1.0;
				}
			}
			for(unsigned int bin=0; bin<maxbin+1; bin++)
			{
				if(bin==0) 
				 m_file << p2min+bin*dp2 << "  " << layer2[bin]/double(pos_size*m_Nf) << "  " << p4min+bin*dp4 << "  " << layer4[bin]/double(pos_size*m_Nf) << "  " << p2_all/m_Nf << "  " << p4_all/m_Nf << endl;
				else 
				 m_file << p2min+bin*dp2 << "  " << layer2[bin]/double(pos_size*m_Nf) << "  " << p4min+bin*dp4 << "  " << layer4[bin]/double(pos_size*m_Nf) << endl;
			}		
			m_file.close();	
			p2_local_all.clear();
			p4_local_all.clear();	
     	}
		void setVoronoi(bool Voronoi)
		{
		 m_Voronoi = Voronoi;
		}
		void setRcut(double rcut)
		{
		 m_rcut=rcut;
		}		
        virtual void compute();
    private:
        std::ofstream m_file;
		bool m_Voronoi;
        unsigned int m_Nf, pos_size, maxbin;
		double m_rcut, p2_all, p2max, p2min, p4_all, p4max, p4min;
		std::vector<std::vector<double> > p2_local_all, p4_local_all;
		std::vector<unsigned int> delta_t;
	};

//---case 42
class CRYSTALLINITY : public Function
	{
    public:
        CRYSTALLINITY(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
            {
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error CRYSTALLINITY dump");
            }			
            m_Nf=0;
			m_rcut=1.5;
			m_Voronoi=true;
			maxbin=100;
			m_refsij=0.6;
			m_refnxi=8;
        }
        virtual ~CRYSTALLINITY()
        {
			std::vector<unsigned int> score_s(maxbin,0);
		    std::vector <double> NUM_s(maxbin,0.0);
			std::vector <double> NXI_s(maxbin,0.0);
			for(unsigned int i=0; i<maxbin; i++)
			{	
				for(unsigned int j=0; j<m_Nf ; j++)
				{
                 score_s[i] += score_s_all[j][i];
				 NUM_s[i] += NUM_s_all[j][i];
				 NXI_s[i] += NXI_s_all[j][i];
				}
				m_file << score_s[i]/m_Nf << "  " << NUM_s[i]/m_Nf << "  " << NXI_s[i]/m_Nf << endl;					 
			}				
			m_file.close();	
			num_all.clear();
			nxi_all.clear();
			NUM_s_all.clear();
			NXI_s_all.clear();			
     	}
		void setVoronoi(bool Voronoi)
		{
		 m_Voronoi = Voronoi;
		}
		void setRcut(double rcut)
		{
		 m_rcut=rcut;
		}
		void setRefsij(double refsij)
		{
		 m_refsij=refsij;
		}
		void setRefnxi(unsigned int refnxi)
		{
		 m_refnxi=refnxi;
		}	
        virtual void compute();
    private:
        std::ofstream m_file;
		bool m_Voronoi;
        unsigned int m_Nf, pos_size, maxbin, m_refnxi;
		double m_rcut, m_refsij;
		std::vector<std::vector<unsigned int> > num_all;
		std::vector<std::vector<unsigned int> > nxi_all;		
		std::vector<unsigned int> delta_t;
		std::vector<std::vector<unsigned int> > score_s_all;
		std::vector<std::vector <double> > NUM_s_all;
		std::vector<std::vector <double> > NXI_s_all;
	};

//--- case 43                                                                                                                                   
class G6_3D : public Function
    {
    public:
        G6_3D(std::string filename): Function()
		{
			m_file.open(filename.c_str(), ios_base::out);
			if(!m_file.good())
			{
			 cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
			 throw runtime_error("Error G6_3D dump");
			}
			m_Nf=0;
			m_rcut=1.647;
			m_Voronoi=true;
			maxbin=1000;
		}
		virtual ~G6_3D() 
		{
			std::vector<double> r_sum(maxbin, 0.0);
			std::vector<double> g_sum(maxbin, 0.0);
			std::vector<double> g3D_sum(maxbin, 0.0);			
			for(unsigned int i=0; i<maxbin; i++)
			{	
				unsigned int count=0;
				for(unsigned int j=0; j<m_Nf ; j++)
				{
	             r_sum[i] += r_all[j][i];
				 g_sum[i] += g_all[j][i];
				 g3D_sum[i] += g3D_all[j][i];			 
				 count += 1;
				}
				m_file << r_sum[i]/count << "  "  << g_sum[i]/count << "  "  << g3D_sum[i]/count << "  " << g3D_sum[i]/g_sum[i] << endl;				
			}  			
			cout << "43. Good Luck! Outputting results of the spatial correlation  function of bond orintational order (G6_3D) to 'g6_3D.log'." << endl;		
			m_file.close();	
            r_all.clear();
			g_all.clear();
			g3D_all.clear();			
			
		};
		void setVoronoi(bool Voronoi)
		{
		 m_Voronoi = Voronoi;
		}
		void setRcut(double rcut)
		{
		 m_rcut=rcut;
		}		
		virtual void compute();			
    private:
		std::ofstream m_file;
		bool m_Voronoi;
		unsigned int m_Nf, pos_size, maxbin;
		double m_rcut;		
		std::vector<std::vector<double> > r_all;
		std::vector<std::vector<double> > g_all;
		std::vector<std::vector<double> > g3D_all;	
	};

//---case 44
class W4W6 : public Function
	{
    public:
        W4W6(std::string filename): Function()
        {
            m_file.open(filename.c_str(), ios_base::out);
            if(!m_file.good())
            {
             cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
             throw runtime_error("Error W4W6 dump");
            }			
            m_Nf=0;
			maxbin=100;
			q4_all=0.0;
			w4_all=0.0;
			q4max=-1.0;
			w4max=-1.0;
			q4min=1.0;
			w4min=1.0;
			q6_all=0.0;
			w6_all=0.0;
			q6max=-1.0;
			w6max=-1.0;
			q6min=1.0;
			w6min=1.0;
			m_rcut=1.647;
			m_Voronoi=false;			
        }	
        virtual ~W4W6()
        {
			double dw4=(w4max-w4min)/double(maxbin);			
			double dw6=(w6max-w6min)/double(maxbin);			
			std::vector<double> layer4(maxbin+1, 0.0);
			std::vector<double> layer6(maxbin+1, 0.0);
			for(unsigned int i=0; i<m_Nf; i++)
			{	
				for(unsigned int k=0; k<pos_size; k++)
				{
				 unsigned int mth;
				 mth=int((w4_local_all[i][k]-w4min)/dw4); 
				 layer4[mth] += 1.0;
				 unsigned int lth;
				 lth=int((w6_local_all[i][k]-w6min)/dw6); 
				 layer6[lth] += 1.0;
				}
			}
			for(unsigned int bin=0; bin<maxbin+1; bin++)
			{
				if(bin==0) 
				 m_file << w4min+bin*dw4 << "  " << layer4[bin]/double(pos_size*m_Nf) << "  " << w6min+bin*dw6 << "  " << layer6[bin]/double(pos_size*m_Nf) << "  " << w4_all/m_Nf << "  " << w6_all/m_Nf << endl;
				else 
				 m_file << w4min+bin*dw4 << "  " << layer4[bin]/double(pos_size*m_Nf) << "  " << w6min+bin*dw6 << "  " << layer6[bin]/double(pos_size*m_Nf) << endl;
			}
            cout << "44. Good Luck! Outputting results of bond order parameters (W4W6) and outputting results to 'w4w6.log'." << endl;
			m_file.close();	
			q4_local_all.clear();
			q6_local_all.clear();
     	}
		void setVoronoi(bool Voronoi)
			{
			m_Voronoi = Voronoi;
			}
		void setRcut(double rcut)
			{
			m_rcut=rcut;
			}
        virtual void compute();
    private:
        std::ofstream m_file;
		bool m_Voronoi;
        unsigned int m_Nf, pos_size, maxbin;
		double m_rcut;
		double q4_all, q4max, q4min, q6_all, q6max, q6min;
		double w4_all, w4max, w4min, w6_all, w6max, w6min;		
		std::vector<std::vector<double> > q4_local_all, w4_local_all;
		std::vector<std::vector<double> > q6_local_all, w6_local_all;
	};


//--- case 45
class MolSpt : public Function
    {
    public:
        MolSpt(): Function()
			{
			m_Nf =0;
			m_nout = 5;
			m_all = false;
			}
		virtual ~MolSpt() 
			{
			};
		virtual void compute();	
		void setOutAll(bool oa)
			{
			m_all=oa;
			}	
		void setOutNM(bool onm)
			{
			m_nout=onm;
			}		
    private:
		unsigned int m_nout;
		bool m_all;
		unsigned int m_Nf;
	};

//--- case 46
class PDI : public Function
    {
    public:
        PDI(std::string filename): Function()
			{
			m_Nf =0;
			m_file.open(filename.c_str(), ios::app);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error PDI dump");
				}
			}
		virtual ~PDI() 
			{
			cout <<"average_chain_length  "<<double(m_av_chain_length)/double(m_Nf) <<"  "<<" average_PDI  "<<double(m_av_pdi)/double(m_Nf) <<"\n";
			};
		virtual void compute();
    private:
		std::ofstream m_file;
		double m_av_pdi;
		double m_av_chain_length;
		unsigned int m_Nf;
	};
	
//--- case 47
class LocalVirial : public Function{
public:
    LocalVirial(): Function()
    {
        m_exclude = {};
        m_precision = 0.5;
        m_plane = "x-y";
        m_param = "virial";
        m_vsite = 1;
    }
    void analyse(std::vector<unsigned int>labs,std::vector< std::string > typs,std::vector<std::vector<double>>poss,std::vector<double> virs,BoxSize box,double precision_,int site1,int site2);
    void set_plane(string planev){m_plane=std::move(planev);}
    void set_param(string paramv){m_param=std::move(paramv);}
    void set_vsite(int vsitev){m_vsite=vsitev;}
    void set_precision(double precisionv){m_precision=precisionv;}
    void set_exclude(std::vector<string> excludev){m_exclude=std::move(excludev);}
    virtual ~LocalVirial(){};
    virtual void compute();
private:
    std::ofstream m_file;
    std::string m_param;
    int m_vsite;
    std::string m_plane;
    double m_precision;
    std::vector<string> m_exclude;
};

//--- case 48
class Viscosity : public Function{
public:
    Viscosity(): Function()
    {
        lx = 10.0;
        ly = 10.0;
        lz = 10.0;
        dt = 0.05;
    }
    void set_lx(double lxv){lx=lxv;}
    void set_ly(double lyv){lx=lyv;}
    void set_lz(double lzv){lx=lzv;}
    void set_dt(double dtv){dt=dtv;}
    virtual ~Viscosity(){};
    virtual void compute();
private:
    double lx, ly, lz, dt;
};

//--- case 49
class Delaunay : public Function{
public:
    Delaunay(): Function()
    {
        m_exclude = {};
    }
    void set_exclude(std::vector<string> excludev){m_exclude=std::move(excludev);}
    virtual ~Delaunay(){};
    virtual void compute();
private:
    std::ofstream m_file;
    std::vector<string> m_exclude;
};

#endif

