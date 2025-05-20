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

#include "Functions.h"
using namespace std;
std::string PYGAMD_VERSION = "1";

void output_version_info()
    {
    cout << "galaTackle -- Data Tackle Plug-Ins"<<endl;	
    cout << "PYGAMD " << PYGAMD_VERSION << endl;
	cout << "PYGAMD - Python GPU-Accelerated Molecular Dynamics Software"<<endl;
	cout << "COPYRIGHT" << endl;
	cout << "	PYGAMD Copyright (c) (2021) You-Liang Zhu, Zhong-Yuan Lu" << endl;
	cout << "LICENSE" << endl;
	cout << "	This program is a free software: you can redistribute it and/or " << endl;
	cout << "	modify it under the terms of the GNU General Public License." << endl; 
	cout << "	This program is distributed in the hope that it will be useful, " << endl;
	cout << "	but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
	cout << "	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl; 
	cout << "	See the General Public License v3 for more details." << endl;
	cout << "	You should have received a copy of the GNU General Public License" << endl;
	cout << "	along with this program. If not, see <http://www.gnu.org/licenses/>." << endl;
	cout << "DISCLAIMER" << endl;
	cout << "	The authors of PYGAMD do not guarantee that this program and its " << endl;
	cout << "	derivatives are free from error. In no event shall the copyright " << endl;
	cout << "	holder or contributors be liable for any indirect, incidental," << endl; 
	cout << "	special, exemplary, or consequential loss or damage that results " << endl;
	cout << "	from its use. We also have no responsibility for providing the " << endl;
	cout << "	service of functional extension of this program to general users." << endl;
	cout << "USER OBLIGATION " << endl;
	cout << "	If any results obtained with PYGAMD are published in the scientific " << endl;
	cout << "	literature, the users have an obligation to distribute this program " << endl;
	cout << "	and acknowledge our efforts by citing the paper \"Y.-L. Zhu et al.," << endl;
	cout << "	J. Comput. Chem. 2013, 34, 2197-2211\" in their article." << endl;
	cout << "CORRESPONDENCE" << endl;
	cout << "	Dr. You-Liang Zhu," << endl; 
	cout << "	Email: ylzhu@pygamd.com" << endl;
    }

int main(int argc,char* argv[])
	{
	std::vector<std::string> filename;
	for(unsigned int i=0; i<(unsigned int)argc;i++ )
		{
		filename.push_back(argv[i]);
		}
	output_version_info();
	cout<<"------------------------------------------------------------------"<<endl;
	cout<<"1  Rg^2                2  Ed^2               3  RDF               "<<endl;
	cout<<"4  bond_distri         5  angle_distri       6  dihedral_distri   "<<endl;
	cout<<"7  stress tensor       8  density            9  unwrapping        "<<endl;
	cout<<"10 MSD                 11 RDF-CM             12 MSD-CM            "<<endl;
	cout<<"13 ents                14 strfac             15 domain size       "<<endl;
	cout<<"16 dynamic strfac      17 config check       18 RDF between types "<<endl;
	cout<<"19 File conversion     20 patch/spot display 21 SSF               "<<endl;
	cout<<"22 ADF                 23 CND                24 MSAD              "<< endl;
	cout<<"25 RMSAD               26 ISF                27 OACF              "<< endl;
	cout<<"28 Q4Q6                29 VORONOI            30 NGP               "<< endl;
	cout<<"31 RNGP                32 VHF                33 RVHF              "<< endl;
	cout<<"34 fpSus               35 RfpSus             36 OvlaF             "<< endl;
	cout<<"37 CISF                38 CAGEISF            39 CAGEMSD           "<< endl;
	cout<<"40 RMSD                41 P2P4               42 CRYSTALLINITY     "<< endl;
	cout<<"43 G6_3D               44 W4W6               45 MolSpt            "<< endl;
	cout<<"46 CL&PDI              47 LocalVirial        48 Viscosity         "<< endl;
	cout<<"49 Delaunay                                                       "<< endl;
	cout<<"------------------------------------------------------------------"<<endl;
	cout<<"Input a number or multiple numbers separated by spaces of above list!"<<endl;

	std::vector<std::string > parameter_list;
	parameter_list.resize(50);
	parameter_list[3] =  "Parameters with default value :maxbin=100|gpu=0|rmax=Lx/2|bondex=false|angleex=false|molex=false";
	parameter_list[4] =  "Parameters with default value :npot=2001";	
	parameter_list[5] =  "Parameters with default value :npot=2001";
	parameter_list[6] =  "Parameters with default value :npot=2001";
	parameter_list[7] =  "Parameters with default value :bondex=true|bodyex=true|diameter=true";
	parameter_list[9] =  "Parameters with default value :unwrap_molecule=true|label_free_particle=particle type|molecule_center_in_box=false\n";
	parameter_list[9] += "                               |shiftx=0.0|shifty=0.0|shiftz=0.0|remove_image=false|add_image_to_pos=true\n";
	parameter_list[9] += "                               |remove_bond_cross_box=false|body_keep=false|image_integrate=false";
	parameter_list[10] = "Parameters with default value :direction=XYZ, candidates are X,Y,Z,XY,YZ,XZ,XYZ";
	parameter_list[11] = "Parameters with default value :maxbin=100|gpu=0|rmax=Lx/2";
	parameter_list[12] = "Parameters with default value :direction=XYZ, candidates are X,Y,Z,XY,YZ,XZ,XYZ";
	parameter_list[14] = "Parameters with default value :qmax=160pi/Lmin|gpu=0|deltaq=2pi/Lmin|direction=XYZ|2D=false";
	parameter_list[15] = "Parameters with default value :kmax=20|qc=0.4690|gpu=0";
	parameter_list[16] = "Parameters with default value :kmax=int(L)|q=7.0";
	parameter_list[17] = "Parameters with default value :bondex=true|angleex=true|dihedralex=true|bodyex=true|rcut=2.0";
	parameter_list[18] = "Parameters with default value :maxbin=100|gpu=0|rmax=Lx/2|bondex=false|angleex=false|molex=false";
	parameter_list[19] = "Parameters with default value :lammps=false|gromacs=false|xml=false";
	parameter_list[20] = "Parameters with default value :separatedis=0.1|filtersphere=false|patch_particle_scale=0.97";
	parameter_list[21] = "Parameters with default value :qnummax=40";
	parameter_list[22] = "Parameters with default value :beta=90|rcut=1.0";	
	parameter_list[23] = "Parameters with default value :beta=90|rcut=1.0";
	parameter_list[24] = "Parameters with default value :dt=0.005";
	parameter_list[25] = "Parameters with default value :dt=0.005";
	parameter_list[26] = "Parameters with default value :q=6.02|dt=0.005";
	parameter_list[27] = "Parameters with default value :dt=0.005";
	parameter_list[28] = "Parameters with default value :voronoi=false|rcut=1.647";
	parameter_list[29] = "Parameters with default value :rcut=1.647";
	parameter_list[30] = "Parameters with default value :dt=0.005";
	parameter_list[31] = "Parameters with default value :dt=0.005";
	parameter_list[32] = "Parameters with default value :dt=0.005";
	parameter_list[33] = "Parameters with default value :dt=0.005";
	parameter_list[34] = "Parameters with default value :q=6.02|dt=0.005";
	parameter_list[35] = "Parameters with default value :dt=0.005";
	parameter_list[36] = "Parameters with default value :a=0.30|dt=0.005";
	parameter_list[37] = "Parameters with default value :q=6.02|dt=0.005";
	parameter_list[38] = "Parameters with default value :q=6.02|dt=0.005|voronoi=false|rcut=1.65";	
	parameter_list[39] = "Parameters with default value :dt=0.005|voronoi=false|rcut=1.65";
	parameter_list[40] = "Parameters with default value :dt=0.005";
	parameter_list[41] = "Parameters with default value :voronoi=false|rcut=1.0";
	parameter_list[42] = "Parameters with default value :voronoi=true|rcut=1.5|refsij=0.6|refnxi=8";
	parameter_list[43] = "Parameters with default value :maxbin=1000|voronoi=true|rcut=1.647";
	parameter_list[44] = "Parameters with default value :voronoi=false|rcut=1.647";
	parameter_list[45] = "Parameters with default value :nmout=5|outall=false";
	parameter_list[46] = "";
	parameter_list[47] = "Parameters with default value :virial=1|plane=x-y|precision=0.5|exclude={}\n"
                             "Example:                    46:virial_matrix=5|plane=x-z|precision=1.0|exclude=Li-Na-K";
	parameter_list[48] = "Parameters with default value :lx=10.0|ly=10.0|lz=10.0|dt=0.05\n"
                         "    Notice ！  The required execution file is data.log";
	parameter_list[49] = "Parameters with default value :exclude={}\n"
                         "Example:                    :exclude=Li-Na-K";
	std::vector<unsigned int > command;
	std::vector<std::string > command1;
	std::vector<std::string> command2;
	while(command.size()==0)
		{
		std::string line; 
		std::cout<<">>>";
		getline(cin, line); 
		std::string temp, temp1;
		
		for(unsigned int i =0; i< line.size(); i++)
			{
			if(line.at(i)!=','&&line.at(i)!=' '&&line.at(i)!=':')
				temp.push_back(line.at(i));
			else if(line.at(i)==':')
				{
				while(i+1!=line.size()&&line.at(i+1)!=','&&line.at(i+1)!=' ')
					{
					i+=1;				
					temp1.push_back(line.at(i));
					}
				}
			if(line.at(i)==','||line.at(i)==' '||i ==(line.size()-1))
				{		
				if(temp.size()!=0)
					{
					unsigned int id;
					stringstream ss(temp);
					ss>>id;
					command.push_back(id);
					command1.push_back(temp);
					command2.push_back(temp1);
					}
				else if(temp1.size()!=0)
					{
					cout<<"Warning!! ignore the parameters '"<<temp1<<"'"<<endl;
					}
				temp.clear();
				temp1.clear();
				}
			}
	
		for(unsigned int i=0;i<command.size();)
			{
			if(command[i]==0||command[i]>49)
				{
				cout<<"Error!! Please input the numbers in the list! '"<<command1[i]<<"' do not work!"<<endl;
				return 0;
				}
			if(command2[i]=="help"||command2[i]=="h"||command2[i]=="Help"||command2[i]=="H")
				{
				std::cout<<command[i]<<". "<<parameter_list[command[i]]<<std::endl;
				std::vector<unsigned int >::iterator it = command.begin()+i;
				command.erase(it);
				std::vector<std::string >::iterator it1 = command1.begin()+i;
				command1.erase(it1);
				std::vector<std::string >::iterator it2 = command2.begin()+i;	
				command2.erase(it2);
				}
			else
				i++;
			}
		}
	std::cout<<"-----------------------------------------"<<endl;
	std::vector<shared_ptr<Function> > object;
	for(unsigned int i=0;i<command.size();i++)
		{
		switch(command[i])
			{
			case 1:
				{
				shared_ptr<Rg2> rg2 = shared_ptr<Rg2>(new Rg2("rg2.log"));
				object.push_back(rg2);
				cout<<"1.  Computing the square of radius of gyration and output result to 'rg2.log'"<<endl;
				break;
				}
			case 2:
				{
				shared_ptr<Ed2> ed2 = shared_ptr<Ed2>(new Ed2("ed2.log"));
				object.push_back(ed2);
				cout<<"2.  Computing the square of distance of end-to-end and output result to 'ed2.log'"<<endl;
				break;			
				}
			case 3:
				{	
				shared_ptr<RDF> rdf = shared_ptr<RDF>(new RDF("rdf.log"));
				object.push_back(rdf);
				cout<<"3.  Computing radial distribution function and output result to '*.rdf' and 'rdf.log'"<<endl;
				cout<<"    "<<parameter_list[3]<<endl;					
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="maxbin")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<1000000)
								{
								rdf->setMaxbin(val);
								cout<<"Note! Change maxbin="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'maxbin="<<val<<"'"<<endl;
							}
						else if(param=="gpu")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0||val<=16)
								{
								rdf->setGPU(val);
								cout<<"Note! Change gpu="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'gpu="<<val<<"'"<<endl;
							}
						else if(param=="rmax")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								rdf->setRmax(val);
								cout<<"Note! Change rmax="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rmax="<<val<<"'"<<endl;
							}
						else if(param=="bondex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									rdf->setBondEx(true);
								else
									rdf->setBondEx(false);
								cout<<"Note! Change bondex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'bondex="<<value<<"'"<<endl;
							}
						else if(param=="angleex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									rdf->setAngleEx(true);
								else
									rdf->setAngleEx(false);
								cout<<"Note! Change angleex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'angleex="<<value<<"'"<<endl;
							}
						else if(param=="molex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									rdf->setMolEx(true);
								else
									rdf->setMolEx(false);
								cout<<"Note! Change molex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'molex="<<value<<"'"<<endl;
							}								
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;			
				}
			case 4:
				{	
				shared_ptr<Bond_distr> bond_distr = shared_ptr<Bond_distr>(new Bond_distr("bond_distr.log"));
				object.push_back(bond_distr);
				cout<<"4.  Computing bond length distribution and output result to 'bond_distr.log'"<<endl;
				cout<<"    "<<parameter_list[4]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="npot")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0)
								{
								bond_distr->setNpot(val);
								cout<<"Note! Change npot="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;				
				}
			case 5:
				{		
				shared_ptr<Angle_distr> angle_distr = shared_ptr<Angle_distr>(new Angle_distr("angle_distr.log"));
				object.push_back(angle_distr);
				cout<<"5.  Computing angle radian distribution and output result to 'angle_distr.log'"<<endl;
				cout<<"    "<<parameter_list[5]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="npot")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0)
								{
								angle_distr->setNpot(val);
								cout<<"Note! Change npot="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;			
				}
			case 6:
				{		
				shared_ptr<Dihedral_distr> dihedral_distr = shared_ptr<Dihedral_distr>(new Dihedral_distr("dihedral_distr.log"));
				object.push_back(dihedral_distr);
				cout<<"6.  Computing dihedral radian distribution and output result to 'dihedral_distr.log'"<<endl;
				cout<<"    "<<parameter_list[6]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="npot")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0)
								{
								dihedral_distr->setNpot(val);
								cout<<"Note! Change npot="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;			
				}
			case 7:
				{
				shared_ptr<StressTensor> stress_tensor = shared_ptr<StressTensor>(new StressTensor("stress_tensor.log"));
				object.push_back(stress_tensor);
				cout<<"7.  Computing stress tensor and output result to 'stress_tensor.log'"<<endl;
				cout<<"    "<<parameter_list[7]<<endl;	
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="bondex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									stress_tensor->setBondEx(true);
								else
									stress_tensor->setBondEx(false);
								cout<<"Note! Change bondex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'bondex="<<value<<"'"<<endl;
							}
						else if(param=="bodyex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									stress_tensor->setBodyEx(true);
								else
									stress_tensor->setBodyEx(false);
								cout<<"Note! Change bodyex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'bodyex="<<value<<"'"<<endl;
							}
						else if(param=="diameter")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									stress_tensor->setDiameterConsider(true);
								else
									stress_tensor->setDiameterConsider(false);
								cout<<"Note! Change bodyex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'diameter="<<value<<"'"<<endl;
							}
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;	
				}	
			case 8:
				{
				shared_ptr<Density> density = shared_ptr<Density>(new Density("density.log"));
				object.push_back(density);
				cout<<"8.  Computing density and output result to 'density.log'"<<endl;
				break;	
				}
			case 9:
				{		
				shared_ptr<Reimage> re_image = shared_ptr<Reimage>(new Reimage());
				object.push_back(re_image);
				cout<<"9.  Computing the image for unfolding molecules and output to '*.reimage.xml(or .mst)' files"<<endl;
				cout<<"    "<<parameter_list[9]<<endl;						
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="shiftx")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							re_image->setShiftX(val);
							cout<<"Note! Change shiftx="<<val<<" successfully."<<endl;
							}
						else if(param=="shifty")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							re_image->setShiftY(val);
							cout<<"Note! Change shifty="<<val<<" successfully."<<endl;
							}
						else if(param=="shiftz")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							re_image->setShiftZ(val);
							cout<<"Note! Change shiftz="<<val<<" successfully."<<endl;
							}
						else if(param=="unwrap_molecule")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									re_image->setUnwrapMolecule(true);
								else
									re_image->setUnwrapMolecule(false);
								cout<<"Note! Change unwrap_molecule="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'unwrap_molecule="<<value<<"'"<<endl;
							}
						else if(param=="molecule_center_in_box")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									re_image->setMoleculeCenterInBox(true);
								else
									re_image->setMoleculeCenterInBox(false);
								cout<<"Note! Change molecule_center_in_box="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'molecule_center_in_box="<<value<<"'"<<endl;
							}
						else if(param=="label_free_particle")
							{
							re_image->setLabelFreeParticle(value);
							cout<<"Note! Change label_free_particle="<<value<<" successfully."<<endl;
							}
						else if(param=="remove_image")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									re_image->setRemoveImage(true);
								else
									re_image->setRemoveImage(false);
								cout<<"Note! Change remove_image="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'remove_image="<<value<<"'"<<endl;
							}
						else if(param=="add_image_to_pos")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									re_image->addImageToPos(true);
								else
									re_image->addImageToPos(false);
								cout<<"Note! Change add_image_to_pos="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'add_image_to_pos="<<value<<"'"<<endl;									
							}								
						else if(param=="remove_bond_cross_box")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									re_image->setRemoveBondCrossBox(true);
								else
									re_image->setRemoveBondCrossBox(false);
								cout<<"Note! Change remove_bond_cross_box="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'remove_bond_cross_box="<<value<<"'"<<endl;									
							}
						else if(param=="body_keep")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									re_image->setBodyKeep(true);
								else
									re_image->setBodyKeep(false);
								cout<<"Note! Change body_keep="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'body_keep="<<value<<"'"<<endl;									
							}
						else if(param=="image_integrate")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									re_image->setImageIntegrate(true);
								else
									re_image->setImageIntegrate(false);
								cout<<"Note! Change image_integrate="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'image_integrate="<<value<<"'"<<endl;									
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;	
				}
			case 10:
				{
				shared_ptr<MSD> msd = shared_ptr<MSD>(new MSD("msd.log"));
				object.push_back(msd);
				cout<<"10.  Computing mean square displacement with low resolution and output result to 'msd.log'"<<endl;
				cout<<"    "<<parameter_list[10]<<endl;						
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="direction")
							{
							if(value=="X"||value=="Y"||value=="Z"||value=="XY"||value=="YZ"||value=="XZ"||value=="XYZ")
								{
								msd->setDirection(value);
								cout<<"Note! Change direction="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the wrong value 'direction="<<value<<"'"<<endl;
							}
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}					
				break;		
				}	
			case 11:
				{	
				shared_ptr<RDFCM> rdf_cm = shared_ptr<RDFCM>(new RDFCM("rdf_cm.log"));
				object.push_back(rdf_cm);
				cout<<"11.  Computing the radial distribution function of center of mass and output result to 'rdf_cm.log'"<<endl;
				cout<<"    "<<parameter_list[11]<<endl;						
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="maxbin")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<1000000)
								{
								rdf_cm->setMaxbin(val);
								cout<<"Note! Change maxbin="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'maxbin="<<val<<"'"<<endl;
							}
						else if(param=="gpu")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0||val<=16)
								{
								rdf_cm->setGPU(val);
								cout<<"Note! Change gpu="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'gpu="<<val<<"'"<<endl;
							}
						else if(param=="rmax")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								rdf_cm->setRmax(val);
								cout<<"Note! Change rmax="<<val<<" successfully."<<endl;
								}
							}
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;				
				}
			case 12:
				{	
				shared_ptr<MSDCM> msd_cm = shared_ptr<MSDCM>(new MSDCM("msd_cm.log"));
				object.push_back(msd_cm);
				cout<<"12.  Computing the mean square displacement of center mass and output result to 'msd_cm.log'"<<endl;
				cout<<"    "<<parameter_list[12]<<endl;					
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="direction")
							{
							if(value=="X"||value=="Y"||value=="Z"||value=="XY"||value=="YZ"||value=="XZ"||value=="XYZ")
								{
								msd_cm->setDirection(value);
								cout<<"Note! Change direction="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the wrong value 'direction="<<value<<"'"<<endl;
							}
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}					
				break;				
				}
			case 13:
				{	
				shared_ptr<Entanglement> ents = shared_ptr<Entanglement>(new Entanglement("ents.log"));
				object.push_back(ents);
				cout<<"13.  Computing the number of entanglements and output result to 'ents.log'"<<endl;
				break;				
				}
			case 14:
				{	
				shared_ptr<STRFAC> strf = shared_ptr<STRFAC>(new STRFAC("strf.log"));
				object.push_back(strf);
				cout<<"14.  Computing structure factor and output result to '*.strf' and 'strf.log'"<<endl;
				cout<<"    "<<parameter_list[14]<<endl;						
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="qmax")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<300)
								{
								strf->setQmax(val);
								cout<<"Note! Change qmax="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'qmax="<<val<<"'"<<endl;
							}
						else if(param=="deltaq")
							{
							float val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<300)
								{
								strf->setDeltaq(val);
								cout<<"Note! Change deltaq="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'deltaq="<<val<<"'"<<endl;
							}								
						else if(param=="gpu")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0||val<=16)
								{
								strf->setGPU(val);
								cout<<"Note! Change gpu="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'gpu="<<val<<"'"<<endl;
							}
						else if(param=="direction")
							{
							if(value=="X"||value=="Y"||value=="Z"||value=="XY"||value=="YZ"||value=="XZ"||value=="XYZ")
								{
								strf->setDirection(value);
								cout<<"Note! Change direction="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the wrong value 'direction="<<value<<"'"<<endl;
							}
						else if(param=="2D")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									strf->set2D(true);
								else
									strf->set2D(false);
								cout<<"Note! Change 2D="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value '2D="<<value<<"'"<<endl;
							}								
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;			
				}
			case 15:
				{	
				shared_ptr<DOMAINSIZE> domsize = shared_ptr<DOMAINSIZE>(new DOMAINSIZE("domsize.log"));
				object.push_back(domsize);
				cout<<"15.  Computing domain size and output result to 'domsize.log'"<<endl;
				cout<<"    "<<parameter_list[15]<<endl;						
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="kmax")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<300)
								{
								domsize->setKmax(val);
								cout<<"Note! Change kmax="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'kmax="<<val<<"'"<<endl;
							}
						else if(param=="qc")
							{
							float val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<1000)
								{
								domsize->setQc(val);
								cout<<"Note! Change qc="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'qc="<<val<<"'"<<endl;
							}
						else if(param=="gpu")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0||val<=16)
								{
								domsize->setGPU(val);
								cout<<"Note! Change gpu="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'gpu="<<val<<"'"<<endl;
							}									
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;			
				}
			case 16:
				{	
				shared_ptr<DSTRFAC> dstrf = shared_ptr<DSTRFAC>(new DSTRFAC("dstrf.log"));
				object.push_back(dstrf);
				cout<<"16.  Computing dynamic structure factor and output result to 'dstrf.log'"<<endl;
				cout<<"    "<<parameter_list[16]<<endl;					
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="kmax")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<300)
								{
								dstrf->setKmax(val);
								cout<<"Note! Change kmax="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'kmax="<<val<<"'"<<endl;
							}
						else if(param=="q")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<100000)
								{
								dstrf->setQ(val);
								cout<<"Note! Change q="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'q="<<val<<"'"<<endl;
							}								
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;			
				}
			case 17:
				{
				shared_ptr<ConfigCheck> config_check = shared_ptr<ConfigCheck>(new ConfigCheck("config_check.log"));
				object.push_back(config_check);
				cout<<"17.  Check the configuration and output result to 'config_check.log'"<<endl;
				cout<<"    "<<parameter_list[17]<<endl;	
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="bondex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									config_check->setBondEx(true);
								else
									config_check->setBondEx(false);
								cout<<"Note! Change bondex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'bondex="<<value<<"'"<<endl;
							}
						else if(param=="angleex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									config_check->setAngleEx(true);
								else
									config_check->setAngleEx(false);
								cout<<"Note! Change angleex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'angleex="<<value<<"'"<<endl;
							}
						else if(param=="dihedralex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									config_check->setDihedralEx(true);
								else
									config_check->setDihedralEx(false);
								cout<<"Note! Change dihedralex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dihedralex="<<value<<"'"<<endl;
							}
						else if(param=="bodyex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									config_check->setBodyEx(true);
								else
									config_check->setBodyEx(false);
								cout<<"Note! Change bodyex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'bodyex="<<value<<"'"<<endl;
							}
						else if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0)
								{
								config_check->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;			
				}
			case 18:
				{	
				shared_ptr<RDFBetweenTypes> rdftypes = shared_ptr<RDFBetweenTypes>(new RDFBetweenTypes("rdf_by_type.log"));
				object.push_back(rdftypes);
				cout<<"18.  Computing radial distribution function between types and output result to '*.type.rdf' and 'rdf_by_type.log'"<<endl;
				cout<<"    "<<parameter_list[18]<<endl;						
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="maxbin")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<1000000)
								{
								rdftypes->setMaxbin(val);
								cout<<"Note! Change maxbin="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'maxbin="<<val<<"'"<<endl;
							}
						else if(param=="gpu")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0||val<=16)
								{
								rdftypes->setGPU(val);
								cout<<"Note! Change gpu="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'gpu="<<val<<"'"<<endl;
							}
						else if(param=="rmax")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								rdftypes->setRmax(val);
								cout<<"Note! Change rmax="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rmax="<<val<<"'"<<endl;
							}
						else if(param=="bondex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									rdftypes->setBondEx(true);
								else
									rdftypes->setBondEx(false);
								cout<<"Note! Change bondex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'bondex="<<value<<"'"<<endl;
							}
						else if(param=="angleex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									rdftypes->setAngleEx(true);
								else
									rdftypes->setAngleEx(false);
								cout<<"Note! Change angleex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'angleex="<<value<<"'"<<endl;
							}
						else if(param=="molex")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									rdftypes->setMolEx(true);
								else
									rdftypes->setMolEx(false);
								cout<<"Note! Change molex="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'molex="<<value<<"'"<<endl;
							}								
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;			
				}
			case 19:
				{
				shared_ptr<FileConversion> mst_file_conversion = shared_ptr<FileConversion>(new FileConversion());
				object.push_back(mst_file_conversion);
				cout<<"19.  File conversion"<<endl;
				cout<<"    "<<parameter_list[19]<<endl;					
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="lammps")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									mst_file_conversion->setLammps(true);
								else
									mst_file_conversion->setLammps(false);
								cout<<"Note! Change lammps="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'lammps="<<value<<"'"<<endl;
							}
						else if(param=="gromacs")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									mst_file_conversion->setGromacs(true);
								else
									mst_file_conversion->setGromacs(false);
								cout<<"Note! Change gromacs="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'gromacs="<<value<<"'"<<endl;
							}
						else if(param=="xml")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									mst_file_conversion->setXML(true);
								else
									mst_file_conversion->setXML(false);
								cout<<"Note! Change xml="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'xml="<<value<<"'"<<endl;
							}									
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;
				}
			case 20:
				{	
				shared_ptr<PatchToParticle> patch = shared_ptr<PatchToParticle>(new PatchToParticle());
				object.push_back(patch);
				cout<<"20.  Convert patches to particles for display"<<endl;
				cout<<"    "<<parameter_list[20]<<endl;						
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="separatedis")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<1000000)
								{
								patch->setSeparateDis(val);
								cout<<"Note! Change separatedis="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'separatedis="<<val<<"'"<<endl;
							}
						else if(param=="filtersphere")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									patch->setFilterSphere(true);
								else
									patch->setFilterSphere(false);
								cout<<"Note! Change filtersphere="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'filtersphere="<<value<<"'"<<endl;
							}
						else if(param=="patch_particle_scale")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<1000000)
								{
								patch->setPatchParticleScale(val);
								cout<<"Note! Change patch_particle_scale="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'patch_particle_scale="<<val<<"'"<<endl;
							}								
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;			
				}
			case 21:
				{
				shared_ptr<SSF> ssf = shared_ptr<SSF>(new SSF("ssf.log"));
				object.push_back(ssf);
				cout<<"21. Computing the static structure factor (SSF) and outputting results to 'ssf.log'." << endl;
				cout<<"    "<<parameter_list[21]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="qnummax")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0&&val<200)
								{
								ssf->setqnummax(val);
								cout<<"Note! Change qnummax="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'qnummax="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 22:                                                                                                                         
				{	
				shared_ptr<ADF> adf=shared_ptr<ADF>(new ADF("adf.log"));
				object.push_back(adf);
				cout<<"22. Computing the angular distribution function (ADF) and outputting results to 'adf.log'." << endl;
				cout<<"    "<<parameter_list[22]<<endl;				
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="beta")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0.0&&val<180.0)
								{
								adf->setBeta(val);
								cout<<"Note! Change beta="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'beta="<<val<<"'"<<endl;
							}
						else if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								adf->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}			 
				break;				
				}				
			case 23:                                                                                                                         
				{	
				shared_ptr<CND> cnd=shared_ptr<CND>(new CND("cnd.log"));
				object.push_back(cnd);
				cout<<"23. Computing the contact number distribution (CND) and outputting results to 'cnd.log'." << endl;
				cout<<"    "<<parameter_list[23]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="beta")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>=0.0&&val<180.0)
								{
								cnd->setBeta(val);
								cout<<"Note! Change beta="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'beta="<<val<<"'"<<endl;
							}
						else if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								cnd->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;				
				}
			case 24:
				{
				shared_ptr<MSAD> msad=shared_ptr<MSAD>(new MSAD("msad.log"));
				object.push_back(msad);
				cout<<"24. Computing the mean square angular displacement (MSAD) and outputting results to 'msad.log'." << endl;
				cout<<"    "<<parameter_list[24]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								msad->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}						
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 25:
				{
				shared_ptr<RMSAD> rmsad=shared_ptr<RMSAD>(new RMSAD("rmsad.log"));
				object.push_back(rmsad);
				cout<<"25. Computing the reorientaional mean square angular displacement (RMSAD) and outputting results to 'rmsad.log'."<<endl;
				cout<<"    "<<parameter_list[25]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								rmsad->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}						
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}					 				 
				break;
				}
			case 26:
				{
				shared_ptr<ISF> isf = shared_ptr<ISF>(new ISF("isf.log"));
				object.push_back(isf);
				cout<<"26. Computing the self part intermediate scattering function (ISF) and outputting results to 'isf.log'." << endl;
				cout<<"    "<<parameter_list[26]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								isf->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}
						else if(param=="q")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								isf->setQ(val);
								cout<<"Note! Change q="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'q="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 27:
				{
				shared_ptr<OACF> oacf=shared_ptr<OACF>(new OACF("oacf.log"));
				object.push_back(oacf);
				cout<<"27. Computing the orientational autocorrelation function (OACF) and outputting results to 'oacf.log'." << endl;
				cout<<"    "<<parameter_list[27]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								oacf->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 28:
				{
				shared_ptr<Q4Q6> q4q6 = shared_ptr<Q4Q6>(new Q4Q6("q4q6.log"));
				object.push_back(q4q6);
				cout<<"28. Computing bond order parameters (Q4Q6) and outputting results to 'q4q6.log'." << endl;
				cout<<"    "<<parameter_list[28]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="voronoi")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									q4q6->setVoronoi(true);
								else
									q4q6->setVoronoi(false);
								cout<<"Note! Change voronoi="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'voronoi="<<value<<"'"<<endl;
							}							
						else if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								q4q6->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 				 
				break;
				}
			case 29:
				{
				shared_ptr<VORONOI> voronoi = shared_ptr<VORONOI>(new VORONOI("voronoi.log"));
				object.push_back(voronoi);
				cout<<"29. Computing the volume of the Voronoi cells (VORONOI) and outputting results to 'voronoi.log'." << endl;
				cout<<"    "<<parameter_list[29]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								voronoi->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}						
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 30:
				{
				shared_ptr<nonGauPar> nongaupar = shared_ptr<nonGauPar>(new nonGauPar("nongaupar.log"));
				object.push_back(nongaupar);
				cout<<"30. Computing non-Gaussian parameter (NGP) and outputting results to 'nongaupar.log'" << endl;
				cout<<"    "<<parameter_list[30]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								nongaupar->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}						
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 31:
				{
				shared_ptr<RnonGauPar> rnongaupar = shared_ptr<RnonGauPar>(new RnonGauPar("rnongaupar.log"));
				object.push_back(rnongaupar);
				cout<<"31. Computing the rotational non-Gaussian parameter (RNGP) and outputting results to 'rnongaupar.log'."<<endl;
				cout<<"    "<<parameter_list[31]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								rnongaupar->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}						
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 32:
				{
				shared_ptr<SVH> svh = shared_ptr<SVH>(new SVH("selfvhf.log"));
				object.push_back(svh);
				cout<<"32. Computing the self part of van Hove fucntion(VHF) and outputting results to 'selfvhf.log'." << endl;
				cout<<"    "<<parameter_list[32]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								svh->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}						
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 33:
				{
				shared_ptr<RSVH> rsvh = shared_ptr<RSVH>(new RSVH("rselfvhf.log"));
				object.push_back(rsvh);
				cout<<"33. Computing the self part of rotational van Hove fucntion (RVHF) and outputting results to 'rselfvhf.log'." << endl;
				cout<<"    "<<parameter_list[33]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								rsvh->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}						
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 34:
				{
				shared_ptr<fpSus> fpsus = shared_ptr<fpSus>(new fpSus("fpsus.log"));
				object.push_back(fpsus);
				cout<<"34. Computing the four-point susceptibility and outputting results to 'fpsus.log'"<<endl;
				cout<<"    "<<parameter_list[34]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								fpsus->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}
						else if(param=="q")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								fpsus->setQ(val);
								cout<<"Note! Change q="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'q="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 35:
				{
				shared_ptr<RfpSus> rfpsus = shared_ptr<RfpSus>(new RfpSus("rfpsus.log"));
				object.push_back(rfpsus);
				cout<<"35. Computing the rotational four-point susceptibility and outputting results to 'rfpsus.log'"<<endl;
				cout<<"    "<<parameter_list[35]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								rfpsus->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}				
			case 36:
				{
				shared_ptr<OVLAF> ovlaf = shared_ptr<OVLAF>(new OVLAF("ovlaf.log"));
				object.push_back(ovlaf);
				cout<<"36. Computing self part of the overlap function (OVLAF) and outputting results to 'ovlaf.log'." << endl;
				cout<<"    "<<parameter_list[36]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								ovlaf->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}
						else if(param=="a")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								ovlaf->setA(val);
								cout<<"Note! Change a="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'a="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 37:
				{
				shared_ptr<CISF> cisf = shared_ptr<CISF>(new CISF("cisf.log"));
				object.push_back(cisf);
				cout<<"37. Computing the coherent intermediate scattering function (CISF) and outputting results to 'cisf.log'." << endl;
				cout<<"    "<<parameter_list[37]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								cisf->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}
						else if(param=="q")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								cisf->setQ(val);
								cout<<"Note! Change q="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'q="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 38:
				{
				shared_ptr<CAGEISF> cageisf = shared_ptr<CAGEISF>(new CAGEISF("cageisf.log"));
				object.push_back(cageisf);
				cout<<"38. Computing the cage-relative self part intermediate scattering function (CAGEISF) and outputting results to 'cageisf.log'." << endl;
				cout<<"    "<<parameter_list[38]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								cageisf->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}
						else if(param=="q")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								cageisf->setQ(val);
								cout<<"Note! Change q="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'q="<<val<<"'"<<endl;
							}							
						else if(param=="voronoi")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									cageisf->setVoronoi(true);
								else
									cageisf->setVoronoi(false);
								cout<<"Note! Change voronoi="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'voronoi="<<value<<"'"<<endl;
							}							
						else if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								cageisf->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 39:
				{
				shared_ptr<CAGEMSD> cagemsd = shared_ptr<CAGEMSD>(new CAGEMSD("cagemsd.log"));
				object.push_back(cagemsd);
				cout<<"39. Computing the cage-relative mean square displacement (CAGEMSD) and outputting results to 'cagemsd.log'." << endl;
				cout<<"    "<<parameter_list[39]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								cagemsd->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}
						else if(param=="voronoi")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									cagemsd->setVoronoi(true);
								else
									cagemsd->setVoronoi(false);
								cout<<"Note! Change voronoi="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'voronoi="<<value<<"'"<<endl;
							}							
						else if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								cagemsd->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}					 
				break;
				}
			case 40:
				{
				shared_ptr<RMSD> rmsd=shared_ptr<RMSD>(new RMSD("rmsd.log"));
				object.push_back(rmsd);
				cout<<"40. Computing the rotaional mean square displacement (RMSD) and outputting results to 'rmsd.log'." << endl;
				cout<<"    "<<parameter_list[40]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="dt")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								rmsd->setDt(val);
								cout<<"Note! Change dt="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'dt="<<val<<"'"<<endl;
							}						
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 
				break;
				}
			case 41:
				{
				shared_ptr<P2P4> p2p4 = shared_ptr<P2P4>(new P2P4("p2p4.log"));
				object.push_back(p2p4);
				cout<<"41. Computing orientational order parameters (P2P4) and outputting results to 'p2p4.log'." << endl;
				cout<<"    "<<parameter_list[41]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="voronoi")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									p2p4->setVoronoi(true);
								else
									p2p4->setVoronoi(false);
								cout<<"Note! Change voronoi="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'voronoi="<<value<<"'"<<endl;
							}							
						else if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								p2p4->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 				 
				break;
				}
			case 42:
				{
				shared_ptr<CRYSTALLINITY> crystallinity = shared_ptr<CRYSTALLINITY>(new CRYSTALLINITY("crystallinity.log"));
				object.push_back(crystallinity);
				cout<<"42. Computing the crystallinity (CRYSTALLINITY) and outputting results to 'crystallinity.log'." << endl;
				cout<<"    "<<parameter_list[42]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="voronoi")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									crystallinity->setVoronoi(true);
								else
									crystallinity->setVoronoi(false);
								cout<<"Note! Change voronoi="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'voronoi="<<value<<"'"<<endl;
							}							
						else if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								crystallinity->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}
						else if(param=="refsij")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								crystallinity->setRefsij(val);
								cout<<"Note! Change refsij="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'refsij="<<val<<"'"<<endl;
							}
						else if(param=="refnxi")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								crystallinity->setRefnxi(val);
								cout<<"Note! Change refnxi="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'refnxi="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 				 
				break;
				}
			case 43:
				{	
				shared_ptr<G6_3D> g6_3D = shared_ptr<G6_3D>(new G6_3D("g6_3D.log"));		   
				object.push_back(g6_3D);
				cout<<"43. Computing the spatial correlation  function of bond orintational order (G6_3D) and outputting the results to 'g6_3D.log'." << endl;
				cout<<"    "<<parameter_list[43]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="maxbin")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0)
								{
								g6_3D->setRcut(val);
								cout<<"Note! Change maxbin="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'maxbin="<<val<<"'"<<endl;
							}
						else if(param=="voronoi")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									g6_3D->setVoronoi(true);
								else
									g6_3D->setVoronoi(false);
								cout<<"Note! Change voronoi="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'voronoi="<<value<<"'"<<endl;
							}							
						else if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								g6_3D->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}					 
				break;
				}
			case 44:
				{
				shared_ptr<W4W6> w4w6 = shared_ptr<W4W6>(new W4W6("w4w6.log"));
				object.push_back(w4w6);
				cout<<"44. Computing bond order parameters (W4W6) and outputting results to 'w4w6.log'." << endl;
				cout<<"    "<<parameter_list[44]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="voronoi")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									w4w6->setVoronoi(true);
								else
									w4w6->setVoronoi(false);
								cout<<"Note! Change voronoi="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'voronoi="<<value<<"'"<<endl;
							}							
						else if(param=="rcut")
							{
							double val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0.0)
								{
								w4w6->setRcut(val);
								cout<<"Note! Change rcut="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'rcut="<<val<<"'"<<endl;
							}							
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}				 				 
				break;
				}
			case 45:
				{
				shared_ptr<MolSpt> mos = shared_ptr<MolSpt>(new MolSpt());
				object.push_back(mos);
				cout<<"45. Separate molecules into files." << endl;
				cout<<"    "<<parameter_list[45]<<endl;
				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|'))
					{
					string::size_type pos=sub_str.find("=",0);
					if(pos!=sub_str.npos)
						{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if(param=="outall")
							{
							if(value=="true"||value=="false")
								{
								if(value=="true")
									mos->setOutAll(true);
								else
									mos->setOutAll(false);
								cout<<"Note! Change outall="<<value<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'outall="<<value<<"'"<<endl;
							}							
						else if(param=="nmout")
							{
							unsigned int val;
							stringstream ss2(value);
							ss2>>val;
							if(val>0 && val <1000000)
								{
								mos->setOutNM(val);
								cout<<"Note! Change nmout="<<val<<" successfully."<<endl;
								}
							else
								cout<<"Warning!! Ignore the unreasonable value 'nmout="<<val<<"'"<<endl;
							}						
						else
							cout<<"Warning!! Ignore the unknown parameter '"<<param<<"'"<<endl;
						}
					else
						cout<<"Warning!! Ignore the nonstandard input '"<<sub_str<<"'"<<endl;
					}
				break;
				}
			case 46:
				{
				shared_ptr<PDI> pdi = shared_ptr<PDI>(new PDI("pdi.log"));
				object.push_back(pdi);
				cout<<"46.  Computing the average length and PDI of polymer chains and output result to 'pdi.log'"<<endl;
				break;
				}
			case 47:
				{
                shared_ptr<LocalVirial> localvirial = shared_ptr<LocalVirial>(new LocalVirial());
                object.push_back(localvirial);
                cout<<"47.  Computing the local virial and output result to file(localvirial.log or localvirial_matrix.log)" << endl;
                cout<<"     "<<parameter_list[47]<<endl;

                std::string text = command2[i];
                stringstream ss1(text);
                std::string sub_str;
                while(getline(ss1,sub_str,'|')){
                    string::size_type pos=sub_str.find('=',0);
                    if(pos!=std::string::npos)
						{
							string param = sub_str.substr(0, pos);
							string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
							if (param=="virial_matrix"){
								int val;
								stringstream ss2(value);
								ss2>>val;
								localvirial->set_param(param);
								if (val <= 0) val=1;
								if (val > 6) val=6;
								localvirial->set_vsite(val);
								cout<<"Note! " + param +"="<< val <<" successfully."<<endl;
							}else if(param=="plane"){
								localvirial->set_plane(value);
								cout<<"Note! Change sel="<<value<<" successfully."<<endl;
							} else if (param=="precision"){
								double val;
								stringstream ss2(value);
								ss2>>val;
								localvirial->set_precision(val);
								cout<<"Note! Change precision="<<val<<" successfully."<<endl;
							}else if (param=="exclude"){
								std::vector<string> val;
								stringstream ss2(value);
								std::string sub_str2;
								std::string s2 = "{";
								while (getline(ss2,sub_str2,'-')){
									val.push_back(sub_str2);
									s2 += sub_str2 + ", ";
								}
								s2 = s2.erase(s2.size()-2, s2.size());
								s2 += "}";
								localvirial->set_exclude(val);
								cout<<"Note! Change exclude="<< s2 <<" successfully."<<endl;
							}
						}
					}
                break;
				}
			case 48:
				{
				shared_ptr<Viscosity> viscosity = std::make_shared<Viscosity>();
				object.push_back(viscosity);
				cout<<"48.  Computing the viscosity and output result to file(viscosity-thermo.log and phi-thermo.log)" << endl;
				cout<<"     "<<parameter_list[48]<<endl;

				std::string text = command2[i];
				stringstream ss1(text);
				std::string sub_str;
				while(getline(ss1,sub_str,'|')){
					string::size_type pos=sub_str.find('=',0);
					if(pos!=std::string::npos)
					{
						string param = sub_str.substr(0, pos);
						string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
						if (param=="lx"){
							double val;
							stringstream ss2(value);
							ss2>>val;
							viscosity->set_lx(val);
							cout<<"Note! " + param +"="<< val <<" successfully."<<endl;
						}else if(param=="ly"){
							double val;
							stringstream ss2(value);
							ss2>>val;
							viscosity->set_ly(val);
							cout<<"Note! " + param +"="<< val <<" successfully."<<endl;
						}else if(param=="lz"){
							double val;
							stringstream ss2(value);
							ss2>>val;
							viscosity->set_lz(val);
							cout<<"Note! " + param +"="<< val <<" successfully."<<endl;
						}else if(param=="dt"){
							double val;
							stringstream ss2(value);
							ss2>>val;
							viscosity->set_dt(val);
							cout<<"Note! " + param +"="<< val <<" successfully."<<endl;
						}
					}
					}
					break;
				}
			case 49:
				{
					shared_ptr<Delaunay> delaunay = std::make_shared<Delaunay>();
					object.push_back(delaunay);
					cout<<"49.  Computing the local volume and output result to file localvolume.log" << endl;
					cout<<"     "<<parameter_list[49]<<endl;

					std::string text = command2[i];
					stringstream ss1(text);
					std::string sub_str;
					while(getline(ss1,sub_str,'|')){
						string::size_type pos=sub_str.find('=',0);
						if(pos!=std::string::npos)
						{
							string param = sub_str.substr(0, pos);
							string value = sub_str.substr(pos+1, sub_str.size()-pos-1);
							if (param=="exclude"){
								std::vector<string> val;
								stringstream ss2(value);
								std::string sub_str2;
								std::string s2 = "{";
								while (getline(ss2,sub_str2,'-')){
									val.push_back(sub_str2);
									s2 += sub_str2 + ", ";
								}
								s2 = s2.erase(s2.size()-2, s2.size());
								s2 += "}";
								delaunay->set_exclude(val);
								cout<<"Note! Change exclude="<< s2 <<" successfully."<<endl;
							}
						}
					}
					break;
				}
			default:
				return 0;
			} 
		}
	cout<<"-----------------------------------------"<<endl;
	bool mst_dcd = false;
	bool xml_dcd = false;
	if(filename.size()==1)
		cout<<"Error! Please specify configuration files!"<<endl;
	if (filename.size()==3)
		{
		std::string file1 = filename[1].c_str();
		std::string file2 = filename[2].c_str();
		std::size_t f1mst = file1.find(".mst");
		std::size_t f2mst = file2.find(".mst");
		
		std::size_t f1xml = file1.find(".xml");
		std::size_t f2xml = file2.find(".xml");
		
		std::size_t f1dcd = file1.find(".dcd");
		std::size_t f2dcd = file2.find(".dcd");
		
		if(f1mst!=std::string::npos&&f2dcd!=std::string::npos)
			mst_dcd = true;
		else if(f2mst!=std::string::npos&&f1dcd!=std::string::npos)
			{
			std::string temp= file1;
			file1 = file2;
			file2 = temp;
			mst_dcd = true;	
			}
		else if(f1xml!=std::string::npos&&f2dcd!=std::string::npos)
			xml_dcd = true;
		else if(f2xml!=std::string::npos&&f1dcd!=std::string::npos)
			{
			std::string temp= file1;
			file1 = file2;
			file2 = temp;
			xml_dcd = true;	
			}
			
		if(mst_dcd)
			{
			string mst_open = file1;
			string dcd_open = file2;
			cout<< "Reading " << mst_open.c_str() << "..." << endl;			
			DCDBuilder build;
			MolInfo mol(build);
			build.readDataFromMST(mst_open.c_str());
			mol.initialize();
			
			for(unsigned int j=0;j<object.size(); j++)
				{
				object[j]->add(&build, &mol);
				}
				
			build.initiate(dcd_open.c_str());
			unsigned int nframe = build.getNframes();
			cout<<"Computing... "<<endl;
			for(unsigned int i=0; i<nframe; i++)
				{
				build.updateDataFromDCD();
				mol.updatePosition0();
				for(unsigned int j=0;j<object.size(); j++)
					{
					object[j]->compute();
					}
				}
			build.outPutInfo();
			mol.outPutInfo();
			}
		else if(xml_dcd)
			{
			string xml_open = file1;
			string dcd_open = file2;
			cout<< "INFO : Reading " << xml_open.c_str() << "..." << endl;			
			DCDBuilder build;
			MolInfo mol(build);
			build.readDataFromXML(xml_open.c_str());
			mol.initialize();
			
			for(unsigned int j=0;j<object.size(); j++)
				{
				object[j]->add(&build, &mol);
				}
				
			build.initiate(dcd_open.c_str());
			unsigned int nframe = build.getNframes();
			cout<<"Computing... "<<endl;
			for(unsigned int i=0; i<nframe; i++)
				{
				build.updateDataFromDCD();
				mol.updatePosition0();
				for(unsigned int j=0;j<object.size(); j++)
					{
					object[j]->compute();
					}
				}
			build.outPutInfo();
			mol.outPutInfo();
			}
		}
		
	//--- mst files in default 
	if(!mst_dcd&&!xml_dcd&&filename.size()>1)
		{
		std::string file1 = filename[1].c_str();
		std::size_t f1mst = file1.find(".mst");
		std::size_t f1xml = file1.find(".xml");
		std::size_t f1log = file1.find(".log");
		if (f1mst!=std::string::npos)
			{
			MSTReader build;
			MolInfo mol(build);
			for(unsigned int j=0;j<object.size(); j++)
				{
				object[j]->add(&build, &mol);
				}
			
			for(unsigned int i=1;i<filename.size(); i++)
				{
				std::string mst_open = filename[i];
				while(build.readDataFromMST(mst_open.c_str()))
					{
					cout<< "Reading " << mst_open.c_str() << "..." << endl;
					if (build.ifchangedNp())
						{
						mol.initialize();
						if (i>1)
							cout<<"The number of particles is changed from "<<build.getLastNParticles()<<" to "<< build.getNParticles()<< endl;
						}
					else if(build.ifchangedNb())
						{
						mol.initialize();
						if (i>1)
							cout<<"The number of bonds is changed from "<<build.getLastNBonds()<<" to "<< build.getBond().size()<< endl;
						}
						
					mol.updatePosition0();
					for(unsigned int j=0;j<object.size(); j++)
						{
						cout<<"Computing... "<<endl<<endl;
						object[j]->compute();
						}				
					}
				}
			build.outPutInfo();
			mol.outPutInfo();
			}
		else if (f1xml!=std::string::npos) 
			{
			XMLBuilder build;
			MolInfo mol(build);
			for(unsigned int j=0;j<object.size(); j++)
				{
				object[j]->add(&build, &mol);
				}
			
			for(unsigned int i=1;i<filename.size(); i++)
				{
				std::string xml_open = filename[i];
				build.readDataFromXML(xml_open.c_str());
			
				cout<< "Reading " << xml_open.c_str() << "..." << endl;
				if (build.ifchangedNp())
					{
					mol.initialize();
					if (i>1)
						cout<<"The number of particles is changed from "<<build.getLastNParticles()<<" to "<< build.getNParticles()<< endl;
					}
				else if(build.ifchangedNb())
					{
					mol.initialize();
					if (i>1)
						cout<<"The number of bonds is changed from "<<build.getLastNBonds()<<" to "<< build.getBond().size()<< endl;
					}
					
				mol.updatePosition0();
				for(unsigned int j=0;j<object.size(); j++)
					{
					cout<<"Computing... "<<endl<<endl;
					object[j]->compute();
					}
				}
			build.outPutInfo();
			mol.outPutInfo();
			}
		else if (f1log!=std::string::npos)
			{
			cout<<"Computing... "<<endl;
			for(unsigned int j=0;j<object.size(); j++)
				{
				object[j]->add_datalog(filename[1]);
				object[j]->compute();
				}
			cout<<"Complete!"<<endl;
			}
		}
		
return 0;	
	} 

