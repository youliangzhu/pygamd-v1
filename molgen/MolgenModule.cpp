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

#include "Molgen.h"


std::string PYGAMD_VERSION = "1";
void output_version_info()
    {
    cout << "MOLGEN -- Molecule Generator Plug-Ins"<<endl;
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

PYBIND11_MODULE(molgen, m)
    {

//	register_info();
    output_version_info();

	export_Molecule(m);
	export_DNAchain(m);
	export_Protein(m);
	export_Object(m);
	export_Generators(m);
    }

