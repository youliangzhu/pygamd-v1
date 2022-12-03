

#include "Xmlbuilder.h"



/*! \param fname File name with the data to load
    The file will be read and parsed fully during the constructor call.
*/

BoxDim::BoxDim()
    {
    xlo = xhi = ylo = yhi = zlo = zhi = 0.0;
    }

/*! \param Len Length of one side of the box
    \post Box ranges from \c -Len/2 to \c +Len/2 in all 3 dimensions
 */
BoxDim::BoxDim(float Len)
    {
    // sanity check
    assert(Len > 0);
    
    // assign values
    xlo = ylo = zlo = -Len/float(2.0);
    xhi = zhi = yhi = Len/float(2.0);
    }

/*! \param Len_x Length of the x dimension of the box
    \param Len_y Length of the x dimension of the box
    \param Len_z Length of the x dimension of the box
 */
BoxDim::BoxDim(float Len_x, float Len_y, float Len_z)
    {
    // sanity check
    assert(Len_x > 0 && Len_y > 0 && Len_z > 0);
    
    // assign values
    xlo = -Len_x/float(2.0);
    xhi = Len_x/float(2.0);
    
    ylo = -Len_y/float(2.0);
    yhi = Len_y/float(2.0);
    
    zlo = -Len_z/float(2.0);
    zhi = Len_z/float(2.0);
    }



Xmlbuilder::Xmlbuilder(const std::string &fname)
    {
    // initialize member variables
    m_timestep = 0;
    m_box_read = false;
    m_num_dimensions = 3;
    
    // initialize the parser map
    m_parser_map["box"] = bind(&Xmlbuilder::parseBoxNode, this, _1);
    m_parser_map["position"] = bind(&Xmlbuilder::parsePositionNode, this, _1);
    m_parser_map["image"] = bind(&Xmlbuilder::parseImageNode, this, _1);
    m_parser_map["velocity"] = bind(&Xmlbuilder::parseVelocityNode, this, _1);
    m_parser_map["mass"] = bind(&Xmlbuilder::parseMassNode, this, _1);
    m_parser_map["diameter"] = bind(&Xmlbuilder::parseDiameterNode, this, _1);
    m_parser_map["type"] = bind(&Xmlbuilder::parseTypeNode, this, _1);
    m_parser_map["body"] = bind(&Xmlbuilder::parseBodyNode, this, _1);	
    m_parser_map["bond"] = bind(&Xmlbuilder::parseBondNode, this, _1);
    m_parser_map["angle"] = bind(&Xmlbuilder::parseAngleNode, this, _1);
    m_parser_map["dihedral"] = bind(&Xmlbuilder::parseDihedralNode, this, _1);
    m_parser_map["charge"] = bind(&Xmlbuilder::parseChargeNode, this, _1);
    m_parser_map["molecule"] = bind(&Xmlbuilder::parseMoleIdNode, this, _1);
    m_parser_map["orientation"] = bind(&Xmlbuilder::parseOrientationNode, this, _1); 
    m_parser_map["quaternion"] = bind(&Xmlbuilder::parseQuaternionNode, this, _1);  
    // read in the file
    readFile(fname);
    }


unsigned int Xmlbuilder::getNumDimensions() const
    {
    return (unsigned int)m_num_dimensions;
    }


unsigned int Xmlbuilder::getNumParticles() const
    {
    assert(m_pos_array.size() > 0);
    return (unsigned int)m_pos_array.size();
    }


unsigned int Xmlbuilder::getNumParticleTypes() const
    {
    assert(m_type_mapping.size() >= 0);
    return (unsigned int)m_type_mapping.size();
    }


BoxDim Xmlbuilder::getBox() const
    {
    return m_box;
    }


unsigned int Xmlbuilder::getTimeStep() const
    {
    return m_timestep;
    }


void Xmlbuilder::setTimeStep(unsigned int ts)
    {
    m_timestep = ts;
    }



void Xmlbuilder::readFile(const string &fname)
    {

    XMLNode root_node;
	bool node_existed = false;
    string node_name[3] = {"hoomd_xml","polymer_xml","galamost_xml"};
    cout<< "INFO : Reading " << fname << "..." << endl;
	XMLResults results;	
	for(unsigned int i =0; i<3; i++)
		{

		root_node = XMLNode::parseFile(fname.c_str(),node_name[i].c_str(), &results);
		if(results.error == eXMLErrorNone)
			{
			node_existed=true;
			cout<< "INFO : Parsing "<< node_name[i]<<" node!"<< endl;		
			break;
			}
		}
    // handle errors
    if (!node_existed)
        {
        // create message
        if (results.error==eXMLErrorFirstTagNotFound)
            {
            cerr << endl << "***Error! Root node of " << fname << " can not be parsed!" << endl << endl;
            throw runtime_error("Error reading xml file");
            }
        ostringstream error_message;
        error_message << XMLNode::getError(results.error) << " in file "
        << fname << " at line " << results.nLine << " col "
        << results.nColumn;
        cerr << endl << "***Error! " << error_message.str() << endl << endl;
        throw runtime_error("Error reading xml file");
        }
 
        
    string xml_version;
    if (root_node.isAttributeSet("version"))
        {
        xml_version = root_node.getAttribute("version");
        }
    else
        {
        cout << "Notice: No version specified in polymer_xml root node: assuming 1.0" << endl;
        xml_version = string("1.0");
        }
        
    // right now, the version tag doesn't do anything: just warn if it is not a valid version
    vector<string> valid_versions;
    valid_versions.push_back("1.0");
    valid_versions.push_back("1.1");
    valid_versions.push_back("1.2");
    valid_versions.push_back("1.3");
    valid_versions.push_back("1.4");	
    bool valid = false;
    vector<string>::iterator i;
    for (i = valid_versions.begin(); i != valid_versions.end(); ++i)
        {
        if (xml_version == *i)
            {
            valid = true;
            break;
            }
        }
    if (!valid)
        cout << endl
             << "***Warning! polymer_xml file with version not in the range 1.0-1.2  specified,"
             << " I don't know how to read this. Continuing anyways." << endl << endl;
             
    // the file was parsed successfully by the XML reader. Extract the information now
    // start by checking the number of configurations in the file
    int num_configurations = root_node.nChildNode("configuration");
    if (num_configurations == 0)
        {
        cerr << endl << "***Error! No <configuration> specified in the XML file" << endl << endl;
        throw runtime_error("Error reading xml file");
        }
    if (num_configurations > 1)
        {
        cerr << endl << "***Error! Sorry, the input XML file must have only one configuration" << endl << endl;
        throw runtime_error("Error reading xml file");
        }
        
    // extract the only configuration node
    XMLNode configuration_node = root_node.getChildNode("configuration");

    if (configuration_node.isAttributeSet("time_step"))
        {
        m_timestep = atoi(configuration_node.getAttribute("time_step"));
        }
    
    // extract the number of dimensions, or default to 3
    if (configuration_node.isAttributeSet("dimensions"))
        {
        m_num_dimensions = atoi(configuration_node.getAttribute("dimensions"));
        }
    else
        m_num_dimensions = 3;
        
    // loop through all child nodes of the configuration
    for (int cur_node=0; cur_node < configuration_node.nChildNode(); cur_node++)
        {
        // extract the name and call the appropriate node parser, if it exists
        XMLNode node = configuration_node.getChildNode(cur_node);
        string name = node.getName();
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        
        std::map< std::string, boost::function< void (const XMLNode&) > >::iterator parser;
        parser = m_parser_map.find(name);
        if (parser != m_parser_map.end())
            parser->second(node);
        else
            cout << "Notice: Parser for node <" << name << "> not defined, ignoring" << endl;
        }
        
    // check for required items in the file
    if (!m_box_read)
        {
        cerr << endl
             << "***Error! A <box> node is required to define the dimensions of the simulation box"
             << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
    if (m_pos_array.size() == 0)
        {
        cerr << endl << "***Error! No particles defined in <position> node" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
    if (m_type_array.size()!= 0 && m_type_array.size() == 0)
        {
        cerr << endl << "***Error! No particles defined in <type> node" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
		
    if (m_Molecular_id.size() != 0 && m_Molecular_id.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_Molecular_id.size() << " m_Molecular != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }		
        
    // check for potential user errors
    if (m_vel_array.size() != 0 && m_vel_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_vel_array.size() << " velocities != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
    if (m_mass_array.size() != 0 && m_mass_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_mass_array.size() << " masses != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
    if (m_diameter_array.size() != 0 && m_diameter_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_diameter_array.size() << " diameters != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
    if (m_image_array.size() != 0 && m_image_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_image_array.size() << " images != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
    if ( m_type_array.size() != 0 && m_type_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_type_array.size() << " type values != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
    if (m_body_array.size() != 0 && m_body_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_body_array.size() << " body values != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }		
    if (m_charge_array.size() != 0 && m_charge_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_charge_array.size() << " charge values != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
        
    // notify the user of what we have accomplished
    cout <<"INFO : --- polymer_xml file read summary" << endl;
    cout <<"INFO : "<< getNumParticles() << " positions at timestep " << m_timestep << endl;
    if (m_image_array.size() > 0)
        cout <<"INFO : "<< m_image_array.size() << " images" << endl;
    if (m_vel_array.size() > 0)
        cout <<"INFO : "<< m_vel_array.size() << " velocities" << endl;
    if (m_mass_array.size() > 0)
        cout <<"INFO : "<< m_mass_array.size() << " masses" << endl;
    if (m_diameter_array.size() > 0)
        cout <<"INFO : "<< m_diameter_array.size() << " diameters" << endl;
    cout <<"INFO : "<< getNumParticleTypes() <<  " particle types" << endl;
    if (m_body_array.size() > 0)
        cout <<"INFO : "<< m_body_array.size() << " particle body values" << endl; 	
    if (m_bonds.size() > 0)
        cout <<"INFO : "<< m_bonds.size() << " bonds" << endl;
    if (m_angles.size() > 0)
        cout <<"INFO : "<< m_angles.size() << " angles" << endl;
    if (m_dihedrals.size() > 0)
        cout <<"INFO : "<< m_dihedrals.size() << " dihedrals" << endl;
    if (m_charge_array.size() > 0)
        cout <<"INFO : "<< m_charge_array.size() << " charges" << endl;
    if (m_orientation.size() > 0)
        cout <<"INFO : "<< m_orientation.size() << " orientations" << endl;
    if (m_quaternion.size() > 0)
        cout <<"INFO : "<< m_quaternion.size() << " quaternions" << endl;				
    if (m_Molecular_mapping.size() > 0)
        cout <<"INFO : "<< m_Molecular_mapping.size() << " kind Molecular" << endl;
    }


void Xmlbuilder::parseBoxNode(const XMLNode &node)
    {
    // first, verify that this is the box node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("box"));
    
    // temporary values for extracting attributes as floats
    float Lx,Ly,Lz;
    istringstream temp;
    
    // use string streams to extract Lx, Ly, Lz
    // throw exceptions if these attributes are not set
    if (!node.isAttributeSet("lx"))
        {
        cerr << endl << "***Error! lx not set in <box> node" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
    temp.str(node.getAttribute("lx"));
    temp >> Lx;
    temp.clear();
    
    if (!node.isAttributeSet("ly"))
        {
        cerr << endl << "***Error! ly not set in <box> node" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
    temp.str(node.getAttribute("ly"));
    temp >> Ly;
    temp.clear();
    
    if (!node.isAttributeSet("lz"))
        {
        cerr << endl << "***Error! lz not set in <box> node" << endl << endl;
        throw runtime_error("Error extracting data from polymer_xml file");
        }
    temp.str(node.getAttribute("lz"));
    temp >> Lz;
    temp.clear();
    
    // initialize the BoxDim and set the flag telling that we read the <box> node
    m_box = BoxDim(Lx,Ly,Lz);
    m_box_read = true;
    }


void Xmlbuilder::parsePositionNode(const XMLNode &node)
    {
    // check that this is actually a position node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("position"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");

    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        float x,y,z;
        parser >> x >> y >> z;
        if (parser.good())
            m_pos_array.push_back(vec(x,y,z));
        }
    }


void Xmlbuilder::parseImageNode(const XMLNode& node)
    {
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("image"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        int x,y,z;
        parser >> x >> y >> z;
        if (parser.good())
            m_image_array.push_back(vec_int(x,y,z));
        }
    }


void Xmlbuilder::parseVelocityNode(const XMLNode &node)
    {
    // check that this is actually a velocity node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("velocity"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        float x,y,z;
        parser >> x >> y >> z;
        if (parser.good())
            m_vel_array.push_back(vec(x,y,z));
        }
    }


void Xmlbuilder::parseMassNode(const XMLNode &node)
    {
    // check that this is actually a velocity node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("mass"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        float mass;
        parser >> mass;
        if (parser.good())
            m_mass_array.push_back(mass);
        }
    }


void Xmlbuilder::parseDiameterNode(const XMLNode &node)
    {
    // check that this is actually a velocity node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("diameter"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        float diameter;
        parser >> diameter;
        if (parser.good())
            m_diameter_array.push_back(diameter);
        }
    }


void Xmlbuilder::parseMoleIdNode(const XMLNode &node)
    {
    // check that this is actually a molecule node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("molecule"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        // dynamically determine the particle types
        string molecule;
        parser >> molecule;
        if (parser.good())
            m_Molecular_id.push_back(getMoleculerId(molecule));
        }
    }


void Xmlbuilder::parseTypeNode(const XMLNode &node)
    {
    // check that this is actually a type node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("type"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        // dynamically determine the particle types
        string type;
        parser >> type;
        if (parser.good())
            m_type_array.push_back(getTypeId(type));
        }
    }
	
void Xmlbuilder::parseBodyNode(const XMLNode &node)
    {
    // check that this is actually a type node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("body"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");

    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        // handle -1 as NO_BODY
        int body;
        parser >> body;
        
        if (parser.good())
            {
            if (body == -1)
                m_body_array.push_back(NO_BODY);
            else
                m_body_array.push_back(body);
            }
        }
    }

void Xmlbuilder::parseBondNode(const XMLNode &node)
    {
    // check that this is actually a bond node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("bond"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        string type_name;
        unsigned int a, b;
        parser >> type_name >> a >> b;
        if (parser.good())
            m_bonds.push_back(Bond(getBondTypeId(type_name), a, b));
        }
    }

void Xmlbuilder::parseAngleNode(const XMLNode &node)
    {
    // check that this is actually a angle node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("angle"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        string type_name;
        unsigned int a, b, c;
        parser >> type_name >> a >> b >> c;
        if (parser.good())
            m_angles.push_back(Angle(getAngleTypeId(type_name), a, b, c));
        }
    }


void Xmlbuilder::parseDihedralNode(const XMLNode &node)
    {
    // check that this is actually a dihedral node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("dihedral"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        string type_name;
        unsigned int a, b, c, d;
        parser >> type_name >> a >> b >> c >> d;
        if (parser.good())
            m_dihedrals.push_back(Dihedral(getDihedralTypeId(type_name), a, b, c, d));
        }
    }



void Xmlbuilder::parseChargeNode(const XMLNode &node)
    {
    // check that this is actually a charge node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("charge"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        float charge;
        parser >> charge;
        if (parser.good())
            m_charge_array.push_back(charge);
        }
    }

void Xmlbuilder::parseOrientationNode(const XMLNode &node)
    {
    // check that this is actually a orientation node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("orientation"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        float x,y,z;
        parser >> x >> y >> z;	
        if (parser.good())
            m_orientation.push_back(vec(x,y,z));
        }
    }	

void Xmlbuilder::parseQuaternionNode(const XMLNode &node)
    {
    // check that this is actually a orientation node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("quaternion"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        float x,y,z,w;
        parser >> x >> y >> z >> w;	
        if (parser.good())
            m_quaternion.push_back(vec4(x,y,z,w));
        }
    }	
	
unsigned int Xmlbuilder::getMoleculerId(const std::string& name) 
    {
    // search for the type mapping
    for (unsigned int i = 0; i < m_Molecular_mapping.size(); i++)
        {
        if (m_Molecular_mapping[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_Molecular_mapping.push_back(name);
    return (unsigned int)m_Molecular_mapping.size()-1;
    }


unsigned int Xmlbuilder::getTypeId(const std::string& name) 
    {
    // search for the type mapping
    for (unsigned int i = 0; i < m_type_mapping.size(); i++)
        {
        if (m_type_mapping[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_type_mapping.push_back(name);
    return (unsigned int)m_type_mapping.size()-1;
    }


unsigned int Xmlbuilder::getBondTypeId(const std::string& name)
    {
    // search for the type mapping
    for (unsigned int i = 0; i < m_bond_type_mapping.size(); i++)
        {
        if (m_bond_type_mapping[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_bond_type_mapping.push_back(name);
    return (unsigned int)m_bond_type_mapping.size()-1;
    }


unsigned int Xmlbuilder::getAngleTypeId(const std::string& name)
    {
    // search for the type mapping
    for (unsigned int i = 0; i < m_angle_type_mapping.size(); i++)
        {
        if (m_angle_type_mapping[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_angle_type_mapping.push_back(name);
    return (unsigned int)m_angle_type_mapping.size()-1;
    }


unsigned int Xmlbuilder::getDihedralTypeId(const std::string& name)
    {
    // search for the type mapping
    for (unsigned int i = 0; i < m_dihedral_type_mapping.size(); i++)
        {
        if (m_dihedral_type_mapping[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_dihedral_type_mapping.push_back(name);
    return (unsigned int)m_dihedral_type_mapping.size()-1;
    }





unsigned int Xmlbuilder::getNumBondTypes()  const
    {
    return (unsigned int)m_bond_type_mapping.size();
    }


unsigned int Xmlbuilder::getNumAngleTypes() const 
    {
    return (unsigned int)m_angle_type_mapping.size();
    }

/*! \return Number of dihedral types determined from the XML file
*/
unsigned int Xmlbuilder::getNumDihedralTypes() const
    {
    return (unsigned int)m_dihedral_type_mapping.size();
    }








