#include "XMLBuilder.h"


XMLBuilder::XMLBuilder()
    {
    // initialize member variables
    m_timestep = 0;
    m_box_read = false;
    m_ndimension = 3;
	m_box = BoxSize(0.0, 0.0, 0.0);
	m_num_particles = 0;
	m_num_bonds = 0;
	m_last_np = 0;
	m_last_nb = 0;
	m_if_trajectory = false;
    // initialize the parser map
    m_parser_map["box"] = bind(&XMLBuilder::parseBoxNode, this, std::placeholders::_1);
    m_parser_map["position"] = bind(&XMLBuilder::parsePositionNode, this, std::placeholders::_1);
    m_parser_map["image"] = bind(&XMLBuilder::parseImageNode, this, std::placeholders::_1);
    m_parser_map["velocity"] = bind(&XMLBuilder::parseVelocityNode, this, std::placeholders::_1);
    m_parser_map["mass"] = bind(&XMLBuilder::parseMassNode, this, std::placeholders::_1);
    m_parser_map["diameter"] = bind(&XMLBuilder::parseDiameterNode, this, std::placeholders::_1);
    m_parser_map["type"] = bind(&XMLBuilder::parseTypeNode, this, std::placeholders::_1);
    m_parser_map["body"] = bind(&XMLBuilder::parseBodyNode, this, std::placeholders::_1);
    m_parser_map["bond"] = bind(&XMLBuilder::parseBondNode, this, std::placeholders::_1);
    m_parser_map["angle"] = bind(&XMLBuilder::parseAngleNode, this, std::placeholders::_1);
    m_parser_map["dihedral"] = bind(&XMLBuilder::parseDihedralNode, this, std::placeholders::_1);
    m_parser_map["charge"] = bind(&XMLBuilder::parseChargeNode, this, std::placeholders::_1);
    m_parser_map["inert"] = bind(&XMLBuilder::parseInertNode, this, std::placeholders::_1);
    m_parser_map["rotation"] = bind(&XMLBuilder::parseRotationNode, this, std::placeholders::_1); 
    m_parser_map["rotangle"] = bind(&XMLBuilder::parseRotangleNode, this, std::placeholders::_1); 
    m_parser_map["h_init"] = bind(&XMLBuilder::parseInitNode, this, std::placeholders::_1);
    m_parser_map["h_cris"] = bind(&XMLBuilder::parseCrisNode, this, std::placeholders::_1);
    m_parser_map["molecule"] = bind(&XMLBuilder::parseMoleculeNode, this, std::placeholders::_1);
    m_parser_map["orientation"] = bind(&XMLBuilder::parseOrientationNode, this, std::placeholders::_1);
    m_parser_map["quaternion"] = bind(&XMLBuilder::parseQuaternionNode, this, std::placeholders::_1);
    m_parser_map["constraint"] = bind(&XMLBuilder::parseConstraintNode, this, std::placeholders::_1);
    m_parser_map["vsite"] = bind(&XMLBuilder::parseVsiteNode, this, std::placeholders::_1);
    m_parser_map["Aspheres"] = bind(&XMLBuilder::parseAsphereNode, this, std::placeholders::_1);
    m_parser_map["Patches"] = bind(&XMLBuilder::parsePatchNode, this, std::placeholders::_1);
    m_parser_map["virial"] = bind(&XMLBuilder::parseVirialNode, this, std::placeholders::_1);
    m_parser_map["virial_matrix"] = bind(&XMLBuilder::parseVirialMatrixNode, this, std::placeholders::_1);
    m_object_name = "XMLBuilder";
    // read in the file
    
    }

XMLBuilder::~XMLBuilder()
    {

    }

void XMLBuilder::clearDataAll()
	{
    m_pos.clear();
    m_image.clear();
    m_vel.clear();
    m_mass.clear();
    m_diameter.clear();
    m_type.clear();
    m_body.clear();
    m_charge.clear();
    m_init.clear();
    m_cris.clear();
    m_force.clear();
    m_virial.clear();
    m_virial_matrix.clear();

    m_bonds.clear();
    m_angles.clear();
    m_dihedrals.clear();
    m_vsites.clear();

    m_orientation.clear();
    m_quaternion.clear();
    m_inert.clear();
    m_rotangle.clear();
    m_rotation.clear();
    m_asphere.clear();
    m_patch.clear();
    m_patch_num.clear();
    m_molecule.clear();
    m_type_map.clear();
    m_bond_type_map.clear();
    m_angle_type_map.clear();
    m_dihedral_type_map.clear();
    m_vsite_type_map.clear();
	}

void XMLBuilder::readDataFromXML(const string &fname)
    {
	clearDataAll();
	readFile(fname);
	m_fname=fname;
    }
	
unsigned int XMLBuilder::getLastNParticles() const
    {
    return m_last_np;
    }
	
unsigned int XMLBuilder::getLastNBonds() const
    {
    return m_last_nb;
    }

unsigned int XMLBuilder::getNDimensions() const
    {
    return (unsigned int)m_ndimension;
    }


unsigned int XMLBuilder::getNParticles() const
    {
    assert(m_pos.size() > 0);
    return (unsigned int)m_pos.size();
    }


unsigned int XMLBuilder::getNParticleTypes() const
    {
    assert(m_type_map.size() >= 0);
    return (unsigned int)m_type_map.size();
    }


BoxSize XMLBuilder::getBox() const
    {
    return m_box;
    }


unsigned int XMLBuilder::getTimeStep() const
    {
    return m_timestep;
    }


void XMLBuilder::setTimeStep(unsigned int ts)
    {
    m_timestep = ts;
    }

void XMLBuilder::readFile(const string &fname)
    {
	m_if_changed_np = false;
	m_if_changed_nb = false;
    XMLNode root_node;
	bool node_existed = false;
    string node_name[3] = {"hoomd_xml","polymer_xml","galamost_xml"};
	XMLResults results;	
	for(unsigned int i =0; i<3; i++)
		{

		root_node = XMLNode::parseFile(fname.c_str(),node_name[i].c_str(), &results);
		if(results.error == eXMLErrorNone)
			{
			node_existed=true;
			m_node_name = node_name[i];		
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
        cout << "Notice: No version specified in galamost_xml root node: assuming 1.0" << endl;
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
             << "***Warning! galamost_xml file with version not in the range 1.0-1.2  specified,"
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
        m_ndimension = atoi(configuration_node.getAttribute("dimensions"));
        }
    else
        m_ndimension = 3;
        
    // loop through all child nodes of the configuration
    for (int cur_node=0; cur_node < configuration_node.nChildNode(); cur_node++)
        {
        // extract the name and call the appropriate node parser, if it exists
        XMLNode node = configuration_node.getChildNode(cur_node);
        string name = node.getName();
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        
        std::map< std::string, std::function< void (const XMLNode&) > >::iterator parser;
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
        throw runtime_error("Error extracting data from galamost_xml file");
        }
		
	double Lz = double (m_box.lz);	
	if(m_ndimension==2&&Lz>0.0)
		{
        cerr << "***Error! two dimensions of the simulation box should be with Lz = 0.0, the Lz = " << Lz <<" in xml files"<< endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
		}
	
	if(m_ndimension==3&&Lz<1.0e-6)
		{
        cerr << "***Error! Lz = 0.0 should be with two dimensions of the simulation, three dimensions defind in xml files "<< endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
		}
	
    if (m_pos.size() == 0)
        {
        cerr << endl << "***Error! No particles defined in <position> node" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_type.size() == 0)
        {
        cerr << endl << "***Error! No particles defined in <type> node" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
		
    if (m_molecule.size() != 0 && m_molecule.size() != m_pos.size())
        {
        cerr << endl << "***Error! " << m_molecule.size() << " molecule != " << m_pos.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_vel.size() != 0 && m_vel.size() != m_pos.size())
        {
        cerr << endl << "***Error! " << m_vel.size() << " velocities != " << m_pos.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_mass.size() != 0 && m_mass.size() != m_pos.size())
        {
        cerr << endl << "***Error! " << m_mass.size() << " masses != " << m_pos.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_diameter.size() != 0 && m_diameter.size() != m_pos.size())
        {
        cerr << endl << "***Error! " << m_diameter.size() << " diameters != " << m_pos.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_image.size() != 0 && m_image.size() != m_pos.size())
        {
        cerr << endl << "***Error! " << m_image.size() << " images != " << m_pos.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if ( m_type.size() != 0 && m_type.size() != m_pos.size())
        {
        cerr << endl << "***Error! " << m_type.size() << " type values != " << m_pos.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_body.size() != 0 && m_body.size() != m_pos.size())
        {
        cerr << endl << "***Error! " << m_body.size() << " body values != " << m_pos.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }		
    if (m_charge.size() != 0 && m_charge.size() != m_pos.size())
        {
        cerr << endl << "***Error! " << m_charge.size() << " charge values != " << m_pos.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
		
	//-check bonds, angles and dihedrals
	unsigned int np = m_pos.size();
	if (np!=m_num_particles)
		{
		m_if_changed_np = true;
		m_last_np = m_num_particles;
		}
	m_num_particles = np;
	
	unsigned int nb = m_bonds.size();
	if (nb!=m_num_bonds)
		{
		m_if_changed_nb = true;
		m_last_nb = m_num_bonds;
		}
	m_num_bonds = nb;

	for (unsigned int i=0; i<m_bonds.size();i++)
		{
		Bond bi = m_bonds[i];
		if(bi.a>=np||bi.b>=np)
			{
			cerr << endl << "***Error! bond '" << bi.type <<" "<<bi.a<<" "<<bi.b<< "' with the particle index not in the range [0, N-1], with N = "<< np<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from galamost_xml file");				
			}
		}
		
	for (unsigned int i=0; i<m_angles.size();i++)
		{
		Angle ai = m_angles[i];
		if(ai.a>=np||ai.b>=np||ai.c>=np)
			{
			cerr << endl << "***Error! angle '" << ai.type <<" "<<ai.a<<" "<<ai.b<<" "<<ai.c<< "' with the particle index not in the range [0, N-1], with N = "<< np<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from galamost_xml file");				
			}
		}	
		
	for (unsigned int i=0; i<m_dihedrals.size();i++)
		{
		Dihedral di = m_dihedrals[i];
		if(di.a>=np||di.b>=np||di.c>=np||di.d>=np)
			{
			cerr << endl << "***Error! dihedral '" << di.type <<" "<<di.a<<" "<<di.b<<" "<<di.c<<" "<<di.d<< "' with the particle index not in the range [0, N-1], with N = "<< np<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from galamost_xml file");				
			}
		}

	for (unsigned int i=0; i<m_vsites.size();i++)
		{
		Dihedral di = m_vsites[i];
		if(di.a>=np||di.b>=np||di.c>=np||di.d>=np)
			{
			cerr << endl << "***Error! vsite '" << di.type <<" "<<di.a<<" "<<di.b<<" "<<di.c<<" "<<di.d<< "' with the particle index not in the range [0, N-1], with N = "<< np<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from galamost_xml file");				
			}
		}			
    }
	
void XMLBuilder::outPutInfo()
	{
    // notify the user of what we have accomplished
    cout <<"----------------------------------------- "<< endl;	
    cout <<"INFO : --- XML file read summary" << endl;
	cout<< "INFO : Parsing "<< m_node_name <<" node!"<< endl;	
    cout <<"INFO : "<< getNParticles() << " positions at timestep " << m_timestep << endl;
    if (m_image.size() > 0){
        cout <<"INFO : "<< m_image.size() << " images" << endl;}
    if (m_vel.size() > 0){
        cout <<"INFO : "<< m_vel.size() << " velocities" << endl;}
    if (m_mass.size() > 0){
        cout <<"INFO : "<< m_mass.size() << " masses" << endl;}
    if (m_diameter.size() > 0){
        cout <<"INFO : "<< m_diameter.size() << " diameters" << endl;}
    cout <<"INFO : "<< getNParticleTypes() <<  " particle types" << endl;
    if (m_body.size() > 0){
        cout <<"INFO : "<< m_body.size() << " particle body values" << endl;}
    if (m_bonds.size() > 0){
        cout <<"INFO : "<< m_bonds.size() << " bonds" << endl;}
    if (m_angles.size() > 0){
        cout <<"INFO : "<< m_angles.size() << " angles" << endl;}
    if (m_dihedrals.size() > 0){
        cout <<"INFO : "<< m_dihedrals.size() << " dihedrals" << endl;}
    if (m_charge.size() > 0){
        cout <<"INFO : "<< m_charge.size() << " charges" << endl;}
    if (m_orientation.size() > 0){
        cout <<"INFO : "<< m_orientation.size() << " orientations" << endl;}
    if (m_quaternion.size() > 0){
        cout <<"INFO : "<< m_quaternion.size() << " quaternions" << endl;}
    if (m_molecule.size() > 0){
		cout <<"INFO : "<< m_molecule.size() << " molecules" << endl;}
		
    if (m_rotangle.size() > 0)
        cout <<" "<< m_rotangle.size() << " rotangles" << endl;
    if (m_force.size() > 0)
        cout <<" "<< m_force.size() << " forces" << endl;
    if (m_virial.size() > 0)
        cout <<" "<< m_virial.size() << " virials" << endl;
    if (m_virial_matrix.size() > 0)
        cout <<"INFO : "<< m_virial_matrix.size() << " virial_matrix" << endl;
    if (m_rotation.size() > 0)
        cout <<" "<< m_rotation.size() << " rotations" << endl;
	}
	

void XMLBuilder::parseBoxNode(const XMLNode &node)
    {
    // first, verify that this is the box node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("box"));
    
    // temporary values for extracting attributes as doubles
    double Lx,Ly,Lz;
    istringstream temp;
    
    // use string streams to extract Lx, Ly, Lz
    // throw exceptions if these attributes are not set
    if (!node.isAttributeSet("lx"))
        {
        cerr << endl << "***Error! lx not set in <box> node" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    temp.str(node.getAttribute("lx"));
    temp >> Lx;
    temp.clear();
    
    if (!node.isAttributeSet("ly"))
        {
        cerr << endl << "***Error! ly not set in <box> node" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    temp.str(node.getAttribute("ly"));
    temp >> Ly;
    temp.clear();
    
    if (!node.isAttributeSet("lz"))
        {
        cerr << endl << "***Error! lz not set in <box> node" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    temp.str(node.getAttribute("lz"));
    temp >> Lz;
    temp.clear();
    
    // initialize the BoxSize and set the flag telling that we read the <box> node
    m_box = BoxSize(Lx,Ly,Lz);
    m_box_read = true;
    }


void XMLBuilder::parsePositionNode(const XMLNode &node)
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
        double x,y,z;
        parser >> x >> y >> z;
        if (parser.good())
            m_pos.push_back(vec(x,y,z));
        }
    }


void XMLBuilder::parseImageNode(const XMLNode& node)
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
            m_image.push_back(vec_int(x,y,z));
        }
    }


void XMLBuilder::parseVelocityNode(const XMLNode &node)
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
        double x,y,z;
        parser >> x >> y >> z;
        if (parser.good())
            m_vel.push_back(vec(x,y,z));
        }
    }
	

void XMLBuilder::parseMassNode(const XMLNode &node)
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
        double mass;
        parser >> mass;
        if (parser.good())
            m_mass.push_back(mass);
        }
    }


void XMLBuilder::parseDiameterNode(const XMLNode &node)
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
        double diameter;
        parser >> diameter;
        if (parser.good())
            m_diameter.push_back(diameter);
        }
    }


void XMLBuilder::parseMoleculeNode(const XMLNode &node)
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
        int mol;
        parser >> mol;
        if (parser.good())
			{
			if (mol == -1)
                m_molecule.push_back(NO_INDEX);
            else
                m_molecule.push_back(mol);
			}
        }
    }


void XMLBuilder::parseTypeNode(const XMLNode &node)
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
            m_type.push_back(getTypeId(type));
        }
    }
	
void XMLBuilder::parseBodyNode(const XMLNode &node)
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
                m_body.push_back(NO_BODY);
            else
                m_body.push_back(body);
            }
        }
    }

void XMLBuilder::parseBondNode(const XMLNode &node)
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
			{
			unsigned int bond_id = getBondTypeId(type_name);
			m_bonds.push_back(Bond(type_name, a, b, bond_id));	
			}  
        }
    }

void XMLBuilder::parseConstraintNode(const XMLNode &node)
    {
    // check that this is actually a constraint node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("constraint"));
    
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
			{
			unsigned int bond_id = getBondTypeId(type_name);
			m_bonds.push_back(Bond(type_name, a, b, bond_id, "c"));	
			}
        }
    }		
	
void XMLBuilder::parseAngleNode(const XMLNode &node)
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
			{
			unsigned int angle_id = getAngleTypeId(type_name);
			m_angles.push_back(Angle(type_name, a, b, c, angle_id));
			}
        }
    }

void XMLBuilder::parseDihedralNode(const XMLNode &node)
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
			{
			unsigned int dihedral_id = getDihedralTypeId(type_name);
			m_dihedrals.push_back(Dihedral(type_name, a, b, c, d, dihedral_id));
			}
        }
    }
	
void XMLBuilder::parseVsiteNode(const XMLNode &node)
    {
    // check that this is actually a dihedral node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("vsite"));
    
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
			{
			unsigned int vsite_id = getVsiteTypeId(type_name);
			m_vsites.push_back(Dihedral(type_name, a, b, c, d, vsite_id));
			}
        }
    }	

void XMLBuilder::parseChargeNode(const XMLNode &node)
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
        double charge;
        parser >> charge;
        if (parser.good())
            m_charge.push_back(charge);
        }
    }

void XMLBuilder::parseOrientationNode(const XMLNode &node)
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
        double x,y,z;
        parser >> x >> y >> z;	
        if (parser.good())
            m_orientation.push_back(vec(x,y,z));
        }
    }
	
	
void XMLBuilder::parseInertNode(const XMLNode &node)
    {
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("inert"));
    
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        double x,y,z;
        parser >> x >> y >> z;
			
        if (parser.good())
            m_inert.push_back(vec(x,y,z));
        }
    }

void XMLBuilder::parseRotangleNode(const XMLNode &node)
    {
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("rotangle"));
    
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        double x,y,z;
        parser >> x >> y >> z;
			
        if (parser.good())
            m_rotangle.push_back(vec(x,y,z));
        }
    }
	
void XMLBuilder::parseRotationNode(const XMLNode &node)
    {
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("rotation"));
    
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        double x,y,z;
        parser >> x >> y >> z;
			
        if (parser.good())
            m_rotation.push_back(vec(x,y,z));
        }
    }

void XMLBuilder::parseInitNode(const XMLNode &node)
    {
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("h_init"));
    
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        unsigned int x;
        parser >> x;
        if (parser.good())
            m_init.push_back(x);
        }
    }

void XMLBuilder::parseCrisNode(const XMLNode &node)
    {
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("h_cris"));
    
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        unsigned int x;
        parser >> x;
        if (parser.good())
            m_cris.push_back(x);
        }
    }
	
void XMLBuilder::parseQuaternionNode(const XMLNode &node)
    {
    // check that this is actually a quaternion node
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
        double x,y,z,w;
        parser >> x >> y >> z >> w;	
        if (parser.good())
            m_quaternion.push_back(vec4(x,y,z,w));
        }
    }
	
	
void XMLBuilder::parseAsphereNode(const XMLNode &node)
    {
    // check that this is actually a Aspheres node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("Aspheres"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");
    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
		string na;
        double x,y,z,w,m,n;
        parser >> na >> x >> y >> z >> w >> m >> n;	
        if (parser.good())
            m_asphere.push_back(str_vec6(na, x, y, z, w, m, n));
        }
    }
	
void XMLBuilder::parsePatchNode(const XMLNode &node)
    {
    // check that this is actually a Patches node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("Patches"));
    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");

    istringstream parser;
    parser.str(all_text);
	unsigned int count=0;
    while (parser.good())
        {
		string na;
		unsigned int num;
		parser >> na >> num;
		if (parser.good())
			{
			m_patch_num.push_back(str_vec6(na, double(num), double(count), double(count+num), 0.0, 0.0, 0.0));
			count += num;				
			for (unsigned int i =0; i<num;i++)
				{
				string pa;
				double x,y,z,w;
				parser >> pa >> x >> y >> z >> w;
				if (parser.good())
					m_patch.push_back(str_vec6(pa, x, y, z, w, 0.0, 0.0));			
				}			
			}
        }
    }	

void XMLBuilder::parseVirialNode(const XMLNode &node)
    {
        // check that this is actually a position node
        string name = node.getName();
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        assert(name == string("virial"));

        // extract the data from the node
        string all_text;
        for (int i = 0; i < node.nText(); i++)
            all_text += string(node.getText(i)) + string("\n");

        istringstream parser;
        parser.str(all_text);
        while (parser.good())
        {
            double x;
            parser >> x;
            if (parser.good())
                m_virial.push_back(x);
        }
    }

void XMLBuilder::parseVirialMatrixNode(const XMLNode &node)
{
    // check that this is actually a position node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("virial_matrix"));

    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++)
        all_text += string(node.getText(i)) + string("\n");

    istringstream parser;
    parser.str(all_text);
    while (parser.good())
    {
        double x,y,z,a,b,c;
        parser >> x >> y >> z >> a >> b >> c;
        if (parser.good())
            m_virial_matrix.push_back({x,y,z,a,b,c});
    }
}

unsigned int XMLBuilder::getTypeId(const std::string& name) 
    {
    // search for the type map
    for (unsigned int i = 0; i < m_type_map.size(); i++)
        {
        if (m_type_map[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_type_map.push_back(name);
    return (unsigned int)m_type_map.size()-1;
    }


unsigned int XMLBuilder::getBondTypeId(const std::string& name)
    {
    // search for the type map
    for (unsigned int i = 0; i < m_bond_type_map.size(); i++)
        {
        if (m_bond_type_map[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_bond_type_map.push_back(name);
    return (unsigned int)m_bond_type_map.size()-1;
    }


unsigned int XMLBuilder::getAngleTypeId(const std::string& name)
    {
    // search for the type map
    for (unsigned int i = 0; i < m_angle_type_map.size(); i++)
        {
        if (m_angle_type_map[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_angle_type_map.push_back(name);
    return (unsigned int)m_angle_type_map.size()-1;
    }


unsigned int XMLBuilder::getDihedralTypeId(const std::string& name)
    {
    // search for the type map
    for (unsigned int i = 0; i < m_dihedral_type_map.size(); i++)
        {
        if (m_dihedral_type_map[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_dihedral_type_map.push_back(name);
    return (unsigned int)m_dihedral_type_map.size()-1;
    }

unsigned int XMLBuilder::getVsiteTypeId(const std::string& name)
    {
    // search for the type map
    for (unsigned int i = 0; i < m_vsite_type_map.size(); i++)
        {
        if (m_vsite_type_map[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_vsite_type_map.push_back(name);
    return (unsigned int)m_vsite_type_map.size()-1;
    }

unsigned int XMLBuilder::getNBondTypes()  const
    {
    return (unsigned int)m_bond_type_map.size();
    }


unsigned int XMLBuilder::getNAngleTypes() const 
    {
    return (unsigned int)m_angle_type_map.size();
    }
	

/*! \return Number of dihedral types determined from the XML file
*/
unsigned int XMLBuilder::getNDihedralTypes() const
    {
    return (unsigned int)m_dihedral_type_map.size();
    }

/*! \return Number of dihedral types determined from the XML file
*/
unsigned int XMLBuilder::getNVsiteTypes() const
    {
    return (unsigned int)m_vsite_type_map.size();
    }






