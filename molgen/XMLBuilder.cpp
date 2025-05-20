#include "XMLBuilder.h"


XMLBuilder::XMLBuilder(const std::string &fname) : m_fname(fname)
    {
    // initialize member variables
    m_timestep = 0;
    m_box_read = false;
    m_num_dimensions = 3;
    
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
    m_parser_map["h_init"] = bind(&XMLBuilder::parseInitNode, this, std::placeholders::_1);
    m_parser_map["h_cris"] = bind(&XMLBuilder::parseCrisNode, this, std::placeholders::_1);		
    m_parser_map["molecule"] = bind(&XMLBuilder::parseMoleculeNode, this, std::placeholders::_1);
    m_parser_map["orientation"] = bind(&XMLBuilder::parseOrientationNode, this, std::placeholders::_1);
    m_parser_map["quaternion"] = bind(&XMLBuilder::parseQuaternionNode, this, std::placeholders::_1);	
    m_parser_map["constraint"] = bind(&XMLBuilder::parseConstraintNode, this, std::placeholders::_1);
    m_parser_map["vsite"] = bind(&XMLBuilder::parseVsiteNode, this, std::placeholders::_1);
    m_parser_map["aspheres"] = bind(&XMLBuilder::parseAsphereNode, this, std::placeholders::_1);
    m_parser_map["patches"] = bind(&XMLBuilder::parsePatchNode, this, std::placeholders::_1);		
    m_object_name = "XMLBuilder";
    // read in the file
    readFile(fname);
    }

XMLBuilder::~XMLBuilder()
    {

    }
	
unsigned int XMLBuilder::getNDimensions() const
    {
    return (unsigned int)m_num_dimensions;
    }


unsigned int XMLBuilder::getNParticles() const
    {
    assert(m_pos_array.size() > 0);
    return (unsigned int)m_pos_array.size();
    }


unsigned int XMLBuilder::getNParticleTypes() const
    {
    assert(m_type_mapping.size() >= 0);
    return (unsigned int)m_type_mapping.size();
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
	if(m_num_dimensions==2&&Lz>0.0)
		{
        cerr << "***Error! two dimensions of the simulation box should be with Lz = 0.0, the Lz = " << Lz <<" in xml files"<< endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
		}
	
	if(m_num_dimensions==3&&Lz<1.0e-6)
		{
        cerr << "***Error! Lz = 0.0 should be with two dimensions of the simulation, three dimensions defind in xml files "<< endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
		}
	
    if (m_pos_array.size() == 0)
        {
        cerr << endl << "***Error! No particles defined in <position> node" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_type_array.size() == 0)
        {
        cerr << endl << "***Error! No particles defined in <type> node" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
		
    if (m_molecule_array.size() != 0 && m_molecule_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_molecule_array.size() << " molecule != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_vel_array.size() != 0 && m_vel_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_vel_array.size() << " velocities != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_mass_array.size() != 0 && m_mass_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_mass_array.size() << " masses != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_diameter_array.size() != 0 && m_diameter_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_diameter_array.size() << " diameters != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_image_array.size() != 0 && m_image_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_image_array.size() << " images != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if ( m_type_array.size() != 0 && m_type_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_type_array.size() << " type values != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
    if (m_body_array.size() != 0 && m_body_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_body_array.size() << " body values != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }		
    if (m_charge_array.size() != 0 && m_charge_array.size() != m_pos_array.size())
        {
        cerr << endl << "***Error! " << m_charge_array.size() << " charge values != " << m_pos_array.size()
             << " positions" << endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
        }
		
	//-check bonds, angles and dihedrals
	unsigned int np = m_pos_array.size();
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
    cout <<"INFO : --- galamost_xml file read summary" << endl;
	cout<< "INFO : Parsing "<< m_node_name <<" node!"<< endl;	
    cout <<"INFO : "<< getNParticles() << " positions at timestep " << m_timestep << endl;
    if (m_image_array.size() > 0){
        cout <<"INFO : "<< m_image_array.size() << " images" << endl;}
    if (m_vel_array.size() > 0){
        cout <<"INFO : "<< m_vel_array.size() << " velocities" << endl;}
    if (m_mass_array.size() > 0){
        cout <<"INFO : "<< m_mass_array.size() << " masses" << endl;}
    if (m_diameter_array.size() > 0){
        cout <<"INFO : "<< m_diameter_array.size() << " diameters" << endl;}
    cout <<"INFO : "<< getNParticleTypes() <<  " particle types" << endl;
    if (m_body_array.size() > 0){
        cout <<"INFO : "<< m_body_array.size() << " particle body values" << endl;}
    if (m_bonds.size() > 0){
        cout <<"INFO : "<< m_bonds.size() << " bonds" << endl;}
    if (m_angles.size() > 0){
        cout <<"INFO : "<< m_angles.size() << " angles" << endl;}
    if (m_dihedrals.size() > 0){
        cout <<"INFO : "<< m_dihedrals.size() << " dihedrals" << endl;}
    if (m_charge_array.size() > 0){
        cout <<"INFO : "<< m_charge_array.size() << " charges" << endl;}
    if (m_orientation.size() > 0){
        cout <<"INFO : "<< m_orientation.size() << " orientations" << endl;}
    if (m_quaternion.size() > 0){
        cout <<"INFO : "<< m_quaternion.size() << " quaternions" << endl;}
    if (m_molecule_array.size() > 0){
		cout <<"INFO : "<< m_molecule_array.size() << " molecules" << endl;}
	if (m_mass_array.size()==0)
		{
		m_mass_array.resize(m_pos_array.size());
		for(unsigned int i=0; i<m_mass_array.size();i++)
			m_mass_array[i] = 1.0;
        cout <<"INFO : "<<" set mass to be 1.0 by default!" << endl;			
		}		
	}

bool XMLBuilder::updatePositionFromXML(const string &fname)
    {
	m_fname = fname;
	m_box_read = false;
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
	
	
    if (configuration_node.isAttributeSet("natoms"))
        {
        unsigned int natoms = atoi(configuration_node.getAttribute("natoms"));
		if(natoms!=m_pos_array.size())
			return false;
        }	

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
 
        if(name=="position")
            updatePositionNode(node);
        if(name=="image"&&m_image_array.size()==m_pos_array.size())
            updateImageNode(node);		
        if(name=="box")
            parseBoxNode(node);
		if(name=="velocity"&&m_vel_array.size()==m_pos_array.size())
            updateVelocityNode(node);	
		if(name=="orientation"&&m_orientation.size()==m_pos_array.size())
            updateOrientationNode(node);
		if(name=="quaternion"&&m_quaternion.size()==m_pos_array.size())
            updateQuaternionNode(node);		
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
	if(m_num_dimensions==2&&Lz>0.0)
		{
        cerr << "***Error! two dimensions of the simulation box should be with Lz = 0.0, the Lz = " << Lz <<" in xml files"<< endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
		}
	
	if(m_num_dimensions==3&&Lz<1.0e-6)
		{
        cerr << "***Error! Lz = 0.0 should be with two dimensions of the simulation, three dimensions defind in xml files "<< endl << endl;
        throw runtime_error("Error extracting data from galamost_xml file");
		}
	return true;
    }
	
void XMLBuilder::updatePositionNode(const XMLNode &node)
    {
    // check that this is actually a position node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("position"));
    unsigned int num =0;
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++){
		all_text += string(node.getText(i)) + string("\n");}
	unsigned int np = m_pos_array.size();
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {			
        double x,y,z;
        parser >> x >> y >> z;
        if (parser.good())
			{
			if (num<np)
				m_pos_array[num] = vec(x,y,z);
			num += 1;
			}
        }
	if(num!=np)
		{
		cerr << "***Error! The number of particles for updating position "<< num <<" is larger than the one of the stored frame "<< np <<" !"<< endl << endl;
		throw runtime_error("Error XMLBuilder::updatePositionNode");			
		}			
    }
	
void XMLBuilder::updateImageNode(const XMLNode& node)
    {
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("image"));
    unsigned int num =0;    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++){
		all_text += string(node.getText(i)) + string("\n");}
	unsigned int np = m_image_array.size();    
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        int x,y,z;
        parser >> x >> y >> z;
        if (parser.good())
			{
			if (num<np)
				m_image_array[num] = vec_int(x,y,z);
			num += 1;
			}
        }
	if(num!=np)
		{
		cerr << "***Error! The number of particles for updating image "<< num <<" is larger than the one of the stored frame "<< np <<" !"<< endl << endl;
		throw runtime_error("Error XMLBuilder::updateImageNode");			
		}			
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
            m_pos_array.push_back(vec(x,y,z));
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
            m_image_array.push_back(vec_int(x,y,z));
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
            m_vel_array.push_back(vec(x,y,z));
        }
    }
	
void XMLBuilder::updateVelocityNode(const XMLNode &node)
    {
    // check that this is actually a velocity node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("velocity"));
    unsigned int num =0;    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++){
		all_text += string(node.getText(i)) + string("\n");}
 	unsigned int np = m_pos_array.size();   
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        double x,y,z;
        parser >> x >> y >> z;
        if (parser.good())
			{
			if (num<np)
				m_vel_array[num]=vec(x,y,z);
			num += 1;		
			}
        }
	if(num!=np)
		{
		cerr << "***Error! The number of particles for updating velocity "<< num <<" is larger than the one of the stored frame "<< np <<" !"<< endl << endl;
		throw runtime_error("Error XMLBuilder::updateVelocityNode");			
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
            m_mass_array.push_back(mass);
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
            m_diameter_array.push_back(diameter);
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
                m_molecule_array.push_back(NO_INDEX);
            else
                m_molecule_array.push_back(mol);
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
            m_type_array.push_back(getTypeId(type));
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
                m_body_array.push_back(NO_BODY);
            else
                m_body_array.push_back(body);
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
            m_charge_array.push_back(charge);
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
	
void XMLBuilder::updateOrientationNode(const XMLNode &node)
    {
    // check that this is actually a orientation node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("orientation"));
    unsigned int num =0;    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++){
		all_text += string(node.getText(i)) + string("\n");}
 	unsigned int np = m_pos_array.size();   
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        double x,y,z;
        parser >> x >> y >> z;
        if (parser.good())
			{
			if (num<np)
				m_orientation[num]=vec(x,y,z);
			num += 1;		
			}
        }
	if(num!=np)
		{
		cerr << "***Error! The number of particles for updating orientation "<< num <<" is larger than the one of the stored frame "<< np <<" !"<< endl << endl;
		throw runtime_error("Error XMLBuilder::updateOrientationNode");			
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
	
void XMLBuilder::updateQuaternionNode(const XMLNode &node)
    {
    // check that this is actually a quaternion node
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    assert(name == string("quaternion"));
    unsigned int num =0;    
    // extract the data from the node
    string all_text;
    for (int i = 0; i < node.nText(); i++){
		all_text += string(node.getText(i)) + string("\n");}
 	unsigned int np = m_pos_array.size();   
    istringstream parser;
    parser.str(all_text);
    while (parser.good())
        {
        double x,y,z,w;
        parser >> x >> y >> z >> w;
        if (parser.good())
			{
			if (num<np)
				m_quaternion[num]=vec4(x,y,z,w);
			num += 1;		
			}
        }
	if(num!=np)
		{
		cerr << "***Error! The number of particles for updating quaternion "<< num <<" is larger than the one of the stored frame "<< np <<" !"<< endl << endl;
		throw runtime_error("Error XMLBuilder::updateQuaternionNode");			
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

unsigned int XMLBuilder::getTypeId(const std::string& name) 
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


unsigned int XMLBuilder::getBondTypeId(const std::string& name)
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


unsigned int XMLBuilder::getAngleTypeId(const std::string& name)
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


unsigned int XMLBuilder::getDihedralTypeId(const std::string& name)
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

unsigned int XMLBuilder::getVsiteTypeId(const std::string& name)
    {
    // search for the type mapping
    for (unsigned int i = 0; i < m_vsite_type_mapping.size(); i++)
        {
        if (m_vsite_type_mapping[i] == name)
            return i;
        }
    // add a new one if it is not found
    m_vsite_type_mapping.push_back(name);
    return (unsigned int)m_vsite_type_mapping.size()-1;
    }

unsigned int XMLBuilder::getNBondTypes()  const
    {
    return (unsigned int)m_bond_type_mapping.size();
    }


unsigned int XMLBuilder::getNAngleTypes() const 
    {
    return (unsigned int)m_angle_type_mapping.size();
    }
	

/*! \return Number of dihedral types determined from the XML file
*/
unsigned int XMLBuilder::getNDihedralTypes() const
    {
    return (unsigned int)m_dihedral_type_mapping.size();
    }

/*! \return Number of dihedral types determined from the XML file
*/
unsigned int XMLBuilder::getNVsiteTypes() const
    {
    return (unsigned int)m_vsite_type_mapping.size();
    }






