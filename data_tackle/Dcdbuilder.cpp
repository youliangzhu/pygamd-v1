#include "Dcdbuilder.h"
using namespace std;
#include "time.h"


static unsigned int read_int(fstream &file)
    {
    unsigned int val;
    file.read((char *)&val, sizeof(unsigned int));
    return val;
    }

Dcdbuilder::Dcdbuilder() : mst_reader()
	{
	m_object_name = "Dcdbuilder";
	}
	
void Dcdbuilder::initiate(const string &fname)
	{
    cout<< "INFO : Reading " << fname << "..." << endl;		
    file.open(fname.c_str(), ios::in | ios::binary);
    if (!file.good())
        {
        cerr << endl << "***Error! Error openning DCD" << endl << endl;
        throw runtime_error("Error Dcdbuilder");
        }

	m_frames_count = 0;
	m_nparticles = 0;
	m_num_frames_written = 0;
    // initialize member variables
    m_last_written_step = 0; 
	m_error = false;
	read_file_header();		
	}
	
	
void Dcdbuilder::read_file_header()
    {

    unsigned int val = read_int(file);
    if(val != 84)
		m_error = true;

    // the next 4 bytes in the file must be "CORD"
    char cord_data[5];
    file.read(cord_data, 4);
    m_num_frames_written = read_int(file);     // Number of frames in file, none written yet
    m_start_timestep= read_int(file);;    // Starting timestep
    m_period = read_int(file);;     // Timesteps between frames written to the file
    m_last_written_step = read_int(file);      // Number of timesteps in simulation
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);         // timestep (unused)
    val = read_int(file);        // include unit cell
    if(val != 1)
		m_error = true;
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);
    val = read_int(file);   // Pretend to be CHARMM version 24
    if(val != 24)
		m_error = true;
    val = read_int(file); 
    if(val != 84)
		m_error = true;
    val = read_int(file); 
    if(val != 164)
		m_error = true;
    val = read_int(file); 
    if(val != 2)
		m_error = true;
		
    char title_string[81];
    file.read(title_string, 80);
    
    char time_str[81];
    file.read(time_str, 80);
    
    val = read_int(file); 
    if(val != 164)
		m_error = true;
    val = read_int(file); 
    if(val != 4)
		m_error = true;
    m_nparticles = read_int(file); 
    val = read_int(file); 
    if(val != 4)
		m_error = true;

    if (m_error)
        {
        cerr << endl << "***Error! Error reading wrong data DCD header" << endl << endl;
        throw runtime_error("Error reading DCD file");
        }
    
    // check for errors
    if (!file.good())
        {
        cerr << endl << "***Error! Error reading DCD header" << endl << endl;
        throw runtime_error("Error reading DCD file");
        }
    }


void Dcdbuilder::read_frame_header()
    {
    double unitcell[6];
    // set box dimensions
    unitcell[0] = 0.0;
    unitcell[2] = 0.0;
    unitcell[5] = 0.0;
    // box angles are 90 degrees
    unitcell[1] = 0.0;
    unitcell[3] = 0.0;
    unitcell[4] = 0.0;
    
    unsigned int val = read_int(file); 
    if(val != 48)
			m_error = true;	
    file.read((char *)unitcell, 48);
    val = read_int(file); 
    if(val != 48)
			m_error = true;	
	
    m_box = BoxSize(unitcell[0],unitcell[2],unitcell[5]);
    // check for errors
    if (!file.good()||m_error)
        {
        cerr << endl << "***Error! Error reading DCD frame header" << endl << endl;
        throw runtime_error("Error reading DCD file");
        }
    }


void Dcdbuilder::read_frame_data()
    {
    // write x coords
    m_nparticles  = read_int(file)/sizeof(float) ;
	m_x_array.resize(m_nparticles);
    file.read((char *)&m_x_array[0], m_nparticles * sizeof(float));
    m_nparticles  = read_int(file)/sizeof(float) ;

    m_nparticles  = read_int(file)/sizeof(float) ;
	m_y_array.resize(m_nparticles);
    file.read((char *)&m_y_array[0], m_nparticles * sizeof(float));
    m_nparticles  = read_int(file)/sizeof(float) ;
	

    // write z coords
    m_nparticles  = read_int(file)/sizeof(float) ;
	m_z_array.resize(m_nparticles);
    file.read((char *)&m_z_array[0], m_nparticles * sizeof(float));
    m_nparticles  = read_int(file)/sizeof(float) ;
 
    // check for errors
    if (!file.good())
        {
        cerr << endl << "***Error! Error reading DCD frame data" << endl << endl;
        throw runtime_error("Error reading DCD file");
        }
		
	unsigned int np = m_pos.size();
	if(m_nparticles!=np)
		{
		cerr << "***Error! The number of particles for updating position "<< m_nparticles <<" is not equal to the one of the stored frame "<< np <<" !"<< endl << endl;
		throw runtime_error("Error Dcdbuilder::read_frame_data");			
		}	

	double Lx = double (m_box.lx);
	double Ly = double (m_box.ly);		 
	double Lz = double (m_box.lz);
	
	double Lxinv = 0.0;
	double Lyinv = 0.0;
	double Lzinv = 0.0;	
	if(Lx!=0.0)
		Lxinv = 1.0/Lx;
	if(Ly!=0.0)
		Lyinv = 1.0/Ly;
	if(Lz!=0.0)
		Lzinv = 1.0/Lz;	
	if(m_image.size()!=np)
		m_image.resize(np);	
	for (unsigned int i=0; i<np; i++)
		{
		double x = (double) m_x_array[i];
		double y = (double) m_y_array[i];
		double z = (double) m_z_array[i];
		
		int ix = (int) rint(x*Lxinv);
		int iy = (int) rint(y*Lyinv);	
		int iz = (int) rint(z*Lzinv);
		
		x -= double(ix)*Lx;
		y -= double(iy)*Ly;
		z -= double(iz)*Lz;

		m_pos[i] = vec(x,y,z);
		m_image[i] = vec_int(ix, iy, iz);
		}
    }


unsigned int Dcdbuilder::getNframes() const
    {
    return m_num_frames_written;
    }


unsigned int Dcdbuilder::getTimeStep() const
    {
    return (m_frames_count-1)*m_period + m_start_timestep;
    }
	
std::string Dcdbuilder::getFilename()
	{		
	unsigned int timestep = getTimeStep();
	std::ostringstream oss;
	oss<<"timestep." << setfill('0') << setw(10) << timestep;	
	return oss.str();
	}

void Dcdbuilder::updateDataFromDCD()
    {
	if(m_frames_count<m_num_frames_written)
		{
		read_frame_header();
		read_frame_data();
		}
	m_frames_count += 1;
	if(m_frames_count>=m_num_frames_written)
		file.close();
	}

