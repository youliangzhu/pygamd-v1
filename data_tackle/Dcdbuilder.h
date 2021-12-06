#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>


#include "mst_reader.h"

using namespace std;

#ifndef __DCD_BULDER_H__
#define __DCD_BULDER_H__


class Dcdbuilder : public mst_reader
    {
    public:
        //! Loads in the file and parses the info
        Dcdbuilder();
		
        void initiate(const string &fname); 	
         //! Helper function to read the input file
        void read_frame_header();
		
        void read_frame_data();
		
        void read_file_header();  	
		
        void updateDataFromDCD(); 
		
        virtual unsigned int getNframes() const;
		
        virtual unsigned int getTimeStep() const;
		
		virtual std::string getFilename();	
    private:
        fstream file;          

		unsigned int m_num_frames_written;
		unsigned int m_frames_count;
		unsigned int m_nparticles;
		unsigned int m_last_written_step;
		unsigned int m_period;
		unsigned int m_start_timestep;
		bool m_error;
        std::vector< float > m_x_array;            //!< x position of all particles loaded
        std::vector< float > m_y_array;            //!< y position of all particles loaded
        std::vector< float > m_z_array;            //!< z position of all particles loaded
    };



#endif



