
#include "Xmlbuilder.h"

 
 int main(int argc,char* argv[])
   {
 	std::vector<std::string> filename;
	for(unsigned int i=0; i<(unsigned int)argc;i++ )
		{
		filename.push_back(argv[i]);
		}     
	for (unsigned int i=1;i<filename.size(); i++)
		{
		// Generate a filename with the timestep padded to ten zeros
		string xml_open = filename[i];
		Xmlbuilder build(xml_open.c_str());
		unsigned int N = build.getNumParticles();

		
		BoxDim box = build.getBox();
		std::vector<vec> pos = build.getPos();
		std::vector<vec> si = build.getOrientation();
		std::vector<unsigned int> type = build.getType();

		float Lx = box.xhi - box.xlo;
		float Ly = box.yhi - box.ylo;
		float Lz = box.zhi - box.zlo;
		
		for (unsigned int j=0; j<N; j++)
			{
			pos[j].x += Lx/2.0;
			pos[j].y += Ly/2.0;
			pos[j].z += Lz/2.0;			
			}

	
		string::size_type xp = xml_open.find("xml");
		string outs;
		if (xp!=xml_open.npos)
			outs = xml_open.replace(xp,xp+3, "pdb");
		else
			outs = xml_open+".pdb";
	
		ofstream s(outs.c_str());
		
		if (!s.good())
			{
			cerr << endl << "***Error! Unable to open dump file for writing: " << outs << endl << endl;
			throw runtime_error("Error writting pdb dump file");
			}	
		
		for (unsigned int i = 0; i< N; i++)
		{
		unsigned int typ = type[i]; 
			if(typ == 0)
					{
					if(i==0)
					{
					s <<"CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1"<<"\n";
					}
					s <<"ATOM"<<setfill(' ')<<setw(7)<<3*i+1<<"  "<<"N"<<"       "<<"X"<<"   "<<"0"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(3)<<setw(12)<<pos[i].x<<setw(8)<<pos[i].y<<setw(8)<<pos[i].z<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<setw(6)<<0.0<<setw(6)<<0.0<<"\n";
					s <<"ATOM"<<setfill(' ')<<setw(7)<<3*i+2<<"  "<<"O"<<"       "<<"X"<<"   "<<"0"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(3)<<setw(12)<<pos[i].x+si[i].x*0.01<<setw(8)<<pos[i].y+si[i].y*0.01<<setw(8)<<pos[i].z+si[i].z*0.01<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<setw(6)<<0.0<<setw(6)<<0.0<<"\n";
	//                s <<"ATOM"<<setfill(' ')<<setw(7)<<3*i+3<<"  "<<"O"<<"       "<<"X"<<"   "<<"0"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(3)<<setw(12)<<pos[i].x-si[i].x*0.01<<setw(8)<<pos[i].y-si[i].y*0.01<<setw(8)<<pos[i].z-si[i].z*0.01<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<setw(6)<<0.0<<setw(6)<<0.0<<"\n"; 
			}
			}

		


		if (!s.good())
			{
			cerr << endl << "***Error! Unexpected error writing polymer dump file" << endl << endl;
			throw runtime_error("Error writting polymer dump file");
			}
			
		s.close();		
		
	}
}

  

