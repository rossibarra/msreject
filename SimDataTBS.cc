#include <SimDataTBS.hpp>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace Sequence {

	SimDataTBS::SimDataTBS(): SimData() {}

	std::ostringstream is;
	int SimDataTBS::fromfile( FILE * openedfile ){
		char ch; 				
		int rv; //return value from fscanf   	
		while(1)	{
			rv = fscanf(openedfile,"%c",&ch);
			if (rv == EOF) return rv;	
			if ( ch=='/') {
				int token;
				float tempval;
				bool done = false;
		
				while (! done)	{
					token = fgetc(openedfile);
					if (  (!isspace(token)) && char(token) != '\n' && char(token) != '/' ) {
						ungetc(token,openedfile);
						fscanf(openedfile,"%f",&tempval);                
						is << tempval << ' ';
					}
					else if ( char(token) == '\n' ){
						done = true;
					}
				}
			}	
			
			if (ch == ':') break;
		}
		__tbsline = is.str(); is.str(""); 	
		unsigned ss;
		rv = fscanf(openedfile,"%u",&ss);
		if (rv == EOF) return rv;
	
		std::vector<double> _positions;
		std::vector<std::string> _data;
				
		if (ss > 0) {
			_positions.resize(ss);
			
			while(1){
				rv=fscanf(openedfile,"%c",&ch);
				if (rv == EOF) return rv;	    
				if (ch == ':') break;
			}
			for (unsigned i = 0; i < ss; ++i){
				rv=fscanf(openedfile,"%lf",&_positions[i]);
				if (rv == EOF) return rv;
			}
			char *seq = new char[ss+2];
			while(1) {
				rv=fscanf(openedfile,"%s",seq);
				
				if (rv == EOF) {
					rv=1;
					break;
				}
				else if ( strcmp(seq,"//") == 0) { 		
					int token;
					float tempval;
					bool done = false;
	
					while (! done) {
						token = fgetc(openedfile);
						if (  (!isspace(token)) && char(token) != '\n' && char(token) != '/' ) {
							ungetc(token,openedfile);
							fscanf(openedfile,"%f",&tempval);                
							is << tempval << ' ';
						}
						else if ( char(token) == '\n' ){
							done = true;
						}
					}
					break;
				}
				else {
					_data.push_back( std::string(seq) );  
				}
			}
			delete [] seq;
		}
		else if (ss == 0){
			_positions.resize(0);
			
			char *seq = new char[ss+2];
			while(1) {
				rv=fscanf(openedfile,"%s",seq);
				if (rv == EOF) {
					rv=1;
					break;
				}
				if ( strcmp(seq,"//") == 0) { 		
					int token;
					float tempval;
					bool done = false;
	
					while (! done) {
						token = fgetc(openedfile);
						if (  (!isspace(token)) && char(token) != '\n' && char(token) != '/' ) {
							ungetc(token,openedfile);
							fscanf(openedfile,"%f",&tempval);                
							is << tempval << ' ';
						}
						else if ( char(token) == '\n' ){
							done = true;
						}
					}
					break;
				}
			}
			delete [] seq;
			_data.resize(0);
		}
		//assign data into base class
		PolyTable::assign(&_positions[0],ss,&_data[0],_data.size());
		return rv;
	}
}
