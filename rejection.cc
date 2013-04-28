/******************************************************************************
 *  Copyright (C) 2007 Jeff Ross-Ibarra <rossibarra@gmail.com>			       	*
 *  Portions of the code (C) 2002-2007 K. Thornton and (C) 2007 R. Cartwright	*
 *                                                                     	   	*
 *  This program is free software: you can redistribute it and/or modify  		*
 *  it under the terms of the GNU General Public License as published by   	*
 *  the Free Software Foundation, either version 3 of the License, or      	*
 *  (at your option) any later version.                                    	*
 *                                                                         	*
 *  This program is distributed in the hope that it will be useful,        	*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         	*
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          	*
 *  GNU General Public License for more details.                           	*
 *                                                                         	*
 *  You should have received a copy of the GNU General Public License      	*
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  	*
 *****************************************************************************/

/* CURRENTLY ONLY WORKS FOR TWO POPS BECAUSE OF THE PRIVATE/SHARED CODE
	Todo: allow for multiple arbitrary pops and different summary stats
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <utility>
#include <numeric>
#include <math.h>
#include <SimDataTBS.hpp>
#include <Sequence/FST.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/PolyTableFunctions.hpp>

using namespace std;
using namespace Sequence;

void get_tbs( const SimDataTBS d, double * avtheta, double * avrho, double * avnm1, double * avnm2, std::vector<unsigned> & config);

struct sumstat {
	double mean, var;
};

bool is_accepted ( double cutoff, int dataloci, vector<sumstat> * sumdata, vector<double> observed, bool verb, bool meanonly );

void calcstat( double statistic, struct sumstat *statinfo, int num );

int main(int argc, char *argv[])
{
  	int opt=1, errors=0, maxloci=0, model=99;
  	bool verbose=false, meanonly=false, acceptall=false, rm_singletons=false;
  	vector< double > observed;
  	string line, temp;
	const char * ofn = "crap";
	const char * dfn = "crap";
	std::vector<unsigned int> config;
	double cutoff=-99;
	
	// get variables
	while( opt<argc ){
		switch( argv[opt][1] ){
			case 'd': //observed data file name
				opt++;
		    	dfn = argv[opt++];   
		  		break;
		  	case 'c': //observed data file name
				opt++;
		    	cutoff = atof(argv[opt++]);   
		  		break;
			case 'A': // accepta all sims regardless of cutoff
				opt++;
				acceptall=true;
		  		break;
			case 'f': //file name
				opt++;
		    	ofn = argv[opt++];   
		  		break;
			case 'l': // loci
				opt++;
				maxloci = atoi(argv[opt++]);
				break;
			case 'm': // loci
				opt++;
				model = atoi(argv[opt++]);
				break;
			case 'V': // verbose with fabs and data
				opt++;
				verbose=true;
				break;
			case 'M': // rejects based only on mean not variance
				opt++;
				meanonly=true;
				break;
			case 'q': // remove singletons
				opt++;
				rm_singletons=true;
				break;
			default : 
				++errors; 
				break;
		}
	} 
	
	// check files works
	ofstream o(ofn);                         
  	if (!o){                                                                       
  		cerr << "need an outfile\n";                                            
    	exit(1);                                                              
  	}                                                                     
  	o.close();
  	
   // check datafile and read in data
	ifstream data(dfn);                         
  	if (!data){                                                                       
  		cerr << "need an infile\n";                                            
    	exit(1);                                                              
  	}         
  	while  ( getline(data,line) ){
    	stringstream parser (line); 
    	float tempstat;
    	while(parser >> temp){
      	tempstat=atof(temp.c_str());
			observed.push_back(tempstat);	  		
      	parser >> std::ws;
    	}
	}
	data.close();
  	
	// check all variables present
	if( cutoff<0 || maxloci==0 || errors ){
		cerr << "\nOptions are:\n"
		<< "\n-l loci\n-c cutoff\n-A accepts all sims"
		<< "-d observed_data_file\n-f output_filename" 
		<< "\n-q remove singletons from simulated data before comparing\n-M rejects on mean only (not variances)\nand -V verbose data output\n\n";
		exit(1);
	}
  
  	// declare some libsequence objects.
  	SimDataTBS d;
  
  	// um, yeah.
	#if __GNUG__ && __GNUC__ >= 3
  		std::ios_base::sync_with_stdio(true);
	#endif
		
	// vectors etc. 				
 	int rv, nelements=0, loci=0, dataloci=0;
	vector< std::string > tbsparams(100,"");
	vector< double > meantheta(100,0.0), meanrho(100,0.0);
	double avrho=0.0, avtheta=0.0;
	vector< double > meannm1(100,0.0), meannm2(100,0.0);
	double avnm1=0.0, avnm2=0.0;
	vector < sumstat > sumdata(observed.size()/2);
	for( unsigned i=0; i<sumdata.size(); i++){
		sumdata[i].mean=0; sumdata[i].var=0;
	}
	assert(observed.size()/sumdata.size()==2);
	
	// read and calculate
	while( (rv=d.fromfile(stdin)) != EOF ){
      if(nelements == 100){
         	o.open(ofn,std::ios::app);
 	       	for(unsigned i=0;i<100;i++){
 	       		o << meantheta[i] << " " << meanrho[i] << " " 
 	       		<< meannm1[i] << " " << meannm2[i] << " " 
 	       		<< tbsparams[i] << endl;
 	      	}
         	o.close();
         	nelements=0;
      }      
	
      //we get the tbs line, calculate avg theta, rho, nm1, nm2, and get config[]
 		get_tbs( d, &avtheta, &avrho, &avnm1, &avnm2, config );

		// prevents segfaults from empty d, counts the locus but not for means/variance
		if(d.empty()){ 
			loci++; 
			if(loci==maxloci){
				if( is_accepted( cutoff, dataloci, &sumdata, observed, verbose, meanonly ) == true || acceptall == true ){ 
					tbsparams[nelements]=d.tbsline();
					meanrho[nelements]=avrho/maxloci;
					meantheta[nelements]=avtheta/maxloci;
					meannm1[nelements]=avnm1/maxloci;
					meannm2[nelements]=avnm2/maxloci;
					nelements++; 
				} 
				for( unsigned i=0; i<sumdata.size(); i++){
					sumdata[i].mean=0; sumdata[i].var=0;
				}
				loci=0; avtheta=0.0; avrho=0.0; avnm1=0; avnm2=0; dataloci=0;
			}			
			continue; 
		}
		dataloci++; loci++;
		unsigned total = std::accumulate(config.begin(),config.end(),0,plus<int>());       
	  	FST fst(&d, config.size(), &config[0]);
		
		//calculate temporary vector of FSTish stats	
	  // calcstat( fst.HBK(),  &sumdata[0], dataloci );
		calcstat( fst.shared(0,1).size(), &sumdata[0], dataloci );
		calcstat( fst.fixed(0,1).size(), &sumdata[1], dataloci );
		std::pair< std::set< double >, std::set< double > > temppriv = fst.Private(0,1);
		
		calcstat( temppriv.first.size(), &sumdata[2], dataloci );
		calcstat( temppriv.second.size(), &sumdata[3], dataloci );
		
	  	SimDataTBS d2;
	  	if(d.size() != total){
	  		std::cerr << "oh crap\n";
	    	exit(1);
	  	}
	  	int sum = 0;
	  	for(unsigned int i = 0 ; i < config.size() ; ++i)
	  	{	
	  		d2.assign(&*d.pbegin(),d.numsites(),&d[sum],config[i]);
				RemoveInvariantColumns(&d2);
				PolySIM P2(&d2); 
//				calcstat ( P2.DandVK(), &sumdata[4+i], dataloci);
//				calcstat ( P2.ThetaPi(), &sumdata[5+2*i], dataloci);
//				cout << "pi: " << P2.ThetaPi() << endl; // << sumdata[5+2*i] << endl;
	   }
		
	 	// check to see how well sims match observed	 	
		if(loci==maxloci){
	 		// if all elements meet criteria, keepme = size of vector observed
			if( is_accepted( cutoff, dataloci, &sumdata, observed, verbose, meanonly ) == true || acceptall == true ){ 
				tbsparams[nelements]=d.tbsline();
				meanrho[nelements]=avrho/maxloci;
				meantheta[nelements]=avtheta/maxloci;
				meannm1[nelements]=avnm1/maxloci;
				meannm2[nelements]=avnm2/maxloci;
				nelements++; 
			} 
			for( unsigned i=0; i<sumdata.size(); i++){
				sumdata[i].mean=0; sumdata[i].var=0;
			}
			loci=0; avtheta=0.0; avrho=0.0; avnm1=0; avnm2=0; dataloci=0;
		}
	}
	
	//write last runs to file if needed
	if(nelements>0)
   {
         	o.open(ofn,std::ios::app);
 	       	for(int i=0;i<nelements;i++){
 	       		o << model << " " << meantheta[i] << " " << meanrho[i] << " " 
 	       		<< meannm1[i] << " " << meannm2[i] << " " 
 	       		<< tbsparams[i] << endl;
       		}
         	o.close();
         	nelements=0;
	}
}

// get config and average theta, rho from tbsline.  see SimDataTBS for tbsline stuff
void get_tbs( const SimDataTBS d, double * avtheta, double * avrho, double * avnm1, double * avnm2, std::vector<unsigned> & config){
	stringstream tparams (d.tbsline());
	string tbsdata;
	int countdata=0;
	config.clear();
	double temprho=0, temptheta=0, tempbp=0;
	while(tparams >> tbsdata){
		if(countdata==1) temptheta=atof(tbsdata.c_str());
		else if(countdata==2) temprho=atof(tbsdata.c_str());
		else if(countdata==3) tempbp=atof(tbsdata.c_str());
		else if(countdata==4 || countdata==5) config.push_back(atoi(tbsdata.c_str()));
		else if(countdata==6) *avnm1+=atof(tbsdata.c_str());
		else if(countdata==7) *avnm2+=atof(tbsdata.c_str());
		tparams >> std::ws;
		countdata++;	
	}
	*avtheta+=temptheta/tempbp;
	*avrho+=temprho/tempbp;
	assert(countdata>6);
}

void calcstat( double statistic, struct sumstat *statinfo, int num )
/* variance algorithm from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
n = 0
mean = 0
S = 0
foreach x in data:
  n = n + 1
  delta = x - mean
  mean = mean + delta/n
  S = S + delta*(x - mean)      // This expression uses the new value of mean
end for
variance = S/(n - 1) */
{
	double delta = statistic - statinfo->mean;
	statinfo->mean+=(delta)/(num);
	statinfo->var+= delta*( statistic - statinfo->mean );
}

// if sumdata[].mean and sumdata[].var are within cutoff of observed[] for each sumstat, returns true, otherwise false;
bool is_accepted ( double cutoff, int dataloci, vector<sumstat> * sumdata, vector<double> observed, bool verb, bool meanonly  ) {
	unsigned keepme=0; double cut=0;
	for(unsigned i = 0 ; i < (*sumdata).size() ; ++i){
		//final step of variance calc 
		(*sumdata)[i].var=(*sumdata)[i].var/(dataloci-1);
		
		//sets cutoff to 10% or 0.1 or cutoff whichever is bigger.
		if(observed[i*2]==0) cut = cutoff < 0.1 ? 0.1 : cutoff;
		else cut = (cutoff*observed[i*2]) < 0.1 ? 0.1 : cutoff*observed[i*2];
		keepme += (fabs( ((*sumdata)[i].mean - observed[i*2]) )<=cut);
//			cout << "mean: " << observed[i*2] << " " << cutoff << " " << cut << " " << (*sumdata)[i].mean << " " << fabs( ((*sumdata)[i].mean - observed[i*2]) ) << " " << keepme << endl;

		//only rejects from variance if meanonly=false
		if(meanonly==false){
			if(observed[i*2+1]==0) cut = cutoff < 0.1 ? 0.1 : cutoff;
			else cut = (cutoff*observed[i*2+1]) < 0.1 ? 0.1 : cutoff*observed[i*2+1];
			keepme += (fabs( ((*sumdata)[i].var - observed[i*2+1]) )<=cut);
//			cout << "var: " << observed[i*2+1] << " " << cutoff << " " << cut << " " << (*sumdata)[i].var << " " << fabs( ((*sumdata)[i].var - observed[i*2+1]) ) << " " << keepme << endl;
		}
	}
	//prints fabs and simdata to stdout if -V is used
	if( verb==true){
		
		 cout << "\nfabs: ";
 		for(unsigned i = 0 ; i < (*sumdata).size() ; ++i){
			if(observed[i*2]==0) cout << fabs( ((*sumdata)[i].mean - observed[i*2]) ) << ' ' ;
	// 	else cout << fabs( ((*sumdata)[i].mean - observed[i*2]) ) << ' ';
			else cout << fabs( ((*sumdata)[i].mean - observed[i*2])/observed[i*2] ) << ' ';
	//		else cout << "\nsimdat: " << (*sumdata)[i].mean << "\nobs: " << observed[i*2] << "\nfabs: " << fabs( ((*sumdata)[i].mean - observed[i*2])/observed[i*2] ) << ' ';
		
			if(meanonly==false){
				if(observed[i*2+1]==0) cout << (*sumdata)[i].var << fabs( ((*sumdata)[i].var - observed[i*2+1]) ) << ' ';
	//			else cout << fabs( ((*sumdata)[i].var - observed[i*2+1]) ) << ' '; 
				else cout << fabs( ((*sumdata)[i].var - observed[i*2+1])/observed[i*2+1] ) << ' '; 
	//			else cout << "\nsimdat: " << (*sumdata)[i].var << "\nobs: " << observed[i*2+1] << "\nfabs: " << fabs( ((*sumdata)[i].var - observed[i*2+1])/observed[i*2+1] )  << ' ';
			}
		}
		
		cout << keepme << "\nsimdata: " ;
		for(unsigned i = 0 ; i < (*sumdata).size() ; ++i){
			cout << (*sumdata)[i].mean << ' ' ;
			if(meanonly==false) cout << (*sumdata)[i].var << ' ';
		}
			
/*		cout << "\nobserved: ";
 		for(unsigned i = 0 ; i < observed.size() ; ++i){
			cout << observed[i] << ' '; 
 		}*/
 		cout << endl;
 	}
 	if( meanonly==true ){ 
 		if( keepme==observed.size()/2 ) return true;
 		else return false;
 	}
 	else if( keepme==observed.size() ) return true;
	else return false;
}
