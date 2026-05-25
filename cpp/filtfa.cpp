#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>

#include <seqan/seq_io.h>


using namespace seqan;


/* LAZYPIPE PROJECT: C++ CODE FOR NGS PIPELINE (2022)
 *
 * Select (or filter) fasta seqs from a large DB such as NCBI nt based on a list of accessions.
 *
 * Author: Ilya Plyusnin, University of Helsinki (2022)
 * Creadit: Lazypipe project, https://doi.org/10.1093/ve/veaa091
 */


void print_usage(char* name){
      std::cerr << "\nUSAGE: "<< name <<" -i fasta.in -o fasta.out -f flt  [-m mode] [-v]\n\n"
		<<"\n"

		<<"-i  [fasta]      : input fasta sequences\n"
		<<"-o  [fasta]      : filtered/selected fasta sequences\n"
		<<"-f  [file]       : new-line separated list of sequence ids to filter/select\n"
		<<"-m  [str]        : operation mode (optional):\n"
		<<"                      filter: remove seqs listed in -f flt [default]\n"
		<<"                      select: select seqs listed in -f flt\n\n"
		<<"Credit: Lazypipe project, https://doi.org/10.1093/ve/veaa091\n\n"
		<<"\n";
}



/*
 * Returns string prefix up to the first occurance of delimiter
 */
inline std::string get_prefix(const std::string& str, const char delim){
	
	unsigned int pos = str.find_first_of(delim);
	if(pos > str.length()){
		return str;
	}
	else{
		return str.substr(0,pos);
	}
}


int main (int argc, char** argv) {

	if(argc<5){ print_usage(argv[0]); exit(1);}

	char* prog_name 	= argv[0];
	CharString fileNameIn	= "";
	CharString fileNameOut	= "";
	CharString fileNameFlt	= "";
	bool filter		= true;
	bool verbal		= false;
	CharString tmp;
	
	int option= 0;
	while ((option = getopt(argc, argv,"i:o:f:m:v")) != -1) {
        switch (option) {						
		case 'i':
			fileNameIn = CharString(optarg);
			break;
		case 'o':
			fileNameOut = CharString(optarg);
			break;
		case 'f':
			fileNameFlt = CharString(optarg);
			break;
		case 'm':
			tmp = CharString(optarg);
			if(tmp == "select"){
				filter = false;
			}
			else if(tmp == "filter"){
				filter = true;
			}
			else{
				std::cerr << "ERROR: invalid arg -m "<<tmp<<"\n";
				exit(1);
			}
			break;
		case 'v':
			verbal = true;
			break;
             	default:
	     		print_usage(prog_name); 
                	exit(1);
        	}
    	}
	
	
	/* READING SEQ ACCESSIONS TO SELECT/FILTER */
	if(verbal){
		std::cerr << "\treading "<<fileNameFlt;
	}
	std::unordered_map<std::string,bool> flt_map;	
	std::ifstream in( toCString(fileNameFlt), std::ios::in);
	if(!in.is_open()){
		std::cerr << "ERROR: failed to open \'"<< fileNameFlt <<"\'\n"; exit(1);
	}
	char buffer[1000];
	while(in.getline(buffer,1000)){
		std::string seqid = std::string(buffer);
		flt_map[seqid] = true;	
	}
	in.close();
	
	if(verbal){
		std::cerr << ": found "<< flt_map.size() << " uniq names\n";
	}
	
	
	/* READING FASTA IN + SELECTING/FILTERING FASTA OUT */
	std::cerr << "# reading "<< fileNameIn <<"\n";
	
	SeqFileIn fileIn( toCString(fileNameIn) );
	SeqFileOut fileOut( toCString(fileNameOut) );
	StringSet<CharString> ids;
	StringSet<IupacString> seqs;
	unsigned int batch_size	= 1000;
	unsigned int report_batch= batch_size*1000;
	unsigned int seq_read	= 0;
	unsigned int seq_flt	= 0;
	unsigned int seq_sel	= 0;
	
	while(!atEnd(fileIn) ){
		try{
			readRecords(ids,seqs,fileIn,batch_size);
			seq_read += length(seqs);
			
			//writeRecords(fileOut,ids,seqs);
			for(unsigned int ind=0; ind<length(ids); ind++){
				CharString id 	= ids[ind];
				std::string id2	= std::string( toCString(id) );
				
				if(length(id) == 0){
					std::cerr << "WARNING: skipping seq with an empty id\n";
					seq_flt++;
					continue;
				}
				id2 = get_prefix(id2,' ');
				// DEBUG
				//std::cerr << "id: '"<<id << "'\nid2: '"<< id2 <<"'\n"; exit(1);
				
				
				if( (filter && flt_map.count(id2)==0) || (!filter && flt_map.count(id2)>0) ){
					writeRecord(fileOut,ids[ind],seqs[ind]);
					seq_sel++;
				}
				else{
					seq_flt++;
				}
			}
			clear(ids);
			clear(seqs);
			
			if( verbal && (seq_read%report_batch == 0) ){
				std::cerr << "\t"<<(seq_read/1000)<<"k seqs read\n";
			}		
		}
		catch (IOError const & e){
			std::cerr << "ERROR: IOError:\n"<<e.what()<<"\n";
			return 1;}
		catch(ParseError const &e){
			std::cerr << "ERROR: badly formatted record:\n"<<e.what()<<"\n";
			return 1;
		}
	}
	
	if(verbal){
		if(filter){
			fprintf(stderr,"# removed %u/%u (%2.2f%%) reads from %s\n",seq_flt,seq_read,((seq_flt+0.0001)/seq_read)*100.0, toCString(fileNameIn) );
		}
		else{
			fprintf(stderr,"# selected %u/%u (%2.2f%%) reads from %s\n",seq_sel,seq_read,((seq_sel+0.0001)/seq_read)*100.0, toCString(fileNameIn) );
		}	
	}
	
	return 0;
}
	
