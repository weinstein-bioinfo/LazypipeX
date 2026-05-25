#include <fstream>
#include <getopt.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>

#include <seqan/seq_io.h>

using namespace seqan;


// LAZYPIPE PROJECT: C++ CODE FOR NGS PIPELINE (2022)
//
// FILTERING A SET OF READS FROM *.FASTQ FILE
// Author: Ilya Plyusnin, University of Helsinki (2022)
// Creadit: Lazypipe project, https://doi.org/10.1093/ve/veaa091





void print_usage(const std::string name){
      std::cerr << "\nUSAGE: "<< name <<" -1 read1 [-2 read2] -o read1 [-O read2] [-s suffix] -f flt  [-m mode]\n"
		<<"\n"

		<<"-1  [fastq]      : fastq-file with reads, can be forward reads for paired-end data\n"
		<<"-2  [fastq]      : fastq-file with reads, can be reverce reads for paired-end data, can be ommited\n"
		<<"-o  [fastq]      : print filtered -1 reads to this file (takes precidence over -s option)\n"
		<<"-O  [fastq]      : print filtered -2 reads to this file\n"
		<<"-s  [str]        : print filtered reads to read1_suffix [and read2_suffix]\n"			
		<<"-f  [file]       : new-line separated list of read names to filter. For example a list of refgen reads.\n"
		<<"-m  [str]        : optional str parameter, possible values:\n"
		<<"                      filter: filter reads in the flt file [default]\n"
		<<"                      select: select reads in the flt file\n\n"
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

int main(int argc, char** argv) {
	
	std::string progname	= std::string(argv[0]);
	if(argc < 4){
		print_usage(progname);
		exit(1);
	}	
	
	// paramenters
	std::string in_read1		= "";
	std::string in_read2		= "";
	std::string out_read1		= "";
	std::string out_read2		= "";
	std::string flt				= "";
	std::string suffix			= "";
	std::string mode			= "";
	bool filter					= true;
	bool verbal					= false;
	
	int option= 0;
	while ((option = getopt(argc, argv,"1:2:o:O:s:f:m:v")) != -1) {
        switch (option) {						
		case '1':
			in_read1 = std::string(optarg);
			break;
		case '2':
			in_read2 = std::string(optarg);
			break;
		case 'f':
			flt  = std::string(optarg);
			break;
		case 's':
			suffix = std::string(optarg);
			break;
		case 'o':
			out_read1 = std::string(optarg);
			break;
		case 'O':
			out_read2 = std::string(optarg);
			break;
		case 'm':
			mode = std::string(optarg);
			if(mode == "select"){
				filter = false;
			}
			else if(mode == "filter"){
				filter = true;
			}
			else{
				std::cerr << "ERROR: invalid arg -m "<<mode<<"\n";
				exit(1);
			}
			break;
		case 'v':
			verbal = true;
			break;
             	default:
	     		print_usage(progname); 
                	exit(1);
        	}
    	}
	
	// CHECK ARGUMENTS
	if(out_read1=="" && suffix!=""){
		out_read1 = in_read1+"_"+suffix; 
		out_read2 = in_read2+"_"+suffix;
	}
	if(in_read1 == ""){
		std::cerr << "ERROR: missing input reads\n\n";
		print_usage(progname);
		exit(1);
	}
	if(out_read1 == ""){
		std::cerr << "ERROR: missing arguments: -o read1 [-O read2] OR -s suffix\n\n";
		print_usage(progname);
		exit(1);
	}

	
	/* READING SEQ ACCESSIONS TO SELECT/FILTER */
	if(verbal){
		std::cerr << "\treading "<<flt;
	}
	std::unordered_map<std::string,bool> flt_map;	
	std::ifstream in( flt, std::ios::in);
	if(!in.is_open()){
		std::cerr << "ERROR: failed to open \'"<<flt<<"\'\n"; exit(1);
	}
	char buffer[1000];
	while(in.getline(buffer,1000)){
		std::string seqid 	= std::string(buffer);
		flt_map[seqid] 		= true;	
	}
	in.close();
	if(verbal){
		std::cerr << ": found "<< flt_map.size() << " uniq names\n";	
	}


	// COLLECT FASTQ FILES TO PROCESS
	std::vector<std::string> reads_in;
	std::vector<std::string> reads_out;
	if(in_read1 !=""){
		reads_in.push_back( in_read1 );
		reads_out.push_back( out_read1 );
	}
	if(in_read2 !=""){
		reads_in.push_back( in_read2 );
		reads_out.push_back( out_read2 );
	}
	
	/* READING FASTQ IN + SELECTING/FILTERING FASTQ OUT */ 
	for(unsigned int k=0; k<reads_in.size(); k++){
		if(verbal){
			std::cerr << "\tprocessing "<<reads_in[k]<<"\n";
		}
	
		SeqFileIn fileIn( toCString(reads_in[k]) );
		SeqFileOut fileOut( toCString(reads_out[k]) );
		StringSet<CharString> ids;
		StringSet<IupacString> seqs;
		StringSet<CharString> quals;
		unsigned int batch_size	= 1000;
		unsigned int report_batch= batch_size*1000;
		unsigned int seq_read	= 0;
		unsigned int seq_flt	= 0;
		unsigned int seq_sel	= 0;	

		while(!atEnd(fileIn) ){
			try{
				readRecords(ids,seqs,quals,fileIn,batch_size);
				seq_read += length(seqs);

				//writeRecords(fileOut,ids,seqs);
				for(unsigned int ind=0; ind<length(ids); ind++){
					CharString id_char		= ids[ind];
					std::string id			= std::string( toCString(id_char) );

					if(length(id) == 0){
						std::cerr << "\tWARNING: skipping seq with an empty id\n";
						seq_flt++;
						continue;
					}
					id = get_prefix(id,' ');
					if(id.length()== 0){
						fprintf(stderr,"\tWARNING: skipping/removing read %u: empty/malformed @name field: '%s'\n",seq_read,toCString(id_char));
						seq_flt++;
						continue;
					}
					if( id[0]=='@' ){
						id = id.substr(1);
					}
					if( id.length()== 0){
						fprintf(stderr,"\tWARNING: skipping/removing read %u: empty/malformed @name field: '%s'\n",seq_read,toCString(id_char));
						seq_flt++;
						continue;
					}			
					if( id.length()>2 && (id.substr(id.length()-2,id.length()) == "/1" || id.substr(id.length()-2,id.length()) == "/2")){
						id = id.substr(0,id.length()-2);
					}


					if( (filter && flt_map.count(id)==0) || (!filter && flt_map.count(id)>0) ){
						writeRecord(fileOut,ids[ind],seqs[ind],quals[ind]);
						seq_sel++;
					}
					else{
						seq_flt++;
					}
				}
				clear(ids);
				clear(seqs);
				clear(quals);

				if( verbal && (seq_read%report_batch == 0) ){
					std::cerr << "\t"<<(seq_read/1000)<<"k seqs read\n";
				}		
			}
			catch (IOError const & e){
				std::cerr << "ERROR: IOError:\n"<<e.what()<<"\n";
				return 1;
			}
			catch(ParseError const &e){
				std::cerr << "ERROR: badly formatted record:\n"<<e.what()<<"\n";
				return 1;
			}
		}
		if(verbal){
			if(filter){
				fprintf(stderr,"\t removed %u/%u (%2.2f%%) reads from %s\n",seq_flt,seq_read,((seq_flt+0.0001)/seq_read)*100.0, reads_in[k].c_str() );
			}
			else{
				fprintf(stderr,"\t selected %u/%u (%2.2f%%) reads from %s\n",seq_sel,seq_read,((seq_sel+0.0001)/seq_read)*100.0, reads_in[k].c_str() );
			}
		}	
	}
		
	
}




