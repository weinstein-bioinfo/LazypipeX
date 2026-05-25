#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "bioio.h"


/* LAZYPIPE PROJECT: C++ CODE FOR NGS PIPELINE (2022)
 *
 * Retrieve reads for contig or taxid based on contig taxonomy assignments in resdir
 *
 * Author: Ilya Plyusnin, University of Helsinki (2022)
 * Creadit: Lazypipe project, https://doi.org/10.1093/ve/veaa091
 */
using namespace std;

void print_usage(const string name){
	cerr << "USAGE: "<< name <<" -t taxid(s) -c contid -s species [-r dir -a annot -1 read1 -2 read2]\n"
		<<"\n"
		<<"Retrieve reads for contig/species/taxid based on annotations in resdir\n"
		<<"\n\n"
		<<"-c str            : contig id\n"
		<<"-s str            : species name\n"
		<<"-t str            : taxid(s) (csv)\n"
		<<"-r dir            : directory with pipeline results [results]\n"
		<<"                    MUST include: readid_contigid.tsv + annot_table.tsv\n"
		<<"-a str            : annotation file [annot_table.tsv]\n"
		<<"-1 file           : forward reads [read1.trim[.hflt].fq.gz ]\n"
		<<"-2 file           : reverse reads [read2.trim[.hflt].fq.gz ]\n"
		<<"-w dir            : work directory [ . ]\n"
		<<"-p str            : output prefix [contid/species/taxid]\n"
		<<"-v                : verbal mode [false]\n"
		<<"\n\n"
		<<"Credit: Lazypipe project, https://doi.org/10.1093/ve/veaa091\n\n"
		<<"\n";
}

inline std::vector<std::string> split(const std::string& str, const char delim)
{
    std::stringstream ss {str};
    std::string item;
    std::vector<std::string> result {};
    while (std::getline(ss, item, delim)){
    	if( isspace(item.back()) ){
		item.erase(item.length()-1,1);
	}
    	result.emplace_back(item);
	}

    return result;
}


// Test file for accessibility/existance
inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

/*
 * Reads key-value pairs from tsv file to an unordered_map-map structure.
 */
void readto_map_map(const string file, unsigned int keyi, unsigned int vali, unordered_map<string,unordered_map<string,bool> > &map_map){

	ifstream in(file, ios::in);
	if(!in.is_open()){
		cerr << "\nERROR: failed to open \'"<<file<<"\'\n\n"; exit(1);
	}
	char buffer[10000];
	vector<string> sp;
	while(in.getline(buffer,10000)){
		if(buffer[0] == '@'){
			continue;
		}
		sp = split(string(buffer),'\t');
		if(  keyi<sp.size() && vali<sp.size()){
			if( map_map.count(sp[keyi]) == 0){
				unordered_map<string,bool> map;
				map_map[sp[keyi]] = map;
			}
			map_map[sp[keyi]][sp[vali]] = true;
		}	
	}
	in.close();
}
/*
 * This overload reads header and calls basic version.
 */
void readto_map_map(const string file, const string key_header, const string val_header, unordered_map<string,unordered_map<string,bool> > &map_map){

	ifstream in(file, ios::in);
	if(!in.is_open()){
		cerr << "\nERROR: failed to open \'"<<file<<"\'\n"; exit(1);
	}
	char buffer[10000];
	vector<string> sp;

	//read header
	int keyi = -1;
	int vali = -1;
	while(in.getline(buffer,10000)){
		sp = split(string(buffer),'\t');
		for(unsigned int col=0; col<sp.size(); col++){
			if(sp[col] == key_header){
				keyi = col;
			}
			if(sp[col] == val_header){
				vali = col;
			}
			if(keyi>=0 && vali>=0)
				break;
		}
		break;
	}
	in.close();

	if(keyi<0 || vali<0){
		cerr << "\nERROR: missing key/value columns "<<key_header<<"/"<<val_header<<" from file "<<file <<"\n\n";
		exit(1);
	}

	readto_map_map(file,keyi,vali,map_map);
}

int main(int argc, char** argv) {
	
	string prog_name= string(argv[0]);	
	// paramenters
	string taxid		= "";
	string contid		= "";
	string species      = "";
	string resdir		= "results";
	string wrkdir		= ".";
	string prefix		= "";
	string annot_file   = "";
	string read1		= "";
	string read2		= "";
	string read1_out	= "";
	string read2_out	= "";
	bool verbal	= false;
	
	// PARSING COMMANDLINE OPTIONS
	int option= 0; // -t taxid -c contid -r resdir [-1 read1 -2 read2 -w wrkdir -p prefix]\n"
	while ((option = getopt(argc, argv,"a:c:s:t:r:w:p:1:2:v")) != -1) {
        switch (option) {
        case 'a':
        	annot_file = string(optarg);
        	break;
		case 'c':
			contid	= string(optarg);
			break;
		case 's':
			species = string(optarg);
			break;
		case 't':
			taxid	= string(optarg);
			break;
		case 'r':
			resdir	= string(optarg);
			break;
		case '1':
			read1	= string(optarg);
			break;
		case '2':
			read2	= string(optarg);
			break;
		case 'w':
			wrkdir	= string(optarg);
			break;
		case 'p':
			prefix	= string(optarg);
			break;
		case 'v':
			verbal= true;
			break;	
		case '?':
        		cerr << "option -"<<optopt<<" requires an argument\n\n";
        		exit(1);
             	default:
	     		print_usage(prog_name); 
                	exit(1);
        	}
    	}
	if(taxid=="" && contid=="" && species==""){
		cerr<<"ERROR: missing options -t taxid or -c contid or -s species\n\n";
		print_usage(prog_name); exit(1);
	}
	if(resdir==""){
		cerr<<"ERROR: missing options -r resdir\n\n";
		print_usage(prog_name); exit(1);
	}
	if(read1 == ""){
		/* MAYBE IN DISTANT FUTURE: READING GZIPED FILES
		if( exists(resdir + "/read1.trim.hflt.fq.gz") &&  exists(resdir + "/read2.trim.hflt.fq.gz")){
			read1	= resdir + "/read1.trim.hflt.fq.gz";
			read2	= resdir + "/read1.trim.hflt.fq.gz";
		}
		else if( exists(resdir+"/read1.trim.fq.gz") && exists(resdir+"/read1.trim.fq.gz") ){
			read1	= resdir + "/read1.trim.fq.gz";
			read2	= resdir + "/read1.trim.fq.gz";
		}
		*/
		if( exists(resdir + "/reads/read1.trim.hflt.fq") &&  exists(resdir + "/reads/read2.trim.hflt.fq")){
			read1	= resdir + "/reads/read1.trim.hflt.fq";
			read2	= resdir + "/reads/read2.trim.hflt.fq";
		}
		else if( exists(resdir+"/reads/read1.trim.fq") && exists(resdir+"/reads/read2.trim.fq") ){
			read1	= resdir + "/reads/read1.trim.fq";
			read2	= resdir + "/reads/read2.trim.fq";
		}
		else{
			cerr << "\nERROR: missing input reads: "+resdir+"/reads/read*.trim[.hflt].fq"+"\n\n";
			print_usage(prog_name);
			exit(1);
		}
	}
	if(read2 == ""){
		read2	= read1;
		size_t found = read2.find("_R1");
		if(found != string::npos){
			read2.replace(found,3,"_R2");
		}
		else{
			cerr << "\nERROR: could not guess read2 from read1\n\n";
			print_usage(prog_name);
			exit(1);
		}
		if( !exists(read2) ){
			cerr << "\nERROR: missing read2 file: "<< read2 << "\n\n",
			print_usage(prog_name);
			exit(1);
		}
	}
	if( prefix == "" ){
		if( contid!=""){
			prefix = contid;
		}
		else if(species != ""){
			prefix = species;
			replace(prefix.begin(), prefix.end(), ' ', '_');
		}
		else if(taxid !=""){
			prefix = taxid; }
		else{
			prefix = ""; }
	}
	read1_out= resdir +"/reads/"+prefix+"_r1.fq";
	read2_out= resdir +"/reads/"+prefix+"_r2.fq";

	// annot_file
	if( annot_file == "" ){
		annot_file = resdir+"/annot_table.tsv";
	}
	if( !exists(annot_file) ){
		cerr << "\nERROR: missing annot_table: "<<annot_file<< "\n\n";
	}
	// DONE PARSING OPTIONS


	// -c contid overrides -t taxid(s)
	vector<string> contid_list;
	if(contid !=""){
		contid_list.push_back( contid );
	}
	// -t taxid(s) or -s species
	else if( taxid != "" || species !="" ){
		if(verbal){ cerr << "\t# reading file: "<< annot_file << "\n";}
		
		unordered_map<string,unordered_map<string,bool> > taxid_contid_map;
		unordered_map<string,unordered_map<string,bool> > species_contid_map;
		if(taxid != ""){
			readto_map_map(annot_file, "staxid", "contig", taxid_contid_map);
		}
		if(species != ""){
			readto_map_map(annot_file, "species", "contig", species_contid_map);
		}

		vector<string> taxid_list 	= split(taxid,',');
		vector<string> contid_list2;

		if(species != ""){
			if( species_contid_map.count(species) > 0){
				for(auto it = species_contid_map[species].begin(); it != species_contid_map[species].end(); ++it ){
					contid_list2.push_back( it->first );
				}
			}
			if( contid_list2.size() == 0){
				cerr << "\tWARNING: no contigs found for "+species+"\n";
				exit(0);
			}
		}
		else if(taxid_list.size()>0){
			for(unsigned int i=0; i<taxid_list.size(); i++){
				if( taxid_contid_map.count(taxid_list[i]) > 0){
					for(auto it = taxid_contid_map[taxid_list[i]].begin(); it != taxid_contid_map[taxid_list[i]].end(); ++it ){
						contid_list2.push_back( it->first );
					}
				}
			}
			if( contid_list2.size() == 0){
				cerr << "\tWARNING: no contigs found for taxid(s) "+taxid+"\n";
				exit(0);
			}
		}

		contid_list.insert(contid_list.end(), contid_list2.begin(), contid_list2.end());
	}
	if(verbal){ cerr << "\t# contids found: "<<contid_list.size()<<"\n";}
	
	
	// READ contid2readid map:
	//
    string readid_contid_file = resdir+"/readid_contigid.tsv";
	if(verbal){ cerr << "\t# reading "<< readid_contid_file <<"\n";}
	unordered_map<string,unordered_map<string,bool> > contid_readid_map;
	readto_map_map(readid_contid_file, 1, 0, contid_readid_map);

	// CREATE mask for readids to retrieve:
	//
	unordered_map<string,bool> readid_map;
	
	for(unsigned int i=0; i<contid_list.size(); i++){
		if(contid_readid_map.count(contid_list[i]) == 0){
			fprintf(stderr,"\t# skipping contig with no reads: %s\n",contid_list[i].c_str());
			continue;
		}
		unordered_map<string,bool> map_tmp	= contid_readid_map[contid_list[i]];
		for(auto it = map_tmp.begin(); it != map_tmp.end(); ++it)
			readid_map[ it->first ] = true;
	}
	if(verbal){ cerr << "\t# readids found: "<<readid_map.size()<<"\n";}


	// Filtering reads
	vector<string> read_files;
	read_files.push_back(read1);
	read_files.push_back(read2);
	vector<string> read_files_out;
	read_files_out.push_back(read1_out);
	read_files_out.push_back(read2_out);
	
	for(unsigned int k=0; k<read_files.size(); k++){
		if(verbal){ cerr << "\t# processing "<<read_files[k]<<"\n";}
	
		std::ifstream fastq {read_files[k], std::ios::binary};
		std::ofstream fout (read_files_out[k], std::ofstream::out);
		unsigned int seqn		=  bioio::count_fastq_records(fastq);
		unsigned int seqn_corr	= 0;
		unsigned int sel_num	= 0;
		vector<string>	sp;
		
		for(unsigned int i=0; i<seqn; i++,seqn_corr++){
			auto seq = bioio::detail::read_fastq_record<string,string,string>(fastq);
			if(seq.name.length()== 0){
			// count_fastq_records counts "^@" lines in fastq, which can be larger than rec num due to matches in quality string
				break;
			}
			string name = bioio::detail::split(seq.name,' ')[0];
			//cerr << "name:"<<name<<"\n"; exit(1);
			if( name[0]=='@' ){
				name = name.substr(1);
			}
			if( name.substr(name.length()-2,name.length()) == "/1" || name.substr(name.length()-2,name.length()) == "/2"){
				name = name.substr(0,name.length()-2);
			}
			//cerr<< "seqname:\""<<name<<"\"\n";exit(1);
			
			if( readid_map.count(name)>0 ){
				fout << seq.name << "\n";
				fout << seq.seq << "\n";
				fout << "+\n";
				fout << seq.qual<< "\n";
				sel_num++;
			}
		}
		fastq.close();
		fout.close();
		if(verbal){
			fprintf(stderr,"\t# selected %u/%u (%2.2f%%) reads to %s\n",
				sel_num,seqn_corr,((sel_num+0.1)/seqn_corr)*100.0,read_files_out[k].c_str());
		}
	}
		
	
}




