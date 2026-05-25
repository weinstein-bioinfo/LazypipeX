#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
//#include <cstring>

#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <getopt.h>
#include <seqan/seq_io.h>


/* LAZYPIPE PROJECT: C++ CODE FOR NGS PIPELINE (2025)
 *
 * Retrieve contigs mapped to a given taxid based on results in resdir
 *
 * Author: Ilya Weinstein, University of Helsinki (2025)
 * Creadit: Lazypipe project, https://doi.org/10.1093/ve/veaa091
 */
using namespace std;

void print_usage(const string name){
	std::cerr << "\nUSAGE: "<< name <<" -t taxid(s) [-r resdir -v] 1> contigs.fa\n"
		<<"\n"
		<<"Retrieve contigs assigned to a given taxid\n"
		<<"\n"
		<<"-s str            : Species name\n"
		<<"-t str            : staxid id (or CSV of ids) from annot_table.tsv\n"
		<<"-r dir            : Directory with pipeline results. Default: results.\n"
		<<"                    MUST include files: annot_table.tsv + contings.fa\n"
		<<"-v                : Verbal mode. Default: false\n"
		<<"-h                : Print this user manual\n"
		<<"\n"
		<<"output            : Retrieved contigs are printed to STDOUT\n"
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

inline std::string get_suffix(const std::string& str, const char delim){
	unsigned int pos = str.find_first_of(delim);
	if(pos > str.length()){
		return str;
	}
	else{
		return str.substr(pos+1, std::string::npos);
	}
}


// Test file for accessibility/existance
inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

/**
  * Returns column index corresponding to column name.
  *
  * file 		tsv-file with headers
  * colname     colname to search for
  * returns     colindex or -1 if colname was not found
  */
unsigned int get_colname_ind(const string file, const string colname){
	ifstream in(file, ios::in);
	if(!in.is_open()){
		std::cerr << "\nERROR: failed to open \'"<<file<<"\'\n";
		exit(1);
	}
	char buffer[10000];
	vector<string> sp;
	while(in.getline(buffer,10000)){
		sp = split(string(buffer),'\t');
		for(unsigned int i=0 ; i<sp.size(); i++){
			if( sp[i]==colname){
				in.close();
				return i;
			}
		}
		std::cerr << "WARNING: tsv file "<<file<<" has no colname "<<colname<<"\n";
		in.close();
		return -1;
	}
	return -1;
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
	
	string prog_name			= string(argv[0]);	
	// paramenters
	string species				= "";
	string taxid					= "";
	string resdir				= "results";
	string contigs_file			= "";
	string annot_file			= "";
	bool verbal					= false;
	bool help					= false;
	
	// PARSING PARAMETERS
	int option= 0;
	while ((option = getopt(argc, argv,"s:t:c:r:p:1:2:vh")) != -1) {
        switch (option) {
        case 's':
        		species	= string(optarg);
        		break;
		case 't':
			taxid	= string(optarg);
			break;
		case 'r':
			resdir	= string(optarg);
			break;
		case 'v':
			verbal= true;
			break;
		case 'h':
			help = true;
			break;
		case '?':
        		std::cerr << "option -"<<optopt<<" requires an argument\n\n";
        		exit(1);
             	default:
	     		print_usage(prog_name); 
                	exit(1);
        	}
    	}
	if(help){
		print_usage(prog_name);
		exit(1);
	}
	if(species=="" && taxid==""){
		std::cerr<<"ERROR: please specify -t taxid or -s species_name\n"; print_usage(prog_name); exit(1);
	}
	if(resdir==""){ 		
		std::cerr<<"ERROR: please specify -r resdir\n\n"; print_usage(prog_name); exit(1);
	}

	// contid_taxid_[score] file
	if( exists(resdir+"/annot_table.tsv") ){
		annot_file = resdir+"/annot_table.tsv";
	}
	else{
		std::cerr << "\nERROR: missing file: "<< (resdir+"/annot_table.tsv") << "\n\n";
		exit(1);
	}
	// contig file
	if( exists(resdir+"/contigs.fa") ){
		contigs_file = resdir+"/contigs.fa";
	}
	else{
		std::cerr << "\nERROR: missing file: "<< (resdir+"/contigs.fa") << "\n\n";
		exit(1);
	}
	// DONE PARSING OPTIONS
		


	vector<string> contid_list;
	if( taxid != "" || species !="" ){
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

	
	// convert to a set (implemented as a map)
	unordered_map<string,bool> contid_set;
	for (auto it = contid_list.begin(); it != contid_list.end(); ++it) {
		contid_set[*it] = true;
	}
	

	// Retrieving Contigs by id
	
	if(verbal){ 
		std::cerr << "\t# reading "<<contigs_file<<"\n";
	}
	
	seqan::SeqFileIn fileIn( contigs_file.c_str() );
	seqan::StringSet<seqan::CharString> ids;
	seqan::StringSet<seqan::IupacString> seqs;
	unsigned int batch_size	= 1000;
	//unsigned int report_batch= batch_size*1000;
	unsigned int seq_read	= 0;
	unsigned int seq_sel	= 0;
	
	while(!atEnd(fileIn) ){
		try{
			readRecords(ids,seqs,fileIn,batch_size);
			seq_read += length(seqs);
			
			//writeRecords(fileOut,ids,seqs);
			for(unsigned int ind=0; ind<length(ids); ind++){
				seqan::CharString id 	= ids[ind];
				std::string id2	= std::string( toCString(id) );
				
				if(length(id) == 0){
					std::cerr << "WARNING: skipping seq with an empty id\n";
					continue;
				}
				id2 = get_prefix(id2,'_');
				id2 = get_suffix(id2,'=');
				// DEBUG
				//std::cerr << "id: '"<<id << "'\nid2: '"<< id2 <<"'\n"; exit(1);
				
				
				if(  contid_set.count(id2) > 0 ){
					//writeRecord(fileOut,ids[ind],seqs[ind]);
					//std::cerr << "\tfound contig: "<<id<<"\n";
					std::cout << ">"<<id2 <<"\n"
									<<seqs[ind]<<"\n";
					seq_sel++;
				}
			}
			clear(ids);
			clear(seqs);	
		}
		catch (seqan::IOError const & e){
			std::cerr << "ERROR: IOError:\n"<<e.what()<<"\n";
			return 1;}
		catch(seqan::ParseError const &e){
			std::cerr << "ERROR: badly formatted record:\n"<<e.what()<<"\n";
			return 1;
		}
	}	
	
}
