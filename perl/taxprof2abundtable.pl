#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Converts CAMI Taxonomic Profile to Lazypipe abund_table.tsv format
#

my $usage = "\nUSAGE: $0 taxprof.cami [--readn num | --readncol num] 1> abund_table.tsv\n".
			"\n".
			"taxprof.cami   : taxomic profile in CAMI format\n".
			"readn          : number of reads in the metagenome [100]\n".
			"readncol       : column (1-based) in taxprof.cami with readn [-1]. Overrides --readn option\n".
			"abund_table.tsv: abundance table\n".
			"\n";


my $readn_tot		= 100;
my $readn_col		= -1;
my $contign 		= "NA";
my $ignore_strains	= 1;

GetOptions('readn=i' => \$readn_tot, 'readncol=i' => \$readn_col) or die $usage;
if(scalar @ARGV<1) { die $usage; }
my $cami_file		= shift(@ARGV);
if($readn_col > 0){
	$readn_col--;
}
#print "readn_tot: $readn_tot\n";
#print "readn_col: $readn_col\n";exit(1);



print STDERR "# reading $cami_file\n";
open(IN,"<$cami_file") or die "Can\'t open $cami_file: $!\n";
my $l;
my $ln = 0;

print "readn\tcontign\ttaxid\t",
		"species\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\t",
		"species_id\tgenus_id\tfamily_id\torder_id\tclass_id\tphylum_id\tsuperkingdom_id\n";

while($l=<IN>){
	my $ln++;
	#chomp($l);
	if($l =~ m/^#/){ next; }
	if($l =~ m/^@/){ next; }
	
	my @sp = split(/\t/,$l,-1);
	if($sp[1] =~ m/species|strain/i){
	
		my $taxid		= $sp[0];
		my @taxids  	= split(/\|/,$sp[2],-1);
		my @taxnames	= split(/\|/,$sp[3],-1);
		my $pc			= $sp[4];
		my $readn		= 0;
		if($readn_col > -1){
			if(!defined($sp[$readn_col])){
				print STDERR "WARNING: undefined --readncol in line $ln:\n";
				print STDERR "$l\n";
				next;
			}
			$readn		= $sp[$readn_col];
		}
		else{
			$readn		= ($readn_tot*$pc/100.0);
		}
		
		print join("\t",$readn,
						$contign,
						$taxid,
						$taxnames[6],$taxnames[5],$taxnames[4],$taxnames[3],$taxnames[2],$taxnames[1],$taxnames[0],
						$taxids[6],$taxids[5],$taxids[4],$taxids[3],$taxids[2],$taxids[1],$taxids[0]), "\n";
	}
}
close(IN); 


