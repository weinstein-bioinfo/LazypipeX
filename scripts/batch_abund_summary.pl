#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Temp qw(tempdir);
use Sort::Naturally;
use Getopt::Long qw(GetOptions);
use List::Util;


# Creates a metasummary for a batch of Lazypipe(NP) runs/results.
# This version creates metasummaries from abund_table.tsv(s).
# 

my $usage = "\nUSAGE: $0 --maxdir num --ta|target ba|ph|vi  [-v] resdir 1> metasummary.csv\n".
			"\n".
			"maxdir  [int]    : limit summary to this many directories [false]\n".
			"target  [str]    : print metasummary for Bacteria (\"ba\"), Phages (\"ph\") or Viruses (\"vi\") [vi]\n".
			"resdir  [dir]    : root result directory containing results for each sample\n".
			"toptax  [str]    : Taxrank used for selecting Bacteria/Viruses/Phages [division]\n".
			"numth   [int]    : Number of threads\n".
			"wrkdir  [dir]    : Root for temporary directories and files [./wrkdir]\n".
			"v                : verbal mode [false]\n".			
			"\n".
			"Collects bacteria/bphage/virus hits from  abund_tables.tsv(s) in res-dir/sample-dir(s) and prints metasummary\n\n";
		
# PARAMS
my $maxdir				= !1;
my $target				= "vi";
my $toptaxrank			= "division";
my $verbal				= 0;
my $numth				= 8;
my $wrkdir				= './wrkdir';
my $CLEANUP				= 1;
GetOptions(	'maxdir=i'		=> \$maxdir,
			'ta|target=s'	=> \$target,
			'toptax=s'		=> \$toptaxrank,
			'numth=i'		=> \$numth,
			'wrkdir=s'		=> \$wrkdir,
			'v'         		=> \$verbal);
			
# CHECK OPTIONS
if(scalar @ARGV < 1){ die $usage; }
if( !($target eq 'ba' || $target eq 'ph' || $target eq 'vi') ){
	die "ERROR: unknown option: --target $target\n";
}

# START WORKING
system("mkdir -p $wrkdir");
my $tmpdir			= tempdir("batch_summary_XXXXXXXX", DIR => $wrkdir, CLEANUP => $CLEANUP);
my $abund_tsv_ta		= "$tmpdir/abund_tsv.$target.tmp";
my $resdir			= shift(@ARGV);
my @sampledirs 		= get_dirs($resdir);
@sampledirs 			= sort { ncmp($a,$b) } @sampledirs;


# print header
print join(",","Sample","readn","contign","assembly.size","prob","Species","Genus","Family","Order","Class"),"\n";

# iterate sample-dirs
my $dirc			= 0;
foreach my $dir(@sampledirs){
	$dirc++;
	if($maxdir && ($dirc > $maxdir)){
		last;
	}
	# IN:
	my $sample  		= basename($dir);
	my $abund_tsv	= "$dir/abund_table.tsv";
	
	print STDERR "\treading $sample..\n";
	
	# check the $annot_tsv file exists and has annotations
	if(!(-e $abund_tsv) || get_line_count($abund_tsv) < 2){
		print STDERR "\t\tWARNING: $sample: missing/empty abund_table: skipping\n";
		next;
	}
	if($target eq 'ba' ){
		system_call("csvtk filter2 -tj $numth -f '\$$toptaxrank==\"Bacteria\"' $abund_tsv 1> $abund_tsv_ta", $verbal);
	}
	elsif($target eq 'ph'){
		system_call("csvtk filter2 -tj $numth -f '\$$toptaxrank==\"Phages\"' $abund_tsv 1> $abund_tsv_ta", $verbal);
	}
	else{
		system_call("csvtk filter2 -tj $numth -f '\$$toptaxrank==\"Viruses\"' $abund_tsv 1> $abund_tsv_ta", $verbal);
	}
	
	# check there are valid hits for the target
	if( get_line_count($abund_tsv_ta) < 2){
		print STDERR "\t\tWARNING: $sample: no hits for $target: skipping\n";
		next;}
	my %sp2readn		= read_tohash($abund_tsv_ta, "species","readn");
	my %sp2contign	= ();
	eval{	
		%sp2contign	= read_tohash($abund_tsv_ta, "species","contign");
		1;
	} or do{
		print STDERR "\t\tWARNING: $sample: missing contign values: ignoring\n";
	};
	my %sp2assembly	= ();
	my %sp2prob		= ();
	eval{	%sp2assembly		= read_tohash($abund_tsv_ta, "species","assembly.size");
		1;
	} or do {
		print STDERR "\t\tWARNING: $sample: missing assembly.size values: ignoring\n";
	};	
	eval{	%sp2prob		= read_tohash($abund_tsv_ta, "species","prob");
		1;
	} or do {
		print STDERR "\t\tWARNING: $sample: missing prob values: ignoring\n";
	};

	my %sp2genus		= read_tohash_uniq($abund_tsv_ta,"species","genus");
	my %sp2family	= read_tohash_uniq($abund_tsv_ta,"species","family");
	my %sp2order		= read_tohash_uniq($abund_tsv_ta,"species","order");
	my %sp2class		= read_tohash_uniq($abund_tsv_ta,"species","class");
			
	# iterate species in alphanum order and print summaries
	for my $sp(sort keys %sp2readn){
		my @readn		= split(/;/,$sp2readn{$sp});
		my @contign		= split(/;/,$sp2contign{$sp} || '0');
		my @assembly		= split(/;/,$sp2assembly{$sp} || '0');
		my @prob			= split(/;/,$sp2prob{$sp} || '0');
		my $genus		= $sp2genus{$sp};
		my $family		= $sp2family{$sp};
		my $order		= $sp2order{$sp};
		my $class		= $sp2class{$sp};
		my $sp_quoted	= '"'.$sp.'"';
		
		print join(",",$sample,sum(\@readn),sum(\@contign),sum(\@assembly),sum(\@prob),$sp_quoted,$genus,$family,$order,$class),"\n";	
	}
	
	# cleanup
	system("rm -f $abund_tsv_ta");
}




# Returns all directories in a path
sub get_dirs{
	my $path	= shift();
	my @dirs	= ();
	opendir( my $DIR, $path );
	while ( my $entry = readdir $DIR ) {
		next unless -d $path . '/' . $entry;
    	next if $entry eq '.' or $entry eq '..';
    	push(@dirs, "$path/$entry");
	}
	closedir $DIR;
	return @dirs;
}

#
# Read one-row tsv file with headers to a hash
# 
# USAGE:
# my %h = read_named_tsv($tsv_file)
# my $val1 = $h{'val1'}
# my $val2 = $h{'val2'}
#
sub read_named_tsv{
	my $file 	= shift();
	 
	open (IN, "<$file") or die "couldn't open the file $file $!";
	
	my $headers = <IN>;chomp($headers);
	my $values	= <IN>;chomp($values);
	if(!defined($headers)  || !defined($values)){
		die "ERROR: failed to read $file\n";
	}
	my @sp_h= split(/\t/, $headers, -1);
	my @sp_v= split(/\t/, $values, -1);
	if(scalar(@sp_h) != scalar(@sp_v)){
		die "ERROR: header and value columns do not match in: $file\n";
	}
	
	my %hash = ();
	for(my $i=0; $i<scalar(@sp_h); $i++){
		$hash{$sp_h[$i]} = $sp_v[$i];
	}
	return %hash;
}

# Read tsv file to hash. Supports multiple values (catenated with ";")
# 
# USAGE:
# my %hash = read_tohash(my_file.tsv,keycolname,valcolname)
# 
sub read_tohash{
	# INPUT
	my $file= shift;
	my $key = shift;
	my $val = shift;
	my $keyi = -1;
	my $vali = -1;
	my $sep = ';';
	
	# OUT
	my %hash= ();
	
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	
	# READ HEADER
	my $ln=0;
	my $l=<IN>;
	chomp($l);
	my @headers= split(/\t/,$l,-1);
	#print STDERR "headers: ",join(";",@headers),"\n";
	for(my $i=0; $i<scalar(@headers); $i++){
		if($headers[$i] eq $key){
			$keyi = $i;
		}
		if($headers[$i] eq $val){
			$vali = $i;
		}
	}
	if($keyi <0){
		die "ERROR: key col=$key not found in file $file\n";
	}
	if($vali <0){
		die "ERROR: value col=$val not found in file $file\n";
	}
	
	while($l=<IN>){
        	$ln++;
		#if($l =~ m/^[@#!]/){ next; }
		chomp($l);
		my @sp= split(/\t/,$l,-1);
		
		if( !defined($hash{$sp[$keyi]}) ){
        	$hash{$sp[$keyi]} = $sp[$vali];
		}
		else{
			$hash{$sp[$keyi]} = join($sep, $hash{$sp[$keyi]}, $sp[$vali]);
		}
	}
	close(IN);
	return %hash;
}
# Read tsv file to hash. Supports multiple values (catenated with ";"). Removes duplicate values.
# 
# USAGE:
# my %hash = read_tohash_uniq(my_file.tsv,keycolname,valcolname)
# 
sub read_tohash_uniq{
	# INPUT
	my $file= shift;
	my $key = shift;
	my $val = shift;
	my $keyi = -1;
	my $vali = -1;
	my $sep = ';';
	
	# OUT
	my %hash= ();
	
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	
	# READ HEADER
	my $ln=0;
	my $l=<IN>;
	chomp($l);
	my @headers= split(/\t/,$l,-1);
	#print STDERR "headers: ",join(";",@headers),"\n";
	for(my $i=0; $i<scalar(@headers); $i++){
		if($headers[$i] eq $key){
			$keyi = $i;
		}
		if($headers[$i] eq $val){
			$vali = $i;
		}
	}
	if($keyi <0){
		die "ERROR: key col=$key not found in file $file\n";
	}
	if($vali <0){
		die "ERROR: value col=$val not found in file $file\n";
	}
	
	while($l=<IN>){
        	$ln++;
		#if($l =~ m/^[@#!]/){ next; }
		chomp($l);
		my @sp= split(/\t/,$l,-1);
		
		if( !defined($hash{$sp[$keyi]}) ){
			my %tmp = ($sp[$vali] => 1);
        	$hash{$sp[$keyi]} = \%tmp;
		}
		else{
			$hash{$sp[$keyi]}->{$sp[$vali]} = 1;
		}
	}
	close(IN);
	
	# convert hash->hash to hash-> str;str;str
	for my $k(keys %hash){
		$hash{$k} 	= join($sep, sort keys %{$hash{$k}});
	}
	
	return %hash;
}

# SOME SIMPLE STATISTICS

# RETURNS LIST MEDIAN
sub median{
	my @values 	= @{shift(@_)};
	my $mid 	= int @values/2;
	@values 	= sort {$a <=> $b} @values;
	my $median;	
	if (@values % 2) {
    		$median = $values[ $mid ];
	}
	else {
    		$median = ($values[$mid-1] + $values[$mid])/2;
	}
	return $median;
}

sub max{
	my @values 	= @{shift(@_)};
	my $max		= (scalar(@values)>0) ? $values[0] : !1;
	foreach my $val(@values){
		if($val > $max){ $max = $val };
	}
	return $max;
}

sub sum{
	my @values 	= @{shift(@_)};
	my $sum		= 0;
	foreach my $val(@values){
		$sum += $val;
	}
	return $sum;
}
# Returns 1-based column index in TSV/CSV-file
# Usage:
# my $taxidi = colind("my_tsv_file.tsv", "taxid")
# my $taxidi = colind("my_csv_file.csv", "taxid",",")
sub colind{

	my $file	= shift();
	my $colname = shift();
	my $sep 	= (scalar(@_)) ? shift(): "\t";
	
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	my $l= <IN>;
	close(IN);
	chomp($l);
	my @sp= split(/$sep/,$l,-1);
	for(my $i=0; $i<scalar(@sp); $i++){
		if($sp[$i] eq $colname){
			return ($i+1);
		}
	}
	return -1;
}

# system_call($system_call, $verbal=false)
#
# $system_call	String with the system call to excecute
# $verbal	Set true to print the call prior to excecuting the call. Optional, default is false.
# 
sub system_call{
	my $call	= shift;
	my $verbal	= (scalar(@_)>0) ? shift(@_): 0;
	if($verbal){
		print STDERR "\t$call\n";
	}
	my @args= ("bash","-c",$call);
	system(@args) == 0 or die $!;
}

sub get_line_count{
	# in:
	my $file	= shift;
	# out:
	my $linen	= 0;
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	while( my $l=<IN> ){
		$linen++;
	}
	close(IN);
	return $linen;
}

