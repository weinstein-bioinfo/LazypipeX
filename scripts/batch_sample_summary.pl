#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Temp qw(tempdir);
use Sort::Naturally;
use Getopt::Long qw(GetOptions);
use List::Util;
use POSIX qw(floor);

# Creates a metasummary for a batch of Lazypipe(NP) runs/results.
# This version summarises assembly/annotation data by sample, using assembly.stats and abund_table.csv data.
# 

my $usage = "\nUSAGE: $0 --maxdir num [-v] resdir 1> summary.csv\n".
			"\n".
			"maxdir  [int]    : limit summary to this many directories [false]\n".
			"resdir  [dir]    : root result directory containing results for each sample\n".
			"toptax  [str]    : Taxrank used for grouping annotations [division]\n".
			"numth   [int]    : Number of threads\n".
			"wrkdir  [dir]    : Root for temporary directories and files [./wrkdir]\n".
			"v                : verbal mode [false]\n".			
			"\n".
			"Collects data from assembly.stats.tsv(s) & abund_table.tsv(s) in resdir/sample-dir(s) and prints metasummary\n\n";
		
# PARAMS
my $maxdir				= !1;
my $toptaxrank			= "division";
my $verbal				= 0;
my $numth				= 8;
my $wrkdir				= './wrkdir';
my $CLEANUP				= 1;
GetOptions(	'maxdir=i'		=> \$maxdir,
			'toptax=s'		=> \$toptaxrank,
			'numth=i'		=> \$numth,
			'wrkdir=s'		=> \$wrkdir,
			'v'         		=> \$verbal);
			
# CHECK OPTIONS
if(scalar @ARGV < 1){ die $usage; }

# START WORKING
system("mkdir -p $wrkdir");
my $tmpdir			= tempdir("batch_summary_XXXXXXXX", DIR => $wrkdir, CLEANUP => $CLEANUP);
#my $abund_tsv_ta		= "$tmpdir/abund_tsv.$target.tmp";
my $resdir			= shift(@ARGV);
my @sampledirs 		= get_dirs($resdir);
@sampledirs 			= sort { ncmp($a,$b) } @sampledirs;


# TMP
my $assembly_stats_tmp	= "$tmpdir/assembly.stats.tsv";
# OUT
my %assembly_stats;
my %abund_tables;

# iterate sample-dirs and parse
#	assembly.stats.tsv(s) to a collection of regular hashes
#	abund_table.tsv(s) to a collection of array-hashes
my $dirc			= 0;
foreach my $dir(@sampledirs){
	$dirc++;
	if($maxdir && ($dirc > $maxdir)){
		last;
	}
	# IN:
	my $sample  				= basename($dir);
	my $assembly_stats		= "$dir/assembly.stats.tsv";
	my $abund_table			= "$dir/abund_table.tsv";
	
	print STDERR "\treading $sample..\n" if($verbal);
	
	# assembly.stats.tsv(s)
	if(-e $assembly_stats){
		system("echo -e 'key\tvalue' 1> $assembly_stats_tmp");
		system("cat $assembly_stats 1>> $assembly_stats_tmp");
		my %tmp	= read_tohash_uniq($assembly_stats_tmp,"key","value");
		$assembly_stats{$sample}		= \%tmp;
	}
	else{
		print STDERR "\t\tWARNING: $sample/assembly.stats.tsv missing\n";
	}
	
	# abund_table.tsv(s)
	if(-e $abund_table){
		my @tmp	= read_csv2arrayhash($abund_table);
		$abund_tables{$sample}	= \@tmp;
	}
	else{
		print STDERR "\t\tWARNING: $sample/abund_table.tsv missing\n";
	}
}

# SUMMARISE COLLECTED DATA
	# Collect unique toptaxranks called "div" here for simplicity
my %div_readn	= ();
foreach my $sample(keys %abund_tables){
	my %tmp	= map{ $_->{$toptaxrank} => $_->{readn}} @{$abund_tables{$sample}};
	foreach my $tax(keys %tmp){
		$div_readn{$tax} = 0 if (!defined($div_readn{$tax}));
		$div_readn{$tax} += $tmp{$tax};
	}
}
	# optionally filter topranks with few reads: TODO
my @divs				= sort {$div_readn{$a} <=> $div_readn{$b}} keys %div_readn;
my @divs_nospaces	= map{ s/ /_/gr } @divs;

print STDERR "\tprinting summary..\n" if($verbal);
# print header
my @headers_stats			= ("reads","reads.flt","reads.hgflt","reads.contigs","contigs","sum","orfs","N50");
my @headers_reads_bydiv		= map{'reads.' . $_} @divs_nospaces;
my @headers_conts_bydiv		= map{'contigs.' . $_} @divs_nospaces;
my @headers_assembly_bydiv	= map{'assembly.' . $_} @divs_nospaces;
print join(",","Sample",@headers_stats,@headers_reads_bydiv,@headers_conts_bydiv,'reads.mapped','reads.unmapped','contigs.mapped','contigs.unmapped','assembly.mapped.bp','assembly.unmapped.bp'),"\n";
# print data
foreach my $sample(sort keys %assembly_stats){
	# fields from stats
	my %stats		= %{$assembly_stats{$sample}};
	my @vals_stats	= @stats{@headers_stats};
	# fields from abund_table
	my @abund_table	= @{$abund_tables{$sample}};
	my %div2reads	= map{ $_ => 0} @divs;
	my %div2contigs	= map{ $_ => 0} @divs;
	my %div2assembly= map{ $_ => 0} @divs;
	
	foreach my $rec(@abund_table){
		my $div				= $rec->{$toptaxrank};
		$div2reads{$div}		= 0 if (!defined($div2reads{$div}));
		$div2contigs{$div}	= 0 if (!defined($div2contigs{$div}));
		$div2reads{$div}		+= $rec->{readn};
		$div2contigs{$div}	+= $rec->{contign};
		$div2assembly{$div}	+= $rec->{'assembly.size'} || 0;
	}
	# round read-numbers
	$div2reads{$_} 		= floor($div2reads{$_} + 0.5) for keys %div2reads;
	$div2assembly{$_}	= floor($div2assembly{$_} + 0.5) for keys %div2assembly;
	
	# Mapped & unmapped reads
	my $reads_mapped	 	= 0; $reads_mapped += $_ for values %div2reads;
	my $reads_unmapped	= $stats{'reads.contigs'} - $reads_mapped;
	# mapped & unmapped contigs
	my $contigs_mapped	= 0; $contigs_mapped += $_ for values %div2contigs;
	my $contigs_unmapped= $stats{'contigs'} - $contigs_mapped;
	# mapped/unmapped assemlby
	my $assembly_mapped	= 0; $assembly_mapped += $_ for values %div2assembly;
	my $assembly_unmapped= $stats{sum} - $assembly_mapped;
	
	# print values
	print join(",",$sample,@vals_stats,@div2reads{@divs},@div2contigs{@divs},$reads_mapped,$reads_unmapped,$contigs_mapped,$contigs_unmapped,$assembly_mapped,$assembly_unmapped),"\n";
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

# Read csv/tsv file to a table represented by array-of-hashes structure.
# Will guess delimter from the first line.
#
# USAGE:
# my @table = read_csv2arrayhash("data.tsv")
# or
# my @table = read_csv2arrayhash("data.csv")
#
# say $table[10]->{$colname};
# 
sub read_csv2arrayhash{
	my $subid	= "read_csv2arrayhash()";
	
	# INPUT
	my $file		= shift;
	
	# OUT
	my @table	= ();
	
	open(IN,"<$file") or die "ERROR: $subid: Can\'t open $file: $!\n";
	
	# GUESS DELIM
	my $del		= ',';
	my $l		= <IN>;
	my $ln		= 1;
	chomp($l);
	$l =~ s/\r|\n$//;
	my @csv_sp	= split(',',$l,-1);
	my @tsv_sp	= split('\t',$l,-1);
	if(scalar(@csv_sp) >= scalar(@tsv_sp)){
		$del		= ',';
	}
	else{
		$del		= '\t';
	}
	
	# READ HEADER
	my @headers= split(/$del/,$l,-1);
	while($l=<IN>){
        	$ln++;
		chomp($l);
		$l =~ s/\r|\n$//;
		my @sp= split(/$del/,$l,-1);
		if(scalar(@sp)<scalar(@headers)){
			print $STDERR "WARNING: $subid: missing columns at line $ln\n: skipping";
			next;
		}
		my %row = ();
		for(my $i=0; $i<scalar(@sp) && $i<scalar(@headers); $i++){
			$row{$headers[$i]} = $sp[$i];
		}
		push(@table,\%row);
	}
	close(IN);

	return @table;
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

