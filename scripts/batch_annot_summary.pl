#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Temp  qw(tempdir);
use Sort::Naturally;
use Getopt::Long qw(GetOptions);
use List::Util qw(uniqstr);


# PROJECT	: Lazypipe v3 Manuscript
# Date		: 2024 5. Sep
# Desc:		Creates a metasummary for a batch of Lazypipe runs/results using annot_table.tsv
# 			Aiming for format: Sample <\t>Contid.id <\t>Contig.length <\t>Search <\t>Species <\t>Genus <\t>Family <\t>Accession <\t>qcov <\t>pident <\t>scov <\n>

my $usage = "\nUSAGE: $0 --maxdir num --ta|target ba|ph|vi --minbp num --bits num [-v] resdir 1> metasummary.csv\n".
			"\n".
			"maxdir  [int]    : limit summary to this many directories [false]\n".
			"target  [str]    : print metasummary for Bacteria (\"ba\"), Phages (\"ph\") or Viruses (\"vi\") [vi]\n".
			"minbp   [int]    : min contig length in basepairs [1000]\n".
			"bits    [int]    : min bitscore [0]\n".
			"qcov    [num]    : min qcov, number in range [0,1] [0]\n".
			"resdir  [dir]    : root result directory containing results for each sample\n".
			"toptax  [str]    : Taxrank used for selecting Bacteria/Viruses/Phages [division]\n".
			"tophit           : List only tophits for each contig and search\n".
			"                   By default will create ';'-sperated list of all valid hits for each contig\n".
			"numth   [int]    : Number of threads\n".
			"wrkdir  [dir]    : Root for temporary directories and files [./wrkdir]\n".
			"v                : verbal mode [false]\n".
			"\n".
			"Collects bacteria/phage/virus hits from annot_table.tsv(s) in res-dir/sample-dir(s) and prints metasummary\n\n";
		
# PARAMS
my $maxdir				= !1;
my $target				= "vi";
my $minbp				= 1000;
my $bitscore				= 0;
my $qcov					= 0;
my $toptaxrank			= "division";
my $verbal				= 0;
my $tophit				= 0;
my $numth				= 8;
my $wrkdir				= "./wrkdir";

GetOptions(	'maxdir=i'		=> \$maxdir,
			'ta|target=s'	=> \$target,
			'minbp=i'		=> \$minbp,
			'bits=i'			=> \$bitscore,
			'qcov=f'			=> \$qcov,
			'toptax=s'		=> \$toptaxrank,
			'tophit'			=> \$tophit,
			'numth=i'		=> \$numth,
			'wrkdir=s'		=> \$wrkdir,
			'v'         		=> \$verbal);
		
if(scalar @ARGV < 1){ die $usage; }
if( !($target eq 'ba' || $target eq 'ph' || $target eq 'vi') ){
	die "ERROR: unknown option: --target $target\n";
}

# CREATE tmpdir
system("mkdir -p $wrkdir");
my $tmpdir		= tempdir("batch_summary_XXXXXXXX", DIR => $wrkdir, CLEANUP => 1);
my $resdir		= shift(@ARGV);
my @sampledirs 	= get_dirs($resdir);
@sampledirs 		= sort { ncmp($a,$b) } @sampledirs; 
my @meta_table;
my $dirc			= 0;

# print header
print join(",","Sample","Contig.id","Contig.length","Search","Species","Genus","Family","Order","Class","Accession","Qcov.top","Qcov.max","Scov.max","Bits.max","Pident.max","host.source","genome.composition"),"\n";

# iterate sample-dirs
foreach my $dir(@sampledirs){
	$dirc++;
	if($maxdir && ($dirc > $maxdir)){
		last;
	}
	# in:
	my $sample  			= basename($dir);
	my $annot_tsv		= "$dir/annot_table.tsv";
	# tmp:
	my $annot_tsv_ta			= "$tmpdir/annot_tsv.$target.tmp";

	print STDERR "\treading $sample..\n";
	if($verbal){
		print STDERR "\treading $annot_tsv..\n";}
	
	if( !(-e $annot_tsv) || get_line_count($annot_tsv) < 2){
		print STDERR "\t\tannotation file is missing or empty: skipping $dir\n";
		next;
	}	
	
	if($target eq 'ba' ){
		system_call("csvtk filter2 -tj $numth -f '\$$toptaxrank==\"Bacteria\" && \$clen>=$minbp && \$bitscore>=$bitscore && \$qcov>=$qcov' $annot_tsv 1> $annot_tsv_ta", $verbal);
	}
	elsif($target eq 'ph'){
		system_call("csvtk filter2 -tj $numth -f '\$$toptaxrank==\"Phages\" && \$clen>=$minbp && \$bitscore>=$bitscore && \$qcov>=$qcov' $annot_tsv 1> $annot_tsv_ta", $verbal);
	}
	else{
		system_call("csvtk filter2 -tj $numth -f '\$$toptaxrank==\"Viruses\" && \$clen>=$minbp && \$bitscore>=$bitscore && \$qcov>=$qcov' $annot_tsv 1> $annot_tsv_ta", $verbal);
	}
	
	# check there are valid hits for the target
	if( get_line_count($annot_tsv_ta) < 2){
		print STDERR "\t\tno hits for $target: skipping $dir\n";
		next;
	}
	
	# get the list of search-values that have valid contig annotations
	my %search_list = read_tohash($annot_tsv_ta,"search","search");
	
	foreach my $search (sort keys %search_list){
		
		my $annot_tsv_ta_search	= "$tmpdir/annot_tsv.$target.$search.tmp";
		
		system_call("csvtk filter2 -tj $numth -f '\$search == \"$search\"' $annot_tsv_ta 1> $annot_tsv_ta_search", $verbal);
		
		if($tophit){
			system_call("cat $annot_tsv_ta_search | ".
				"csvtk sort -tj $numth -k contig:N -k bitscore:nr  | ".
				"csvtk uniq -tj $numth -f contig 1> $annot_tsv_ta_search.tmp",$verbal );
			system("mv $annot_tsv_ta_search.tmp $annot_tsv_ta_search");
		}
		
		# Old code: parse data in multiple passes
		#my %cont2clen	= read_tohash_uniq( $annot_tsv_ta_search, "contig","clen");
		#my %cont2sp		= read_tohash_uniq( $annot_tsv_ta_search, "contig","species");
		#my %cont2ge		= read_tohash_uniq( $annot_tsv_ta_search, "contig","genus");
		#my %cont2fa		= read_tohash_uniq( $annot_tsv_ta_search, "contig","family");
		#my %cont2acc		= read_tohash_uniq( $annot_tsv_ta_search, "contig","sseqid");
		#my %cont2qcov	= read_tohash( $annot_tsv_ta_search, "contig","qcov");
		#my %cont2pident	= read_tohash( $annot_tsv_ta_search, "contig","pident");
		#my %cont2scov	= read_tohash( $annot_tsv_ta_search, "contig","scov");
		#my %cont2bits	= read_tohash( $annot_tsv_ta_search, "contig","bitscore");
		#my %cont2host	= read_tohash_uniq( $annot_tsv_ta_search, "contig","host.source");
		#my %cont2comp	= read_tohash_uniq( $annot_tsv_ta_search, "contig","genome.composition");
		
		# New code:
		#	- read CSV to an array-hash structure in a single pass
		#	- restructure array-hash to a collection of contigt_to_value-list hashes
		my @annot		= read_csv2arrayhash($annot_tsv_ta_search);
		my (%cont2clen, %cont2sp, %cont2ge, %cont2fa, %cont2or, %cont2cl,%cont2acc, %cont2qcov, %cont2pident, %cont2scov, %cont2bits, %cont2host, %cont2comp);
		foreach my $rec(@annot){
			my $cont				= $rec->{contig} || die "ERROR: missing field in $annot_tsv_ta_search: contig\n";
			# init arrays on the go
			$cont2sp{$cont}		= [] if(!defined($cont2sp{$cont}));
			$cont2ge{$cont}		= [] if(!defined($cont2ge{$cont}));
			$cont2fa{$cont}		= [] if(!defined($cont2fa{$cont}));
			$cont2or{$cont}		= [] if(!defined($cont2or{$cont}));
			$cont2cl{$cont}		= [] if(!defined($cont2cl{$cont}));
			$cont2acc{$cont}		= [] if(!defined($cont2acc{$cont}));
			$cont2qcov{$cont}	= [] if(!defined($cont2qcov{$cont}));
			$cont2pident{$cont}	= [] if(!defined($cont2pident{$cont}));
			$cont2scov{$cont}	= [] if(!defined($cont2scov{$cont}));
			$cont2bits{$cont}	= [] if(!defined($cont2bits{$cont}));
			$cont2host{$cont}	= [] if(!defined($cont2host{$cont}));
			$cont2comp{$cont}	= [] if(!defined($cont2comp{$cont}));
			
			# check rec
			my @ks_must_have 	= ("species","genus","family","sseqid","qcov","pident","scov","bitscore");
			foreach my $k(@ks_must_have){
				die "ERROR: missing field: $k\n" if(!defined($rec->{$k}));
			}
			
			# fill the arrays
			$cont2clen{$cont}	= $rec->{clen} || 'NA';
			push(@{$cont2sp{$cont}}, $rec->{species});
			push(@{$cont2ge{$cont}}, $rec->{genus});
			push(@{$cont2fa{$cont}}, $rec->{family});
			push(@{$cont2or{$cont}}, $rec->{order} || 'NA');
			push(@{$cont2cl{$cont}}, $rec->{class} || 'NA');
			push(@{$cont2acc{$cont}}, $rec->{sseqid});
			push(@{$cont2qcov{$cont}}, $rec->{qcov});
			push(@{$cont2pident{$cont}}, $rec->{pident});
			push(@{$cont2scov{$cont}}, $rec->{scov});
			push(@{$cont2bits{$cont}}, $rec->{bitscore});
			push(@{$cont2host{$cont}}, $rec->{"host.source"} || 'NA');
			push(@{$cont2comp{$cont}}, $rec->{"genome.composition"} || 'NA');
		}
		
		# print
		for my $cont(sort{$cont2clen{$b} <=> $cont2clen{$a}} keys %cont2clen ){
			my $clen			= $cont2clen{$cont};
			my $sp			= '"'. join(';', sort {$a cmp $b} uniqstr @{$cont2sp{$cont}}) .'"';
			my $ge			= '"'. join(';', sort {$a cmp $b} uniqstr @{$cont2ge{$cont}}) .'"';
			my $fa			= '"'. join(';', sort {$a cmp $b} uniqstr @{$cont2fa{$cont}}) .'"';
			my $or			= '"'. join(';', sort {$a cmp $b} uniqstr @{$cont2or{$cont}}) .'"';
			my $cl			= '"'. join(';', sort {$a cmp $b} uniqstr @{$cont2cl{$cont}}) .'"';
			my $acc			= '"'. join(';', sort {$a cmp $b} uniqstr @{$cont2acc{$cont}}) .'"';
			my @qcov			= @{$cont2qcov{$cont}};
			my $qcov			= $qcov[0];
			my $qcov_max		= max(\@qcov);
			my @pident		= @{$cont2pident{$cont}};
			my $pident_max	= max(\@pident);
			my @scov			= @{$cont2scov{$cont}};
			my $scov_max		= max(\@scov);
			my @bits			= @{$cont2bits{$cont}};
			my $bits_max		= max(\@bits);
			
			my $host			= '"'. join(';', sort {$a cmp $b} uniqstr @{$cont2host{$cont}}) .'"';
			my $comp			= '"'. join(';', sort {$a cmp $b} uniqstr @{$cont2comp{$cont}}) .'"';
			
			print join(",",$sample, $cont, $clen, $search, $sp, $ge, $fa, $or, $cl, $acc, $qcov,$qcov_max,$scov_max,$bits_max,$pident_max, $host, $comp),"\n";	
		}	
	}
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
		if($val eq 'NA'){ next; }
		if($val > $max){ $max = $val };
	}
	return $max;
}

sub sum{
	my @values 	= @{shift(@_)};
	my $sum		= 0;
	foreach my $val(@values){
		if($val eq 'NA'){ next; }
		$sum += $val;
	}
	return $sum;
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

