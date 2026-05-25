package Lazypipe::Utils;
use strict;
use warnings;
use Exporter;
use MIME::Base64;
use File::Basename;
our @ISA			= qw( Exporter );
our @EXPORT		= qw( read_tsv2hash read_tsv2hash_noheaders read_tsv2kvahash read_tsv2kuvahash read_tsv2hashtable read_tsv2arraytable read_file2array write_array2file write_hash2file write_table_html ncol nlines colind mcolind median max sum system_call );
our @EXPORT_OK	= qw( format_int filebin2uri );
our $VERBAL		= !1;

###
# UTIL FUNCTIONS FOR LAZYPIPE PROJECT
#
# Credit:
#
# Plyusnin, I., Vapalahti, O., Sironen, T., Kant, R., & Smura, T. (2023).
# Enhanced Viral Metagenomics with Lazypipe 2. Viruses, 15(2), 431.
#
# Plyusnin,I., Kant,R., Jaaskelainen,A.J., Sironen,T., Holm,L., Vapalahti,O. and Smura,T. (2020) 
# Novel NGS Pipeline for Virus Discovery from a Wide Spectrum of Hosts and Sample Types. Virus Evolution, veaa091
#
# Contact: grp-lazypipe@helsinki.fi
###



# Read tsv file to a simple key-value hash. 
# For each key returns the last value encountered.
# 
# USAGE:
# my %hash = read_tsv2hash(my_file.tsv,keycolname,valcolname)
# 
sub read_tsv2hash{
	# INPUT
	my $file= shift;
	my $key = shift;
	my $val = shift;
	my $keyi = -1;
	my $vali = -1;
	
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
	if($keyi <0 || $keyi>=(scalar @headers)){
		die "ERROR: key col=$key not found in file $file\n";
	}
	if($vali <0 || $vali>=(scalar @headers)){
		die "ERROR: value col=$val not found in file $file\n";
	}
	
	while($l=<IN>){
        	$ln++;
		#if($l =~ m/^[@#!]/){ next; }
		chomp($l);
		my @sp= split(/\t/,$l,-1);
		
		$hash{$sp[$keyi]} = $sp[$vali];
	}
	close(IN);
	return %hash;
}

# Read tsv file to a simple key-value hash. 
# For each key returns the last value encountered.
# 
# USAGE:
# my %hash = read_tsv2hash_noheaders(my_file.tsv,keycol,valcol)
# 
sub read_tsv2hash_noheaders{
	# INPUT
	my $file= shift;
	my $keyi = shift;
	my $vali = shift;
	
	# OUT
	my %hash= ();
	
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	
	# READ TABLE
	my $ln=0;
	while(my $l=<IN>){
        	$ln++;
		#if($l =~ m/^[@#!]/){ next; }
		chomp($l);
		my @sp= split(/\t/,$l,-1);
		
		$hash{$sp[$keyi]} = $sp[$vali];
	}
	close(IN);
	return %hash;
}

# Read tsv file to a key-value-array hash. 
# For each key returns an array of all associated values. These values are not unique.
# 
# USAGE:
# my %hash = read_tsv2kvahash(my_file.tsv,keycolname,valcolname)
# 
sub read_tsv2kvahash{
	# INPUT
	my $file= shift;
	my $key = shift;
	my $val = shift;
	my $keyi = -1;
	my $vali = -1;
	
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
			my @tmp	= ($sp[$vali]);
        		$hash{$sp[$keyi]} = \@tmp;
		}
		else{
			push(@{$hash{$sp[$keyi]}}, $sp[$vali]);
		}
	}
	close(IN);
	return %hash;
}


# Read tsv file to a key-unique-value-array hash.
# For each key returns an array of unique associated values. Each value in the value-array is unique.
# 
# USAGE:
# my %hash = read_tsv2kuvahash(my_file.tsv,keycolname,valcolname)
# 
sub read_tsv2kuvahash{
	# INPUT
	my $file= shift;
	my $key = shift;
	my $val = shift;
	my $keyi = -1;
	my $vali = -1;
	
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
	
	# convert hash->\%hash to hash->\@array
	for my $k(keys %hash){
		my @tmp		= sort keys %{$hash{$k}};
		$hash{$k} 	= \@tmp;
	}
	return %hash;
}

# Read tsv file to a hash of rows indexed by keys.
#
# USAGE:
# my %table = read_tsv2hashtable(data.tsv,keycolname)
# $table{$rowkey}->{$colname};
# 
sub read_tsv2hashtable{
	# INPUT
	my $file= shift;
	my $key = shift;
	my $keyi = -1;
	
	# OUT
	my %hash= ();
	
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	
	# READ HEADER
	my $ln=0;
	my $l=<IN>;
	chomp($l);
	my @headers= split(/\t/,$l,-1);
	for(my $i=0; $i<scalar(@headers); $i++){
		if($headers[$i] eq $key){
			$keyi = $i;
		}
	}
	if($keyi <0){
		die "ERROR: key col=$key not found in file $file\n";
	}
	
	while($l=<IN>){
        	$ln++;
		chomp($l);
		my @sp= split(/\t/,$l,-1);
		if(scalar(@sp)<scalar(@headers)){
			print $STDERR "WARNING: read_tohashtable: missing columns at line $ln\n: skipping";
			next;
		}
		my %row = ();
		for(my $i=0; $i<scalar(@sp) && $i<scalar(@headers); $i++){
			$row{$headers[$i]} = $sp[$i];
		}
		$hash{$sp[$keyi]} = \%row;
	}
	close(IN);
	
	return %hash;
}

# Reads lines in a text file to an array.
#
# USAGE:
# my @array = read_file2array($file)
#
sub read_file2array{
	# INPUT
	my $file= shift;
	# OUT
	my @array= ();
	
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	while(my $l=<IN>){
        chomp($l);
		push(@array,$l);
	}
	close(IN);
	return @array;
}

# Read tsv file to an array of row arrays.
# Checks that all rows have the same number of columns as the first row
#
# USAGE:
# my @table = read_tsv2arraytable('mydata.tsv');
# # iterate and print all rows:
# my rowi = 1;
# foreach my $rowp(@table){
#	say "$rowi  :",paste(',', @{$rowp});	
# }
# # value on row 3, col 3
# say "[3,3]: ",$table[3]->[3];
# 
sub read_tsv2arraytable{
	my $subid	= "read_tsv2arraytable()";
	# INPUT
	my $file		= shift;
	
	# OUT
	my @table	= ();
	
	open(my $IN,"<$file") or die "$subid: Can\'t open $file: $!\n";
	my ($l,$ln,$coln) = 0;
	
	# First row
	if(defined($IN) && ($l=<$IN>)){
		$ln++;
		chomp($l);
		my @sp	= split(/\t/,$l,-1);
		$coln	= scalar(@sp);
		push(@table,\@sp);
	}	
	
	while($l=<$IN>){
        	$ln++;
		chomp($l);
		my @sp	= split(/\t/,$l,-1);
		if(scalar(@sp) != $coln){
			print STDERR "WARNING: $subid: invalid column number on line $ln: skipping\n";
			next;
		}
		push(@table, \@sp);
	}
	close($IN);
	
	return @table;
}

# Writes array to a text file.
#
# USAGE:
# write_array2file($file, \@array)
#
sub write_array2file{
	my $file 	= shift;
	my @array	= @{shift(@_)};
	open(OUT,">$file") or die "Can\'t write to $file: $!\n";
	foreach my $el(@array){
		print OUT "$el\n";
	}
	close(OUT);	
}

# Writes hash to a text file. 
#
# USAGE:
# write_hash2file(file,\%hash,sep="\t")
#
sub write_hash2file{
	my $file		= shift;
	my $hashp	= shift;
	my $sep		= scalar(@_)>0 ? shift(@_): "\t";
	open(OUT,">$file") or die "Can\'t open $file: $!\n";
	foreach my $k(sort keys %{$hashp}){
		print OUT "$k",$sep,$hashp->{$k},"\n";
	}
	close(OUT);
}

# Prints table (array of array-refs) to html-table
# USAGE:
# write_table_html(table=> \@summary_table, table_attrs=>$table_attrs, th_attrs=>\@hd_attrs, fh=>*OUT);
sub write_table_html{
	my (%args)			= @_;
	my @table	 		= @{$args{table}};
	my $table_attrs		= defined($args{table_attrs}) ? $args{table_attrs} : "";
	my @th_attrs			= @{$args{th_attrs}};
	my $fh				= $args{fh};
	
	print $fh "<table $table_attrs >\n";
	
	# print table header
	print $fh "\t<thead><tr>";
	my @row			= @{shift(@table)};
	for(my $coli = 0; $coli<scalar(@row); $coli++){
		# add <span>-element needed by sorted-tables (placeholder for arrow-tag)
		print $fh "<th $th_attrs[$coli]>$row[$coli] &emsp;<span></span></th>\n";
	}
	print $fh "\t</tr></thead>\n";

	# print tbody
	print $fh "\t<tbody>\n";
	foreach my $rowp(@table){
		print $fh "\t<tr>";
		my @row = @{$rowp};
		foreach my $col(@row){
			print $fh "<td>$col</td>";
		}
		print $fh "</tr>\n";
	}
	print $fh "</tbody>\n";
	
	print $fh "</table>\n";
}


# SOME SIMPLE STATISTICS

# USAGE: median(\@values)
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
# USAGE: max(\@values)
sub max{
	my @values 	= @{shift(@_)};
	my $max		= (scalar(@values)>0) ? $values[0] : !1;
	foreach my $val(@values){
		if($val > $max){ $max = $val };
	}
	return $max;
}
# USAGE: sum(\@values)
sub sum{
	my @values 	= @{shift(@_)};
	my $sum		= 0;
	foreach my $val(@values){
		$sum += $val;
	}
	return $sum;
}

# Adds thousand separator to int
# Usage: format_int(1234,' ')
sub format_int{
	my $a 	= shift;
	my $sep	= shift;
	my $b 	= reverse $a;
	my @c 	= unpack("(A3)*", $b);
	my $d 	= join $sep, @c;
	my $e	= reverse($d);
	return($e);
}

# Encode file as a binary-base64 BLOB competible with IGV. File can be ASCII/UTF-8/binary.
# 
# From IGV documentation data-uri format:
# 	data:application/gzip;base64,<uuencoded string>
# DEPENDENCIES: igv_reports/create_datauri python script
#
# USAGE: file2uri_igvreports($file)
#
sub file2uri_igvreports{
	my $data_file	= shift;
	my $uri 	= `create_datauri $data_file`;
	chomp($uri);
	return $uri;
}

# Encode base64-encoded uri compatible with IGV. Assumes input file is binary.
# 
# From IGV documentation data-uri format:
# 	data:application/gzip;base64,<uuencoded string>
#
# USAGE:	
# my $uri	= filebin2uri($file_binary)
# 
sub filebin2uri{
	my $data_file	= shift;
	open my $in, '<', $data_file or die;
	binmode $in;
	my $data 	= '';
	my $length	= 1000;
	while (1) {
    		my $success = read $in, $data, $length, length($data);
    		die $! if not defined $success;
    		last if not $success;
	}
	close $in;

	#my $encoded	= encode_base64url($data,'');
	my $encoded	= encode_base64($data,'');
	return "data:application/gzip;base64,". $encoded;
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
# Returns 1-based column index in TSV/CSV-file that matches a given regular expression.
# Returns -1 if no match was found.
# Usage:
# my $score_col = colind("dbhits.tsv","AS:i:[0-9]+")
# 
sub mcolind{
	my $file	= shift();
	my $regex	= shift();
	my $sep 	= (scalar(@_)) ? shift(): "\t";
	
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	my $l= <IN>;
	close(IN);	
	chomp($l);
	my @sp= split(/$sep/,$l,-1);
	for(my $i=0; $i<scalar(@sp); $i++){
		if($sp[$i] =~ m/$regex/g ){
			return ($i+1);
		}
	}
	return -1;
}
# Returns column number in TSV/CSV-file
# Usage:
# my $ncols = ncol("dbhits.tsv")
# my $ncols = ncol("dbhits.csv",",");
sub ncol{
	my $file	= shift();
	my $sep 	= (scalar(@_)) ? shift(): "\t";
	
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	my $l= <IN>;
	close(IN);	
	chomp($l);
	my @sp= split(/$sep/,$l,-1);
	return scalar(@sp);
}
# Returns number of lines in a file
sub nlines{
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

# system_call($system_call, $verbal)
#
# $system_call	String with the system call to excecute
# $v	erbal		Defaults to $NGSlib::VERBAL
# 
sub system_call{
	my $call		= shift;
	my $verbal	= (scalar(@_)>0) ? shift(@_): ( defined($VERBAL) ? $VERBAL : 0 );
	if($verbal){
		print STDERR "\t$call\n";
	}
	my @args= ("bash","-c",$call);
	system(@args) == 0 or die $!;
}

1;


