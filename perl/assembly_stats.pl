#!/usr/bin/perl
use strict;
use warnings;

# 
# Prints stats for an assembly
my $usage= "\nUSAGE: $0 contigs.fa metrics format > stats.tsv\n".
	   "\n".
	   "contigs.fa  \t: fasta file with assembly contig sequences\n".
	   "metrics [str]\t: comma-separed list of metrics contnum,min,max,mean,bp,Nxx,Lbpxx,LNxx\n".
	   "\tcontigs \t: number of contigs\n".
	   "\tmin     \t: min contig length\n".
	   "\tmax     \t: max contig length\n".
	   "\tmean    \t: mean contig length\n".
	   "\tsum     \t: total contig length\n".
	   "\tNxx     \t: N50, N30 or any other Nxx metric\n".
	   "\tLNxx    \t: number of contigs >= xx, i.e. LN1000 returns number of contigs >=1000\n".
	   "\tLbpxx   \t: total length of contigs >= xx\n".
	   "format [str]\t: col|row[,names]\n".
	   "\tcol     \t: metrics printed in a column, ordered as in \'metrics\' param\n".
	   "\trow     \t: metrics printed in a row, ordered as in \'metrics\' param\n".
	   "\tnames   \t: will print metric names along with the values\n".
	   "\n";
if(scalar(@ARGV)<3) { 
	print "$usage";
	exit(1);
}


my $contigs_fasta	= shift(@ARGV);
my @metrics			= split(/,/,shift(@ARGV));
my $format			= shift(@ARGV);

# chech params
for(my $i=0; $i<scalar(@metrics);$i++){
	if(	$metrics[$i] =~ m/^(contigs)$/i){ 	$metrics[$i]= 'contigs';}
	elsif( $metrics[$i] =~ m/(min)/i ){ 		$metrics[$i]= 'min';}
	elsif( $metrics[$i] =~ m/(max)/i ){ 		$metrics[$i]= 'max';}
	elsif( $metrics[$i] =~ m/(mean)/i ){ 		$metrics[$i]= 'mean';}
	elsif( $metrics[$i] =~ m/(sum)/i ){ 		$metrics[$i]= 'sum';}
	elsif( $metrics[$i] =~ m/^N([0-9]{2})$/g ){	$metrics[$i]= "N,$1";}
	elsif( $metrics[$i] =~ m/^LN([0-9]+)$/g ){	$metrics[$i]= "LN,$1";}
	elsif( $metrics[$i] =~ m/^Lbp([0-9]+)$/g ){	$metrics[$i]= "Lbp,$1";}
	else{
		die "ERROR: unknown metric: $metrics[$i]\n";
	}
}
if( !($format =~m/(row|col)/ig) ){
	die "ERROR: invalid format: $format\n";
}

my %contigs 			= %{readfasta($contigs_fasta)};
my @contig_len			= ();

foreach my $seqid( keys %contigs){
	push(@contig_len, length($contigs{$seqid}));
}
#@contig_len 			= sort {$b <=> $a} @contig_len;
my $contig_len_sorted 	= 0;


# CALCULATE METRICS
my %metric_values	= ();
foreach my $m(@metrics){

	my @pair= split(/,/,$m);
	
	if($pair[0] eq 'contigs'){
		$metric_values{$m}	= scalar(@contig_len);
	}
	if($pair[0] eq 'min'){
		$metric_values{$m}	= min(\@contig_len);
	}
	if($pair[0] eq 'max'){
		$metric_values{$m}	= max(\@contig_len);
	}
	if($pair[0] eq 'mean'){
		$metric_values{$m}	= mean(\@contig_len);
	}
	if($pair[0] eq 'sum'){
		$metric_values{$m}	= sum(\@contig_len);
	}
	if($pair[0] eq 'N'){
		if(!$contig_len_sorted){
			@contig_len 		= sort {$b <=> $a} @contig_len;
			$contig_len_sorted  = 1;
		}
		$metric_values{$m} 	= get_NXX(\@contig_len,$pair[1]);
	}
	if($pair[0] eq 'LN'){
		$metric_values{$m}	= get_LN(\@contig_len,$pair[1]);
	}
	if($pair[0] eq 'Lbp'){
		$metric_values{$m}	= get_Lbp(\@contig_len,$pair[1]);
	}			
}

# PRINT RESULTS
if( ($format =~ m/(col)/i) ){
	if(  ($format=~m/(names)/i)  ){
		foreach my $metric( @metrics ){
			my $name	= $metric;
			$name 		=~ s/,//g;
			print "$name\t$metric_values{$metric}\n";
		}
	}
	else{
		foreach my $metric( @metrics ){
			print "$metric_values{$metric}\n";
		}
	}
}
else{
	if($format=~m/(names)/i){
		my @name_list = ();
		foreach my $metric( @metrics ){
			my $name	= $metric;
			$name 		=~ s/,//g;
			push(@name_list,$name);
		}
		print "",join("\t",@name_list),"\n";
	}
	my @values = ();
	foreach my $metric( @metrics){
		push(@values,$metric_values{$metric});
	}
	print "",join("\t",@values),"\n";
}


sub readfasta{
  	my $file		= shift(@_);
	my %sequence;
	my $header;
	my $temp_seq;
	 
	open (IN, "<$file") or die "couldn't open the file $file $!";
	
	while (<IN>){	
		chop;
		next if /^\s*$/; #skip empty line 
		if ($_ =~ s/^>//)  #when see head line
		{	
			$header= $_;
			if ($sequence{$header}){print colored("#CAUTION: SAME FASTA HAS BEEN READ MULTIPLE TIMES.\n#CAUTION: PLEASE CHECK FASTA SEQUENCE:$header\n","red")};
			if ($temp_seq) {$temp_seq=""} # If there is alreay sequence in temp_seq, empty the sequence file
			
		}
		else # when see the sequence line 
		{
		   s/\s+//g;
		   $temp_seq .= $_;
		   $sequence{$header}=$temp_seq; #update the contents
		}
	
	}
	
	return \%sequence;
}


# Returns N50, N90 or NX statistic
# NOTE: the array must be sorted 
sub get_NXX{
	my @array	= @{shift(@_)};
	my $N		= shift(@_);
	my $sum		= sum(\@array);
	my $sumN	= $sum*$N/100.0;
	
	my $sumi= 0;
	for(my $i=0; $i<scalar(@array); $i++){
		$sumi+= $array[$i];
		if($sumi>$sumN){
			return $array[$i];
		}
	}
}

#
# Returns LNbp statistic: number of baisepairs in contigs >= L
#
# get_Lbp(contig_len,L)
# contig_len [list]	: array of contig lengths
# L [num]		: min contig length to count
#
# example: get_Lbp(contig_len,1000) return number of bps in contigs >= 1000 bp
#
sub get_Lbp{
	my @array	= @{shift(@_)};
	my $L		= shift(@_);
	my $Lbp	= 0;
	for(my $i=0; $i<scalar(@array); $i++){
		if($array[$i] >= $L){
		$Lbp += $array[$i];}
	}
	return $Lbp;
}
#
# Returns LN statistic: number of contigs >= L
#
sub get_LN{
	my @array	= @{shift(@_)};
	my $L		= shift(@_);
	my $LN		= 0;
	for(my $i=0; $i<scalar(@array); $i++){
		if($array[$i] >= $L){
		$LN++;}
	}
	return $LN;
}


## BASIC SET OPERATIONS
sub sum{
	my @array	= @{shift(@_)};
	if(scalar(@array)<1){
		return undef;
	}
	my $sum		= 0;
	foreach(@array){
		$sum += $_;
	}
	return $sum;
}
sub mean{
	my @array	= @{shift(@_)};
	if(scalar(@array)<1){
		return undef;
	}
	return sum(\@array)/scalar(@array);
}

sub min{
	my @array	= @{shift(@_)};
	if(scalar(@array)<1){
		return undef;
	}
	my $min = $array[0];
	foreach(@array){
		if($_ < $min){
			$min = $_;}
	}
	return $min;
}
sub max{
	my @array	= @{shift(@_)};
	if(scalar(@array)<1){
		return undef;
	}
	my $max = $array[0];
	foreach(@array){
		if($_ > $max){
			$max = $_;}
	}
	return $max;
}
