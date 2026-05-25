#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd;


# Converts abundance table to text file that can be converted to Krona graph with ktImportText
#
# input: abundance table
#	Table MUST be in tab-deleted format with a header row.
#	Table MUST include the following headers, any order is accepted, other headers are ignored:
#		taxid
#		readn
#		species
#		genus
#		family
#		order
#		class
#		phylum
#		$toptaxrank
#	--sample str	: name of the biological sample (optional)
#	--tail num	: number in percents [0,100], taxa with least abundance summing to this amount will be ignored
#   --toptaxrank str: top taxonomy rank to use [division]
#
# output: taxonomy profile in CAMI format (see https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd)


my $usage	= "USAGE: $0 abund_table [-s|sample id -t|tail num] --toptaxrank division 1> krona.txt\n".
		"\n".
		"Converts abundance table to Krona text format.\n".
		"Krona text can be converted to Krona graph with ktImportText\n\n";
	
my $help			= !1;
my $sampleid		= cwd();
my $tail			= 0;
my $toptaxrank	= 'division';

GetOptions(	'help|h'			=> \$help, 
			'sample|s=s'		=>\$sampleid, 
			'toptaxrank=s'	=> \$toptaxrank,
			'tail|t=f'		=>\$tail) or die $usage;
if((scalar @ARGV < 1) || $help){ 
	die $usage;
}


# ABUND TABLE <HEADER>
my @headers_taxnames		= ($toptaxrank,'phylum','class','order','family','genus','species');
my @headers_required		= (@headers_taxnames,'readn','taxid');
my %headeri				= ();
my $file 				= shift(@ARGV);
open(IN,"<$file") or die "Can\'t open $file: $!\n";
my $l=<IN>;
chomp($l);
$l =~ s/^#//;
my @headers= split(/\t/,$l,-1);
for my $h(@headers_required){
	my $i = array_search(\@headers,$h);
	if($i < 0){ print STDERR "ERROR: $0: missing header: $h\n"; die $usage;}
	$headeri{$h} = $i;
}

# TAXA IN THE TAIL
my %tail_taxa		= ();
if($tail > 0){
my %species_readn	= ();
my $readn_tot		= 0;
while($l=<IN>){
	if($l =~ m/^[#!]/){next;}
	chomp($l);
	my @sp			= split(/\t/,$l,-1);
	my $species		= $sp[$headeri{species}]; 
	my $readn		= $sp[$headeri{readn}];
	$readn_tot		+= $readn;
	$species_readn{$species}		= 0 if(!defined($species_readn{$species}));
	$species_readn{$species}		+= $readn;
}

my $cumsum = 0;
for my $k(sort {$species_readn{$a} <=> $species_readn{$b}} keys %species_readn){
	$cumsum += $species_readn{$k};
	if($cumsum > ($readn_tot*$tail/100)){
		last;
	}
	$tail_taxa{$k}= 1;
}
close(IN);


# Read abund_table again, line-by-line and print krona.data to STDOUT
open(IN,"<$file") or die "Can\'t open $file: $!\n";
	my $l=<IN>;# header
}
while($l=<IN>){
	if($l =~ m/^[#!]/){next;}
	chomp($l);
	my @sp = split(/\t/,$l,-1);
	# reorder data by header
	my %data;
	for my $h(@headers_required){
		$data{$h}	= $sp[$headeri{$h}];
	}
	next if($tail>0 && defined($tail_taxa{$data{species}}));	
	
	my @taxpathsn_line;
	for my $h(@headers_taxnames){
		if(defined($data{$h}) && $data{$h} ne ''){
			push(@taxpathsn_line,$data{$h});
		}
	}
	
	# READN TAB-DELIMITED TAXPATH
	printf("%u\t%s\n", $data{readn}, join("\t",@taxpathsn_line));
}
close(IN);


sub array_search {
    my ($arr, $elem) = @_;
    my $idx= -1;
    for my $i (0..$#$arr) {
        if ($arr->[$i] eq $elem) {
            $idx = $i;
            last;
        }
    }
    return $idx;            
}






