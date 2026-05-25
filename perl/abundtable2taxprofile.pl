#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);


# Converts abundance table to taxonomy profile
#
# input: abundance table
#	Table MUST be in tab-deleted format with a header row.
#	Table MUST include the following headers, any order is accepted, other headers are ignored:
#		taxid
#		readn
#		species
#		species_id
#		genus
#		genus_id
#		family
#		family_id
#		order
#		order_id 
#		class
#		class_id
#		phylum
#		phylum_id
#		superkingdom
#		superkingdom_id
#	--sample str	: name of the biological sample (optional)
#	--tail num	: number in percents [0,100], taxa with least abundance summing to this amount will be ignored
#	--hgtaxid	: taxid for the host genome. When given, this will be excluded from the tailed data.
#
# output: taxonomy profile in CAMI format (see https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd)
#
# credit: Ilya Plyusnin, University of Helsink (Ilja.Pljusnin@helsinki.fi)


my $usage=      "USAGE: $0 abund_table -s|sample id [-t|tail num --hgtaxid] 1> taxonomy_profile\n".
		"\n".
		"Converts abundance table to CAMI taxonomy profile\n\n";
		
my $help 	= !1;
my $sampleid= "";
my $tail	= 0;
my $hgtaxid	= 0;
GetOptions('help|h' 	=> \$help,
	'sample|s=s'	=> \$sampleid,
	'tail|t=f'	=> \$tail,
	'hgtaxid=i'	=> \$hgtaxid) or die $usage;
if((scalar @ARGV < 1) || $help){ die $usage; }


# ABUND TABLE <HEADER>
my @headers_taxnames = ('superkingdom','phylum','class','order','family','genus','species');
my @headers_taxids = ();
for my $h(@headers_taxnames){ push(@headers_taxids,join('_',$h,'id'))};
my @headers_required;
push(@headers_required,@headers_taxnames,@headers_taxids,'readn','taxid');
#my @headers_optional = ('strain');
my %headeri;
my $file = shift(@ARGV);
open(IN,"<$file") or die "Can\'t open $file: $!\n";
my $l=<IN>;
chomp($l);
$l =~ s/^#//;
my @headers= split(/\t/,$l,-1);
#DEBUG
#print "headers:$headers[0]\n";
for my $h(@headers_required){
	my $i = array_search(\@headers,$h);
	if($i < 0){ print STDERR "ERROR: $0: missing header: $h\n"; die $usage;}
	$headeri{$h} = $i;
}


# TAXA IN THE TAIL
my %tail_taxa;
if($tail > 0){
my %abund_species; 	# species_taxid > readn
my $readn_tot	= 0;
my $readn_host	= 0;
while($l=<IN>){
	if($l =~ m/^[#!]/){next;}
	chomp($l);
	my @sp		= split(/\t/,$l,-1);
	my $species_id  = $sp[$headeri{'species_id'}]; 
	my $abund	= $sp[$headeri{'readn'}];
	if($hgtaxid &&  ($hgtaxid eq $species_id)){
			$readn_host = $abund;
	}	
	$readn_tot	+= $sp[$headeri{'readn'}];
	if(!defined($abund_species{$species_id})){ $abund_species{$species_id} = 0;}
	$abund_species{$species_id} += $sp[$headeri{'readn'}];
}

my $cumsum = 0;
for my $k(sort {$abund_species{$a} <=> $abund_species{$b}} keys %abund_species){
	$cumsum += $abund_species{$k};
	if($cumsum > (($readn_tot-$readn_host)*$tail/100)){	# host reads excluded from tail cut-off
		last;
	}
	$tail_taxa{$k}= $abund_species{$k};
}
# DEBUG
#my $tmp=0;
#for my $k(sort {$tail_taxa{$a} <=> $tail_taxa{$b}} keys %tail_taxa){
#	$tmp+= $tail_taxa{$k};
#	print "$k\t",$tail_taxa{$k},"\t$tmp\t",($tmp/$readn_tot*100),"\n";
#}exit(1);
close(IN);
open(IN,"<$file") or die "Can\'t open $file: $!\n";
my $l=<IN>;
}



# ABUND TABLE <DATA>
my %abund; 	# taxid > readn
my %rank;	# taxid > rank
my %taxpath;	# taxid > taxpath	: 1224 > 2|1224
my %taxpathsn;	# taxid > taxpathsn	: 1224 > Bacteria|Proteobacteria
my $readn_tot =0;

while($l=<IN>){
	if($l =~ m/^[#!]/){next;}
	chomp($l);
	my @sp = split(/\t/,$l,-1);
	# reorder data by header
	my %data;
	for my $h(@headers_required){
		$data{$h}	= $sp[$headeri{$h}];
	}
	my (@taxpath_line,@taxpathsn_line);
	for my $h(@headers_taxids){
		push(@taxpath_line,$data{$h});
	}
	for my $h(@headers_taxnames){
		push(@taxpathsn_line,$data{$h});
	}
	my $species_id = $data{'species_id'};
	
	if( $tail > 0){	# smalles taxa are ignored, e.g. the smallest 1%
		if(defined($tail_taxa{$species_id})){
			next;
		}
	}
	
	for(my $i=0; $i<scalar(@headers_taxnames); $i++){
		my $h_taxname	= $headers_taxnames[$i];
		my $h_taxid 	= $headers_taxids[$i];
		my $taxid 	= $data{$h_taxid};
		my $taxname	= $data{$h_taxname};
		
		if($taxid eq ""){ # skipping fields that have no taxid
			next;
		}
		
		if( !defined($abund{$taxid})  ){ $abund{$taxid} 		= 0;}
		if( !defined($rank{$taxid}) ){ $rank{$taxid} 		= $h_taxname; }
		if( !defined($taxpath{$taxid}) ){ $taxpath{$taxid} 	= join('|', @taxpath_line[0..$i]); }
		if( !defined($taxpathsn{$taxid}) ){ $taxpathsn{$taxid} 	= join('|', @taxpathsn_line[0..$i]); }
		
		$abund{$taxid} += $data{'readn'};
		
		#DEBUG
		#if($taxid eq ""){
		#	print STDERR $data{'readn'},"\n";
		#}
	}
	$readn_tot += $data{'readn'};
}
close(IN);


# TAXPROFILE <HEADER>
print '@',"Version:0.9.1\n";
print '@',"SampleID:$sampleid\n";
print '@',"Ranks:",join('|',@headers_taxnames),"\n";
print '@',"TaxonomyID:NCBI taxonomy\n";
print '@@',"TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_LAZYPIPE_READN\n";

# TAXPROFILE <OUTPUT>
my %rank_order;	# taxon rank str > 0..10
for(my $i=0; $i<scalar @headers_taxnames; $i++){
	$rank_order{$headers_taxnames[$i]} = $i;
}
for my $taxid( sort {$rank_order{$rank{$a}} <=> $rank_order{$rank{$b}} or $abund{$b} <=> $abund{$a}} keys %abund){
	#print '',join('\t',$taxid,$rank{$taxid}, $taxpath{$taxid}, $taxpathsn{$taxpath},$abund{$taxid}),"\n";
	print $taxid;
	print "\t",$rank{$taxid};
	print "\t",$taxpath{$taxid};
	print "\t",$taxpathsn{$taxid};
	printf("\t%.6f", ($abund{$taxid}/$readn_tot*100));
	print "\t",$abund{$taxid},"\n";
}


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






