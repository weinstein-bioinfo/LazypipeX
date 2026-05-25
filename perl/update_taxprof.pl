#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);


# Updates taxonomic profile using the supplied NCBI taxonomy dump
#
# Credit:
# Plyusnin,I., Kant,R., Jaaskelainen,A.J., Sironen,T., Holm,L., Vapalahti,O. and Smura,T. (2020) 
# Novel NGS Pipeline for Virus Discovery from a Wide Spectrum of Hosts and Sample Types. Virus Evolution, veaa091
#
# Contact: grp-lazypipe@helsinki.fi
#


my $usage=      "\nUSAGE: $0 -t|taxonomy path taxprof.txt 1> taxprof.updated.txt\n".
				"\n".
				"t|taxonomy     : path : Path to NCBI taxonomy. MUST contain files \$taxonomy/delnodes.dmp and \$taxonomy/merged.dmp\n".
				"taxprof.txt    : file : CAMI taxonomic profile\n".
				"h|help                : print this help message\n".
				"\n".
				"Updates taxonomic profile according to supplied NCBI taxonomy:\n".
				"\t- delete taxa listed in \$taxonomy/delnodes.dmp\n".
				"\t- merge taxa listed in \$taxonomy/merged.dmp\n\n".
				"Credit:\n".
				"Plyusnin,I., Kant,R., Jaaskelainen,A.J., Sironen,T., Holm,L., Vapalahti,O. and Smura,T. (2020)\n".
				"Novel NGS Pipeline for Virus Discovery from a Wide Spectrum of Hosts and Sample Types. Virus Evolution, veaa091\n".
				"Contact: grp-lazypipe\@helsinki.fi\n\n";
my $help 		= !1;
my $taxonomy	= "\$data/taxonomy";
my $v			= 1;

GetOptions('t|taxonomy=s' 	=> \$taxonomy,
			'h|help'		=> \$help) or die $usage;
$taxonomy		=~ s/\$(\w+)/$ENV{$1}/g;	# expand env vars

if((scalar @ARGV < 1) || $help){ die $usage; }
my $taxprof 	= shift(@ARGV);
my $delnodes	= "$taxonomy/delnodes.dmp";
my $merged		= "$taxonomy/merged.dmp";
my $taxnames	= "$taxonomy/names.dmp";
my $taxnodes	= "$taxonomy/nodes.dmp";
if(!(-e $delnodes)){	die "ERROR: expected file missing: $delnodes\n";	}
if(!(-e $merged)){		die "ERROR: expected file missing: $merged\n";  	}
if(!(-e $taxnames)){	die "ERROR: expected file missing: $taxnames\n";	}
if(!(-e $taxnodes)){	die "ERROR: expected file missing: $taxnodes\n";	}

print STDERR "\treading $delnodes\n";
my %delnodes_h	= ();
read_tohash($delnodes,0,0,\%delnodes_h);
print STDERR "\treading $merged\n";
my %merged_h	= ();
read_tohash($merged,0,2,\%merged_h);
print STDERR "\treading $taxnames\n";
my %taxnames_h	= read_taxnames($taxnames);
print STDERR "\treading $taxnodes\n";
my %taxnodes_h	= read_taxnodes($taxnodes);

# DEBUG
#my $id = 216572;
#my @node = split(/\t/,$taxnodes_h{$id},-1);
#print STDERR "taxid=$id \t",join(":",@node),"\n"; exit(1);

# READING AND UPDATING TAXPROFILE

my %hi; 			# header_label > column
my @required		= ("taxid","taxpath","taxpathsn","percentage");
my $nodes_deleted 	= 0;
my $nodes_merged	= 0;
my $parsed_headers	= 0;
my $parsed_ranks	= 0;
my @ranks			= ();	

open(IN,"<$taxprof") or die "Can\'t open $taxprof: $!\n";
my $ln = 0;
while(my $l=<IN>){
	$ln++;
	chomp($l);
	
	if($l =~ m/^\@Ranks:/i){
		my $str 	= $l;
		$str		=~ s/^\@Ranks://i;
		$str  		=~ s/^\s+|\s+$//g;
		@ranks		= split(/\|/,$str,-1);
		print "$l\n";
		$parsed_ranks = 1;
		#print STDERR "ranks:\n";foreach my $rank(@ranks){print STDERR "\t$rank\n";}; exit(1);
		next;
	}
	
	if($l =~ s/^@@//){
		my @sp		= split(/\t/,$l,-1);
		for(my $i=0; $i<scalar(@sp); $i++){
			$hi{lc($sp[$i])} = $i;
		}
		foreach my $h(@required){
			if( !defined($hi{$h})){  die "ERROR: missing header \"$h\" in $taxprof\n";}
			#print STDERR "\t$h  > $hi{$h}\n";
		}
		print "\@\@$l\n";
		$parsed_headers = 1;
		next;
	}
	
	if($l =~ m/^[@#!]/){
		print $l,"\n";
		next;
	}
	
	if(!$parsed_headers){
		print STDERR "ERROR: missing \@\@header line in $taxprof\n";
	}
	if(!$parsed_ranks){
		print STDERR "ERROR: missing \@Ranks line in $taxprof\n";
	}
	
	my @sp			= split(/\t/,$l,-1);
	my $taxid		= $sp[$hi{'taxid'}];
	my @taxpath 	= split(/\|/,$sp[$hi{'taxpath'}],-1);
	my @taxpathsn	= split(/\|/,$sp[$hi{'taxpathsn'}],-1);
	#my $percentage	= $sp[$hi{'percentage'}];
	
	# updating taxids + taxnames
	if( defined($delnodes_h{$taxid})){
		$nodes_deleted++;
		print STDERR "\tdeleting $taxid\n";
		next;
	}
	if( defined($merged_h{$taxid})){
		print STDERR "\tmerging $taxid to $merged_h{$taxid}\n";
		
		my $taxid_new		= $merged_h{$taxid};
		my @taxpath_new 	= get_parent_taxids(\%taxnodes_h,\@ranks,$taxid_new);
		my @taxpathsn_new	= ();
		
		foreach my $t(@taxpath_new){
			if(defined($taxnames_h{$t})){
				push(@taxpathsn_new, $taxnames_h{$t});
			}
			else{
				push(@taxpathsn_new, "\t");
			}
		}
		
		$sp[$hi{'taxid'}] 		= $taxid_new;
		$sp[$hi{'taxpath'}] 	= join("|",@taxpath_new);
		$sp[$hi{'taxpathsn'}] 	= join("|",@taxpathsn_new);
		$l  					= join("\t",@sp);
		$nodes_merged++;
	}
	print $l,"\n";
}
close(IN);


# print stats
print STDERR "\n";
print STDERR "\tdeleted: $nodes_deleted nodes\n";
print STDERR "\tmerged:  $nodes_merged nodes\n\n";




# read_tohash(file,key_ind,val_ind,\%hash)
# 
# Does not support multiple values for a single key. The last value encountered is saved.
sub read_tohash{
	my $file= shift;
	my $keyi= shift;
	my $vali= shift;
	my $hashp= shift;
	#print STDERR "# read_tohash from $file\n";
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	my $ln=0;
	while(my $l=<IN>){
        	$ln++;
		if($l =~ m/^[@#!]/){ next; }
		chomp($l);
		my @sp= split(/\t/,$l,-1);
        	$hashp->{$sp[$keyi]} = $sp[$vali];
	}
	close(IN);
}

# read taxonomy nodes
#
# USAGE:
# my %taxnodes = read_taxnodes("$taxonomy/nodes.dmp");
# my %taxnames = read_taxnames("$taxonomy/names.dmp");
#
# my $id= 1656;
# my @sp= split(/\t/, $taxnodes{$id});
# print "id         :",$sp[0],"\n";
# print "parent id  :",$sp[1],"\n";
# print "rank       :",$sp[2],"\n";
# print "name       :",$taxnames{$sp[0]},"\n";
#
sub read_taxnodes{
	my $file = shift;
	my %taxnodes;
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	my $ln=0;
	while(my $l=<IN>){
		$ln++;
		chomp($l);
		my @sp= split(/\t\|\t/,$l);
		if(scalar(@sp)<3){
			die "ERROR: missing data at $ln:$file\n";
		}
		$taxnodes{$sp[0]}= join("\t",$sp[0],$sp[1],$sp[2]);
	}
	return %taxnodes;
}

# read taxonomy node names
sub read_taxnames{
	my $file = shift;
	my %taxnames;
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	$ln=0;
	while(my $l=<IN>){
		$ln++;
		chomp($l);
		my @sp= split(/\t\|\t/,$l);
		for(my $j=0; $j<scalar(@sp); $j++){
			$sp[$j] =~ s/^(\|\t)|(\t\|)$//g; # removing leading and trailing \t|\t characters
		}
		if(scalar(@sp)<4){
			die "ERROR: missing data at $ln:$file\n";
		}
		if( !defined($taxnames{$sp[0]}) ){ # if not defined we add any name available
			$taxnames{$sp[0]}= $sp[1];
		}
		elsif( lc($sp[3]) eq 'scientific name' ){ # else we replace saved name with scientific
			$taxnames{$sp[0]}= $sp[1];
		}
	}
	return %taxnames;
}


# Returns parent tax ids, including the node itself		
sub get_parent_taxids{
	# in:
	my $taxnode_ref	= shift(@_);
	my @ranks		= @{shift(@_)};
	my $taxid		= shift(@_);
	
	# out:
	my @parent_taxids;
	
	my %ranks_h		= ();
	foreach my $rank(@ranks){ $ranks_h{$rank} = 1; }
			
	while( defined($taxnode_ref->{$taxid}) ){
		my $node	= $taxnode_ref->{$taxid};
		my @sp		= split(/\t/,$node,-1);
			
		if( defined($ranks_h{$sp[2]}) ){
			push(@parent_taxids,$taxid);
		}
		if($sp[0] == $sp[1]){ 	# node is it's own parent
			last;
		}
		if($taxid == 1){		# root node
			last;
		}
		$taxid		= $sp[1];	# parent taxid
	}
	return reverse(@parent_taxids);
}
