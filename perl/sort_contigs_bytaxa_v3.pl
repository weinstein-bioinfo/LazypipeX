#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# 1) Sorts contig sequences first into main classes
#	1) division			: all contig hits belong to divsion X, e.g. "Viruses"
#	2) divisionX-divisionY	: all contig hits belong to division X or Y, e.g. "Bacteria-Viruses"
#	3) mixture			: contig hits are from a mixture of divisions
#	4) unknown			: no hits for this contig
#	
# 2) In each main class sort into subclasses
#	1) familyX			: all hits belong to the same family
#	2) familyX-familyY	: hits are a mixture of two families
#	3) mixture			: hits are a mixture of tree or more families
#	4) unknown			: all hits belong to an unknown family
# 
# 3) In each subclass sort by genus
#	1) genusX
#	2) genusX-genusY
#	3) mixture
#	4) unknown

my $usage=  "\nUSAGE: $0 -c contigs.fa -a annot.tsv -r dir --toptaxrank division  [-v]\n".
			"\n".
			"-c|contings [fasta]: contigs fasta file\n".
			"-a|annot    [tsv]  : tsv file with contig annotations. This MUST include these columns:\n".
			"                       contig\n".
			"                       \$toptaxrank\n".
			"                       order\n".
			"                       family\n".
			"                       genus\n".
			"                       species\n".
			"-r|res [dir]       : print sorted contigs to this directory\n".
			"--toptaxrank str: top taxonomy rank to use [division]\n".
			"-v      : verbal\n\n";

my $contigs_file 		= !1;
my $contigs_annot_file  = !1;
my $contigs_dir			= !1;
my $verb 				= !1;
my $toptaxrank			= 'division';
GetOptions(	'c|contigs=s' 	=> \$contigs_file,
			'a|annot=s'		=> \$contigs_annot_file,
			'r|res=s'		=> \$contigs_dir,
			'toptaxrank=s'	=> \$toptaxrank,
    			'v'				=> \$verb) or die $usage;
if(!$contigs_file || !$contigs_annot_file || !$contigs_dir){
	die "ERROR: missing arguments\n\n$usage";
}

if( !(-e $contigs_file)){
	die "ERROR: missing file $contigs_file\n";
}

if( !(-e $contigs_annot_file)){
	die "ERROR: missing file $contigs_annot_file";
}


# READING BLAST AND TAXONOMY TO ARRAY
# each contig is classified to main class and a number of subclasses based on taxonomy linked to blast hits
# classification is saved in "seqid" -> "main\tsubclass" hash
my %class_hash;

if($verb){ print STDERR "\t# reading $contigs_annot_file\n";}

open(IN,"<$contigs_annot_file") or die "Can\'t open $contigs_annot_file: $!\n";

my %header; # header_label > column
my @required= ("contig","$toptaxrank","order","family","genus","species");
my $l=<IN>; chomp($l);
my @sp		= split(/\t/,$l,-1);
for(my $i=0; $i<scalar(@sp); $i++){
	$header{lc($sp[$i])} = $i;
}
foreach my $h(@required){
	if( !defined($header{$h})){  die "# ERROR: missing header \"$h\"\n";}
}
$l=<IN>; chomp($l);
@sp = split(/\t/,$l,-1);
my $contig_id	= $sp[$header{'contig'}];
my $new_id;
my (%tax_di,%tax_or,%tax_fa,%tax_ge,%tax_sp);

if(defined($header{bphage}) && ($sp[$header{bphage}] eq 'yes') ){
	$tax_di{Phages}		= 1;
}
else{
	$tax_di{ $sp[$header{$toptaxrank}] }	= 1;
}
$tax_or{ $sp[$header{order}] }		= 1;
$tax_fa{ $sp[$header{family}] } 		= 1;
$tax_ge{ $sp[$header{genus}] }		= 1;
$tax_sp{ $sp[$header{species}] }		= 1;

while($l=<IN>){
	chomp($l);
	@sp 		= split(/\t/,$l,-1);
	$new_id = $sp[$header{'contig'}];
	if( $contig_id ne $new_id ){
		my $class 			= classify_contig(\%tax_di,\%tax_or,\%tax_fa,\%tax_ge,\%tax_sp);
		$class_hash{$contig_id} = $class;
		%tax_di 		= ();
		%tax_or		= ();
		%tax_fa		= ();
		%tax_ge		= ();
		%tax_sp		= ();
	}

	$contig_id 		= $new_id;
	if(defined($header{'bphage'}) && ($sp[$header{'bphage'}] eq 'yes') ){
		$tax_di{Phages}		= 1;
	}
	else{
		$tax_di{$sp[$header{$toptaxrank}]} ? ($tax_di{$sp[$header{$toptaxrank}]}++) : ($tax_di{$sp[$header{$toptaxrank}]}=1);
	}
	$tax_or{$sp[$header{'order'}]} ? ($tax_or{$sp[$header{'order'}]}++) : ($tax_or{$sp[$header{'order'}]}=1);
	$tax_fa{$sp[$header{'family'}]} ? ($tax_fa{$sp[$header{'family'}]}++) : ($tax_fa{$sp[$header{'family'}]}=1);
	$tax_ge{$sp[$header{'genus'}]}  ?  ($tax_ge{$sp[$header{'genus'}]}++)  : ($tax_ge{$sp[$header{'genus'}]}=1);
	$tax_sp{$sp[$header{'species'}]} ? ($tax_sp{$sp[$header{'species'}]}++) : ($tax_sp{$sp[$header{'species'}]}=1);

}
close(IN);
# LAST RECORD
if(scalar(keys %tax_di) > 0){
	my $class = classify_contig(\%tax_di,\%tax_or,\%tax_fa,\%tax_ge,\%tax_sp);
	$class_hash{$contig_id} = $class;
}


# CREATE SORTED MAP STRUCTURES FOR FASTER CONTIG SORTING:
#
# map{division}->{family}->list_of_contigids
# map{division}->{family}->{genus}->list_of_contigids
# map{division}->{family}->{genus}->{species}->list_of_contigids
my %di_fa_ids = ();
my %di_fa_ge_ids = ();
my %di_fa_ge_sp_ids = ();
foreach my $contid( keys %class_hash){
	my $class 					= $class_hash{$contid};
	my($di,$or,$fa,$ge,$sp) 		= split(/\t/,$class,-1);
	
	# init map{$di}
	if( !defined($di_fa_ids{$di}) ){
		my %tmp 		= ();
		$di_fa_ids{$di} = \%tmp;
	}
	if( !defined($di_fa_ge_ids{$di})){
		my %tmp 			= ();
		$di_fa_ge_ids{$di}  = \%tmp;
	}
	if( !defined($di_fa_ge_sp_ids{$di})){
		my %tmp 			= ();
		$di_fa_ge_sp_ids{$di} = \%tmp;
	}
	# init map{$di}->{$fa}
	if( !defined($di_fa_ids{$di}->{$fa}) ){
		my @tmp 		= ();
		$di_fa_ids{$di}->{$fa} = \@tmp;
	}
	if( !defined($di_fa_ge_ids{$di}->{$fa})){
		my %tmp 			= ();
		$di_fa_ge_ids{$di}->{$fa}  = \%tmp;
	}
	if( !defined($di_fa_ge_sp_ids{$di}->{$fa})){
		my %tmp 			= ();
		$di_fa_ge_sp_ids{$di}->{$fa} = \%tmp;
	}
	# init map{$di}->{$fa}->{$ge}
	if( !defined($di_fa_ge_ids{$di}->{$fa}->{$ge})){
		my @tmp 			= ();
		$di_fa_ge_ids{$di}->{$fa}->{$ge}  = \@tmp;
	}
	if( !defined($di_fa_ge_sp_ids{$di}->{$fa}->{$ge})){
		my %tmp 			= ();
		$di_fa_ge_sp_ids{$di}->{$fa}->{$ge} = \%tmp;
	}
	# init map{$di}->{$fa}->{$ge}->{$sp}
	if( !defined($di_fa_ge_sp_ids{$di}->{$fa}->{$ge}->{$sp})){
		my @tmp 			= ();
		$di_fa_ge_sp_ids{$di}->{$fa}->{$ge}->{$sp} = \@tmp;
	}
	
	# put the contid to appropriate collections
	push(@{$di_fa_ids{$di}->{$fa}}, $contid);
	push(@{$di_fa_ge_ids{$di}->{$fa}->{$ge}}, $contid);
	push(@{$di_fa_ge_sp_ids{$di}->{$fa}->{$ge}->{$sp}}, $contid);
}


# READ CONTIG SEQS AND SORT THEM ACCORDING TO TAXONOMY

if($verb){ print STDERR "\t# reading $contigs_file\n";}
my %contigs = 	%{readfasta($contigs_file)};

if($verb){ print STDERR "\t# sorting contigs by taxonomy\n";}
system("rm -fR $contigs_dir");
system("mkdir -p $contigs_dir");

foreach my $di(sort keys %di_fa_ge_sp_ids){
	
	my $di_dir	= "$contigs_dir/$di";
	system("mkdir -p \"$di_dir\"");
	
	foreach my $fa(sort keys %{$di_fa_ge_sp_ids{$di}} ){
		
		# Print collective fasta for Family contigs
		my $fal  	 = $fa; 		# label used in dir structure
		$fal		 =~ s/\s|\/|\\|\(|\)/_/g;
		my $fa_fasta = "$contigs_dir/$di/$fal".'.fa';
		my $fa_dir	 = "$contigs_dir/$di/$fal";
			
		open(OUT,">$fa_fasta") or die "Can\'t open $fa_fasta: $!\n";
		foreach my $seqid(sort @{$di_fa_ids{$di}->{$fa}} ){
			if( !defined($contigs{$seqid})){
				print STDERR "WARNING: no seqid=$seqid found in $contigs_file: SKIPPING\n";
				next;
			}
			print OUT ">$seqid\n";
			print OUT "$contigs{$seqid}\n";
		}
		close(OUT);
		system("mkdir -p \"$fa_dir\"");
		
		foreach my $ge(sort keys %{$di_fa_ge_sp_ids{$di}->{$fa}} ){ 
		
			# Print collective fasta(s) for Genus Contigs
			
			my $gel 	= $ge;
			$gel		=~s/\s|\/|\\|\(|\)/_/g;
			my $ge_fasta= "$contigs_dir/$di/$fal/$gel".'.fa';
			my $ge_dir	= "$contigs_dir/$di/$fal/$gel";
			
			open(OUT, ">$ge_fasta") or die "Can\'t open $ge_fasta: $!\n";
			foreach my $seqid(sort @{$di_fa_ge_ids{$di}->{$fa}->{$ge}} ){
				
				if( !defined($contigs{$seqid}) ){
				
					print STDERR "WARNING: no seqid=$seqid found in $contigs_file: SKIPPING\n";
					next;
				}
			
				print OUT ">$seqid\n";
				print OUT "$contigs{$seqid}\n";
			}
			close(OUT);
			
			system("mkdir -p \"$ge_dir\"");
			
			foreach my $sp(sort keys %{$di_fa_ge_sp_ids{$di}->{$fa}->{$ge}} ){
				
				# Print collective fasta for Species contigs
				my $spl 	= $sp;
				$spl		=~ s/\s|\/|\\|\(|\)/_/g;
				my $sp_fasta= "$contigs_dir/$di/$fal/$gel/$spl".'.fa';
			
				open(OUT, ">$sp_fasta") or die "Can\'t open $sp_fasta: $!\n";
				foreach my $seqid(sort @{$di_fa_ge_sp_ids{$di}->{$fa}->{$ge}->{$sp}} ){
				
					if( !defined($contigs{$seqid})){
						print STDERR "WARNING: no seqid=$seqid found in $contigs_file: SKIPPING\n";
						next;
					}
				
					print OUT ">$seqid\n";
					print OUT "$contigs{$seqid}\n";
				}
				close(OUT);	
			}
		}
	}
}

#
# Parameters: references to hashes of the form  'taxonomy_unit_name' -> number_of_occurences_in_the_contig
# 1: \%tax_division
# 2: \%tax_order	# used for phage classification
# 3: \%tax_family
# 4: \%tax_genus
# 5: \%tax_sciname	# used for phage classification
#
# Returns string:
# DivisionX[-DivisionY]|mixture|unknown tab OrderX[-OrderY]|mixture|unknown tab FamilyX-[FamilyY]|mixture|unknown tab GenusX[-GenusY]|mixture|unknown
#
sub classify_contig{
	my %tax_division		= %{shift(@_)};
	my %tax_order		= %{shift(@_)};
	my %tax_family		= %{shift(@_)};
	my %tax_genus		= %{shift(@_)};
	my %tax_species		= %{shift(@_)};
	
	my @hash_list = (\%tax_division,\%tax_order,\%tax_family,\%tax_genus,\%tax_species);
	my $class = '';
	
	for(my $i=0; $i<5; $i++){
		my %h = %{shift(@hash_list)};
		if($i>0){ $class .= "\t";}
		if($h{''}){ delete $h{''};}
		my @k = sort keys %h;
		
		if(scalar(@k)==0){
			$class .= 'unknown';
		}
		elsif(scalar(@k)==1){
			$class .= $k[0];
		}
		elsif( scalar(@k)==2 ){
			$class .= join('+',@k);
		}
		else{
			$class .= 'mixture';
		}
	}
	
	my($division,$order,$family,$genus,$species) = split(/\t/,$class,-1);

	return join("\t",$division,$order,$family,$genus,$species);
}
	
sub readfasta{
  	my $file		= shift(@_);
	my %sequence;
	my $header;
	my $temp_seq;
	
	#suppose fasta files contains multiple sequences;
	 
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
