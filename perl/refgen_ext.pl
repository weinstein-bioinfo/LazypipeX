#!/usr/bin/perl
use strict;
use warnings;
use Spreadsheet::ParseExcel::Workbook;
use Spreadsheet::ParseExcel::Worksheet;
use Spreadsheet::ParseXLSX;
use File::Basename;
use Getopt::Long qw(GetOptions);
use Config::General;
use Cwd;
#use Bio::DB::GenBank;
#use Bio::DB::Query::GenBank;
#use Bio::DB::NCBIHelper;
use Bio::SeqIO;
use Env;



# Extend viral contigs with AlignGraph and a reference genome. Reference genome is searched from a list of sources in this order:
#	1. local RefSeq+GeneBank collection of viral sequences (loaded from www.genome.jp/virushostdb) based on viral hits
#	2. NCBI complete genomes based on viral hits.
#
# Input: Lazypipe resdir with files: summary[1,2]*.xlsx, contigs.fa, contigs folder and read files
#
# Output: *_refgenext.fa file (in the resdir): extended+notextended contigs
#
# NOTE: attempts to extend contigs for genum if:  genum has csumq score==1 + enough reads + at least 2 contigs + reference genome is found
# 	these thresholds can be controlled by parameters --th_readn --th_contign --th_csumq

my $usage	= "USAGE: $0 resdir [-w|wrk dir --refgen dir --numth int --th_readn int --th_contign int --th_csumq int -v]\n\n";


# OPTIONS FROM CONFIG FILE
my $config_default_file	= "pipeline.default.config";
my $config_my_file	= "pipeline.my.config";	
if( !(-e($config_my_file)) && !(-e($config_default_file))){
	die "ERROR: expecting file: $config_default_file\n";
}
my $config_file	= ( -e($config_my_file)) ? $config_my_file: $config_default_file;
my $conf = Config::General->new(
   -ConfigFile      => "$config_file",
   -InterPolateVars => 1,
   -LowerCaseNames  => 1,
   -AutoTrue        => 1);
my %opt = $conf->getall();
my $wrkdir	= $opt{'wrkdir'};
system("mkdir -p $wrkdir") == 0 or die;
my $refgendir	= $opt{'genomes_viral'}."/fasta";
if( !(-e $refgendir) ){ die "ERROR: expecting viral refgenomes here: $refgendir\n";}


# DEFAULT OPTIONS
my $numth	= 16;
my $verbal	= !1;
my $th_csumq	= 1;
my $th_readn	= 250;
my $th_contign	= 2;

# COMMAND LINE OPTIONS
GetOptions('wrk|w=s'	=> \$wrkdir,
	   'refgen=s'	=> \$refgendir,
	   'numth=i'	=> \$numth,
	   'th_readn=i'	=> \$th_readn,
	   'th_contign=i'=>\$th_contign,
	   'th_csumq=i'	=> \$th_csumq,
	   'v'		=> \$verbal) or die $usage;
	   
if(scalar(@ARGV) < 1){ die $usage; }	
my $resdir 		= shift(@ARGV);
my @ls = <$resdir/summary1*.xlsx>;	my $summary1 = $ls[0];
@ls = <$resdir/summary2*.xlsx>;		my $summary2 = $ls[0];
my $contigs_fasta	= "$resdir/contigs.fa";
my %contigs_hash 	= ();
my $call;


# READ CONTIGS to hash: tested
if($verbal){ print STDERR "# reading $contigs_fasta..\n"; }

my $in = Bio::SeqIO->new(-file => "$contigs_fasta");
while(my $seq = $in->next_seq() ) {
	my $id = $seq->id();	# assuming format "contig=id_\w+"
	my @sp = split('_',$id,-1);
	$id = $sp[0];
	@sp = split('=',$id,-1);
	if(scalar @sp>1){
		$id = $sp[1];
	}
	$seq->id($id);
	$contigs_hash{$id} = $seq;
}


# Hash with lists of hit records ordered by viral taxid: vi_taxid > \@vi_record_list, each record is a hash ref
# Only includes records that pass a filter: readn >= th_readn && csumq <= th_csumq  && contign >= th_contign
my %vitaxid_record = ();

# Taxonomy infor for each viral taxid:  taxid => hash, where hash has keys (species,genus,family)
my %vitaxid_taxonomy = ();

# PARSING SUMMARY1
if($verbal){print STDERR "# parsing $summary1..\n";}

	my $parser 	= Spreadsheet::ParseXLSX->new;
	my $workbook 	= $parser->parse($summary1);
	my $worksheet_name = "not_found";
	my %colname_ind = ();
	my @colname_must_list = ("readn","readn_pc","csumq","contign","family","genus","species","species_id");
	
	foreach my $ws($workbook->worksheets()){
		my $name = $ws->get_name();
		if(($name =~ /vi/gi)  && ($name =~ /species/gi)){
			$worksheet_name= $name;
			last;
		}
	}
	if($worksheet_name eq "not_found"){ print STDERR "# No viral hits found: exiting\n"; exit(0);}

	my $ws	= $workbook->worksheet($worksheet_name);
	my ($col_min,$col_max) =  $ws->col_range();
	my ($row_min,$row_max) = $ws->row_range();
	for(my $i=$col_min; $i<=$col_max;$i++){
		# mapping colnames to col index
		$colname_ind{$ws->get_cell(0,$i)->unformatted()} = $i;
	}
	foreach my $name(@colname_must_list){
		if( !defined($colname_ind{$name})){
			print STDERR "ERROR: $summary1: $worksheet_name: missing required column: $name\n";
			exit(1);
		}
	}
	
	for(my $j=$row_min+1; $j<=$row_max; $j++){
		my $readn	= $ws->get_cell($j,$colname_ind{'readn'})->unformatted();
		my $readn_pc	= $ws->get_cell($j,$colname_ind{'readn_pc'})->unformatted();
		my $csumq	= $ws->get_cell($j,$colname_ind{'csumq'})->unformatted();
		my $contign	= $ws->get_cell($j,$colname_ind{'contign'})->unformatted();
		my $species_id	= $ws->get_cell($j,$colname_ind{'species_id'})->unformatted();
		my $species	= $ws->get_cell($j,$colname_ind{'species'})->unformatted();
		my $genus	= $ws->get_cell($j,$colname_ind{'genus'})->unformatted();
		my $family	= $ws->get_cell($j,$colname_ind{'family'})->unformatted();
		
		# FILTER:
		if($csumq <= $th_csumq  && $readn >= $th_readn  && $contign >= $th_contign){
			my @record_list_tmp	= ();
			$vitaxid_record{$species_id} = \@record_list_tmp;
			my %taxonomy		= ('species'=>$species, 'genus'=>$genus, 'family'=>$family);
			$vitaxid_taxonomy{$species_id} = \%taxonomy;
		}
	}
	if(scalar keys %vitaxid_record < 1){ print STDERR "# no viral taxa with contigs to extend: exiting\n"; exit(0);}
	
	
# PARSING SUMMARY 2
if($verbal){print STDERR "# parsing $summary2..\n";}

	$parser 	= Spreadsheet::ParseXLSX->new;
	$workbook 	= $parser->parse($summary2);
	$worksheet_name = "not_found";
	%colname_ind = ();
		# taxid is from db hits, species_id is assigned by the link_taxonomy.pl
	#@colname_must_list = ("contig","sid","score","hitlength","querylength","taxid","species","species_id");
	@colname_must_list = ("contig","score","taxid","species","species_id"); 
	
	foreach my $ws($workbook->worksheets()){
		my $name = $ws->get_name();
		if(($name =~ /viruses/gi) ){
			$worksheet_name= $name;
			last;
		}
	}
	if($worksheet_name eq "not_found"){ print STDERR "# No viral hits found: exiting\n"; exit(0);}

	$ws	= $workbook->worksheet($worksheet_name);
	($col_min,$col_max) =  $ws->col_range();
	($row_min,$row_max) = $ws->row_range();
	for(my $i=$col_min; $i<=$col_max;$i++){
		# mapping colnames to col index, we use regex matching for hitlength/querylength to ensure we find variations of these
		my $tmp = $ws->get_cell(0,$i)->unformatted();
		if( $tmp =~ m/hitlength/gi ){
			$colname_ind{'hitlength'} = $i;
		}
		elsif( $tmp =~ m/querylength/gi){
			$colname_ind{'querylength'} = $i;
		}
		elsif( $tmp =~ m/bits|score/gi ){ # bits from SANSparallel, score from centrifuge
			$colname_ind{'score'}	= $i;
		}
		else{
			$colname_ind{$ws->get_cell(0,$i)->unformatted()} = $i;
		}
	}
	foreach my $name(@colname_must_list){
		if( !defined($colname_ind{$name})){
			print STDERR "ERROR: $summary2: $worksheet_name: missing required column: $name\n";
			exit(1);
		}
	}
	
	# Parsing rows
	for(my $j=$row_min+1; $j<=$row_max; $j++){
		my $contig	= $ws->get_cell($j,$colname_ind{'contig'})->unformatted();
		my $score	= $ws->get_cell($j,$colname_ind{'score'})->unformatted();
		my $taxid	= $ws->get_cell($j,$colname_ind{'taxid'})->unformatted();	# taxid is the more specific id from the quired DB
		my $species	= $ws->get_cell($j,$colname_ind{'species'})->unformatted();
		my $species_id	= $ws->get_cell($j,$colname_ind{'species_id'})->unformatted();	# species_id is the species-level id up the taxpath
		#my $sid	= $ws->get_cell($j,$colname_ind{'sid'})->unformatted();
		#my $hitlength	= $ws->get_cell($j,$colname_ind{'hitlength'})->unformatted();
		#my $querylength= $ws->get_cell($j,$colname_ind{'querylength'})->unformatted();		
				
		
		if( defined($vitaxid_record{$species_id}) ){
			#if($sid =~ m/species|genus|family|order|no rank/gi){
			#  next; # centrifuge repors scores for taxon-groups
			#}
			my %tmp = ('contig'=>$contig,'score'=>$score, 'taxid'=>$taxid, 'species'=>$species,'species_id'=>$species_id);
			push(@{$vitaxid_record{$species_id}},\%tmp);
		}		
	}
	

# FOR EACH VIRAL TAXON ITERATE HITS IN THE ORDER OF DECREASING ALIGNMENT SCORE
#my $gb = Bio::DB::GenBank->new();
my $seq;
if($verbal){ print STDERR "# search for refgenome(s)\n";}
my ($taxid_refgen_fa,$taxid_contigs_fa,$taxid_r1_fq,$taxid_r2_fq,$taxid_r1_fa,$taxid_r2_fa,$taxid_contigs_ext,$taxid_contigs_noext);

for my $taxid(sort keys %vitaxid_record){ # here taxid is actually species_id assigned by the link_taxonomy.pl

	$taxid_refgen_fa	= "$wrkdir/$taxid"."_refgen.fa";
	$taxid_contigs_fa	= "$wrkdir/$taxid"."_contigs.fa";
	$taxid_contigs_ext	= "$wrkdir/$taxid"."_contigs_ext.fa";
	$taxid_contigs_noext	= "$wrkdir/$taxid"."_contigs_noext.fa";
	$taxid_r1_fq		= "$wrkdir/$taxid"."_R1.fq";
	$taxid_r2_fq		= "$wrkdir/$taxid"."_R2.fq";
	$taxid_r1_fa		= "$wrkdir/$taxid"."_R1.fa";
	$taxid_r2_fa		= "$wrkdir/$taxid"."_R2.fa";


	if($verbal){ print STDERR "\tprocess taxid=$taxid\n"; }
	
	my @record_list		= sort {$b->{'score'} <=> $a->{'score'}} @{$vitaxid_record{$taxid}};
	my $refgen_found	= !1;
	for my $rec(@record_list){
		#if( ($rec->{'sid'}) =~ m/species|genus|family|order|no rank/gi){
		#	next; # centrifuge repors scores for taxon-groups
		#}
		
		# 1: Search for refgen from local database
		my $str		= "$refgendir/".($rec->{'taxid'}).".fa*";
		my @file_list 	= <"$str">;
		if( scalar(@file_list) > 0){
			$taxid_refgen_fa	= $file_list[0];
			print STDERR "\trefgen found in local DB: $taxid_refgen_fa\n";
			$refgen_found 	= 1;
		}
		# 2: search for refgen from GeneBank
		#else{
		#    eval{
		#	$seq = $gb->get_Seq_by_version($rec->{'sid'});
		#	my $seq_desc  = $seq->description();
		#	my $match_str = 'complete genome';
		#	if( $seq_desc =~ m/$match_str/gi ){
		#		print STDERR "\trefgen loaded from GeneBank: ",$rec->{'sid'},"\n";
		#		# write refgen to wrkdir
		#		my $seqio_out = Bio::SeqIO->new(-file => ">$taxid_refgen_fa",-format => 'Fasta');
		#		$seqio_out->write_seq($seq);
		#		$refgen_found	= 1;
		#	}
		#    } or do{
		# 	my $err = $@ || 'unknown'; 
		#  	print STDERR "WARNING: failed to fetch acc=",($rec->{'sid'}),": $err\n";
		#    };
		#}
		if( $refgen_found ){
			last;
		}
	}
	if( !$refgen_found ){
		print STDERR "\tno refgen found for taxid=$taxid: skipping\n";
		next;
	}
	
	# extract reads to wrkdir
	my %taxid_hash = ();
	$taxid_hash{$taxid} = 1;
	for my $rec(@record_list){ $taxid_hash{ $rec->{'taxid'} }= 1; }
	my $taxid_list = join(',',sort keys %taxid_hash);
	
	system_call("cpp/retrieve_reads -t $taxid_list -r $resdir -w $wrkdir -p $taxid -v");
	system("mv $resdir/$taxid"."_R1.fq $taxid_r1_fq") == 0 or die;
	system("mv $resdir/$taxid"."_R2.fq $taxid_r2_fq") == 0 or die;
	if( (-s $taxid_r1_fq) == 0 || (-s $taxid_r2_fq) == 0){
		print STDERR "\tno reads found for taxid=$taxid: skipping\n";
		next;
	}
	# converting reads to fasta
	my $in = Bio::SeqIO->new(-file => "$taxid_r1_fq", -format => 'Fastq');
	my $out= Bio::SeqIO->new(-file => ">$taxid_r1_fa",-format => 'Fasta', -width => 32000); # Bio::SeqIO splits lines 
	while(my $seq = $in->next_seq() ) {
		$out->write_seq($seq);
	}
	$in = Bio::SeqIO->new(-file => "$taxid_r2_fq", -format => 'Fastq');
	$out= Bio::SeqIO->new(-file => ">$taxid_r2_fa",-format => 'Fasta', -width => 32000);
	while(my $seq = $in->next_seq() ) {
		$out->write_seq($seq);
	}

		
	# extract contigs to wrkdir
	my $seqio_out = Bio::SeqIO->new(-file => ">$taxid_contigs_fa", -format => 'Fasta');
	my %contid_hash = ();
	for my $rec(@record_list){
		my $contid = $rec->{'contig'};
		if( defined($contid_hash{$contid}) ){ next;}
		else{ $contid_hash{$contid} = 1; }
		
		if( !defined($contigs_hash{$contid}) ){
			print STDERR "WARNING: no contig $contid in contig.fa: skipping\n";
		}
		my $seq = $contigs_hash{$contid};
		$seqio_out->write_seq($seq);
	}
	
	# Run AlignGraph:
	system_call("cpp/AlignGraph --numth $numth --read1 $taxid_r1_fa --read2 $taxid_r2_fa --contig $taxid_contigs_fa --genome $taxid_refgen_fa --distanceLow 150 --distanceHigh 1300 --extendedContig $taxid_contigs_ext --remainingContig $taxid_contigs_noext --ratioCheck --coverage 5");
	if( (-s $taxid_contigs_ext)==0 ){
		print STDERR "\tno contigs extended for taxid=$taxid: skipping\n";
		next;
	}
	
	# Read extended contigs, add unextended contigs and create contig file with both groups included to sorted-contigs folder
	my %taxonomy = %{$vitaxid_taxonomy{$taxid}};
	my $species = $taxonomy{'species'};
	$species =~ s/ /_/g;
	my $contigs_ext_final = "$resdir/contigs/Viruses/". $taxonomy{'family'} ."/". $species ."_rgext.fa";	
	my @contigs_ext_list;
	$in = Bio::SeqIO->new(-file => "$taxid_contigs_ext");
	while(my $seq = $in->next_seq() ) {
		push(@contigs_ext_list,$seq);
	}
	
	$in = Bio::SeqIO->new(-file => "$taxid_contigs_noext");
	while(my $seq = $in->next_seq() ) {
		push(@contigs_ext_list,$seq);
	}
	
	$out= Bio::SeqIO->new(-file => ">$contigs_ext_final",-format => 'Fasta', -width => 100);
	for my $seq(@contigs_ext_list) {
		$out->write_seq($seq);
	}
	if($verbal){ print STDERR "\tclustered contigs writen to $contigs_ext_final\n"; }
	if($verbal){ print STDERR "\tclustering stats: from ",(scalar keys %contid_hash)," to ",(scalar @contigs_ext_list),"\n"; }
}




# Parses sample label from summary1*.xlsx filename
sub get_label{
	my $fin = shift;
	$fin 	= basename($fin);
	$fin =~ s/summary1-//gi;
	$fin =~ s/.xlsx//i;
	return $fin;
}

sub system_call{
	my $call= shift;
	print STDERR "$call\n";
	my @args= ("bash","-c",$call);
	system(@args) == 0 or die "";
}
		