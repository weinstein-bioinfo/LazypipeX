#! /usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);
use YAML::Tiny;
use File::Temp  qw(tempdir);
use POSIX qw(strftime);
my $install_dir;
BEGIN{ $install_dir	= defined($ENV{'LAZYPIPE_INSTALL_DIR'}) ? $ENV{'LAZYPIPE_INSTALL_DIR'} : dirname(__FILE__) };
use lib "$install_dir/perl";	# load from perl-subdir
use Lazypipe::Utils;
use Lazypipe::Utils qw(filebin2uri format_int);
use Lazypipe::SeqAn;

#
# LAZYPIPE: NGS PIPELINE FOR VIRUS DISCOVERY AND METAGENOMICS
# 
# PERL INTERFACE
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
#


my $perl_scripts 			= "$install_dir/perl";
my $R_scripts 				= "$install_dir/R";
my $config_file				= (-e "./config.yaml")? "config.yaml" : "$install_dir/config.yaml";
my $BPHAGE_FILTER			= "$perl_scripts/ICTV.bphage.filter.tsv";
my $VIRUS_FAMILY_HOST		= "$perl_scripts/ICTV.virus.family.host.tsv";
my $VIRUS_GENUS_HOST			= "$perl_scripts/ICTV.virus.genus.host.tsv";
my $VIRUS_FAMILY_GENCOMP		= "$perl_scripts/ICTV.virus.family.gencomp.tsv";

# GLOBAL CONSTANTS
my $PIPELINE_NAME		= "Lazypipe";
my $PIPELINE_VERSION		= "3.1";


my $usage= 	"\nUSAGE: $0 -1 file [-2 file] -r|res dir -s|sample str -p main\n".
		"\n".
		"Input:\n".
		"-1|read1 file    : PE forward reads in fastq\n".
		"-2|read2 file    : PE reverse reads in fastq [guess from --read1]\n".
		"--se             : Reads are SE-reads. Any --read2 file will be ignored [false]\n".
		"--hostgen str    : List of host/contaminant genome fasta or database keys in config.yaml\n".
		"                   e.g. --hostgen Homo_sapiens,Ixodes_scapularis\n".
		"--config file    : Configuration file [$config_file]\n".
		"\n".
		"Output:\n".
		"--logs   dir     : Root directory for logs [logs]. Logs will be printed to logs-dir/sample/\n".		
		"-r|res   dir     : Root directory for results [results]. Results will be printed to res-dir/sample/\n".
		"-s|sample str    : Sample label [read1 filename]\n".
		"--tmpdir dir     : Root for temporary directory\n".
		"\n".
		"Parameters:\n".
		"-p|pipe str      : List of steps to perform, eg --pipe pre,flt,ass,ann1 [main]\n".
		"     pre|preprocess : Preprocess reads, i.e. filter low quality reads\n".
		"     flt|filter     : Filter host reads with --hostgen fasta or database\n".
		"     ass|assemble   : Assemble reads to contigs\n".
		"     rea|realign    : Realign reads to contigs\n".
		"     ann1|annot1    : Run 1st round annotation\n".
		"     ann2|annot2    : Run 2nd round annotaiton\n".
		"     rep|report     : Create reports\n".
		"     rgrep|rgreport : Create reference genome reports\n".
		"     sta|stats      : Create assembly stats + QC plots\n".
		"     pack           : Pack results to a tarball\n".
		"     clean          : Clean up intermediate/temporary files.\n".
		"     main           : Run main steps: pre,flt,ass,rea,ann1,ann2,rep,sta,pack,clean [default]\n".
		"     all            : Run all steps: pre,flt,ass,rea,ann1,ann2,rep,rgrep,sta,pack,clean\n".		
		"--ann1 key       : List of database keys defining 1st round annotation\n".
		"                   For each key their MUST be a database defined in config.yaml\n".
		"                   e.g. --ann1 minimap.nt.abv\n".
		"--ann2 target:key: List of target-key pairs defining 2nd round annotations\n".
		"                   \$target  : valid targets are 'ab' (archaea and bacteria), 'ph' (bacteriophages), 'vi' (viruses), 'un' (unmapped)\n".
		"                   \$key     : database key in config.yaml\n".
		"                   e.g. --ann2 vi:blastn.nt.vi,ab:blastn.nt.ab\n".
		"--anns key       : Apply annotation-strategy defined in config.yaml under the supplied key. Overrides any --ann1/ann2 options\n".
		"                   e.g. --anns vi.nt\n".
		"--append         : Append annotation to the existing annot1.tsv or annot2.tsv [false]\n".
		"--ass str        : Assembler: megahit|spades [megahit]\n".
		"--gen str        : Gene prediction: mga|prod [mga]\n".
		"--pre str        : Use fastp|trimm|none to preprocess reads [fastp]\n".
		"--clean          : Delete intermediate files after each step [false]\n".
		"--pack_reads     : Pack trimmed and background filtered reads to *.tar.gz\n".
		"-t|numth num     : Number of threads\n".
		"-w|wmodel str    : Weighting model for abundance estimation: taxacount|bitscore|bitscore2 [bitscore]\n".
		"--databases      : List installed reference databases \n".
		"--filters        : List installed background filters\n".
		"-v               : Verbal mode [false]\n".
		"-h|help          : Print this manual\n".
		"\n".
		"NOTE: command line options take precedence over $config_file options\n".
		"\n".
		"CREDIT:\n".
		"Plyusnin,I., Kant,R., Jaaskelainen,A.J., Sironen,T., Holm,L., Vapalahti,O. and Smura,T. (2020)\n".
		"Novel NGS Pipeline for Virus Discovery from a Wide Spectrum of Hosts and Sample Types. Virus Evolution, veaa091\n\n".
		"CONTACT:\n".
		"grp-lazypipe\@helsinki.fi\n\n"; 


# Read options from config.yaml
	# read config.yaml from command line, if specified
Getopt::Long::Configure("pass_through");
GetOptions('config=s' => \$config_file);
my $yaml 			= YAML::Tiny->read( $config_file );
my %opt 				= %{$yaml->[0]};
my $commandline 		= join " ", $0, @ARGV;
my $time 			= strftime "%Y/%m/%d %H:%M:%S", localtime;
GetOptions(\%opt, 'read1|1=s','read2|2=s','se','res|r=s','logs=s','tmpdir=s','sample|s=s','numth|t=i',
			'pre=s','ass=s','gen=s','ann1|annot1=s','ann2|annot2=s','anns|annstrat=s','append','hostgen=s','hgtaxid=i','wmodel|w=s',
			'min_read2hostgen_score=i',
			'min_read2hostgen_mapq=i',
			'min_read2contig_score=i',
			'min_read2contig_mapq=i',
			'min_orf_length=i',
			'min_sans_bits=i',
			'min_blastp_bits=i',
			'min_blastn_bits=i',
			'min_minimap_DPpeak_score=i',
			'min_psearch_bits_abund=i',
			'min_nsearch_bits_abund=i',
			'min_qcov_abund=f',
			'min_qcov_annot=f',
			'refgenrep_min_clen_sum=i',
			'refgenrep_min_alen_sum=i',
			'refgenrep_max_refgen=i',
			'min_contig2hostgen_score=i',
			'min_contig_length=i',
			'tail=i',
			'tail_contig=i',
			'taxgroups=s',
			'toptaxrank=s',
			'pack_reads',
			'trimm_sample_name',
			'databases',
			'filters',
			'pipe|p=s','help|h','v','clean') or die $usage;
%opt					= options_format(\%opt);
my $VERBAL			= $opt{v};
$Lazypipe::SeqAn::VERBAL = $VERBAL;
$Lazypipe::Utils::VERBAL	 = $VERBAL;

# Choosing not to save complete config.yaml but to print just the commandline to History.log
# options_save2yaml(\%opt);
system("echo $time\t:$commandline >> $opt{res}/History.log");
my %pipe 			= %{$opt{pipe}};

# RUN PIPELINE STEP-BY-STEP
	## PREPROCESSING
pipe_preprocess(\%opt) if( $pipe{prepro} );
pipe_filter_hostreads(\%opt) if($pipe{filter});
pipe_assemble(\%opt) if($pipe{assemble});
pipe_filter_hostcont(\%opt) if($pipe{filter} && $pipe{assemble});
realign_reads_contigs(\%opt) if( $pipe{realign} );
	## ANNOTATION
pipe_annotation_round1(\%opt) if($pipe{ann1} && defined($opt{ann1}) && $opt{ann1} && scalar(@{$opt{ann1}})>0 );
pipe_annotation_round2(\%opt) if($pipe{ann2} && defined($opt{ann2}) && $opt{ann2} && scalar(@{$opt{ann2}})>0 );
pipe_annotation_ictv(\%opt) if( $pipe{ictv});
	## REPORTING, STATISTICS, PACKING AND CLEANUP
generate_reports(\%opt) if( $pipe{report} );
generate_refgen_reports(\%opt) if( $pipe{rgreport});
generate_stats(\%opt) if( $pipe{stats} );
pack_files(\%opt) if( $pipe{pack} );
clean(\%opt) if( $pipe{clean} );
# PIPELINE CALLS END HERE


sub pipe_annotation_round1{
	print STDERR "# ANNOTATION ROUND1\n";
	
	my $subid			= "pipe_annotation_round1()";
	
	# IN
	my %opt				= %{shift()};
	my $contigs 			= "$opt{res}/contigs.fa";
	my $contigs_info		= "$contigs.info.tsv";
	
	# OUT
	my $contigs_ann1_ab		= "$opt{res}/contigs.ann1.ab.fa";
	my $contigs_ann1_vi		= "$opt{res}/contigs.ann1.vi.fa";
	my $contigs_ann1_ph		= "$opt{res}/contigs.ann1.ph.fa";
	my $contigs_ann1_un		= "$opt{res}/contigs.ann1.un.fa";
	my $annot1				= "$opt{res}/annot1.tsv";
	my $dbhits_tmp		 	= "$opt{res}/dbhits.tsv.tmp";
	my $dbhits_blastn		= "$opt{res}/dbhits.blastn.tsv";
	my $dbhits_blastp		= "$opt{res}/dbhits.blastp.tsv";
	my $dbhits_blastx		= "$opt{res}/dbhits.blastx.tsv";
	my $dbhits_diamondx		= "$opt{res}/dbhits.diamondx.tsv";
	my $dbhits_diamondp		= "$opt{res}/dbhits.diamondp.tsv";
	my $dbhits_hmmscan		= "$opt{res}/dbhits.hmmscan.tsv";
	my $dbhits_minimap		= "$opt{res}/dbhits.minimap.paf";
	my $dbhits_sans			= "$opt{res}/dbhits.sans.tsv";
	my $log_fltcont 			= "$opt{logs}/$opt{sample}/filter.contigs.log";
	my $log_blastn			= "$opt{logs}/$opt{sample}/annot.blastn.log";
	my $log_blastp			= "$opt{logs}/$opt{sample}/annot.blastp.log";
	my $log_blastx			= "$opt{logs}/$opt{sample}/annot.blastx.log";
	my $log_diamondx			= "$opt{logs}/$opt{sample}/annot.diamondx.log";
	my $log_diamondp			= "$opt{logs}/$opt{sample}/annot.diamondp.log";
	my $log_hmmscan			= "$opt{logs}/$opt{sample}/annot.hmmscan.log";
	my $log_minimap 			= "$opt{logs}/$opt{sample}/annot.minimap.log";
	my $log_sans				= "$opt{logs}/$opt{sample}/annot.sans.log";
	my $log_report			= "$opt{logs}/$opt{sample}/generate_reports.log";
	my $log_rgreport			= "$opt{logs}/$opt{sample}/generate_rgreports.log";
	
	# PARAMS
	my $cond_ab		= "\'\$division==\"Bacteria\"'";
	my $cond_vi		= "\'\$division==\"Viruses\" && \$bphage==\"no\"\'";
	my $cond_ph		= "\'\$division==\"Phages\" || \$bphage==\"yes\"\'";
	my @anns 		= @{$opt{ann1}};
	my $append 		= $opt{append};

	# UPDATE TAXONOMY
	if( $opt{taxonomy}->{update} ){
		update_taxonomy(\%opt);
	}
	if(nlines($contigs)<2){
		print STDERR "\n\tWARNING: $subid: empty contig file: no annotations will be produced\n\n";
		return;
	}
	
	foreach my $ann(@anns){
		my $contigs_ta 	= ($ann->{target} eq 'un') ? $contigs_ann1_un : $contigs;	# just two options in ann1
		
		if($ann->{search} eq "minimap"){
			annotate_minimap(seqs		=>$contigs_ta,
							db			=>$ann->{db},
							taxonomy		=>$opt{taxonomy}->{db},
							dbhits		=>$dbhits_minimap,
							annot		=>$annot1,
							append		=>$append,
							log			=>$log_minimap,
							min_bits		=>$opt{min_minimap_DPpeak_score},
							numth		=>$opt{numth},
							retain_ties	=>0);
							
			if($opt{clear}){		system("rm -f $dbhits_minimap");}
		}
		elsif($ann->{search} eq "sans"){
			annotate_sans(	seqs			=>$contigs_ta,
							seqinfo		=>$contigs_info,
							taxonomy		=>$opt{taxonomy}->{db},
							dbhits		=>$dbhits_sans,
							annot		=>$annot1,
							append		=>$append,
							log			=>$log_sans,
							min_bits		=>$opt{min_sans_bits},
							numth		=>$opt{numth},
							orf_finder	=>$opt{gen},
							min_orf_length=>$opt{min_orf_length},
							retain_ties =>0);
							
			if($opt{clear}){		system("rm -f $dbhits_sans");}
		}
		elsif($ann->{search} eq "blastn"){
			annotate_blastn(	seqs			=> $contigs_ta,
							db			=> $ann->{db},
							dbhits		=> $dbhits_blastn,
							annot		=> $annot1,
							append		=> $append,
							log			=> $log_blastn,
							numth		=> $opt{numth},
							min_bits		=> $opt{min_blastn_bits},
							retain_ties => 0);
							
			if($opt{clear}){		system("rm -f $dbhits_blastn");}
		}
		elsif($ann->{search} eq "blastp"){
			annotate_blastp(seqs			=>$contigs_ta,
							seqinfo		=>$contigs_info,
							db			=>$ann->{db},
							dbhits		=>$dbhits_blastp,
							annot		=>$annot1,
							append		=>$append,
							log			=>$log_blastp, 
							numth			=>$opt{numth},
							min_bits			=>$opt{min_blastp_bits},
							orf_finder		=>$opt{gen},
							min_orf_length	=>$opt{min_orf_length},
							retain_ties		=>0);
							
			if($opt{clear}){		system("rm -f $dbhits_blastp");}
		}
		elsif($ann->{search} eq "blastx"){
			annotate_blastp(seqs			=>$contigs_ta,
							seqinfo		=>$contigs_info,
							db			=>$ann->{db},
							dbhits		=>$dbhits_blastx,
							annot		=>$annot1,
							append		=>$append,
							log			=>$log_blastx, 
							numth			=>$opt{numth},
							min_bits			=>$opt{min_blastp_bits},
							retain_ties		=>0);
							
			if($opt{clear}){		system("rm -f $dbhits_blastp");}
		}		
		elsif($ann->{search} eq "diamondx"){
			annotate_diamondx(seqs		=>$contigs_ta,
							  db			=>$ann->{db},
							  annot		=>$annot1,
							  dbhits		=>$dbhits_diamondx,
							  append		=>$append,
							  filter_tophits=> 1,
							  retain_ties	=> 0,
							  log		=>$log_diamondx, 
							  numth		=>$opt{numth},
							  min_score	=>$opt{min_diamond_bits},
							  min_orf	=>$opt{min_orf_length},
							  max_target_seqs	=> 5,
							  sensitivity		=> 'very-sensitive',
							  tmpdir		=>$opt{tmpdir});
			if($opt{clear}){		system("rm -f $dbhits_diamondx");}
		}
		elsif($ann->{search} eq "diamondp"){
			annotate_diamondp(seqs		=>$contigs_ta,
							  seqinfo	=>$contigs_info,
							  db			=>$ann->{db},
							  annot		=>$annot1,
							  append		=>$append,
							  dbhits		=>$dbhits_diamondp,
							  log		=>$log_diamondp, 
							  filter_tophits=> 1,
							  retain_ties	=> 0,
							  numth		=>$opt{numth},
							  orf_finder	=>$opt{gen},
							  min_score	=>$opt{min_diamond_bits},
							  min_orf	=>$opt{min_orf_length},
							  max_target_seqs	=> 5,
							  sensitivity		=> 'very-sensitive',
							  tmpdir		=>$opt{tmpdir});
			if($opt{clear}){		system("rm -f $dbhits_diamondp");}
		}
		
		elsif($ann->{search} eq "hmmscan"){
			annotate_hmmscan(annot		=>$annot1,
							append		=>$append,
							db			=>$ann->{db},
							dbtaxid		=>10239,
							dbhits		=>$dbhits_hmmscan,
							seqs			=>$contigs_ta,
							seqinfo		=>$contigs_info,
							log			=>$log_hmmscan,
							min_bits		=>$opt{min_hmmscan_bits},
							max_eval		=>$opt{max_hmmscan_eval},
							numth		=>$opt{numth},
							orf_finder	=>$opt{gen},
							min_orf_length=>$opt{min_orf_length}
							);
			if($opt{clear}){		system("rm -f $dbhits_hmmscan");}
		}
		else{
			print "\tWARNING: unknown \$search option in: --ann1 $ann->{search}: skipping\n";
			next;
		}
		filter_unseqs(seqs=>$contigs, seqs_un=>$contigs_ann1_un, annot=>$annot1, seqidh=>'qseqid', log=>$log_fltcont);
		
		# in consecutive steps append
		$append		= 1;
	}
	# after all steps extract contig-classes:
	add_division_field(annot=>$annot1, taxonomy=>$opt{taxonomy}->{db},numth=>$opt{numth}, overwrite=>1);
	add_bphage_field(annot=>$annot1, phfilter=>$BPHAGE_FILTER,log=>$log_fltcont, taxonomy=>$opt{taxonomy}->{db},numth=>$opt{numth}, overwrite=>1);
	filter_seqs(seqs=>$contigs,seqs_flt=>$contigs_ann1_vi, annot=>$annot1, seqidh=>'qseqid', condition=>$cond_vi, log=>$log_fltcont);
	filter_seqs(seqs=>$contigs,seqs_flt=>$contigs_ann1_ph, annot=>$annot1, seqidh=>'qseqid', condition=>$cond_ph, log=>$log_fltcont);
	filter_seqs(seqs=>$contigs,seqs_flt=>$contigs_ann1_ab, annot=>$annot1, seqidh=>'qseqid', condition=>$cond_ab, log=>$log_fltcont);
	filter_unseqs(seqs=>$contigs,seqs_un=>$contigs_ann1_un, annot=>$annot1, seqidh=>'qseqid', log=>$log_fltcont);
	return;
}
sub pipe_annotation_round2{
	print STDERR "# ANNOTATION ROUND2\n";
	
	my $subid			= "pipe_annotation_round2()";
	# IN
	my %opt				= %{shift()};
	my $contigs 			= "$opt{res}/contigs.fa";
	my $contigs_info		= "$contigs.info.tsv";
	
	# OUT
	my $contigs_ann2_ab		= "$opt{res}/contigs.ann2.ab.fa";
	my $contigs_ann2_vi		= "$opt{res}/contigs.ann2.vi.fa";
	my $contigs_ann2_ph		= "$opt{res}/contigs.ann2.ph.fa";
	my $contigs_ann2_un		= "$opt{res}/contigs.ann2.un.fa";
	my $annot2				= "$opt{res}/annot2.tsv";
	my $dbhits_tmp		 	= "$opt{res}/dbhits.tsv.tmp";
	my $dbhits_blastn		= "$opt{res}/dbhits.blastn.tsv";
	my $dbhits_blastp		= "$opt{res}/dbhits.blastp.tsv";
	my $dbhits_blastx		= "$opt{res}/dbhits.blastx.tsv";
	my $dbhits_diamondx		= "$opt{res}/dbhits.diamondx.tsv";
	my $dbhits_diamondp		= "$opt{res}/dbhits.diamondp.tsv";
	my $dbhits_hmmscan		= "$opt{res}/dbhits.hmmscan.tsv";
	my $dbhits_minimap		= "$opt{res}/dbhits.minimap.paf";
	my $dbhits_sans			= "$opt{res}/dbhits.sans.tsv";
	my $log_fltcont 			= "$opt{logs}/$opt{sample}/filter.contigs.log";
	my $log_blastn			= "$opt{logs}/$opt{sample}/annot.blastn.log";
	my $log_blastp			= "$opt{logs}/$opt{sample}/annot.blastp.log";
	my $log_blastx			= "$opt{logs}/$opt{sample}/annot.blastx.log";
	my $log_diamondx			= "$opt{logs}/$opt{sample}/annot.diamondx.log";
	my $log_diamondp			= "$opt{logs}/$opt{sample}/annot.diamondp.log";
	my $log_hmmscan			= "$opt{logs}/$opt{sample}/annot.hmmscan.log";
	my $log_minimap 			= "$opt{logs}/$opt{sample}/annot.minimap.log";
	my $log_sans				= "$opt{logs}/$opt{sample}/annot.sans.log";
	my $log_report			= "$opt{logs}/$opt{sample}/generate_reports.log";
	my $log_rgreport			= "$opt{logs}/$opt{sample}/generate_rgreports.log";
	
	# PARAMS
	my $cond_ab		= "\'\$division==\"Bacteria\"'";
	my $cond_vi		= "\'\$division==\"Viruses\" && \$bphage==\"no\"\'";
	my $cond_ph		= "\'\$division==\"Phages\" || \$bphage==\"yes\"\'";
	my @anns 		= @{$opt{ann2}};
	my $append 		= $opt{append};
	my %targets 		= ();
	
	
	foreach my $ann(@anns){
		
		my $contigs_ta 				= sprintf("$opt{res}/contigs.ann1.%s.fa", $ann->{target});
		$targets{$ann->{target}} 	= 1;
		
		if(!(-e $contigs_ta) || nlines($contigs_ta)<2){
			print STDERR "WARNING: $subid: empty contigs fasta: skipping\n";
			next;
		}
		
		if($ann->{search} eq 'blastn'){
			annotate_blastn(	seqs		=> $contigs_ta,
							db		=> $ann->{db},
							dbhits	=> $dbhits_blastn,
							annot	=> $annot2,
							append	=> $append,
							log		=> $log_blastn,
							numth	=> $opt{numth},
							min_bits=> $opt{min_blastn_bits});
							
			if($opt{clear}){		system("rm -f $dbhits_blastn");}	
		}
		elsif($ann->{search} eq "blastp"){
			annotate_blastp(seqs			=>$contigs_ta,
							seqinfo		=>$contigs_info,
							db			=>$ann->{db},
							dbhits		=>$dbhits_blastp,							
							annot		=>$annot2,
							append		=>$append,
							log			=>$log_blastp, 
							numth			=>$opt{numth},
							min_bits			=>$opt{min_blastp_bits},
							orf_finder		=>$opt{gen},
							min_orf_length	=>$opt{min_orf_length});
							
			if($opt{clear}){		system("rm -f $dbhits_blastp");}
		}
		elsif($ann->{search} eq "blastx"){
			annotate_blastp(seqs			=>$contigs_ta,
							seqinfo		=>$contigs_info,
							db			=>$ann->{db},
							dbhits		=>$dbhits_blastx,
							annot		=>$annot2,
							append		=>$append,
							log			=>$log_blastx, 
							numth			=>$opt{numth},
							min_bits			=>$opt{min_blastp_bits},
							retain_ties		=>0);
							
			if($opt{clear}){		system("rm -f $dbhits_blastp");}
		}	
		elsif($ann->{search} eq "minimap"){
			annotate_minimap(seqs		=>$contigs_ta,
							db			=>$ann->{db},
							dbhits		=>$dbhits_minimap,
							annot		=>$annot2,
							append		=>$append,
							log			=>$log_minimap,
							min_bits		=>$opt{min_minimap_DPpeak_score},
							numth		=>$opt{numth},
							taxonomy		=>$opt{taxonomy}->{db});
							
			if($opt{clear}){		system("rm -f $dbhits_minimap");}
		}
		elsif($ann->{search} eq "sans"){
			annotate_sans(	seqs			=>$contigs_ta,
							seqinfo		=>$contigs_info,
							dbhits		=>$dbhits_sans,
							annot		=>$annot2,
							append		=>$append,
							log			=>$log_sans,
							min_bits		=>$opt{min_sans_bits},
							numth		=>$opt{numth},
							orf_finder	=>$opt{gen},
							min_orf_length=>$opt{min_orf_length});
							
			if($opt{clear}){		system("rm -f $dbhits_sans");}
		}
		elsif($ann->{search} eq "hmmscan"){
			annotate_hmmscan(annot		=>$annot2,
							append		=>$append,
							db			=>$ann->{db},
							dbtaxid		=>10239,
							dbhits		=>$dbhits_hmmscan,
							seqs			=>$contigs_ta,
							seqinfo		=>$contigs_info,
							log			=>$log_hmmscan,
							min_bits		=>$opt{min_hmmscan_bits},
							max_eval		=>$opt{max_hmmscan_eval},
							numth		=>$opt{numth},
							orf_finder	=>$opt{gen},
							min_orf_length=>$opt{min_orf_length}
							);
			if($opt{clear}){		system("rm -f $dbhits_hmmscan");}
		}
		else{
			print "\tWARNING: unknown \$search option in: --ann2 $ann->{str}: skipping\n";
			next;
		}
		# in consecutive steps append
		$append = 1;
	}

	# after annot2 extract confirmed ab/ph/vi/un contigs:
	add_division_field(annot=>$annot2, taxonomy=>$opt{taxonomy}->{db},numth=>$opt{numth}, overwrite=>1);
	add_bphage_field(annot=>$annot2, phfilter=>$BPHAGE_FILTER,log=>$log_fltcont, taxonomy=>$opt{taxonomy}->{db},numth=>$opt{numth}, overwrite=>1);
	if( defined($targets{ab}) ){
		filter_seqs(seqs=>$contigs,seqs_flt=>$contigs_ann2_ab,annot=>$annot2,seqidh=>'qseqid',condition=>$cond_ab, log=>$log_fltcont);
	}
	if( defined($targets{ph})){
		filter_seqs(seqs=>$contigs,seqs_flt=>$contigs_ann2_ph,annot=>$annot2,seqidh=>'qseqid',condition=>$cond_ph, log=>$log_fltcont);
	}
	if( defined($targets{vi}) ){
		filter_seqs(seqs=>$contigs,seqs_flt=>$contigs_ann2_vi,annot=>$annot2,seqidh=>'qseqid',condition=>$cond_vi, log=>$log_fltcont);
	}	
	if( defined($targets{ab}) && defined($targets{vi}) ){
		filter_unseqs(seqs=>$contigs,seqs_un=>$contigs_ann2_un,annot=>$annot2,seqidh=>'qseqid', log=>$log_fltcont);
	}
	return;
}


sub pipe_preprocess{
	print STDERR "\n\n# PREPROCESS READS\n\n";
	my $subid		= "pipe_preprocess()";
	
	# in:
	my %opt 			= %{shift() };
	my $r1 			= $opt{'read1'};
	my $r2			= $opt{'read2'} || undef;
	
	# out:
	my $p1			= "$opt{res}/reads/read1.trim.fq.gz";
	my $up1			= "$opt{res}/reads/read1.trim.unpaired.fq.gz";
	my $p2			= "$opt{res}/reads/read2.trim.fq.gz";
	my $up2			= "$opt{res}/reads/read2.trim.unpaired.fq.gz";
	my $log  		= "$opt{logs}/$opt{sample}/prepro_reads.log";
	my $fastp_html	= "$opt{res}/reports/fastp.report.html";
	my $fastp_json	= "$opt{res}/reports/fastp.json";
	
	system("rm -f $log");
	system("touch $log");
	system("mkdir -p $opt{res}/reads");
	
	if( $opt{'se'} ){ # SE-reads
		if( $opt{'pre'} eq 'fastp' ){
			system_call("fastp --thread $opt{'numth'} -j $fastp_json -h $fastp_html -i $r1 -o $p1 $opt{'par_fastp'}  2>> $log", $opt{'v'});			
		}
		elsif( $opt{'pre'} eq 'trimm' ){
			system_call("trimmomatic SE -threads $opt{'numth'} $r1 $p1 $opt{'par_trimm'} &>> $log", $opt{'v'});
		}
		elsif( $opt{'pre'} eq 'none' ){
			print STDERR "\tno preprocessing\n";
			system_call("cp $r1 $p1", $opt{'v'});
		}
	}
	else{		# PE-reads
		if( $opt{'pre'} eq 'fastp' ){
			system_call(
				"fastp --thread $opt{'numth'} -j $fastp_json -h $fastp_html -i $r1 -I $r2 -o $p1 -O $p2 --unpaired1 $up1 --unpaired2 $up2 $opt{'par_fastp'} 2>> $log", $opt{'v'});
		}
		elsif( $opt{'pre'} eq 'trimm' ){
			system_call(
				"trimmomatic PE -threads $opt{'numth'} $r1 $r2 $p1 $up1 $p2 $up2 $opt{'par_trimm'} &>> $log", $opt{'v'});
		}
		elsif( $opt{'pre'} eq 'none' ){
			print STDERR "\tno preprocessing\n";
			system_call("cp $r1 $p1", $opt{'v'});
			system_call("cp $r2 $p2", $opt{'v'});
		}
	}
	
	# rm temp
	#system("rm -f $fastp_json");
}

sub pipe_filter_hostreads{
	my $subid 		= "pipe_filter_hostreads()";
	
	# IN:
	my %opt			= %{shift()};
	my $r1			= "$opt{res}/reads/read1.trim.fq.gz";
	my $r2	        = ($opt{se})? undef: "$opt{res}/reads/read2.trim.fq.gz";
	
	foreach my $hostdb( @{$opt{hostdb}} ){
		
		my $hostpref		= $hostdb->{latinName} || $hostdb->{commonName} || $hostdb->{name} || $hostdb->{db};
		$hostpref		=~ s/\..*$//g;
		$hostpref		=~ s/\s/_/g;
		
		my $r1_pass		= "$opt{res}/reads/read1.trim.hflt.fq.gz";
		my $r1_pass_tmp	= "$opt{res}/reads/read1.trim.hflt.tmp.fq.gz";	# in case $r1 and $r1_pass are the same file
		my $r1_flt		= "$opt{res}/reads/$hostpref.read1.fq.gz";
		my $r2_pass		= ($opt{se})? undef: "$opt{res}/reads/read2.trim.hflt.fq.gz";
		my $r2_pass_tmp	= ($opt{se})? undef: "$opt{res}/reads/read2.trim.hflt.tmp.fq.gz";
		my $r2_flt		= ($opt{se})? undef: "$opt{res}/reads/$hostpref.read2.fq.gz";
			
		filter_host_reads(
			r1 			=> $r1,
			r1_pass		=> $r1_pass_tmp,
			r1_flt		=> $r1_flt,
			r2 			=> $r2,
			r2_pass		=> $r2_pass_tmp,
			r2_flt		=> $r2_flt,
			hostdb 		=> $hostdb->{db},
			res 			=> $opt{res},
			log 			=> "$opt{logs}/$opt{sample}/filter.host.reads.log",
			numth		=> $opt{numth},
			bitscore 	=> $opt{min_read2hostgen_score},
			mapq			=> $opt{min_read2hostgen_mapq},
			tmpdir 		=> $opt{tmpdir} || undef,
			gzip			=> $opt{gzip});
		
		if(-e $r1_pass_tmp){
			system("mv $r1_pass_tmp $r1_pass");		}
		if(-e $r2_pass_tmp){
			system("mv $r2_pass_tmp $r2_pass");		}
		
		# channel $r1/2_pass to the next host-filter input
		$r1		= $r1_pass;
		$r2		= $r2_pass || undef;
	}
}

sub pipe_filter_hostcont{
	my $subid		= "pipe_filter_hostcont()";
	
	# IN:
	my %opt			= %{shift()};
	my $contigs		= "$opt{res}/contigs.fa";
	
	if(nlines($contigs)<2){
		return;
	}
	
	if( defined($opt{hostdb}) )	{
		foreach my $hostdb( @{$opt{hostdb}} ){
		
			my $hostpref		= $hostdb->{latinName} || $hostdb->{commonName} || $hostdb->{name} || $hostdb->{db};
			$hostpref		=~ s/\..*$//g;
			$hostpref		=~ s/\s/_/g;
			
			my $contigs_pass	= "$opt{res}/contigs.pass.fa";
			my $contigs_flt	= "$opt{res}/$hostpref.contigs.fa";
		
			filter_host_contigs(
				contigs 		=> $contigs,
				contigs_pass=> $contigs_pass,
				contigs_flt	=> $contigs_flt,
				hostdb  		=> $hostdb->{db},
				res 			=> $opt{res},
				log 			=> "$opt{logs}/$opt{sample}/filter.host.contigs.log",
				numth		=> $opt{numth},
				bitscore		=> $opt{min_contig2hostgen_score},
				tmpdir 		=> $opt{tmpdir} || undef);
			
			system("mv $contigs_pass $contigs");
		}
	}
}



sub pipe_assemble{
	print STDERR "\n\n# ASSEMBLE\n\n";	
	
	my $subid	= "pipe_assemble()";
	# in:
	my %opt 		= %{shift()};
	my $r1 		= (-e "$opt{res}/reads/read1.trim.hflt.fq.gz") ? "$opt{res}/reads/read1.trim.hflt.fq.gz" : "$opt{res}/reads/read1.trim.fq.gz";
	my $r2		= (-e "$opt{res}/reads/read2.trim.hflt.fq.gz") ? "$opt{res}/reads/read2.trim.hflt.fq.gz" : "$opt{res}/reads/read2.trim.fq.gz";
	
	# out:
	my $contigs 			= "$opt{res}/contigs.fa";
	my $contigs_host		= "$opt{res}/contigs.host.fa";
	my $contigs_info		= "$contigs.info.tsv";
	my $orfs_nt			= "$opt{res}/contigs.orfs.nt.fa";
	my $orfs_aa			= "$opt{res}/contigs.orfs.aa.fa";
	my $assembler_out 	= "$opt{res}/assembler_out";
	
	my $log 		= "$opt{logs}/$opt{sample}/assemble.log";
	
	# params:
	my $par_megahit = defined($opt{par_megahit})? $opt{par_megahit} : "";
	my $par_spades  = defined($opt{par_spades}) ? $opt{par_spades}  : "";
	
	system("rm -f $log");
	system("touch $log");	
	
	if($opt{'ass'} eq 'megahit'){
		system_call("rm -fR $assembler_out", $opt{'v'}); # Megahit will complain if that dir exists
		if( $opt{'se'} ){
		system_call("megahit -t $opt{'numth'} $par_megahit --read $r1 --out-dir $assembler_out &>> $log", $opt{'v'});
		}
		else{
		system_call("megahit -t $opt{'numth'} $par_megahit -1 $r1 -2 $r2 --out-dir $assembler_out &>> $log", $opt{'v'});
		}
		
		system_call("seqkit seq -n $assembler_out/final.contigs.fa | ".
					"cut -d' ' --output-delimiter=\$'\\t' -f1,3,4 | ".
					"sed 's/multi=\\|len=//g' | ".
					"sed 's/_/./g' | ".
					"csvtk add-header -tn seqid,coverage,length | ".
					"csvtk sort -t -k seqid:N 1> $contigs_info 2>> $log", $opt{'v'});
		
		system_call("seqkit seq --only-id $assembler_out/final.contigs.fa | ".
					"seqkit replace -p '_' -r '.' | ".
					"seqkit sort -Nw0 1> $contigs 2>> $log", $opt{'v'});
	
	}
	elsif($opt{'ass'} eq 'spades'){
		system_call("rm -fR $assembler_out", $opt{'v'});
		if( $opt{'se'} ){
		system_call("spades.py -t $opt{'numth'} $par_spades -1 $r1 -o $assembler_out &>> $log", $opt{'v'});
		}
		else{
		system_call("spades.py -t $opt{'numth'} $par_spades -1 $r1 -2 $r2 -o $assembler_out &>> $log", $opt{'v'});
		}
		
		system_call("seqkit seq -n $assembler_out/scaffolds.fasta | ".
					"sed 's/length_\\|cov_//gi' | ".
					"sed 's/NODE_/scaffold./' | ".
					"cut -d'_' --output-delimiter=\$'\\t' -f1,2,3 | ".
					"csvtk add-header -tn seqid,length,coverage | ".
					"csvtk sort -t -k seqid:N 1> $contigs_info 2>> $log", $opt{'v'});
		
		system_call("seqkit replace -p '^[A-Za-z]+_([0-9]+).*' -r 'scaffold.\$1' $assembler_out/scaffolds.fasta | ".
					"seqkit sort -Nw0 1> $contigs 2>> $log", $opt{'v'});
	}
	else{
		die "ERROR: invalid --ass $opt{'ass'}";
	}
	
	if( defined($opt{'min_contig_length'})){
		system_call("cat $contigs | seqkit seq -w90 -m $opt{'min_contig_length'} 1> $contigs.tmp", $opt{'v'});
		system_call("mv $contigs.tmp $contigs", $opt{'v'});
	}
	
	if( $opt{'clean'} ){
		system_call("rm -fr $assembler_out", $opt{'v'});
	}
	
	if( nlines($contigs)<2){
		print STDERR "\tWARNING: $subid: exporting empty assembly\n\n";
		return;
	}
	
	detect_orfs(seqs=>$contigs, orfs_nt=>$orfs_nt, orfs_aa=>$orfs_aa, 
				orf_finder=>$opt{gen}, min_orf_length=>$opt{min_orf_length});
}

sub realign_reads_contigs{
	print STDERR "\n# REALIGN READS TO CONTIGS\n\n";
	my $subid	= "realign_reads_contigs()";
	
	# in:
	my %opt 		= %{shift()};
	my $r1 		= (-e "$opt{res}/reads/read1.trim.hflt.fq.gz") ? "$opt{res}/reads/read1.trim.hflt.fq.gz" : "$opt{res}/reads/read1.trim.fq.gz";
	my $r2		= (-e "$opt{res}/reads/read2.trim.hflt.fq.gz") ? "$opt{res}/reads/read2.trim.hflt.fq.gz" : "$opt{res}/reads/read2.trim.fq.gz";	
	my $contigs = "$opt{'res'}/contigs.fa";

	# tmp:
	my $sam 		= "$opt{res}/contigs.bwa.sam";
	my $bam 		= "$opt{res}/contigs.top.bam";
	
	# out:
	my $idxstats= "$opt{res}/contigs.idxstats";
	my $idmap   = "$opt{res}/readid_contigid.tsv";
	my $log 		= "$opt{logs}/$opt{sample}/realign_reads.log";
	
	# check input
	die "ERROR: $subid: missing forward reads\n" if(!(-e $r1));
	die "ERROR: $subid: missing forward reads\n" if(!(-e $r2));
	die "ERROR: $subid: missing contigs\n" if(!(-e $contigs));
	if(nlines($contigs)<2){
		print STDERR "\tWARNING: $subid: empty contigs file: exporting empty idxstats/idmap files\n\n";
		system("touch $idxstats $idmap");
		return;
	}
	
	
	system("rm -f $log; touch $log");
	system_call("bwa index $contigs &>> $log");
	if( $opt{se} ){
		system_call("bwa mem -t $opt{numth} $contigs $r1 1> $sam 2>> $log");
	}
	else{
		system_call("bwa mem -t $opt{numth} $contigs $r1 $r2 1> $sam 2>> $log");
	}
	
	system_call("sambamba view -t $opt{numth} -S -h -F \"not(unmapped) and mapping_quality>=$opt{min_read2contig_mapq} and [AS]>=$opt{min_read2contig_score}\" $sam 1> $sam.tmp 2>> $log");
	system_call("mv $sam.tmp $sam");
	system_call("samtools sort -@ $opt{numth} -n -T $opt{tmpdir} $sam | ".					# sort by read name
				"samtools view -F4 -h 1> $sam.tmp 2>> $log");
	filter_tophits_SAM(in=>"$sam.tmp", out=>$sam);
	system("rm -f $sam.tmp");
	system_call("samtools sort -@ $opt{numth} -T $opt{tmpdir} -o $bam $sam 2>> $log");		# sort by chromosome posit
	system_call("samtools index -@ $opt{numth} $bam 2>> $log");
	system_call("samtools idxstats $bam 1> $idxstats 2>> $log");
	system_call("samtools view  $bam | cut -f1,3 1> $idmap 2>> $log");
	
	if($opt{clean}){
		system_call("rm -f $sam $bam $contigs.amb $contigs.ann $contigs.bwt $contigs.pac $contigs.sa");
	}
}


# USAGE:
# generate_reports( annot=>$annot, opt=>\%opt);
#
sub generate_reports{
	print STDERR "\n# GENERATE REPORTS\n\n";
	my $subid					= "generate_reports()";
	
	# in:
	my %opt						= %{ shift() };
	my $annot					= (-e "$opt{res}/annot2.tsv" && nlines("$opt{res}/annot2.tsv")>1) ? "$opt{res}/annot2.tsv" : "$opt{res}/annot1.tsv";
	my $r1_flt					= "$opt{res}/reads/read1.trim.fq.gz";
	my $r1_hgflt					= "$opt{res}/reads/read1.trim.hflt.fq.gz";
	my $contigs 					= "$opt{res}/contigs.fa";
	my $contigs_stats			= "$opt{res}/contigs.idxstats";
	# dependencies in global variables:
	# $VIRUS_FAMILY_HOST
	# $VIRUS_GENUS_HOST
	# $VIRUS_FAMILY_GENCOMP
	
	# tmp:
	my $annot_lca				= "$annot.lca.tmp";
	my $annot_nucl				= "$annot.nucl.tmp";
	my $annot_prot 				= "$annot.prot.tmp";
	my $annot_hmm				= "$annot.hmm.tmp";
	my $annot_union				= "$annot.union.tmp";
	my $bphage_map				= "$annot.bphage.tmp";
	my $readn_taxid 				= "$opt{res}/readn_taxid.tmp";
	
	# out:
	
	my $annot_table				= "$opt{res}/annot_table.tsv";
	my $annot_excel 				= "$opt{res}/annot_table.xlsx";
	my $abund_table 				= "$opt{res}/abund_table.tsv";
	my $abund_excel 				= "$opt{res}/abund_table.xlsx";		
	my $taxprofile				= "$opt{res}/taxprofile.txt";
	my $krona_data				= "$opt{res}/reports/krona.data.txt";
	my $krona_graph 				= "$opt{res}/reports/krona.report.html";
	my $contigs_dir 				= "$opt{res}/contigs";
	my $log 						= "$opt{logs}/$opt{sample}/generate_reports.log";
	
	# params:
	my $threads 					= ($opt{numth} < 8) ? $opt{numth} : 8; # min(8,numth)
	my $par_abund_table 			= "-h -w $opt{wmodel} --conttail $opt{tail_contig}";
	my $par_taxprofile			= "--sample $opt{sample} --tail $opt{tail}";
	my $par_kronagraph			= "-s $opt{sample} -t $opt{tail}";
	my $par_hgtaxid 				= ($opt{hgtaxid}) ? $opt{hgtaxid} : "";
	my $min_qcov_annot			= defined($opt{min_qcov_annot})? $opt{min_qcov_annot} : 0;
	my $min_qcov_abund			= defined($opt{min_qcov_abund})? $opt{min_qcov_abund} : 0;
	my $min_psearch_bits_abund	= defined($opt{min_psearch_bits_abund}) ? $opt{min_psearch_bits_abund} : 0;
	my $min_nsearch_bits_abund	= defined($opt{min_nsearch_bits_abund}) ? $opt{min_nsearch_bits_abund} : 0;
	my $toptaxrank				= defined($opt{toptaxrank}) ? $opt{toptaxrank}: "division"; # use division after 2025/Marck NCBI taxonomy update
	my $taxranks					= "species,genus,family,$toptaxrank";
	my $taxgroups				= $opt{taxgroups} || "all";

	
	
	if( $opt{'hgtaxid'} && (-e $r1_flt && -e $r1_hgflt ) ){
		my $flt 					= (`$opt{'gzip'} -d -c  $r1_flt  | wc -l`)/4;
		my $hgflt				= (`$opt{'gzip'} -d -c  $r1_hgflt | wc -l`)/4;
		my $hg_readn 			= $flt - $hgflt;
		
		$par_abund_table			= "$par_abund_table --hgabund $hg_readn --hgtaxid $opt{hgtaxid}";
		$par_taxprofile 			= "$par_taxprofile --hgtaxid $opt{hgtaxid}";
	}
	
	system("rm -f $log");
	system("touch $log");

	# UPDATE TAXONOMY
	if( $opt{taxonomy}->{update} ){
		update_taxonomy(\%opt);
	}
	
	# INPUT CHECKS
	if( !(-e $annot) || nlines($annot) < 2){
		print STDERR "\tWARNING: $subid: No reporting due to missing/empty annotation file: $annot\n";
		return();
	}
	
	# CONVERT POSSIBLE LISTS OF STAXIDS TO LCA
	my $taxids_col 	= colind($annot,'staxid');
	system_call("head -n1 $annot | csvtk rename -tf staxid -n staxids | tr '\\n' '\\t' 1> $annot_lca");
	system("echo 'staxid' 1>> $annot_lca");
	system_call("csvtk filter2 -tj $opt{numth} -f '\$staxid!=\"\"' $annot | ".
				"csvtk del-header -t | ".
				"taxonkit lca --data-dir $opt{taxonomy}->{db} -j $opt{numth} -i $taxids_col -s \";\" 1>> $annot_lca 2>> $log");
	
	
	
	# GENERATE CONTIG ANNOTATION TABLES
		# select nucl-annotations and select top score for each uniq contig+staxid
	system_call("cat $annot_lca | ".
				"csvtk filter2 -tj $opt{numth} -f '\$dbtype==\"nucl\" && \$qcov>=$min_qcov_annot' 1> $annot_nucl 2>> $log");
	if(nlines($annot_nucl)>1){
	system_call("cat $annot_nucl | ".
				"csvtk sort -tj $opt{numth} -k qseqid:N -k staxid:N -k bitscore:nr  | ".
				"csvtk uniq -tj $opt{numth} -f qseqid,staxid 1> $annot_nucl.tmp  2>> $log" );
	system_call("mv $annot_nucl.tmp $annot_nucl" );
	}
		
		# select prot-annotations and select top score for each uniq contig+orf+staxid
	system_call("cat $annot_lca | ".
				"csvtk filter2 -tj $opt{numth} -f '\$dbtype==\"prot\" && \$qcov>=$min_qcov_annot' 1> $annot_prot 2>> $log" );
	if(nlines($annot_prot)>1){
		# select topscore for each contig-orf-staxid
	system_call("cat $annot_prot | ".
				"csvtk sort -tj $opt{numth} -k qseqid:N -k orf:N -k staxid:N -k bitscore:nr  | ".
				"csvtk uniq -tj $opt{numth} -f qseqid,orf,staxid 1> $annot_prot.tmp  2>> $log" );
	system_call("mv $annot_prot.tmp $annot_prot" );
	}
	
		# select hmm-annotations, if any, keep all hits
	system_call("cat $annot_lca | ".
				"csvtk filter2 -tj $opt{numth} -f '\$dbtype==\"HMM\" ' 1> $annot_hmm 2>> $log" );
	
		# join nucl, prot and hmm-annotations, sort and print to annot_table
		# rename fields for final table: 'qseqid' > 'contig', 'qseqlen' > 'clen'
	system_call("cat $annot_nucl <(tail -n+2 $annot_prot) <(tail -n+2 $annot_hmm) 1> $annot_table  2>> $log" );
	if(nlines($annot_table) > 1){
	system_call("cat $annot_table | ".
				"csvtk rename -tf qseqid,qseqlen -n contig,clen | ".
				"csvtk sort -tj $opt{numth} -k contig:N -k staxid:N -k dbtype:N 1> $annot_table.tmp 2>> $log" );
	system_call("mv $annot_table.tmp $annot_table" );
	}
		# check that we still have annotations after all the filtering
	if(nlines($annot_table) < 2){
		print STDERR "\n\tWARNING: $subid: No reporting due to empty annotation table: $annot_table\n\n";
		return();
	}	
		# add taxonomy
	my $taxid_col = colind($annot_table,'staxid');
	system("head -n1 $annot_table | tr '\\n' '\\t' 1> $annot_table.tmp");
	system("echo 'species\tgenus\tfamily\torder\tclass\tkingdom' 1>> $annot_table.tmp");
	if(nlines($annot_table) > 1){
	system_call("csvtk filter2 -tj $opt{numth} -f '\$staxid!=\"\"' $annot_table | ".
				"csvtk del-header -t | ".
			  	"taxonkit reformat --data-dir $opt{taxonomy}->{db} -j $opt{numth} -I $taxid_col -f '{s}\\t{g}\\t{f}\\t{o}\\t{c}\\t{k}' -r NA -R NA ".
			  	"1>> $annot_table.tmp 2>> $log");
	}
	system("mv $annot_table.tmp $annot_table");
		# add division-field, in case this has not been added during annotation
	add_division_field(	annot=>$annot_table, 
						taxonomy=>$opt{taxonomy}->{db}, numth=>$opt{numth}, overwrite=>1);
		# ADD bphage-field, redundant if bphage-field already added to $annot, but run this in case it has not
	add_bphage_field(	annot=>$annot_table,
						phfilter=>$BPHAGE_FILTER,
						taxonomy=>$opt{taxonomy}->{db},
						log=>$log,numth=>$opt{numth},overwrite=>1);
						
		# ADD host.source
	my $FIELD_HOST	= "host.source";
	system_call("csvtk join -tj $opt{numth} -f 'genus;genus' -L --na 'NA' $annot_table $VIRUS_GENUS_HOST | ".
				"csvtk rename -tf $FIELD_HOST -n gen_host 1> $annot_table.tmp1 2>> $log");
	system_call("csvtk join -tj $opt{numth} -f 'family;family' -L --na 'NA' $annot_table.tmp1 $VIRUS_FAMILY_HOST | ".
				"csvtk rename -tf $FIELD_HOST -n fam_host | ".
				"csvtk mutate2 -te '\$gen_host!=\"NA\" ? \$gen_host : \$fam_host' -n $FIELD_HOST | ".
				"csvtk cut -tlf '-gen_host,-fam_host'  1> $annot_table.tmp2 2>> $log" );
	system_call("rm -f $annot_table.tmp1");
	system_call("mv $annot_table.tmp2 $annot_table" );
		# ADD genome.composition
	my $FIELD_GENCOMP	= "genome.composition";
	system_call("csvtk join -tj $opt{numth} -f 'family;family' -L --na 'NA' $annot_table $VIRUS_FAMILY_GENCOMP 1> $annot_table.tmp 2>> $log");
	system_call("mv $annot_table.tmp $annot_table");

		# convert $annot_table to excel file
	system_call("$opt{call_R} $R_scripts/print_annot_table.R  $annot_table $annot_excel $toptaxrank $taxgroups 2>> $log");
	
		# SORT CONTIGS TO DIRS USING TAXONOMY CLASSIFICATION
	system_call("perl $perl_scripts/sort_contigs_bytaxa_v3.pl -c $contigs -a $annot_table --res $opt{res}/contigs --toptaxrank $toptaxrank -v &>> $log");	
	
	
	# GENERATE ABUNDANCE TABLES
	if( !(-e $contigs_stats)){
		print STDERR "\n\tWARNING: Unable to estimate abundancies: no $contigs_stats file\n";
		print STDERR "\tYou can generate $contigs_stats by running --pipe rea \n";
		return();
	}
		# filter nucl-annotations by qcov+bitscore, select qseqid,staxid,bitscore
	system_call("cat $annot_lca | ".
				"csvtk filter2 -tj $opt{numth} -f '\$dbtype==\"nucl\" && \$qcov>=$min_qcov_abund && \$bitscore>=$min_nsearch_bits_abund' | ".
				"csvtk cut -tf qseqid,staxid,bitscore,qseqlen  1> $annot_union 2>> $log" );
		# filter prot-annotations by qcov+bitscore, sum bitscore over orfs, select qseqid,staxid,bitscore
	system_call("cat $annot_lca | ".
				"csvtk filter2 -tj $opt{numth} -f '\$dbtype==\"prot\" && \$qcov>=$min_qcov_abund && \$bitscore>=$min_psearch_bits_abund' | ".
				"csvtk summary -i -tj $opt{numth}  -f qseqlen:first -f bitscore:sum -g qseqid,staxid | ".
				"csvtk rename -tf qseqlen:first,bitscore:sum -n qseqlen,bitscore | ".
				"csvtk cut -tf qseqid,staxid,bitscore,qseqlen | ".
				"csvtk del-header 1>> $annot_union 2>> $log");
	if( nlines($annot_union)<2 ){
		print STDERR "\n\tWARNING: $subid: Unable to estimate abundancies: no annotations passing threshold\n";
		system("rm -f $annot_nucl $annot_prot $annot_union $bphage_map $readn_taxid");
		return();
	}
		# select top-scoring annotation for each contig (retain ties)
	system_call("cat $annot_union | ".
				"csvtk sort -tj $opt{numth} -k qseqid:N -k bitscore:nr 1> $annot_union.tmp 2>> $log" );
	system_call("mv $annot_union.tmp $annot_union");
	filter_tophits(dbhits=>$annot_union,dbhits_flt=>"$annot_union.tmp", qcol=>'qseqid', bitscol=>'bitscore', retain_ties=>1);
	system_call("mv $annot_union.tmp $annot_union");
	
		# rename 'qseqid' to 'contig' for final annotation-table
	system_call("csvtk rename -tf qseqid,qseqlen -n contig,contig.len $annot_union 1> $annot_union.tmp 2>> $log");
	system("mv $annot_union.tmp $annot_union");
	
	if( nlines($annot_union)<2 ){
		print STDERR "\n\tWARNING: $subid: Unable to estimate abundancies: no annotations passign threshold\n";
		system("rm -f $annot_nucl $annot_prot $annot_union $bphage_map $readn_taxid");
		return();
	}

		# estimate abundancies:
	system_call("perl $perl_scripts/get_abund_table.pl $par_abund_table $annot_union $contigs_stats 1> $readn_taxid 2> $log" );
	
	if( nlines($readn_taxid)<2 ){
		print STDERR "\n\tWARNING: $subid: printing empty abundance table\n";
		system_call("touch $abund_table");
		system("rm -f $annot_nucl $annot_prot $annot_union $bphage_map $readn_taxid");
		return();
	}	
	
		# add taxonomy to abundancies and print to $abund_table:
	system_call("csvtk filter2 -t -f '\$taxid!=\"\"' $readn_taxid | ".
				"csvtk del-header -t  | ".
			  	"taxonkit reformat --data-dir $opt{taxonomy}->{db} -j $threads -I 4 -t -f '{s}\\t{g}\\t{f}\\t{o}\\t{c}\\t{p}\\t{k}' -r NA -R NA | ".
			  "csvtk add-header -I -t -n readn,contign,assembly.size,taxid,species,genus,family,order,class,phylum,superkingdom,species_id,genus_id,family_id,order_id,class_id,phylum_id,superkingdom_id 1> $abund_table 2>> $log");
		# add division-field
	add_division_field(	annot=>$abund_table, 
						taxonomy=>$opt{taxonomy}->{db}, numth=>$opt{numth}, overwrite=>1);
		# ADD bphage-field
	add_bphage_field(	annot=>$abund_table,
						phfilter=>$BPHAGE_FILTER,
						taxonomy=>$opt{taxonomy}->{db},
						log=>$log,numth=>$opt{numth},overwrite=>1);
				
		# ADD host.source
	system_call("csvtk join -tj $opt{numth} -f 'genus;genus' -L --na 'NA' $abund_table $VIRUS_GENUS_HOST | ".
				"csvtk rename -tf $FIELD_HOST -n gen_host 1> $abund_table.tmp1 2>> $log");
	system_call("csvtk join -tj $opt{numth} -f 'family;family' -L --na 'NA' $abund_table.tmp1 $VIRUS_FAMILY_HOST | ".
				"csvtk rename -tf $FIELD_HOST -n fam_host | ".
				"csvtk mutate2 -te '\$gen_host!=\"NA\" ? \$gen_host : \$fam_host' -n $FIELD_HOST | ".
				"csvtk cut -tlf '-gen_host,-fam_host'  1> $abund_table.tmp2 2>> $log" );
	system_call("rm -f $abund_table.tmp1");
	system_call("mv $abund_table.tmp2 $abund_table" );
		# ADD genome.composition
	system_call("csvtk join -tj $opt{numth} -f 'family;family' -L --na 'NA' $abund_table $VIRUS_FAMILY_GENCOMP 1> $abund_table.tmp 2>> $log");
	system_call("mv $abund_table.tmp $abund_table");
	
		# convert $abund_table to excel file:
	system_call("$opt{'call_R'} $R_scripts/print_abund_table.R ".
				"$abund_table $abund_excel $taxranks $taxgroups $opt{tail} 2>> $log");

	# CREATE KRONA GRAPH	
	system_call("perl $perl_scripts/abundtable2krona.pl $par_kronagraph $abund_table  --toptaxrank $toptaxrank 1> $krona_data 2>> $log");
	system_call("ktImportText $krona_data -o $krona_graph -u \"http://krona.sourceforge.net\" 2>> $log");

	# TAXONOMIC PROFILE
	#system_call("perl $perl_scripts/abundtable2taxprofile.pl $par_taxprofile $abund_table 1> $taxprofile 2>> $log");
	
	# REMOVE TMP
	#system("rm -f $annot_nucl $annot_prot $annot_hmm $annot_union $bphage_map $readn_taxid");

}



# USAGE:
# generate_igv_reports(
#		annot=>$annot_table,
#		contigs=>$contigs_vi,
#		log=>$log_igv,
#		opt=>\%opt, 
#		min_clen_sum => 0, 
#		min_alen_sum => 0,
#		max_refgen => 10);
#
sub generate_refgen_reports{
	print STDERR "\n# GENERATE REFGEN REPORTS\n\n";
	
	# IN:
	my %opt 					= %{shift()};
	my $annot_tsv			= "$opt{res}/annot_table.tsv";
	my $contigs_fa			= (-e "$opt{res}/contigs.ann2.vi.fa")? "$opt{res}/contigs.ann2.vi.fa" : 
								((-e "$opt{res}/contigs.ann1.vi.fa")? "$opt{res}/contigs.ann1.vi.fa" : "$opt{res}/contigs.fa");
	# out:
	my $resdir				= "$opt{res}/reports";
	my $report_html			= "$resdir/refgen.report.html";
	my $log 					= "$resdir/refgen.report.log";
	# tmp:
	my $annot_ta_tsv			= "$opt{tmpdir}/annot.ta.tsv";
	my $contigs_ta_fa		= "$opt{tmpdir}/contigs.ta.fa";
	my $refgens_txt			= "$opt{tmpdir}/refgens.txt";
	my $dataset_zip			= "$opt{tmpdir}/ncbi_dataset.zip";
	my $dataset_genomic_fa	= "$opt{tmpdir}/ncbi_dataset/data/genomic.fna";
	my $refgen_metainfo 		= "$opt{tmpdir}/refgen.metainfo.tsv";
	my $refgen_annot			= "$opt{tmpdir}/refgen.annot.tsv";
	# params:
	my $min_clen_sum 		= $opt{refgenrep_min_clen_sum} || 1000;
	my $min_alen_sum 		= $opt{refgenrep_min_alen_sum} || 0;
	my $max_refgen			= $opt{refgenrep_max_refgen} || 10;
	my $toptaxrank			= defined($opt{toptaxrank})? $opt{toptaxrank}: 'division';
	my $target				= "Viruses";
	my $dbtype				= "nucl";
	my $par_minimap_s		= 10;
	
	# INPUT CHECKS
	if( !(-e $annot_tsv)  || nlines($annot_tsv)<2 ){
		print STDERR "\n\tNo reporting due to missing/empty annotation file: $annot_tsv\n\n";
		return();
	}
	
	system("rm -f $log");
	system("touch $log");
	system("rm -fr $resdir/data");
	
		# Collect vi-annotations
	system_call("csvtk filter2 -tf '\$dbtype==\"$dbtype\" && \$$toptaxrank==\"$target\" && \$bphage==\"no\"' $annot_tsv 1> $annot_ta_tsv  2>> $log");
		# parse annotation 
	my %sp2clen			= read_tsv2kvahash($annot_ta_tsv,"species","clen");
	my %sp2alen			= read_tsv2kvahash($annot_ta_tsv,"species","alen");
	my %sp2contig		= read_tsv2kuvahash($annot_ta_tsv,"species","contig");
	my %sp2refgenid		= read_tsv2kuvahash($annot_ta_tsv,"species","sseqid");
	#my %sp2staxid 		= read_tsv2kuvahash($annot_ta_tsv,"species","staxid");
	my %sp2genus			= read_tsv2hash($annot_ta_tsv,"species","genus");
	#my %sp2family		= read_tsv2hash($annot_ta_tsv,"species","family");
	
		# Filter target species list
	my %splist			= ();
	foreach my $sp(keys %sp2clen){
		my $clen_sum		= sum($sp2clen{$sp});
		my $alen_sum		= sum($sp2alen{$sp});
		if($clen_sum >= $min_clen_sum && $alen_sum >= $min_alen_sum){
			$splist{$sp} = 1;
		}
	}
	# convert "gb/ref/emb|seqid|" returned by blastn to "seqid"
	# limit to $max_refgen
	foreach my $sp(keys %splist){
		my @seqids1		= @{$sp2refgenid{$sp}};
		my @seqids2		= ();
		foreach my $id(@seqids1){
			$id =~ s/\w+\|([A-Za-z0-9_\.]+)\|/$1/g;
			push(@seqids2,$id);
		}
		if(scalar(@seqids2) > $max_refgen){
			@seqids2 = @seqids2[0..($max_refgen-1)];
		}
		$sp2refgenid{$sp} = \@seqids2;
	}
	
	my %contigs 			= %{readfasta($contigs_fa)};
	my @summary_table	= ();
	
	# Create refgen/igv reports
	foreach my $sp( sort { $sp2genus{$a} cmp $sp2genus{$b} } keys %splist){
		
		# Collect target-contigs to fasta
		my @ids		= @{$sp2contig{$sp}};
		open(OUT,">$contigs_ta_fa") or die "Can\'t open $contigs_ta_fa: $!\n";
		foreach my $seqid( @ids ){
			if( !defined($contigs{$seqid})){
				print STDERR "WARNING: no seqid=$seqid found in $contigs_ta_fa: SKIPPING\n";
				next;
			}
			print OUT ">$seqid\n";
			print OUT "$contigs{$seqid}\n";
		}
		close(OUT);
		
		# retrieve genome+annotation data for refgens
		my @refgenids 	= @{$sp2refgenid{$sp}};
		write_array2file($refgens_txt,\@refgenids);
			# try download
		my $call 		= "datasets download virus genome accession ".
						"--include genome,annotation --no-progressbar ".
						"--inputfile $refgens_txt --filename $dataset_zip 2>> $log";
		if($opt{v}){ print STDERR "\t$call\n"}
		if(system($call)!= 0){
			print STDERR "\tWARNING: failed to download annotation for $sp: skipping\n";
			next;
		}				
		system_call("unzip -q -o -d $opt{tmpdir} $dataset_zip 2>> $log");
			# parse metadata
		system_call("dataformat tsv virus-genome --force --package $dataset_zip --fields accession,host-name,geo-location,isolate-lineage 1> $refgen_metainfo 2>> $log");
		my %refgen_minfo = read_tsv2hashtable($refgen_metainfo,"Accession");
		
			# parse gene annotation data to bed file
		system_call("dataformat tsv virus-annotation --force --package $dataset_zip --fields accession,gene-name,gene-genomic-range-start,gene-genomic-range-stop,gene-cds-nuc-fasta-range-start,gene-cds-nuc-fasta-range-stop --elide-header 1>> $refgen_annot 2>> $log");
		
		my $headers_tmp	= "#accession\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb";
		my $spdir 		= $sp; $spdir =~ s/\s+/_/g;
		my $bed			= "$resdir/data/$spdir/refgen.all.bed";
		my @colors		= ('#A202FF','#2072B2','#1260CC','#OOABFF','#00E8FF');	# iterate blue palet
		system("mkdir -p $resdir/data/$spdir");
		
		open(OUT,">$bed") or die "Can\'t open $bed: $!\n";
		print OUT "$headers_tmp\n";
		open(IN, "<$refgen_annot") or die "Can\'t open $refgen_annot: $!\n";
		my $ind	= 0;
		while(my $l=<IN>){
			chomp($l);
			my ($acc,$name,$gene_start,$gene_stop,$cds_start,$cds_stop) = split(/\t/,$l,-1);
			my $start	= ($gene_start ne '')? $gene_start: ($cds_start ne '')? $cds_start : -1;
			my $stop		= ($gene_stop ne '')? $gene_stop: ($cds_stop ne '')? $cds_stop : -1;
			if($start <0 || $stop <0){next;}
			print OUT join("\t",$acc,$start,$stop,$name,0,'.',0,0,$colors[$ind % (scalar(@colors))]),"\n";
			$ind++;
		}
		close(IN);close(OUT);
		
				
		# run minimap > bam > bam.bai
		my $sam 			= "$resdir/data/$spdir/refgen.all.sam";
		my $bam 			= "$resdir/data/$spdir/refgen.all.bam";
		my $coverage_tsv	= "$resdir/data/$spdir/refgen.all.coverage.tsv";
		
		system_call("minimap2 -t $opt{numth} --cs -s $par_minimap_s -a $dataset_genomic_fa $contigs_ta_fa 1> $sam 2>> $log");
		system_call("samtools sort -@ $opt{numth} -o $bam $sam  2>> $log");
		system_call("samtools index $bam 2>> $log");
		system_call("samtools coverage --f 0 $bam 1> $coverage_tsv  2>> $log");
		my %coverage		= read_tsv2hashtable("$coverage_tsv","#rname");
		
		# create the IGV html
		generate_igv_html(
			opt=> \%opt, 
			species=>$sp, 
			refgens=>$dataset_genomic_fa, 
			refgens_metainfo=> \%refgen_minfo,
		 	contigs=>$contigs_ta_fa,
		 	bam=>$bam,
		 	bed=>$bed,
		 	resdir=>"$resdir/data/$spdir", 
		 	split_byrefgen=>1,
		 	use_data_uri=>1, log=>$log);
		
		# Filling @summary_table
		# @headers	= ("Genus","Species","clen.sum","alen.sum","Refgen.len","Refgen.cov","Refgen.igv");
		my $genus		= $sp2genus{$sp};
		my $clen_sum		= format_int(sum( $sp2clen{$sp} ),' ');
		my $alen_sum		= format_int(sum( $sp2alen{$sp} ),' ');
			@refgenids 	= @{$sp2refgenid{$sp}};
		foreach my $id(@refgenids){
			my $refgen_len		= defined($coverage{$id}->{"endpos"})? $coverage{$id}->{"endpos"} : "NA";
			my $refgen_cov		= defined($coverage{$id}->{"coverage"})? $coverage{$id}->{"coverage"} : "NA";
			if($refgen_len ne "NA"){		$refgen_len	= format_int($refgen_len,' '); }
			if($refgen_cov ne "NA"){		$refgen_cov	= sprintf("%.2f",$refgen_cov); }
			my $igv_path_rel		= "data/$spdir/$id.igv.html";
			my $igv_path			= "$resdir/$igv_path_rel";
			my $refgen_igv		= (-e $igv_path)? "<a href='$igv_path_rel'>$id</a>":"$id";
			my @summary_row 		= ($genus, $sp, $clen_sum, $alen_sum,$refgen_len,$refgen_cov,$refgen_igv);
			push(@summary_table,\@summary_row);
		}
	}
	
	# sort @summary_table: by $refgen_cov > by $species
	@summary_table	= 
		sort { ($a->[5] eq "NA" && $b->[5] eq "NA")? 0:
				($a->[5] eq "NA") ? +1 :
				($b->[5] eq "NA") ? -1 :
				($b->[5] <=> $a->[5]) } @summary_table;
	@summary_table	= sort{ $a->[1] cmp $b->[1]} @summary_table; 
	my @headers		= ("Genus ","Species","clen.sum","alen.sum","Refgen.len","Refgen.cov","Refgen.igv");
	unshift(@summary_table, \@headers);
	
	# create refgen-report-html
	open(OUT,">$report_html") or die "ERROR: failed to open: $report_html\n";
	
	if( !$opt{webmode} ){	
		print OUT "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n";
		print OUT "<html>\n";
		print OUT "<head>\n";
		print OUT "<title>Lazypipe RefGen Report</title>\n";
		print OUT "<meta charset=\"utf-8\"/>\n";
		print OUT "<style>\n";
		print_file("$perl_scripts/include/reports.css",\*OUT);
		print OUT "</style>\n";
		print OUT "<script>\n";
		print_file("$perl_scripts/include/reports.js",\*OUT);
		print OUT "</script>\n";
		print OUT "</head>\n";
		print OUT "<body>\n";
	}

	print OUT "<h3> Lazypipe RefGen Report</h3>\n";
	print OUT "<p>Sample: $opt{sample}</p>\n";
	
	# column attributes, printed to hd-tag
	my $table_attrs = "class='sortable'";
	my @th_attrs	= ();
	push(@th_attrs, "class='mixed' type='string'");		# Genus: string
	push(@th_attrs, "class='asc' type='string'");	# Species: string
	push(@th_attrs, "class='mixed' type='int'");		# clen.sum: int
	push(@th_attrs, "class='mixed' type='int'");		# alen.sum: int
	push(@th_attrs, "class='mixed' type='int'");		# refgen.len: int
	push(@th_attrs, "class='desc' type='num'");		# refgen.cov: num
	push(@th_attrs, "class='mixed' type='element'");	# refgen.igv: <a>-element
	
	write_table_html(table=>\@summary_table, table_attrs=>$table_attrs, th_attrs=>\@th_attrs, fh=>\*OUT);
	
	print OUT
	"<script>\n".
	"	init_sortable_tables();\n".
	"</script>\n";
	
	if( !$opt{webmode} ){
		print OUT "</body>\n";
		print OUT "</html>\n";
	}
	close(OUT);
}

sub pipe_annotation_ictv{
	my $subid			= "pipe_annotation_ictv()";
	
	print STDERR "# STARTING $subid\n" if ($VERBAL);
	
	# IN
	my %opt 				= %{shift()};
	my $contigs_ta		= undef;
	$contigs_ta			= "$opt{res}/contigs.ann1.vi.fa" if(-e "$opt{res}/contigs.ann1.vi.fa");
	$contigs_ta			= "$opt{res}/contigs.ann1.vi.fa" if(-e "$opt{res}/contigs.ann2.vi.fa");
	if(!defined($contigs_ta) || nlines($contigs_ta)<2){
		print STDERR "\tWARNING: $subid: missing/empty contigs fasta: no annotations will be produced\n";
		return;
	}
	my $contigs_info		= "$opt{res}/contigs.fa.info.tsv";
		die "ERROR: $subid: missing file: $contigs_info" if(!(-e $contigs_info));
	my $ictv_blastn_db	= $opt{'ann.databases'}->{"blastn.ictv"}->{db} || 
		die "ERROR: $subid: no ICTV BLASTP database. Expecting: ann.databases:blastn.ictv: blastdb\n";
	my $ictv_blastp_db	= $opt{'ann.databases'}->{"blastp.ictv"}->{db} || 
		die "ERROR: $subid: no ICTV BLASTP database. Expecting: ann.databases:blastp.ictv: blastdb\n";
	my $ictv_tsv			= $opt{'ICTV.VMR'}->{db} || 
		die "ERROR: $subid: no ICTV VMR table. Expecting: ICTV.VMR: db: ICTV.VMR.tsv";
	my $contigs_idxstats	= "$opt{res}/contigs.idxstats";
		die "ERROR: $subid: missing file: $contigs_idxstats" if(!(-e $contigs_idxstats));
	
	# OUT
	my $dbhits_blastn		= "$opt{res}/ictv/ictv.dbhits.blastn.tsv";
	my $dbhits_blastp		= "$opt{res}/ictv/ictv.dbhits.blastp.tsv";
	my $annot				= "$opt{res}/ictv/ictv.annot.tsv";
	my $taxid_prob			= "$opt{res}/ictv/ictv.taxid_prob.tsv";
	my $taxid_prob_iter		= "$opt{res}/ictv/ictv.taxid_prob_byiter.tsv";
	my $genus_prob			= "$opt{res}/ictv/ictv.genus_prob.tsv";
	my $genus_prob_iter		= "$opt{res}/ictv/ictv.genus_prob_byiter.tsv";
	# reports
	my $annot_table			= "$opt{res}/ictv/ictv.annot_table.tsv";
	my $annot_excel 			= "$opt{res}/ictv/ictv.annot_table.xlsx";
	my $abund_table 			= "$opt{res}/ictv/ictv.abund_table.tsv";
	my $abund_excel 			= "$opt{res}/ictv/ictv.abund_table.xlsx";		
	#my $krona_data			= "$opt{res}/reports/ictv.krona.data.txt";
	#my $krona_graph			= "$opt{res}/reports/ictv.krona.graph.html";
	my $log					= "$opt{logs}/$opt{sample}/ictv.annot.log";
	
	# TMP
	my $contig_taxid_bits	= "$opt{tmpdir}/ictv.contig_taxid_bits.tsv";
	my $readn_taxid 			= "$opt{tmpdir}/ictv.readn_taxid.tsv";
	my $contigs_ta_un		= "$opt{tmpdir}/ictv.contigs.un.fa";
	my $contigs_ta_un_info	= "$opt{tmpdir}/ictv.contigs.un.info.tsv";
	my $annot_taxonomy		= "$opt{tmpdir}/ictv.annot+taxonomy.tsv";
	my $ictv_tmp				= "$opt{tmpdir}/ictv_vmr_selected.tsv";
	
	# PARAMS
	my $staxid			= "staxid";
	my $staxid_genus		= "";
	my $logLdiff			= 0.01;
	my $maxiter			= 1000;
	my $nospnorep		= 0;
	my $score2prob_temp	= $opt{score2prob_temp} || 10;
	my $numth			= $opt{numth};
	my $estimate_q_prob	= 1;
	my $sp_prob_cutoff	= 0;
	my $ge_prob_cutoff	= 0.001;
	my $par_abund_table 	= "-h -w $opt{wmodel} --conttail $opt{tail_contig}";
	
	# WORK
	system("rm -f $log; touch $log");
	#system("rm -fr $opt{res}/ictv/");
	system("mkdir -p $opt{res}/ictv");
	
	
	# run BLASTN x ICTV isolates
	annotate_blastn(seqs				=> $contigs_ta,
					db				=> $ictv_blastn_db,
					dbhits			=> $dbhits_blastn,
					annot			=> $annot,
					append			=> 0,
					log				=> $log,
					numth			=> $opt{numth},
					min_bits			=> $opt{min_blastn_bits},
					filter_tophits	=> 0,
					max_target_seqs	=> 10);
	# chain to BLASTP
	filter_unseqs(seqs=>$contigs_ta, seqs_un=>$contigs_ta_un, annot=>$annot, seqidh=>'qseqid', log=>$log);
	if( nlines($contigs_ta_un) > 0){
		annotate_blastp(seqs			=>$contigs_ta_un,
					seqinfo			=>$contigs_info,
					db				=>$ictv_blastp_db,
					dbhits			=>$dbhits_blastp,
					annot			=>$annot,
					append			=>1,
					log				=>$log, 
					numth			=>$opt{numth},
					min_bits			=>$opt{min_blastp_bits},
					orf_finder		=>$opt{gen},
					min_orf_length	=>$opt{min_orf_length},
					max_target_seqs	=> 10,
					filter_tophits	=> 0,
					retain_ties		=>0);
	}
	else{
		print STDERR "\t$subid: all contigs annotated with BLASTN\n" if ($VERBAL);
	}
	
	# PARSE ICTV TABLE
	system_call("cat $ictv_tsv | ".
			"csvtk cut -tf 'Isolate ID',Species,Genus,Family,Order,Genome,'Host source' | ".
			"csvtk rename -tf 'Host source' -n 'Host_source' | ".
			"csvtk mutate -tf 'Isolate ID' -n 'Isolate.NID' -p '^VMR([0-9]+)' 1> $ictv_tmp 2>> $log");
	my %ictv_isolatenid2sp= read_tsv2hash($ictv_tmp,"Isolate.NID","Species");
	my %ictv_isolatenid2ge= read_tsv2hash($ictv_tmp,"Isolate.NID","Genus");
		
	
	# RUN EM
		# estimate P(q) from contigs.idxstats + annot-table
	my %contids_vi			= read_tsv2hash($annot,"qseqid","qseqid"); # used for list of viral qseqids
	my $q_num				= scalar(keys %contids_vi)+0.0;
	my %q_prob				= map{ $_ => 1.0 / $q_num } keys %contids_vi;# by default just use P(q) = 1/|Q|
	if($estimate_q_prob){		
		my %contid2readn		= read_tsv2hash_noheaders($contigs_idxstats,0,2);
		my $readn_sum_vi		= 0.0;
		map { $readn_sum_vi += $contid2readn{$_}} keys %contids_vi;
		%q_prob				= map{ $_ => $contid2readn{$_} / $readn_sum_vi } keys %contids_vi;
	}
	my $res		= EM_loop(annot=>$annot, staxid=>$staxid, q_prob=>\%q_prob, logLdiff=>$logLdiff,  maxiter=>$maxiter, nospnorep=>$nospnorep, score2prob_temp=>$score2prob_temp, numth=>$numth);
	my $logL		= $res->{logL};
	my %F		= %{$res->{F}};
	my %Fiter	= %{$res->{Fiter}};
	
	# print F-abundancies (Isolates)
	my @taxids	= sort{ $F{$b} <=> $F{$a} } keys %F;
	my %prob_ge = ();
	foreach my $id(@taxids){
		next if(!defined($ictv_isolatenid2ge{$id}));
		my $ge	= $ictv_isolatenid2ge{$id};
		$prob_ge{$ge} 	= 0 if(!defined($prob_ge{$ge}));
		$prob_ge{$ge}	+= $F{$id};
	}
	
	open(OUT,">$taxid_prob") or die "Can\'t open $taxid_prob: $!\n";
	print OUT "Isolate.NID\tSpecies\tGenus\tprob\tprob.genus\n";
	foreach my $id(@taxids){
		my $sp		= $ictv_isolatenid2sp{$id} || 'NA';
		my $ge		= $ictv_isolatenid2ge{$id} || 'NA';
		my $prob		= $F{$id};
		my $prob_ge	= $prob_ge{$ge} || 0;
		print OUT "$id\t$sp\t$ge\t$prob\t$prob_ge\n";
	}
	close(OUT);
	# F-abundancies for each EM iteration (Isolates)
	my @iters	= sort{ $a <=> $b } keys %Fiter;
	open(OUT,">$taxid_prob_iter") or die "Can\'t open $taxid_prob_iter: $!\n";
	printf OUT "Isolate.NID\tSpecies\t%s\n",join("\t",@iters);
	foreach my $t(sort {$F{$b} <=> $F{$a}} keys %F){
		my @Fiter_t	= map {$Fiter{$_}->{$t}} @iters;
		my $sp		= $ictv_isolatenid2sp{$t} || 'NA';
		printf OUT "$t\t$sp\t%s\n",join("\t",@Fiter_t);
	}
	close(OUT);
	

	# GENERATE CONTIG ANNOTATION TABLES

	if( nlines($annot) < 2){
		print STDERR "\n\tWARNING: $subid: No reporting due to empty annotation file: $annot\n\n";
		return();
	}	
	# start with all annotations, delete division and bphage fields (these are NA from blast)
	system_call("csvtk sort -tj $numth -k qseqlen:rN -k qseqid:N -k staxid:N  $annot | ".
				"csvtk rename -tf staxid -n Isolate.NID 1> $annot_table 2>> $log");
	my $division_col		= colind($annot,'division');
	my $bphage_col		= colind($annot,'bphage');
	system("csvtk cut -tf -$division_col,-$bphage_col $annot_table 1> $annot_table.tmp 2>> $log");
	system("mv $annot_table.tmp $annot_table");
	# add P_taxid
	system("csvtk cut -tf Isolate.NID,prob,prob.genus $taxid_prob 1> $taxid_prob.tmp 2>> $log");
	system_call("csvtk join -tj $numth -f 'Isolate.NID;Isolate.NID' -L --na 'NA' $annot_table $taxid_prob.tmp  | ".
				"csvtk rename -tf 'prob,prob.genus' -n 'Isolate.prob,Genus.prob' 1> $annot_table.tmp 2>> $log");
	system("mv $annot_table.tmp $annot_table");
	system("rm -f $taxid_prob.tmp");
	# filter by P-values
	system_call("csvtk filter -tf \"Genus.prob>=$ge_prob_cutoff\" $annot_table 1> $annot_table.tmp 2>> $log");
	system("mv $annot_table.tmp $annot_table");
	
	# add taxonomy + Genome + Host source from ICTV VMR table
	system_call("csvtk join -tj $numth -f 'Isolate.NID' -L --na 'NA' $annot_table $ictv_tmp 1> $annot_table.tmp 2>> $log");
	system("mv $annot_table.tmp $annot_table");
	# add division-filed(required by R-scripts)
	my $isolatenid_col	= colind($annot_table,'Isolate.NID');
	system_call("csvtk mutate2 -t -n division -e \" 'Viruses' \" $annot_table | ".
				"csvtk cut -tf -$isolatenid_col 1> $annot_table.tmp 2>> $log");
	system("mv $annot_table.tmp $annot_table");
	# add bphage-field
	system_call("csvtk mutate2 -tj $numth -n bphage -e '\$Host_source==\"archaea\" || \$Host_source==\"bacteria\" ? \"yes\":\"no\"' $annot_table 1> $annot_table.tmp 2>> $log");
	system("mv $annot_table.tmp $annot_table");
	

	# convert $annot_table to excel file: does not work: TODO printing to excel
	#system_call("$opt{call_R} $R_scripts/print_annot_table.R  $annot_table $annot_excel domain Viruses 2>> $log");
	system_call("$perl_scripts/write_excel.pl $annot_table 1> $annot_excel 2>> $log");
	
		# SORT CONTIGS TO DIRS USING ICTV CLASSIFICATION
	#system_call("perl $perl_scripts/sort_contigs_bytaxa_v3.pl -c $contigs -a $annot_table --res $opt{res}/contigs --toptaxrank $toptaxrank -v &>> $log");
	
	
	#
	# CREATE ABUNDANCE TABLE
	# select top-scoring annotation for each contig (retain ties)
	system_call("cat $annot | ".
				"csvtk cut -tf qseqid,staxid,bitscore | ".
				"csvtk rename -tf qseqid -n contig | ".
				"csvtk sort -tj $opt{numth} -k contig:N -k bitscore:nr 1> $contig_taxid_bits 2>> $log" );
				
	filter_tophits(dbhits=>$contig_taxid_bits,dbhits_flt=>"$contig_taxid_bits.tmp", qcol=>'contig', bitscol=>'bitscore', retain_ties=>1);
	system_call("mv $contig_taxid_bits.tmp $contig_taxid_bits");
	
	if( nlines($contig_taxid_bits)<2 ){
		print STDERR "\n\tWARNING: $subid: Unable to estimate abundancies: no annotations passign threshold\n";
		return();
	}
		# estimate abundancies:
	system_call("perl $perl_scripts/get_abund_table.pl $par_abund_table $contig_taxid_bits $contigs_idxstats 1> $readn_taxid 2> $log" );
	#system("rm $contig_taxid_bits");
	
	if( nlines($readn_taxid)<2 ){
		print STDERR "\n\tWARNING: $subid: printing empty abundance table\n";
		system_call("touch $abund_table");
		return();
	}	
		# add P_taxid
	system("csvtk cut -tf Isolate.NID,prob,prob.genus $taxid_prob 1> $taxid_prob.tmp 2>> $log");
	system_call("csvtk join -tj $numth -f 'taxid;Isolate.NID' -L --na 'NA' $readn_taxid $taxid_prob.tmp  | ".
				"csvtk rename -tf 'prob,prob.genus' -n 'Isolate.prob,Genus.prob' 1> $abund_table 2>> $log");
	system("rm -f $taxid_prob.tmp");
		# filter by P-values
	system_call("csvtk filter -tf \"Genus.prob>=$ge_prob_cutoff\" $abund_table 1> $abund_table.tmp 2>> $log");
	system("mv $abund_table.tmp $abund_table");
		
		# add taxonomy to abundancies and print to $abund_table:
	system_call("csvtk join -tj $numth -f 'taxid;Isolate.NID' -L --na 'NA' $abund_table $ictv_tmp 1> $abund_table.tmp 2>> $log");
	system("mv $abund_table.tmp $abund_table");
	system("rm -f $readn_taxid $contig_taxid_bits");
	
	system_call("$perl_scripts/write_excel.pl $abund_table 1> $abund_excel 2>> $log");

}


# USAGE
# generate_igv_html(
#	refgens => $refgens_fasta		: Fasta with reference genomes
#	contigs=> $contigs_fasta			: Fasta with contigs
#	bam=> $bam_file					: Contigs alignment to refgens in BAM-format, MUST be sorted
#	resdir=> $resdir					: Directory for results
#	split_byrefgen=>1				: Split BAM by refgens (eg output separate igv-reports for each refgen)
#	log=> $log						: Output log, created on demand
#	$opt=>\%opt						: Global options
#	)
sub generate_igv_html{
	# in:
	my (%args)			= @_;
	my %opt 				= %{$args{opt}};
	my $species			= $args{species};
	my $refgens			= $args{refgens};
	my %refgens_minfo	= %{$args{refgens_metainfo}};
	my $contigs			= $args{contigs};
	my $bam_all			= $args{bam};
	my $bed_all			= $args{bed};
	my $resdir			= $args{resdir};
	
	# params:
	my $split_byrefgen	= defined($opt{split_byrefgen}) ? $opt{split_byrefgen} : 0;
	my $use_data_uri		= defined($args{use_data_uri}) ? $args{use_data_uri}: 0;
	# out:
	my $log 				= $args{log};
	
	if($opt{v}){
		print STDERR "\t# generate_igv_html($species)\n";
	}
	
	if(!(-e $log)){
		system("touch $log");
	}
	if(!(-e $resdir)){
		system("mkdir -p $resdir");
	}
	
	# split_byrefgen:
	my @refseq_ids 	= get_SQ_SN($bam_all);
	
	foreach my $id(@refseq_ids){
		# out:
		my $bam_file			= "$id.bam";
		my $refgen_file		= "$id.fasta";		
		my $bam				= "$resdir/$bam_file";
		my $refgen			= "$resdir/$refgen_file";
		my $contigs_gzip		= "$contigs.gzip";
		my $igv_report		= "$resdir/$id.igv.html";
		my %refgen_minfo		= %{$refgens_minfo{$id}};
		my $bed_file			= "$id.bed";
		my $bed				= "$resdir/$bed_file";
		
		# tmp:
		my $sam				= "$resdir/$id.sam.tmp";
		
		# get bam+bai for this refseq	
		system_call("samtools view -b $bam_all $id 1> $bam 2>> $log");
		system_call("samtools index $bam 2>> $log");
		# get refseq fasta+ fai
		system_call("seqkit grep -p $id -w0 $refgens 1> $refgen 2>> $log");
		system_call("samtools faidx $refgen 2>> $log");
		# get bed file for this refseq
		system_call("head -n1 $bed_all 1> $bed 2>> $log");
		system_call("grep $id $bed_all 1>> $bed 2>> $log || [[ \$? == 1 ]]");	# ignore no match with $? == 1 check
		
		# create gzip/binary files for data_uri-conversion
		system_call("gzip -fc $refgen 1> $refgen.gzip 2>> $log");
		system_call("gzip -fc $refgen.fai 1> $refgen.fai.gzip 2>> $log");
		system_call("gzip -fc $contigs 1> $contigs.gzip 2>> $log");
		system_call("gzip -fc $bed 1> $bed.gzip 2>> $log");
		
		# create @summary_table for contigs
		my @summary_table	= ();
		my @headers	= ("contig.id","contig.len","contig.cov","ali.start","strand");
		push(@summary_table,\@headers);
		# parsing info on contigs/contig-aligments from sam files
			# read sam file: ignore everything after CIGAR-col
		my $sam_headers	= "qname\tflag\trname\tpos\tmapq\tcigar";	
		system("echo \"$sam_headers\" 1> $sam");
		system_call("samtools view $bam 1>> $sam 2>> $log");
		my %stats 	= read_tsv2hashtable($sam, 'qname');
			# calc contig.len + alen + apid
		foreach my $id(sort { $stats{$a}->{'pos'} <=> $stats{$b}->{'pos'}} keys %stats){
			my $length		= format_int( cigar2qlen($stats{$id}->{'cigar'}), ' ');
			my $coverage		= sprintf("%.2f", (cigar2qcov($stats{$id}->{'cigar'}))*100 );
			my $astart		= format_int( $stats{$id}->{'pos'}, ' ');
			my $strand		= ($stats{$id}->{'flag'} & (1 << 4))? '-':'+';	# (1<<4: 0001 0000: qseq is reverse complement)
			my @row			= ($id,$length,$coverage,$astart,$strand);
			push(@summary_table,\@row);
		}
		my $table_attrs = "class='sortable'";
		my @th_attrs		= ();
		push(@th_attrs, "class='mixed' type='string'");	# id
		push(@th_attrs, "class='mixed' type='int'");
		push(@th_attrs, "class='mixed' type='num'");
		push(@th_attrs, "class='asc' type='int'");
		push(@th_attrs, "class='mixed' type='string'");
		
		# create IGV-html
		open(OUT,">$igv_report") or die "ERROR: failed to open: $igv_report\n";
		print OUT "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n";
		print OUT "<html>\n";
		print OUT "<head>\n";
		print OUT "<title>Lazypipe RefGen Report</title>\n";
		print OUT "<style>\n";
			print_file("$perl_scripts/include/reports.css",\*OUT);
		print OUT "</style>\n";
		print OUT "<script>\n";
			print_file("$perl_scripts/include/reports.js",\*OUT);
		print OUT "</script>\n";
		print OUT "<script src=\"https://cdn.jsdelivr.net/npm/igv\@2.13.11/dist/igv.min.js\"></script>\n";
		print OUT "<meta charset=\"utf-8\"/>\n";
		print OUT "</head>\n";
		
		print OUT "<body>\n";
		print OUT "<h3> Lazypipe RefGen Report</h3>\n";
		print OUT "<p>\n";
		print OUT "\tSample: $opt{sample}\n";
		print OUT "\t<br>Species: $species\n";
		print OUT "\t<br>Reference:\n";
		foreach my $k(sort keys %refgen_minfo){
			print OUT "\t<br>&emsp; $k: $refgen_minfo{$k}\n";
		}
		print OUT "</p>\n";
		print OUT "Contigs:\n";
		write_table_html(table=>\@summary_table, table_attrs=>$table_attrs, th_attrs=>\@th_attrs, fh=>\*OUT);
		print OUT "<br>\n";
		
		print OUT "<p>\n";
		print OUT "\t<div id=\"igv_div\" style=\"padding-top:10px;padding-bottom:10px; border:1px solid lightgray\"></div>\n";
		print OUT "</p>\n";
		print OUT "<script type=\"text/javascript\">\n";
		print OUT "    colpalet_blue = ['#A202FF','#2072B2','#1260CC','#OOABFF','#00E8FF'];\n";
		print OUT "    const options = {\n";
        print OUT "    reference:{\n";
        print OUT "        id:		'$id',\n";
        print OUT "        name:	'reference',\n";
        print OUT "        fastaURL:	\"",($use_data_uri)? filebin2uri("$refgen.gzip") : $refgen_file,"\",\n";
        print OUT "        indexURL:	\"",($use_data_uri)? filebin2uri("$refgen.fai.gzip") : "$refgen_file.fai","\",\n";             
        print OUT "        wholeGenomeView: true\n";
        print OUT "    },\n";
        print OUT "		tracks:[\n";
        if( 	nlines($bed) > 1){
        		print OUT "			{	name:	'Genes',\n";
        		print OUT "				type:	'annotation',\n";
			print OUT "				format:	'bed',\n";
			print OUT "				sourceType:	'file',\n";
			print OUT "				url:		\"",($use_data_uri)? filebin2uri("$bed.gzip"): $bed_file,"\",\n";
			print OUT "				indexed:	 false,\n";
			print OUT "				displayMode: 'expanded'\n";
			print OUT "			},\n";
        }
		print OUT "			{   name:		'Contigs',\n"; 
 	    print OUT "            type:		'alignment',\n"; 
        print OUT "            format:     'bam',\n";
 	    print OUT "            showCoverage: false,\n";
 	    print OUT "            showAlignments: true,\n";
 	    print OUT "            pairsSupported: false,\n";
 	    print OUT "            colorBy: 	'strand',\n";
 	    print OUT "            alignmentRowHeight: 14,\n";
        print OUT "            autoHeight: 'true',\n";
        print OUT "            visibilityWindow:50000,\n";
        print OUT "            sourceType: 'file',\n";
	    print OUT "            url:			\"",($use_data_uri)? filebin2uri($bam): $bam_file,"\",\n";
	    print OUT "            indexURL:		\"",($use_data_uri)? filebin2uri("$bam.bai"): "$bam_file.bai","\",\n";
        print OUT "            indexed:    true,\n";
        print OUT "            displayMode: 'expanded'\n";
        print OUT "        }\n";
        print OUT "    ]\n";
        print OUT "    };\n";

        print OUT "	var igvDiv = document.getElementById('igv_div');\n";
        print OUT "	igv.createBrowser(igvDiv, options);\n";
        print OUT "	console.log('Created IGV browser');\n";
        print OUT " 	init_sortable_tables();\n";
        print OUT "</script>\n";
        
        print OUT "</body>\n";
        print OUT "</html>\n";
        close(OUT);
	}
}

# Get @SQ SN-field values from bam file
# USAGE:
# my @SN_list = get_SQ_SN($bam_file)
#
sub get_SQ_SN{
	# in:
	my $bam = shift(@_);
	# tmp:
	my $file = "$bam.h.tmp";
	
	system_call("samtools view -H $bam | grep \"^\@SQ\" 1> $file");
	my @sn_list = ();
	open(IN,"<$file") or die "Can\'t open $file: $!\n";
	while(my $l=<IN>){
		chomp($l);
		my @sp = split(/\t/,$l,-1);
		foreach my $str(@sp){
			$str =~ s/^\s+|\s+$//g;
			if($str =~ m/^SN:([\w\.]+)/g){
				push(@sn_list,$1);
				next;
			}
		}
	}
	close(IN);
	system("rm -f $file");
	return(@sn_list);
}


sub generate_stats{
	
	# in:
	my %opt 			= %{shift()};
	my $r1			= defined( $opt{read1} )? $opt{read1} : 0;	# can be undefined
	my $r2			= defined( $opt{read2} )? $opt{read2} : 0;	# can be undefined
	my $r1_flt		= "$opt{res}/reads/read1.trim.fq.gz";
	my $r1_hgflt		= "$opt{res}/reads/read1.trim.hflt.fq.gz";
	   $r1_hgflt		= (-e $r1_hgflt)? $r1_hgflt : $r1_flt;
	my $contigs 		= "$opt{res}/contigs.fa";
	my $idx 			= "$opt{res}/contigs.idxstats";
	my $orfs_nt		= "$opt{res}/contigs.orfs.nt.fa";
	
	# tmp:
	my $r1_len		= "$opt{res}/reads/read1.len";		# output if -r1 specified
	my $r2_len		= "$opt{res}/reads/read2.len";		# output if -r2 specified
	my $cont_len		= "$opt{res}/contigs.len";
	
	# out:
	my $stats		= "$opt{'res'}/assembly.stats.tsv";
	my $stats_yaml	= "$opt{'res'}/assembly.stats.yaml";
	my $log 			= "$opt{'logs'}/$opt{'sample'}/stats.log";
	# res/figures: contig.hist*.png, read.hist.png
	#my $surv_fig		= "$opt{'res'}/qc.readsurv.jpeg";
	
	
	print STDERR "\n# ASSEMBLY STATS\n\n";

	if((-e $contigs) && (nlines($contigs)>1) && (-e $orfs_nt)){
		system_call("perl $perl_scripts/assembly_stats2.pl --format col,names ".
				  	"--reads $r1 --reads_flt $r1_flt --reads_hgflt $r1_hgflt ".
				  	"--cont $contigs --idx $idx --orfs $orfs_nt 1> $stats 2> $log", $opt{'v'});
		system_call("perl $perl_scripts/assembly_stats.pl $contigs mean,sum,N50,LN500,Lbp500,LN1000,Lbp1000 col,names 1>> $stats 2>> $log", $opt{'v'});
		
		assembly_stats(res_dir=>$opt{res}, numth=>$opt{numth}, yaml=>$stats_yaml);
	}

	print STDERR "\n# QC PLOTS\n\n";
	# extract length
	if($r1 && (-e $r1)){
		system_call("seqkit fx2tab -nil $r1 | cut -f2  1> $r1_len");
	}
	if($r2 && (-e $r2) ){
		system_call("seqkit fx2tab -nil $r2 | cut -f2  1> $r2_len");
	}
	if( (-e $contigs) && nlines($contigs)>1 ){
		system_call("seqkit fx2tab -nil $contigs | cut -f2  1> $cont_len");
	}
	
	# plot histograms
	system("mkdir -p $opt{res}/figures");
	if( (-e $r1_len && -e $r2_len)){
		system_call("$opt{'call_R'} $R_scripts/hist_figures.R read_hist $opt{res}/figures $r1_len $r2_len ");
	}
	elsif( -e $r1_len){
		system_call("$opt{'call_R'} $R_scripts/hist_figures.R read_hist $opt{res}/figures $r1_len ");
	}
	
	if( -e $cont_len){
		system_call("$opt{'call_R'} $R_scripts/hist_figures.R cont_hist $opt{res}/figures $cont_len ");
	}

	# REMOVE TMP
	system("rm -f $r1_len $r2_len $cont_len");
}





sub pack_files{
	print STDERR "\n# PACK FILES FOR SHARING\n\n";
	
	# in:
	my %opt 			= %{shift()};
	my $res			= $opt{res};
	my $dirname		= dirname($opt{res});
	my $basename 	= basename($opt{res});
	
	# out: 
	# $dirname/$basename.tar.gz

	system_call("rm -fR $dirname/$basename.tar" );
	system_call("mkdir -p $dirname/$basename.tar" );
	
	my @files_share = ();
	#push(@files_share, <$opt{'res'}/*.html>);
	#push(@files_share, <$opt{'res'}/*.fa>);
	#push(@files_share, <$opt{'res'}/*.jpeg>);
	#push(@files_share, <$opt{'res'}/*.xlsx>);
	#push(@files_share, <$opt{'res'}/*.tsv>);

		# TSV-files
	push(@files_share, "$res/abund_table.tsv");
	push(@files_share, "$res/annot_table.tsv");
	push(@files_share, "$res/assembly.stats.tsv");
	push(@files_share, "$res/assembly.stats.yaml");
	push(@files_share, "$res/contigs.fa.info.tsv");
	push(@files_share, "$res/contigs.idxstats");
	push(@files_share, "$res/readid_contigid.tsv");
		# XLSX-files
	push(@files_share, "$res/abund_table.xlsx");
	push(@files_share, "$res/annot_table.xlsx");
		# FASTA-files
	push(@files_share, "$res/contigs");
	push(@files_share, "$res/contigs.fa");
	push(@files_share, "$res/contigs.ann1.ab.fa");
	push(@files_share, "$res/contigs.ann1.ph.fa");
	push(@files_share, "$res/contigs.ann1.un.fa");
	push(@files_share, "$res/contigs.ann1.vi.fa");
	push(@files_share, "$res/contigs.ann2.ab.fa");
	push(@files_share, "$res/contigs.ann2.ph.fa");
	push(@files_share, "$res/contigs.ann2.un.fa");
	push(@files_share, "$res/contigs.ann2.vi.fa");
	push(@files_share, "$res/contigs.orfs.aa.fa");
	push(@files_share, "$res/contigs.orfs.nt.fa");
	push(@files_share, "$res/ictv");
	if($opt{pack_reads}){
		push(@files_share, "$res/reads/");
	}
		# IMAGES
	push(@files_share, "$res/figures");
		# REPORTS
	push(@files_share, "$res/reports");
		# VARIOUS
	push(@files_share, "$res/History.log");
	
	# filter existing-files
	my @files_share_flt = ();
	foreach my $file(@files_share){
		if(-e $file){
			push(@files_share_flt,$file);
		}
	}
	my $files_share_str = join(" ",@files_share_flt);
	
	system_call("cp -r $files_share_str $dirname/$basename.tar/" );
	system_call("tar -czf $dirname/$basename.tar.gz -C $dirname/$basename.tar ." );
	system_call("rm -fR $dirname/$basename.tar" );	
}
sub clean{
	my %opt = %{shift()};
	
	# CLEANUP
	system_call("rm -fR $opt{res}/*.bam  $opt{'res'}/*.bai", $opt{'v'});
	system_call("rm -fR $opt{res}/*.sam", $opt{'v'});
	system_call("rm -fR $opt{res}/*Graph*", $opt{'v'});
	system_call("rm -fR $opt{res}/Roadmaps", $opt{'v'});
	system_call("rm -fR $opt{res}/*tmp*", $opt{'v'});
	system_call("rm -fR $opt{res}/assembler_out", $opt{'v'});
	system_call("rm -fR $opt{res}/dbhits.*", $opt{v});
	system_call("rm -fR $opt{res}/contigs.fa.amb $opt{res}/contigs.fa.ann $opt{res}/contigs.fa.bwt $opt{res}/contigs.fa.pac $opt{res}/contigs.fa.sa");
	system_call("rm -fR $opt{res}/hostgen.sam.flt");
}		





#
# Check and format pipeline options
#
sub options_format{
	my $subid	= "options_format()";
	my %opt = %{ shift(@_) };
	
	if($opt{help}){
		print $usage; exit(0);
	}
	
	# List installed databases
	if( $opt{databases} && defined($opt{'ann.databases'}) ){
		print "\n# Listing Installed Reference Databases:\n";
		
		foreach my $k(sort keys %{$opt{'ann.databases'}}){
			my $db			= $opt{'ann.databases'}->{$k};
			my $dbpath		= $db->{db};
			if($dbpath =~ /\$(\w+)/g){
				my $envar	= $1;
				if(defined($ENV{$envar})){
					$dbpath =~ s/\$$envar/$ENV{$envar}/g;
				}
				else{
					print STDERR "WARNING: $subid: undefined environment variable \"$envar\"\n";
				}
			}
			my @dbfiles		= glob("$dbpath*");
			next if( scalar(@dbfiles)<1 );
			print "$k:\n";
			print "\tdb:     $db->{db}\n";
			print "\tname:   $db->{name}\n";
		}
		exit(0);
	}
	
	#List installed background filters
	if( $opt{filters} && defined($opt{'host.databases'})){
		print "\n# Listing Installed Background Filters:\n";
		
		foreach my $k(sort keys %{$opt{'host.databases'}}){
			my $db		= $opt{'host.databases'}->{$k};
			my $dbpath	= $db->{db};
			if($dbpath =~ /\$(\w+)/g){
				my $envar	= $1;
				if(defined($ENV{$envar})){
					$dbpath =~ s/\$$envar/$ENV{$envar}/g;
				}
				else{
					print STDERR "WARNING: $subid: undefined environment variable \"$envar\"\n";
				}
			}
			next if( !(-e "$dbpath.amb") || !(-e "$dbpath.ann") || !(-e "$dbpath.bwt"));
			print "$k:\n";
			print "\tdb:         $db->{db}\n";
			print "\tname (acc): $db->{name} ($db->{accession})\n";
		}
		exit(0);
	}

	# Shift section general.parameters to top-level, honor preference of command-line options
	my %tmp = %{$opt{"general.parameters"}};
	for my $par(keys %tmp){
		if(!defined($opt{$par})){
			$opt{$par}	= $tmp{$par};
		}
	}

	# Expand ENV variables, if any
	foreach my $k1(keys %opt){
			# hash reference
		if( defined(ref($opt{$k1})) 
			&& ref($opt{$k1}) eq ref {}){
			foreach my $k2( sort keys %{$opt{$k1}} ){
					# hash reference
				if( defined(ref($opt{$k1}->{$k2})) && ref($opt{$k1}->{$k2}) eq ref {} ){
					
					foreach my $k3( keys %{$opt{$k1}->{$k2}}){
						
						if($opt{$k1}->{$k2}->{$k3} =~ /\$(\w+)/g){
							my $envar	= $1;
							if(defined($ENV{$envar})){
								$opt{$k1}->{$k2}->{$k3} =~ s/\$$envar/$ENV{$envar}/g;
							}
							else{
								print STDERR "WARNING: $subid: undefined environment variable \"$envar\"\n";
							}
						}
					}
				}
					# hash entry
				else{
					if($opt{$k1}->{$k2} =~ /\$(\w+)/g){
						my $envar	= $1;
						if(defined($ENV{$envar})){
							$opt{$k1}->{$k2} =~ s/\$$envar/$ENV{$envar}/g;
						}
						else{
							print STDERR "WARNING: $subid: undefined environment variable \"$envar\"\n";
						}
					}			
				}	
			}
		}
			# hash entry
		else{
			if($opt{$k1} =~ /\$(\w+)/g){
				my $envar	= $1;
				if(defined($ENV{$envar})){
					$opt{$k1} =~ s/\$$envar/$ENV{$envar}/g;
				}
				else{
					print STDERR "WARNING: $subid: undefined environment variable \"$envar\"\n";
				}
			}		
		}	
	}
	
	# set to lower case
	$opt{pipe}		= lc($opt{pipe});
	$opt{pre} 		= lc($opt{pre});
	$opt{ann1} 		= lc($opt{ann1});
	$opt{ann2}		= lc($opt{ann2});
	$opt{ass} 		= lc($opt{ass});
	$opt{gen}       = lc($opt{gen});
	$opt{wmodel} 	= lc($opt{wmodel});

	# check --pipe
	# PIPELINE STEPS
	if( !$opt{pipe} ){
		die "missing arguments: --pipe <str>\n";
	}
	my %pipeh;
	my @tmp= split(/,/,$opt{pipe},-1);
	foreach my $t(@tmp){
		my $t	= lc($t);
			if(	$t =~ m/^pre/i )		{ $pipeh{'prepro'}		= 1; }
		elsif(	$t =~ m/^flt/i )		{ $pipeh{'filter'}		= 1; }
		elsif(	$t =~ m/^ass/i )		{ $pipeh{'assemble'}		= 1; }
		elsif(	$t =~ m/^rea/i )		{ $pipeh{'realign'}		= 1; }
		elsif(	$t =~ m/^(ann1|annot1)/i ){ $pipeh{'ann1'} 	= 1; }
		elsif(	$t =~ m/^(ann2|annot2)/i ){ $pipeh{'ann2'} 	= 1; }
		elsif(	$t =~ m/^ictv/i )	{ $pipeh{ictv}			= 1;	 }
		elsif(	$t =~ m/^sta/i )		{ $pipeh{'stats'}		= 1; }
		elsif(	$t =~ m/^rep/i )		{ $pipeh{'report'}		= 1; }
		elsif(	$t =~ m/^rgrep/i )	{ $pipeh{'rgreport'}		= 1; }
		elsif(	$t =~ m/^pack/i )		{ $pipeh{'pack'}		= 1; }
		elsif(	$t =~ m/^clean/i )		{ $pipeh{'clean'}	= 1; }
		elsif(  $t =~ m/^main/i ){
			$pipeh{'prepro'} 	= 1;
			$pipeh{'filter'}		= 1;
			$pipeh{'assemble'}	= 1;
			$pipeh{'realign'}	= 1;
			$pipeh{'ann1'}		= 1;
			$pipeh{'ann2'}		= 1;
			$pipeh{'stats'}		= 1;
			$pipeh{'report'}		= 1;
			$pipeh{'pack'}		= 1;
		}
		elsif(  $t =~ m/^all/i ){
			$pipeh{'prepro'} 	= 1;
			$pipeh{'filter'}		= 1;
			$pipeh{'assemble'}	= 1;
			$pipeh{'realign'}	= 1;
			$pipeh{'ann1'}		= 1;
			$pipeh{'ann2'}		= 1;
			$pipeh{'stats'}		= 1;
			$pipeh{'report'}		= 1;
			$pipeh{'rgreport'}	= 1;
			$pipeh{'pack'}		= 1;
			$pipeh{'clean'}		= 1;
		}
		else{
			die "ERROR: invalid argument --pipe $t\n";
		}
	}
		# final check
	$opt{'pipe'} = \%pipeh;
	
	
	# annotation strategies
	if($opt{anns}){
		my $strat	= lc($opt{anns});
		if(defined($opt{'ann.strategies'}) && defined($opt{'ann.strategies'}->{$strat})){
			my $strat_str 	= $opt{'ann.strategies'}->{$strat};
			my ($ann1)	 	= $strat_str =~ m/--ann1\s+([\w\.\,]+)/i;
			my ($ann2)		= $strat_str =~ m/--ann2\s+([\w\.\,\:]+)/i;

			if(!defined($ann1) || $ann1 eq ''){
				die "ERROR: invalid annotation.strategy: $strat_str\n";
			}
			$opt{pipe}->{ann1}	= 1;
			$opt{ann1}			= lc($ann1);
			
			if(!defined($ann2) || $ann2 eq ''){
				$opt{pipe}->{ann2} 	= 0;
				$opt{ann2}			= 0;
			}
			else{
				$opt{pipe}->{ann2} 	= 1;
				$opt{ann2}			= lc($ann2);	
			}		
		}
		else{
			print STDERR "ERROR: undefined annotation strategy: $strat\n";
		}
	}	
	
	
	
	# --read1 --read2: required for --pipe prepro, optional for --pipe stats
		# guess read2
	if( !$opt{se} ){ # PE-reads
		if( !$opt{read2} ){
			$opt{read2}			= $opt{read1};
			my $matched 		= ($opt{read2} =~ s/_R1/_R2/);
			if( !$matched ){
				$matched		= ($opt{read2} =~ s/_r1/_r2/);
			}
			if( !$matched ){
				$matched		= ($opt{read2} =~ s/_f1/_r2/);
			}
			if( !$matched ){
				die "ERROR: could not quess filename for reverse reads. Please specify explicitely --read1 and --read2";
			}
			if( $opt{read1} eq $opt{read2}){
				die "ERROR: could not guess filename for reverse reads. Please specify explicitely --read1 and --read2";
			}
		}
	}
		# check that reads are suppied when running --pipe prepro
	if($pipeh{prepro}){
		if( !$opt{read1} ){
			die "ERROR: missing arguments: --read1\n";
		}
		if( !(-e $opt{read1}) ){
			die "ERROR: check read1 file: $opt{read1} does not exist\n";
		}
		if( !$opt{se} && !$opt{read2} ){
			die "ERROR: missing arguments: --read2\n";
		}
		if( !$opt{se} && !(-e $opt{read2}) ){
			die "ERROR: check read2 file: $opt{read2} does not exist\n";
		}
	}
	
	
	# --pre
	if( $opt{pre} =~ /trimm/gi ){
		$opt{pre} = 'trimm';
	}
	elsif( $opt{pre} =~ /fastp/gi ){
		$opt{pre} 	= 'fastp';
	}
	elsif( $opt{pre} =~ /none/gi ){
		$opt{pre} 	= 'none';
	}
	else{
		print STDERR "\ninvalid option --pre $opt{pre}. Running with --pre fastp\n";
		$opt{pre} = "fastp";
	}
	
	# --hostgen
		# --hostgen <undef>
	if( !defined($opt{hostgen})  ||  !$opt{hostgen} ){
		
	}	# --hostgen <hostgenome fasta|bwa-index>
	elsif( $opt{hostgen} =~ m/\.(fasta|faa|fna|fa)$|\.(fasta|faa|fna|fa)\.gz/gi ){
		my @hostdb_list		= split(',',$opt{hostgen},-1);
		my @hostdb_list2		= ();
		
		foreach my $db(@hostdb_list){
			#if( !(-e "$db.amb") || !(-e "$db.ann") || !(-e "$db.bwt")){
			#	die "ERROR: bwa index files .amb/.ann/.bwt not found for host-database $db\n";
			#}	indexed on demand
			my $name				= basename( $db );
			$name				=~ s/\..*$//g;
			my %hostdb			= ();
			$hostdb{accession}	= undef;
			$hostdb{db}			= $db;
			$hostdb{name}		= $name;
			$hostdb{latinName}	= undef;
			$hostdb{commonName}	= undef;
			$hostdb{taxid}		= undef;
			
			push(@hostdb_list2,\%hostdb);
		}
		$opt{hostdb}		= \@hostdb_list2;
		
	}	# --hostgen <list of host.database keys>
	else{
		if(!defined($opt{'host.databases'})){
			die "ERROR: host.databases undefined (check your config.yaml)\n";
		}
		my @hostdb_keys	= split(',',$opt{hostgen},-1);
		my @hostdb_list	= ();
		foreach my $k(@hostdb_keys){
			if( !defined($opt{'host.databases'}->{$k})){
				die "ERROR: database with key=$k not found in host.databases\n";
			}
			my $hostdb	= $opt{'host.databases'}->{$k};
			if( !(-e "$hostdb->{db}.amb") || !(-e "$hostdb->{db}.ann") || !(-e "$hostdb->{db}.bwt")){
				die "ERROR: bwa index files .amb/.ann./.bwt not found for host-database $k\n";
			}
			push(@hostdb_list, $hostdb);
		}
		$opt{hostdb}		= \@hostdb_list;
	}
	
	# --ass
	if(!(($opt{ass} eq 'megahit') || ($opt{ass} eq 'spades')) ){
			print STDERR "\ninvalid option --ass $opt{ass}. Running with --ass megahit\n";
			$opt{ass} = 'megahit';	
	}
	# --gen
	if( !(($opt{gen} eq 'mga') || ($opt{gen} eq 'prod') || ($opt{gen} eq 'orfipy') )){
		print STDERR "\ninvalid option --gen $opt{gen}. Running with --gen mga\n\n";
		$opt{gen} = 'mga';
	}
	
	# --ann1 <str>: parsing to list-hash structure: $ann1[0]->{search|target|db}
	
	my %valid_searches 	= ('minimap'=>1,'sans'=>1,'blastn'=>1,'blastp'=>1, 'hmmscan'=>1, 'diamondx'=>1);
	my %valid_targets	= ('ab',1,'ph',1,'vi',1,'un',1);	

	
	if($pipeh{ann1} && defined($opt{ann1}) && $opt{ann1}){
		# --ann1 is fixed to start with contigs.all-annotation, and complement that search by searching with contigs.un in consecutive annotations
		my $target			= "all";
		my @ann1 			= ();
		my @ann1_str			= split(/,/,$opt{ann1});		
		
	  foreach my $dbname(@ann1_str){
		if(!defined($opt{'ann.databases'}) || !defined($opt{'ann.databases'}->{$dbname})){
			die "\nERROR: no database in config.yaml: $dbname\n\n";
		}
		my %ann	 			= %{$opt{'ann.databases'}->{$dbname}};
		$ann{target} 		= $target;
		
		if($dbname ne 'sans'){
			# check that database files exist
			my $dbpath		= $ann{db};
				$dbpath		=~ s/\$(\w+)/$ENV{$1}/g;
			my @dbfiles		= glob("$dbpath*");
			if( scalar(@dbfiles)<1){
				die "\nERROR: no database on disk: $ann{db}\n\n";
			}
		}
		push(@ann1, \%ann);
		$target = 'un';
	  }
	  if(scalar(@ann1)>0){
	  	$opt{ann1}	= \@ann1;
	  }
	  else{
	  	$opt{ann1}	= undef;
	  }
	}
	
	# --ann2 <str>: parsing to list-hash structure: $ann2[0]->{search|target|db}
	if($pipeh{ann2} && defined($opt{ann2}) && $opt{ann2}){
		my @ann2 		= ();
		my @ann2_str		= split(/,/,$opt{ann2});
	
	  foreach my $target_dbname(@ann2_str){
	  	
	  	my ($target, $dbname)	= split(/:/,$target_dbname,2);
		$target					= lc($target);
		
		if( !defined($valid_targets{$target}) ){
			die "\nERROR: invalid target option: --ann2 $target:$dbname\n\n";
		}
		
		if(!defined($opt{'ann.databases'}) || !defined($opt{'ann.databases'}->{$dbname})){
			die "\nERROR: no database in config.yaml:  $dbname\n\n";
		}
		
		my %ann 				= %{$opt{'ann.databases'}->{$dbname}};
		$ann{target}			= $target;

		if($dbname ne 'sans'){
			# check that database files exist
			my $dbpath		= $ann{db};
			   $dbpath		=~ s/\$(\w+)/$ENV{$1}/g;
			my @dbfiles		= glob("$dbpath*");
			if( scalar(@dbfiles)<1){
				die "\nERROR: no database on disk: $ann{db}\n\n";
			}
		}
		push(@ann2, \%ann);
	  }
	  if(scalar(@ann2)>0){
	  	$opt{ann2}	= \@ann2;
	  }
	  else{
	  	$opt{ann2}	= undef;
	  }
	}
	
	# --pipe report
	if( !(($opt{wmodel} eq 'taxacount') || ($opt{wmodel} eq 'bitscore') || ($opt{wmodel} eq 'bitscore2')) ){
		print STDERR "\nWARNING: invalid option --wmodel $opt{wmodel}. Running with --wmodel bitscore\n\n";
		$opt{wmodel} = 'bitscore';
	}
	
	# --RES --SAMPLE: RESDIR
	if( !$opt{res} ){
		$opt{res}	= "results";
	}
	if( !$opt{sample} ){
		if( !$opt{read1} ){
			die "ERROR: invalid options, must specify --read1 or --sample\n";
		}
		$opt{sample}	= basename( $opt{read1} );
		$opt{sample}	=~ s/\..*$//g;
		if($opt{trimm_sample_name}){
			$opt{sample}=~ s/_.*//g;
		}
	}	
	$opt{res}		= "$opt{res}/$opt{sample}";
	
	# TAXONOMY
	if( !defined($opt{taxonomy}) || !defined($opt{taxonomy}->{db})  ){
		die "ERROR: missing or invalid 'taxonomy' settings\n";
	}
	
	system_call("mkdir -p $opt{res}");
	system_call("mkdir -p $opt{tmpdir}");
	system_call("mkdir -p $opt{logs}/$opt{sample}");
	system_call("mkdir -p $opt{res}/reads");
	system_call("mkdir -p $opt{res}/figures");
	system_call("mkdir -p $opt{res}/reports/figures");

	# THREADS
	if( !$opt{'numth'} ){
		$opt{'numth'}	= 1;
	}		
		
	# CREATE tmpdir
	if( defined($opt{'tmpdir_keep'}) && $opt{'tmpdir_keep'} ){
		$opt{'tmpdir'} 	 	= tempdir("lazypipe_XXXXXXXX", DIR => $opt{'tmpdir'}, CLEANUP => 0);
	}
	else{
		$opt{'tmpdir'} 	 	= tempdir("lazypipe_XXXXXXXX", DIR => $opt{'tmpdir'}, CLEANUP => 1);
	}	
	
	# R ENVIRONMENT
	$opt{call_R} 		= "Rscript" if( !$opt{call_R} );
	$ENV{TMPDIR} 		= $opt{tmpdir};
	my $renviron_local	= $opt{tmpdir}."/.Renviron";
	system("echo \"TMPDIR=$opt{tmpdir}\" > $renviron_local");
	system("echo \"TMP=$opt{tmpdir}\" >> $renviron_local");
	$ENV{R_ENVIRON_USER}	= $renviron_local;	
	
	# CHECK BINARIES AVAILABLE
	my $pigz = `sh -c 'command -v pigz'`; 
	if($pigz){
		$opt{'gzip'} = "pigz -p $opt{'numth'}";
	}
	else{
		$opt{'gzip'} = "gzip";
	}
	# DEBUG: PRINT OPT
	#foreach my $k(sort keys %opt){ print "$k\t:$opt{$k}\n";};
	#foreach my $k(sort keys %{$opt{'pipe'}}){ 	print "$k\t: ",$opt{'pipe'}->{$k},"\n"; }
	#print "--ann1:\n"; foreach my $ann( @{$opt{ann1}}){ print "$ann->{target}:$ann->{search}:$ann->{db}\n"; }
	#print "--ann2:\n"; foreach my $ann( @{$opt{ann2}}){ print "$ann->{target}:$ann->{search}:$ann->{db}\n"; }
	#print "--hostgen:\n"; foreach my $db( @{$opt{hostdb}}){ print "$db->{name}: $db->{latinName}: $db->{db}\n"; }
	#exit(1);
	return %opt;
}





sub update_taxonomy{
	my $opt 				= shift(@_);
	my $taxonomy			= $opt->{taxonomy};
	my $taxonomy_nodes	= "$taxonomy->{db}/nodes.dmp";
	
	# UPDATING TAXONOMY FILES: this will also load taxonomy on the very first usage
	if( (-e "$taxonomy_nodes") 
			&& ((-M "$taxonomy_nodes") > $taxonomy->{update_time}) ){
				
		system_call("wget -q $taxonomy->{url} -O $taxonomy->{db}/taxdump.tar.gz", $opt{'v'});
		system_call("tar -xzf $taxonomy->{db}/taxdump.tar.gz -C $taxonomy->{db}", $opt{'v'});
	}
	else{
		print STDERR "\ttaxonomy db up to date\n";
	}
}

# USAGE:
# my %fasta_hash = ‰{readfasta($fasta_file)};
#
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
		if ($_ =~ s/^>//){	
			$header= $_;
			if($sequence{$header}){
				print colored("#CAUTION: SAME FASTA HAS BEEN READ MULTIPLE TIMES.\n#CAUTION: PLEASE CHECK FASTA SEQUENCE:$header\n","red");
			}
			if($temp_seq){
				$temp_seq="";
			} # If there is alreay sequence in temp_seq, empty the sequence file
		}
		else{
		   s/\s+//g;
		   $temp_seq .= $_;
		   $sequence{$header}=$temp_seq; #update the contents
		}
	}
	return \%sequence;
}

# Prints file to file-handle
sub print_file{
	my $file		= shift;
	my $fh		= shift;
	open(IN, "<$file") or die "couldn't open the file $file $!";
	while (<IN>){
		print $fh $_;
	}
	close(IN);
}

#
# Gathers assembly stats to a 2D hash and optinally write to one-document YAML
#
# 1st level keys, also available in ordered list $stats{keys}
#	'reads.trimmed'
#	'reads.hostflt'
#	'reads.assembled'
#	'contigs'
#	'contigs500bp'
#	'contigs1000bp'
# 
# USAGE:
# my %stats 		= assembly_stats(res_dir=>$mypath);
# say "Reads after trimming:";
# foreach my $k(@{$stats{'reads.trimmed'}->{keys}}){
#	say "$k:   ". $stats{'reads.trimmed'}->{$k};
# }
# OR
# assembly_stats(res_dir=>$mypath, yaml=>"$mypath/stats.yaml");
# my $yaml = YAML::Tiny->read( "$mypath/stats.yaml" );
#
# foreach my $k1(%{$yaml->[0]}){
#	say "$k1:";
#	foreach my $pair(  @{$yaml->[0]->{$k1}} ){
#		say '   ',(keys %$pair)[0],':',(values %$pair)[0];
#	}
# }
#
sub assembly_stats{
	my $subid		= "assembly_stats()";
	my (%args)		= @_;
	if(!defined($args{res_dir})){
		die "ERROR: $subid: missing argument: 'res_dir'";
	}
		
	# in:
	my $res_dir		= $args{res_dir};
	# $args{yaml}	is optional 
	my $r1_trim		= (-e "$res_dir/reads/read1.trim.fq.gz") ? "$res_dir/reads/read1.trim.fq.gz": "$res_dir/reads/read1.trim.fq";
	my $r1_hgflt		= 
		(-e "$res_dir/reads/read1.trim.hflt.fq.gz") ? "$res_dir/reads/read1.trim.hflt.fq.gz": 
			(-e "$res_dir/reads/read1.trim.hflt.fq") ? "$res_dir/reads/read1.trim.hflt.fq" : $r1_trim;
	my $contigs_fa	= "$res_dir/contigs.fa";
	my $orfs_fa		= "$res_dir/contigs.orfs.nt.fa";
	my $idxstats		= "$res_dir/contigs.idxstats";
	my $readid_contigid = "$res_dir/readid_contigid.tsv";
		# input checks
	my @musthave_files = ($r1_trim,$contigs_fa,$orfs_fa,$idxstats,$readid_contigid);
	foreach my $file(@musthave_files){
		if(!(-e $file)){
			die "ERROR: $subid: missing file: $file";
		}
	}
	# tmp:
	my $r1_contigs	= "$res_dir/read1.cont.fq.gz";
	my $seqkit_tmp	= "$res_dir/seqkit.tmp";
	my $ids_tmp		= "$res_dir/ids.tmp";
	my $ids2_tmp		= "$res_dir/ids2.tmp";
	my $ids3_tmp		= "$res_dir/ids3.tmp";
	
	# out:
	my $stats_yaml	= "$res_dir/assembly.stats.yaml";
	
	# params:
	my $threads 		= $args{numth} || 8;

	
	
	# START WORKING
	#	STATS FOR READS
	my %stats		= ();
	system_call("seqkit grep -j $threads -f <(cut -f1 $readid_contigid) $r1_hgflt 1> $r1_contigs");
	system_call("seqkit stats -j $threads -baT $r1_trim $r1_hgflt $r1_contigs  1> $seqkit_tmp");
		#seqkit stats -a: num_seqs	sum_len	min_len	avg_len	max_len	Q1	Q2	Q3	sum_gap	N50	Q20(%)	Q30(%)	GC(%)
	my %tmp			= read_tsv2hashtable($seqkit_tmp,'file');
	my @read_stats_keys	= ("number_of_reads",'minimum_length','mean_length','median_length','maximum_length','total_length','GC(%)','Q20(%)','Q30(%)');
	my %rename	= 
		(num_seqs	=> "number_of_reads",
		 min_len		=> 'minimum_length',
		 max_len		=> 'maximum_length',
		 avg_len		=> 'mean_length',
		 Q1			=> 'Q1',
		 Q2			=> 'Q2',
		 Q3			=> 'Q3',
		 sum_len		=> 'total_length',
		 N50			=> 'N50',
		 'Q20(%)'	=> 'Q20(%)',
		 'Q30(%)'	=> 'Q30(%)',
		 'GC(%)'		=> 'GC(%)'
		 );
	my %stats_1	= ();
	my $rowp				= $tmp{basename($r1_trim)};
	foreach my $k(sort keys %rename){
		if( defined($rowp->{$k})){
			$stats_1{$rename{$k}}		= $rowp->{$k};
		}}
	my %stats_2	= ();
	$rowp				= $tmp{basename($r1_hgflt)};
	foreach my $k(sort keys %rename){
		if( defined($rowp->{$k})){
			$stats_2{$rename{$k}}		= $rowp->{$k};
		}}
	my %stats_3	= ();
	$rowp				= $tmp{basename($r1_contigs)};
	foreach my $k(sort keys %rename){
		if( defined($rowp->{$k})){
			$stats_3{$rename{$k}}		= $rowp->{$k};
		}}
	
	# STATS FOR CONTIGS
	my @contig_stats_keys	= ("number_of_contigs",'number_of_ORFs','minimum_length','mean_length','median_length','maximum_length','total_length','N50','GC(%)');
	my %renamec	= 
		(num_seqs	=> "number_of_contigs",
		 min_len		=> 'minimum_length',
		 max_len		=> 'maximum_length',
		 avg_len		=> 'mean_length',
		 Q1			=> 'Q1',
		 Q2			=> 'Q2',
		 Q3			=> 'Q3',
		 sum_len		=> 'total_length',
		 N50			=> 'N50',
		 'GC(%)'		=> 'GC(%)'
		 );
	
	# STATS FOR ALL CONTIGS
	system_call("seqkit stats -j $threads -baT $contigs_fa 1> $seqkit_tmp");
	%tmp			= read_tsv2hashtable($seqkit_tmp,'file');
	my %conts_1	= ();
	$rowp		= $tmp{(keys %tmp)[0]};
	foreach my $k(sort keys %renamec){
		if( defined($rowp->{$k})){
			$conts_1{$renamec{$k}}		= $rowp->{$k};
		}
	}
		# number of ORFs
	system_call("seqkit seq -ni $orfs_fa 1> $ids2_tmp");
	$conts_1{'number_of_ORFs'}	= nlines($ids2_tmp);
	
	
	# STATS FOR CONTIGS > 500 BP
	my $min_len	= 500;
	system_call("seqkit seq -j $threads -g --min-len $min_len $contigs_fa | seqkit stats -j $threads -baT 1> $seqkit_tmp");
	%tmp			= read_tsv2hashtable($seqkit_tmp,'file');
	my %conts_2	= ();
	$rowp		= $tmp{(keys %tmp)[0]};
	foreach my $k(sort keys %renamec){
		if( defined($rowp->{$k})){
			$conts_2{$renamec{$k}}		= $rowp->{$k};
		}
	}
		# number of ORFs
	# <new code 18 Sep>
	system_call("seqkit seq -j $threads -nig --min-len $min_len $contigs_fa 1> $ids_tmp");
	system_call("seqkit seq -j $threads -nig $orfs_fa | cut -d'_' -f1 1> $ids2_tmp");
	if(nlines($ids_tmp)>0 && nlines($ids2_tmp)>0){
		system_call("csvtk grep -j $threads -f 1 -P $ids_tmp $ids2_tmp 1> $ids3_tmp");
		$conts_2{number_of_ORFs}	= nlines($ids3_tmp);
	}
	else{
		$conts_2{number_of_ORFs}	= 0;
	}
	# </new code>


	# STATS FOR CONTIGS > 1000 BP
	$min_len	= 1000;
	system_call("seqkit seq -j $threads -g --min-len $min_len $contigs_fa | seqkit stats -j $threads -baT 1> $seqkit_tmp");
	%tmp			= read_tsv2hashtable($seqkit_tmp,'file');
	my %conts_3	= ();
	$rowp		= $tmp{(keys %tmp)[0]};
	foreach my $k(sort keys %renamec){
		if( defined($rowp->{$k})){
			$conts_3{$renamec{$k}}		= $rowp->{$k};
		}
	}
		# number of ORFs
	# <new code 18 Sep>
	system_call("seqkit seq -j $threads -nig --min-len $min_len $contigs_fa 1> $ids_tmp");
	system_call("seqkit seq -j $threads -nig $orfs_fa | cut -d'_' -f1 1> $ids2_tmp");
	if(nlines($ids_tmp)>0 && nlines($ids2_tmp)>0){
		system_call("csvtk grep -j $threads -f 1 -P $ids_tmp $ids2_tmp 1> $ids3_tmp");
		$conts_3{number_of_ORFs}	= nlines($ids3_tmp);
	}
	else{
		$conts_3{number_of_ORFs}	= 0;
	}
	# </new code>	
	
	
	# COLLECT 2nd LEVEL HASHES TO 1st LEVEL HASH
	$stats_1{'median_length'}	= $stats_1{Q2};
	$stats_2{'median_length'}	= $stats_2{Q2};
	$stats_3{'median_length'}	= $stats_3{Q2};
	$stats_1{'keys'}				= \@read_stats_keys;
	$stats_2{'keys'}				= \@read_stats_keys;
	$stats_3{'keys'}				= \@read_stats_keys;
	$conts_1{'median_length'}	= $conts_1{Q2};
	$conts_2{'median_length'}	= $conts_2{Q2};
	$conts_3{'median_length'}	= $conts_3{Q2};
	$conts_1{'keys'}				= \@contig_stats_keys;
	$conts_2{'keys'}				= \@contig_stats_keys;
	$conts_3{'keys'}				= \@contig_stats_keys;
	$stats{'reads.trimmed'} 		= \%stats_1;		# reads after trimming
	$stats{'reads.hostflt'} 		= \%stats_2;		# reads after host filtering 
	$stats{'reads.assembled'} 	= \%stats_3;		# reads after assembling
	$stats{'contigs'} 			= \%conts_1;		# all contigs
	$stats{'contigs.500bp'} 		= \%conts_2;		
	$stats{'contigs.1000bp'} 	= \%conts_3;	
	my @tmp						= ('reads.trimmed','reads.hostflt','reads.assembled','contigs','contigs.500bp','contigs.1000bp');
	$stats{keys}					= \@tmp;
	
	# WRITE YAML
	if(defined($args{yaml})){
		open(OUT, ">$args{yaml}") or die "couldn't open file $args{yaml} $!";
		my $sp				= '   ';
		my ($val,$format)	= '';
		my $field_len		= 20;
		foreach my $k1(@{$stats{keys}}){
			print OUT "$k1:\n";
			foreach my $k2(@{$stats{$k1}->{keys}}){
				$val		= $stats{$k1}->{$k2};
				$format	= "$sp- \%-$field_len"."s\%s\n";
				print OUT sprintf($format, $k2.':'.$sp, $val);
			}
		}
		close(OUT);
	}

	# CLEANUP
	system("rm -f $r1_contigs $seqkit_tmp $ids_tmp $ids2_tmp");

	return %stats;
}

# USAGE:
# add_bphage_field(annot=>$annot, phfilter=>$phfilter, taxonomy=>$taxonomy, log=>$log, numth=>$threads, overwrite=>0)
# annot			: tsv-file with sequence annotation, MUST INCLUDE HEADERS: division,staxid/taxid
# phfilter		: tsv-file listing bphage names/taxids, MUST INCLUDE HEADER: taxid
# taxonomy		: path to NCBI taxonomy dump
# log			: log-file
# [numth]		: threads, [1]
# [overwrite]	: overwrite existing bphage-flag [true]
#
sub add_bphage_field{
	my $SIGNATURE	= "add_bphage_field()";
	if($VERBAL){
		print STDERR "\n\t$SIGNATURE\n";
	}
	
 	# Check input:
  	my (%args)			= @_;	
 	if( !defined($args{annot})){
 		die "ERROR: $SIGNATURE: missing input: annot";}
 	if( !defined($args{phfilter})){
 		die "ERROR: $SIGNATURE: missing input: phfilter";}
 	if( !defined($args{log})){
 		die "ERROR: $SIGNATURE: missing input: log";}
  	if( !defined($args{taxonomy})){
 		die "ERROR: $SIGNATURE: missing input: taxonomy";}
	
	# in/out:
	my $annot 				= $args{annot};
	my $phfilter				= $args{phfilter};
	my $log					= $args{log};
	my $taxonomy				= $args{taxonomy};
	my $numth				= defined($args{numth}) ? $args{numth} : 0;
	my $overwrite			= defined($args{overwrite}) ? $args{overwrite} : 1;
	# tmp:
	my $vi_taxid				= "$annot.vi.taxid.tmp";
	my $vi_taxid_lineage		= "$annot.vi.taxid_lineage.tmp";
	my $annot_tmp			= "$annot.tmp";
	# param:
	my $annot_taxidh			= undef;
	if(colind($annot,'staxid') > 0){
		$annot_taxidh		= 'staxid';
	}elsif(colind($annot,'taxid') > 0){
		$annot_taxidh		= 'taxid';
	}
	else{
		die "ERROR: $SIGNATURE: missing field staxid/taxid in $annot";
	}
	my $annot_taxidi			= colind($annot,$annot_taxidh)-1;


	system_call("cat $annot | ".
				"csvtk filter2 -tf '\$division==\"Viruses\"' | ".
				"csvtk cut -tlf $annot_taxidh | tail -n+2 | uniq 1> $vi_taxid 2>> $log");
	if(nlines($vi_taxid) == 0){
		# no viral hits in $annot > return
		system("rm -f $vi_taxid");
		return;
	}
	system_call("cat $vi_taxid | ".
				"taxonkit lineage -i1 -t --data-dir $taxonomy -j  $numth | ".
				"csvtk cut -tlf1,3 | ".
				"csvtk add-header -tn taxid,lineage 1> $vi_taxid_lineage 2>> $log" );
	
	my %taxid_lineage	= read_tsv2hash($vi_taxid_lineage,"taxid","lineage");
	my %phflt_taxids		= read_tsv2hash($phfilter,"taxid","taxid");
	my %phage_taxids		= ();
	foreach my $taxid(keys %taxid_lineage){
		my @lineage 	= split(/;/,$taxid_lineage{$taxid},-1);
		foreach my $tmp(@lineage){
			if(defined($phflt_taxids{$tmp})){
				$phage_taxids{$taxid}	= 1;
				last;
			}
		}
	}
	
	# adding bphage flag to $annot
	my $divisioni			= colind($annot,'division') -1; # colind returns 1-based index
		# search for an existing bphage-flag:
	my $bphagei				= colind($annot,'bphage') -1;	
	my $has_bphage_header	= ($bphagei >= 0) ? 1 : 0;
	
	open(IN,"<$annot") or die "$SIGNATURE: Can\'t open $annot: $!\n";	
	open(OUT,">$annot_tmp") or die "$SIGNATURE: Can\'t open $annot_tmp: $!\n";
	my $l=<IN>;
	chomp($l);
	if( $has_bphage_header){	
		print OUT $l,"\n";
	}
	else{
		print OUT $l,"\tbphage\n";
	}
	while($l=<IN>){
		chomp($l);
		my @sp 		= split(/\t/,$l,-1);
		my $taxid	= $sp[$annot_taxidi];
		my $bphage	= 'no';
		if($divisioni>0 && lc($sp[$divisioni]) eq 'phages'){
			$bphage		= 'yes';
		}
		elsif(defined($phage_taxids{$taxid})){
			$bphage		= 'yes';
		}
			# preserve any bphage-data in $annot
		if(!$has_bphage_header){
			push(@sp,$bphage);
		}
		elsif( $overwrite ){
			$sp[$bphagei]	= $bphage;
		}
		
		print OUT join("\t",@sp),"\n";
	}
	close(IN);close(OUT);
	system("mv $annot_tmp $annot");
	system("rm -f $vi_taxid $vi_taxid_lineage");
}


# USAGE:
#
# add_division_field(annot=>$annot, taxonomy=>$taxonomy, numth=>$threads)
#
# annot			: tsv-file with sequence annotation, MUST INCLUDE HEADERS: division,staxid/taxid
# taxonomy		: path to NCBI taxonomy dump
# [overwrite]	: overwrite existing division-field [true]
# [numth]		: threads, [1]
#
sub add_division_field{
	my $subid	= "add_division_field()";
	
	# in/out:
	my (%args)				= @_;	
	my $annot 				= $args{annot} || die "ERROR: $subid: missing input: annot";
	my $taxonomy				= $args{taxonomy} || die "ERROR: $subid: missing input: taxonomy";
	my $numth				= $args{numth} || 8;
	my $overwrite			= defined($args{overwrite})? $args{overwrite} : 1;
	my $nodesdmp				= "$taxonomy/nodes.dmp";
	my $mergeddmp			= "$taxonomy/merged.dmp";
	my $divdmp				= "$taxonomy/division.dmp";
	(-e $nodesdmp) || die "ERROR: $subid: missing file: $nodesdmp";
	(-e $mergeddmp) || die "ERROR: $subid: missing file: $mergeddmp";
	(-e $divdmp) || die "ERROR: $subid: missing file: $divdmp";

	# param:
	my $annot_taxidh			= undef;
	if(colind($annot,'staxid') > 0){
		$annot_taxidh		= 'staxid';
	}elsif(colind($annot,'taxid') > 0){
		$annot_taxidh		= 'taxid';
	}
	else{
		die "ERROR: $subid: missing field staxid/taxid in $annot";
	}
	my $annot_taxidi			= colind($annot,$annot_taxidh)-1;

	# tmp:
	my $annot_tmp			= "$annot.tmp";
	
	# read nodes.dmp
	my %taxid_div		= ();
	my $taxidi			= 0;
	my $divi				= 4;
	open(IN,"<$nodesdmp") or die "Can\'t open $nodesdmp: $!\n";
	my $ln				= 0;
	while(my $l=<IN>){
		$ln++;
        chomp($l);
		my @sp= split(/\t\|\t/,$l,-1);
		if(scalar(@sp)<($divi+1) ){
			print STDERR "WARNING: $subid: invalid format on line $ln in $nodesdmp: skipping\n";
			next;
		}
		$taxid_div{$sp[$taxidi]} = $sp[$divi];
	}
	close(IN);
	
	# read merged.dmp
	my %merged		= ();
	open(IN,"<$mergeddmp") or die "Can\'t open $mergeddmp: $!\n";
	$ln				= 0;
	while(my $l=<IN>){
		$ln++;
        chomp($l);
        $l =~ s/\t\|$//;	# remove trailing |tab;
		my @sp= split(/\t\|\t/,$l,2);
		if(scalar(@sp)<2 ){
			print STDERR "WARNING: $subid: invalid format on line $ln in $mergeddmp: skipping\n";
			next;
		}
		$merged{$sp[0]} = $sp[1];
	}
	close(IN);
	
	# read divisions.dmp
	my %divid_name		= ();
	open(IN,"<$divdmp") or die "Can\'t open $divdmp: $!\n";
	while(my $l=<IN>){
		$ln++;
        chomp($l);
		my @sp= split(/\t\|\t/,$l,-1);
		$divid_name{$sp[0]} = $sp[2];
	}
	close(IN);
	
	# link division
	my $divcoli				= colind($annot,'division') -1;	
	my $has_division_header	= ($divcoli >= 0) ? 1 : 0;
	open(IN,"<$annot") or die "ERROR: $subid: Can\'t open $annot: $!\n";	
	open(OUT,">$annot_tmp") or die "ERROR: $subid: Can\'t open $annot_tmp: $!\n";
	my $l=<IN>;
	chomp($l);
	if($has_division_header){
		print OUT $l,"\n";
	}
	else{
		print OUT $l,"\tdivision\n";
	}
	while($l=<IN>){
		chomp($l);
		my @sp 		= split(/\t/,$l,-1);
		my $taxid	= $sp[$annot_taxidi];
		if(defined($merged{$taxid})){
			$taxid		= $merged{$taxid} ;
		}		
		my $division= (defined($taxid_div{$taxid}) && defined($divid_name{$taxid_div{$taxid}}))? $divid_name{$taxid_div{$taxid}}: 'NA';
		if(!$has_division_header){
			push(@sp,$division);
		}
		elsif( $overwrite ){
			$sp[$divcoli]	= $division;
		}
		print OUT join("\t",@sp),"\n";
	}
	close(IN);close(OUT);
	system("mv $annot_tmp $annot");
}
