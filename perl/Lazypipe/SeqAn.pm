package Lazypipe::SeqAn;
use strict;
use warnings;
use Exporter;
use MIME::Base64;
use File::Basename;
#use lib dirname(__FILE__);
#use Utils;
#use Lazypipe::Utils;

our @ISA			= qw( Exporter );
our @EXPORT		= qw( annotate_blastp annotate_blastx annotate_blastn annotate_diamondx annotate_diamondp annotate_minimap annotate_sans annotate_hmmscan detect_orfs EM_loop EM_q_t_logprob filter_seqs filter_unseqs filter_tophits_SAM filter_tophits filter_host_reads filter_host_contigs cigar2alen cigar2pide cigar2qcov cigar2qlen cigar2rlen );
our @EXPORT_OK	= qw( );
our $VERBAL		= !1;
our $VERSION		= 3.0;


###
# SEQUENCE ANALYSIS FUNCTIONS FOR LAZYPIPE PROJECT
#
# Dependencies:
# blastn, blastp
# csvtk
# minimap2
# mga 
# [prod]
# runsanspanz.py
# seqkit
# taxonkit
# 
# Credit:
#
# Plyusnin, I., Vapalahti, O., Sironen, T., Kant, R., & Smura, T. (2023).
# Enhanced Viral Metagenomics with Lazypipe 2. Viruses, 15(2), 431.
#
# Contact: grp-lazypipe@helsinki.fi
###


## GLOBAL CONSTANTS
my @ANNOT_TABLE_FIELDS 	= split(/,/,"search,db,dbtype,qseqid,orf,qseqlen,sseqid,bitscore,alen,pident,qlen,qcov,slen,scov,staxid,sname,bphage,division",-1);


# USAGE:
# annotate_blastn(seqs=>$fasta, db=>$blastn_db, dbhits=>$dbhits, annot=>$annot, append=>!1, log=>$log,
#				numth=>$numth, min_bits=>$bits,taxonomy=>$taxonomy);
# Input params [optional]
# IN:
# seqs		: fasta sequences
# db			: blastn database
# OUT:
# annot		: tsv-file for formatted annotations
# [append]	: append to existing $annot file, [false]
# [dbhits]	: tsv-file for raw annotations, [$seqs.blastn.dbhits.tsv]
# [log]		: log-file, [$seqs.blastn.log]
# PARAMS:
# numth				: threads
# min_bits			: minimum bitscore to keep search result
# [retain_ties]		: include ties when filtering the top scoring hit for each query [true]
# [max_target_seqs] : [5]
# [filter_tophits]	: only report the top scoring hit for each query [true]
#
sub annotate_blastn{
	print STDERR "\n# ANNOTATE WITH BLASTN\n\n";

	my $subid			= "annotate_blastn()";
 	my (%args)			= @_;
 	# Check input:
 	if( !defined($args{seqs})){
 		die "ERROR: $subid: missing input: seqs";}
 	if( !defined($args{db})){
 		die "ERROR: $subid: missing input: db";}
 	if( !defined($args{annot})){
 		die "ERROR: $subid: missing input: annot";}
  	if( !defined($args{numth})){
 		die "ERROR: $subid: missing input: numth";}
 	if( !defined($args{min_bits})){
 		die "ERROR: $subid: missing input: min_bits";}

	# in:
	my $seqs 				= $args{seqs};
	my $db					= $args{db};
	# out:
	my $annot				= $args{annot};
	my $append				= defined($args{append}) ? $args{append}: 0;
	my $dbhits				= (defined($args{dbhits})) ? $args{dbhits}: "$seqs.blastn.dbhits.tsv";
	my $log 					= (defined($args{log})) ? $args{log} : "$seqs.blastn.log";
	# params:
	my $numth				= $args{numth};
	my $min_bits				= $args{min_bits};
	my $retain_ties			= defined($args{retain_ties}) ? $args{retain_ties} : 1;
	my $max_target_seqs		= defined($args{max_target_seqs}) ? $args{max_target_seqs} : 5;
	my $taxids				= defined($args{taxids}) ? "-taxids $args{taxids}" : '';
	my $filter_tophits		= defined($args{filter_tophits}) ? $args{filter_tophits} : 1;
	
	my @fields = split(/,/,"qseqid,sseqid,bitscore,length,pident,qlen,qstart,qend,slen,sstart,send,qcovs,staxid,stitle",-1);
	my $fields_space			= join(' ',@fields);
	my $fields_comma			= join(',',@fields);
	my @fields_lazy 			= @ANNOT_TABLE_FIELDS;
	my $fields_lazy_tab 		= join("\t",@fields_lazy);
	my $fields_lazy_comma	= join(',',@fields_lazy);
	my $params				= "-evalue 10 -max_target_seqs $max_target_seqs $taxids -outfmt '6 $fields_space ' ";
	my $search_name 			= "blastn";
	my $searchdb_name		= basename($db);
	my $dbtype				= "nucl";	
	
	# CHECK INPUT
	if( !defined($db) ){	
		die "ERROR: $subid: undefined blastn db";
	}
	if( !(-e "$db.ndb") ){ 
		die "ERROR: $subid: no blastn db: $db";
	}
	if( !(-e "$seqs")){ 
		die "ERROR: $subid: no query sequences: $seqs";
	}

	# start new annotation
	if(!$append){
		system("echo \"$fields_lazy_tab\" 1> $annot");
		system("rm -f $log");
		system("touch $log");
	}
	# append requested but $annot is empty
	if($append && nlines($annot)<1){
		system("echo \"$fields_lazy_tab\" 1> $annot");
	}
	
	if( nlines($seqs)<1){
		print STDERR "\tWARNING: $subid: empty query sequences: $seqs";
		return();
	}	

		# run blast
	system_call("blastn -num_threads $numth $params -db $db -query $seqs 1> $dbhits 2>> $log");
	if( nlines($dbhits) < 1){
		print STDERR "\n\t$subid: NO HITS FOR $seqs in $db\n\n";
		return();
	} 
	
		# add headers and filter by bitscore
	system_call("cat $dbhits | ".
				"csvtk add-header --ignore-illegal-row -tn $fields_comma | ".
				"csvtk filter -tf 'bitscore>$min_bits' 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");
	if( nlines($dbhits) < 2){
		print STDERR "\n\t$subid: NO HITS FOR $seqs in $db\n\n";
		return();
	}
	
		# sort and filter tophits
	system_call("cat $dbhits | ".
				"csvtk sort -t -k qseqid:N -k bitscore:nr 1> $dbhits.tmp 2>> $log");
	system_call("mv $dbhits.tmp $dbhits");
	if($filter_tophits){
		filter_tophits(dbhits=>$dbhits,dbhits_flt=>"$dbhits.tmp", qcol=>'qseqid', bitscol=>'bitscore', retain_ties=>$retain_ties);
		system("mv $dbhits.tmp $dbhits");
	}	
	
	if( nlines($dbhits) < 2){
		print STDERR "\n\t$subid: NO HITS FOR $seqs in $db\n\n";
		return();
	}
	
		# add search,db,qcov,scov, rename length>alen, stitle>sname, add clen (==qlen for blastn)		
	system_call("cat $dbhits | ".
				"csvtk mutate2 -te '\"$search_name\"' -n 'search' | ".
				"csvtk mutate2 -te '\"$searchdb_name\"' -n 'db' | ".
				"csvtk mutate2 -te '\"$dbtype\"' -n 'dbtype' | ".
				"csvtk mutate2 -te '\"NA\"' -n orf | ".
				"csvtk mutate2 -w5 -te '(\$qend>\$qstart) ? (\$qend-\$qstart)/\$qlen : (\$qstart-\$qend)/\$qlen ' -n qcov | ".
				"csvtk mutate2 -w5 -te '(\$send>\$sstart) ? (\$send-\$sstart)/\$slen : (\$sstart-\$send)/\$slen ' -n scov | ".
				"csvtk rename -tf length,stitle -n alen,sname | ".
				"csvtk mutate  -tf qlen -n qseqlen | ".
				"csvtk mutate2 -te '\"NA\"' -n bphage | ".
				"csvtk mutate2 -te '\"NA\"' -n division 1> $dbhits.tmp 2>> $log");
	system_call("mv $dbhits.tmp $dbhits");
				
									
	system_call("cat $dbhits | ".
				"csvtk cut -tlf $fields_lazy_comma | ".
				"csvtk del-header -t 1>> $annot 2>> $log" );
}

# USAGE:
# annotate_blastp(seqs=>$fasta, seqinfo=>$seqinfo, db=>$db, dbhits=>$dbhits, annot=>$annot, append=>!1, log=>$log,
#				numth=>$numth, min_bits=>$bits, taxonomy=>$taxonomy, orf_finder=>$orf_finder, min_orf_length=>$length);
#
# Input params [optional]:
#
# IN
# seqs		: fasta sequences
# seqinfo	: tsv-file with sequence info, MUST include fields 'seqid','length'
# db			: blastp database
# OUT
# annot		: tsv-file for formatted annotations
# [dbhits]	: tsv-file for raw annotations, [$seqs.blastp.dbhits.tsv]
# [append]	: append to existing $annot file, [false]
# [log]		: log-file, [$seqs.blastp.log]
# PARAMS
# numth		: threads
# min_bits	: minimum bitscore to keep search result
# [filter_tophits]  : keep only tophits for each query [true]
# [retain_ties]		: retain ties when filtering tophits [true]
# [orf_finder]		: software for orf detection [mga]
# [min_orf_length]	: minimum orf length [0]
# [max_target_seqs] : [5]
# [taxids]          : see blastp --taxids [undef]
#
sub annotate_blastp{
	print STDERR "\n# ANNOTATE WITH BLASTP\n\n";

 	my (%args)			= @_;
 	my $subid			= "annotate_blastp()";
 	
 	# Check input:
 	if( !defined($args{seqs})){
 		die "ERROR: $subid: missing input: seqs";}
 	if( !defined($args{seqinfo})){
 		die "ERROR: $subid: missing input: seqinfo";}
 	if( !defined($args{db})){
 		die "ERROR: $subid: missing input: db";}
 	if( !defined($args{annot})){
 		die "ERROR: $subid: missing input: annot";}
  	if( !defined($args{numth})){
 		die "ERROR: $subid: missing input: numth";}
 	if( !defined($args{min_bits})){
 		die "ERROR: $subid: missing input: min_bits";}	
 		
	# in:
	my $seqs 				= $args{seqs};
	my $seqsinfo				= $args{seqinfo};
	my $db					= $args{db};
	# out:
	my $orfs_aa				= "$seqs.orfs.aa.fa";
	my $orfs_nt				= "$seqs.orfs.nt.fa";
	my $annot				= $args{annot};
	my $append				= defined($args{append}) ? $args{append} : 0;	
	my $dbhits				= defined($args{dbhits}) ? $args{dbhits} : "$seqs.blastp.dbhits.tsv";
	my $log 					= defined($args{log}) ? $args{log} : "$seqs.blastp.log";
	# params:
	my $numth				= $args{numth};
	my $min_bits				= $args{min_bits};
	my $orf_finder			= defined($args{orf_finder}) ? $args{orf_finder}: "mga";
	my $min_orf_length		= defined($args{min_orf_length}) ? $args{min_orf_length} : 0;
	my $filter_tophits		= defined($args{filter_tophits}) ? $args{filter_tophits} : 1;
	my $retain_ties			= defined($args{retain_ties}) ? $args{retain_ties} : 1;
	my $taxids				= defined($args{taxids}) ? "-taxids $args{taxids}" : '';
	my $max_target_seqs		= defined($args{max_target_seqs}) ? $args{max_target_seqs} : 5;
	my @fields 				= split(/,/,"qseqid,saccver,bitscore,length,pident,qlen,qstart,qend,slen,sstart,send,qcovs,staxid,stitle",-1);
	my $fields_space			= join(' ',@fields);
	my $fields_comma			= join(',',@fields);
	my @fields_lazy 			= @ANNOT_TABLE_FIELDS;
	my $fields_lazy_tab 		= join("\t",@fields_lazy);
	my $fields_lazy_comma	= join(',',@fields_lazy);
	my $params				= "-evalue 10 -max_target_seqs $max_target_seqs $taxids  -outfmt '6 $fields_space ' ";
	my $search_name 			= "blastp";
	my $searchdb_name		= basename($db);
	my $dbtype				= "prot";

	# CHECK INPUT
	if( !defined($db) ){	
		die "ERROR: $subid: undefined db";
	}
	if( !(-e "$db.pdb") ){ 
		die "ERROR: $subid: no db: $db";
	}
	if( !(-e "$seqs")){ 
		die "ERROR: $subid: no query sequences: $seqs";
	}

	# start new annotation
	if(!$append){
		system("echo \"$fields_lazy_tab\" 1> $annot");
		system("rm -f $log");
		system("touch $log");
	}
	# append requested but $annot is empty
	if($append && nlines($annot)<1){
		system("echo \"$fields_lazy_tab\" 1> $annot");
	}	
	
	if( nlines($seqs)<1){
		print STDERR "\tWARNING: $subid: empty query sequences: $seqs\n";
		return();
	}	

	system("rm -f $orfs_aa $orfs_nt");
	detect_orfs(seqs=>$seqs, orfs_nt=>$orfs_nt, orfs_aa=>$orfs_aa, orf_finder=>$orf_finder, min_orf_length=>$min_orf_length,numth=>$numth);
	if( nlines($orfs_aa)<2 ){
		print STDERR "WARNING: $subid: no orfs found in $seqs\n";
		return();
	}

		# run blast
	system_call("blastp -num_threads  $numth $params -db $db -query $orfs_aa 1> $dbhits 2>> $log" );
	if( nlines($dbhits) < 1){
		print STDERR "\n\tNO HITS FOR $seqs in $db\n\n";
		return();
	}
	
		# add headers and filter by bitscore
	system_call("cat $dbhits | ".
				"csvtk add-header -tn $fields_comma | ".
				"csvtk filter -tf 'bitscore>$min_bits' 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");
	if( nlines($dbhits) < 2){
		print STDERR "\n\tNO HITS FOR $seqs in $db\n\n";
		return();
	}	
		
		# sort and filter tophits with ties
	system_call("cat $dbhits | ".
				"csvtk sort -t -k qseqid:N -k bitscore:nr 1> $dbhits.tmp 2>> $log");
	system("mv $dbhits.tmp $dbhits");
	if($filter_tophits){
		filter_tophits(dbhits=>"$dbhits",dbhits_flt=>"$dbhits.tmp", qcol=>'qseqid', bitscol=>'bitscore', retain_ties=>$retain_ties);
		system("mv $dbhits.tmp $dbhits");
	}
	
	if( nlines($dbhits) < 2){
		print STDERR "\n\tNO HITS FOR $seqs in $db\n\n";
		return();
	}
	
		# add search,db,qcov,scov, rename length>alen, stitle>sname
	system_call("cat $dbhits | ".
				"csvtk mutate2 -te '\"$search_name\"' -n 'search' | ".
				"csvtk mutate2 -te '\"$searchdb_name\"' -n 'db' | ".
				"csvtk mutate2 -te '\"$dbtype\"' -n 'dbtype' | ".
				"csvtk rename -tf qseqid -n qseqid_tmp | ".
				"csvtk mutate -tf qseqid_tmp -p '^([^_]+)' -n qseqid | ".
				"csvtk mutate -tf qseqid_tmp -p 'ORF=([^_]+)\$' -n orf | ".
				"csvtk mutate2 -w5 -te '(\$qend>\$qstart) ? (\$qend-\$qstart)/\$qlen : (\$qstart-\$qend)/\$qlen ' -n qcov | ".
				"csvtk mutate2 -w5 -te '(\$send>\$sstart) ? (\$send-\$sstart)/\$slen : (\$sstart-\$send)/\$slen ' -n scov | ".
				"csvtk rename -tf saccver,length,stitle -n sseqid,alen,sname | ".
				"csvtk mutate2 -te '\"NA\"' -n bphage | ".
				"csvtk mutate2 -te '\"NA\"' -n division 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");

		# add fields: seqlen (length of the nucleotide query)
	system_call("csvtk join -tj $numth -f 'qseqid;seqid' -L --na 'NA' $dbhits $seqsinfo | ".
				"csvtk rename -tf length -n qseqlen 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");
				

	system_call("cat $dbhits | ".
				"csvtk cut -tlf $fields_lazy_comma | ".
				"csvtk del-header -t 1>> $annot 2>> $log" );
}



# USAGE:
# annotate_blastx(seqs=>$fasta, seqinfo=>$seqinfo, db=>$db, dbhits=>$dbhits, annot=>$annot, append=>!1, log=>$log,
#				numth=>$numth, min_bits=>$bits, query_gencode=>1, taxids=10239);
#
# Input params [optional]:
#
# seqs             : nucleotide sequences in fasta format
# seqinfo          : tsv-file with sequence info, MUST include fields 'seqid','length'
# db               : blastx database
# annot            : tsv-file for formatted annotations
# [dbhits]         : tsv-file for raw annotations, [$seqs.blastp.dbhits.tsv]
# [append]         : append to existing $annot file, [false]
# [filter_tophits] : only report the top scoring hit for each query, [true]
# [retain_ties]    : include ties when filtering the top scoring hit for each query, [true]
# [log]            : log-file, [$seqs.blastx.log]
# numth            : threads
# bitscore         : bitscore cutoff
# [evalue]         : see blastx --evalue [10]
# [query_gencode]  : see blastx --query_gencode [1]
# [taxids]         : see blastx --taxids [undef]
# [max_target_seqs]: see blastx --max_target_seqs [5]
#
sub annotate_blastx{
	print STDERR "\n# ANNOTATE WITH BLASTX\n\n";

 	my (%args)			= @_;
 	my $subid			= "annotate_blastx()";
 	
 	# Check input:
 	if( !defined($args{seqs})){
 		die "ERROR: $subid: missing input: seqs";}
 	if( !defined($args{seqinfo})){
 		die "ERROR: $subid: missing input: seqinfo";}
 	if( !defined($args{db})){
 		die "ERROR: $subid: missing input: db";}
 	if( !defined($args{annot})){
 		die "ERROR: $subid: missing input: annot";}
  	if( !defined($args{numth})){
 		die "ERROR: $subid: missing input: numth";}
 	if( !defined($args{bitscore})){
 		die "ERROR: $subid: missing input: bitscore";}
 		
	# in:
	my $seqs 				= $args{seqs};
	my $seqsinfo				= $args{seqinfo};
	my $db					= $args{db};
	# out:
	my $annot				= $args{annot};
	my $dbhits				= defined($args{dbhits})? $args{dbhits} : "$seqs.blastx.dbhits.tsv";
	my $log 					= defined($args{log})? $args{log} : "$seqs.blastx.log";
	# params:
	my $append				= defined($args{append})? $args{append} : 0;
	my $filter_tophits		= defined($args{filter_tophits})? $args{filter_tophits} : 1;
	my $retain_ties			= defined($args{retain_ties})? $args{retain_ties} : 1;
	my $numth				= $args{numth};
	my $bitscore				= $args{bitscore	};
	my $evalue				= defined($args{evalue}) ? $args{evalue} : 10;
	my $query_gencode		= defined($args{query_gencode}) ? $args{query_gencode} : 1;
	my $taxids				= defined($args{taxids}) ? "-taxids $args{taxids}" : '';
	my $max_target_seqs		= defined($args{max_target_seqs}) ? $args{max_target_seqs} : 5;
	
	my @fields 				= split(/,/,"qseqid,saccver,bitscore,length,pident,qlen,qstart,qend,slen,sstart,send,qcovs,staxid,stitle",-1);
	my $fields_space			= join(' ',@fields);
	my $fields_comma			= join(',',@fields);
	my @fields_lazy 			= @ANNOT_TABLE_FIELDS;
	my $fields_lazy_tab 		= join("\t",@fields_lazy);
	my $fields_lazy_comma	= join(',',@fields_lazy);
	my $params				= " -evalue $evalue -query_gencode $query_gencode -max_target_seqs $max_target_seqs $taxids -outfmt '6 $fields_space ' ";
	my $search_name 			= "blastx";
	my $searchdb_name		= basename($db);
	my $dbtype				= "prot";

	# CHECK INPUT
	if( !defined($db) ){	
		die "ERROR: $subid: undefined db";
	}
	if( !(-e "$db.pdb") ){ 
		die "ERROR: $subid: no db: $db";
	}
	if( !(-e "$seqs")){ 
		die "ERROR: $subid: no query sequences: $seqs";
	}

	# start new annotation
	if(!$append){
		system("echo \"$fields_lazy_tab\" 1> $annot");
		system("rm -f $log");
		system("touch $log");
	}
	# append requested but $annot is empty
	if($append && nlines($annot)<1){
		system("echo \"$fields_lazy_tab\" 1> $annot");
	}	
	
	if( nlines($seqs)<1){
		print STDERR "\tWARNING: $subid: empty query sequences: $seqs\n";
		return();
	}	

	# run blast
	system_call("blastx -num_threads $numth $params -db $db -query $seqs 1> $dbhits 2>> $log" );
	if( nlines($dbhits) < 1){
		print STDERR "\n\tNO HITS FOR $seqs in $db\n\n";
		return();
	}
	
		# add headers, sort and filter by bitscore
	system_call("cat $dbhits | ".
				"csvtk add-header -tn $fields_comma | ".
				"csvtk sort -t -k qseqid:N -k bitscore:nr | ".
				"csvtk filter -tf 'bitscore>$bitscore' 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");
	if( nlines($dbhits) < 2){
		print STDERR "\n\tNO HITS FOR $seqs in $db\n\n";
		return();
	}	
		
		# sort and filter tophits with or without ties
	if($filter_tophits){
		filter_tophits(dbhits=>$dbhits,dbhits_flt=>"$dbhits.tmp", qcol=>'qseqid', bitscol=>'bitscore', retain_ties=>$retain_ties);
		system("mv $dbhits.tmp $dbhits");
		if( nlines($dbhits) < 2){
			print STDERR "\n\tNO HITS FOR $seqs in $db\n\n";
			return();
		}
	}

		# add search,db,qcov,scov, rename saccver>sseqid, length>alen, stitle>sname
	system_call("cat $dbhits | ".
				"csvtk mutate2 -te '\"$search_name\"' -n 'search' | ".
				"csvtk mutate2 -te '\"$searchdb_name\"' -n 'db' | ".
				"csvtk mutate2 -te '\"$dbtype\"' -n 'dbtype' | ".
				"csvtk rename -tf qseqid -n qseqid_tmp | ".
				"csvtk mutate -tf qseqid_tmp -p '^([^_]+)' -n qseqid | ".
				"csvtk mutate -tf qseqid_tmp -p 'ORF=([^_]+)\$' -n orf | ".
				"csvtk mutate2 -w5 -te '(\$qend>\$qstart) ? (\$qend-\$qstart)/\$qlen : (\$qstart-\$qend)/\$qlen ' -n qcov | ".
				"csvtk mutate2 -w5 -te '(\$send>\$sstart) ? (\$send-\$sstart)/\$slen : (\$sstart-\$send)/\$slen ' -n scov | ".
				"csvtk rename -tf saccver,length,stitle -n sseqid,alen,sname | ".
				"csvtk mutate2 -te '\"NA\"' -n bphage | ".
				"csvtk mutate2 -te '\"NA\"' -n division 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");

		# add fields: qseqlen (length of the nucleotide query)
	system_call("csvtk join -tj $numth -f 'qseqid;seqid' -L --na 'NA' $dbhits $seqsinfo | ".
				"csvtk rename -tf length -n qseqlen 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");
								

	system_call("cat $dbhits | ".
				"csvtk cut -tlf $fields_lazy_comma | ".
				"csvtk del-header -t 1>> $annot 2>> $log" );
}

# Annotate with diamond blastx
# USAGE:
# annotate_diamondx(seqs=>$fasta, db=>$db, dbhits=>$dbhits, ..)
#
# Input params [optional]:
#	IN:
# seqs				: nucleotide sequences in fasta format
# db					: diamond protein database
#	OUT:
# annot				: tsv-file for formatted annotations
# [append]			: append to existing $annot file, [false]
# [dbhits]			: tsv-file for raw annotations, [$seqs.diamondx.dbhits.tsv]
# [log]				: log-file, [$seqs.diamondx.log]
#	PARAMS:
# [filter_tophits]	: only report the top scoring hit for each query [true]
# [retain_ties]		: include ties when filtering the top scoring hit for each query [true]
# [numth]			: threads [8]
#	diamond blastx:
# [evalue]			: see diamond --evalue [10]
# [min_score]		: see diamond --min-score [undef]
# [query_gencode]	: see diamond --query_gencode [1]
# [min_orf]			: see diamond --min-orf [150]
# [taxonlist]		: see diamond --taxonlist [undef]
# [max_target_seqs]	: see diamond --max-target-seqs [5]
# [sensitivity]		: see diamond --fast/mid-sensitive/sensitive/more-sensitive/very-sensitive/ultra-sensitive [very-sensitive]
# [tmpdir]			: see diamond --tmpdir [undef]
#
sub annotate_diamondx{
	print STDERR "\n# ANNOTATE WITH DIAMOND BLASTX\n\n";

 	my (%args)			= @_;
 	my $subid			= "annotate_diamondx()";

 	# Check input:
 	if( !defined($args{seqs})){
 		die "ERROR: $subid: missing input: seqs";}
	if( !(-e $args{seqs})){ 
		die "ERROR: $subid: no query sequences: $args{seqs}";}
 	if( !defined($args{db})){
 		die "ERROR: $subid: missing input: db";}
 	if( !defined($args{annot})){
 		die "ERROR: $subid: missing input: annot";}
 		
	# in:
	my $seqs 				= $args{seqs};
	my $seqsinfo				= $args{seqinfo};
	my $db					= $args{db};
	# out:
	my $annot				= $args{annot};
	my $dbhits				= defined($args{dbhits}) ? $args{dbhits} : "$seqs.diamondx.dbhits.tsv";
	my $log 					= defined($args{log}) ? $args{log} : "$seqs.dimondx.log";
	# constants
	#my @fields 				= split(/,/,"qseqid,sseqid,bitscore,length,pident,qlen,qstart,qend,slen,sstart,send,qcovs,staxid,stitle",-1);
	my @fields				= split(/,/,"qseqid,sseqid,bitscore,length,pident,qlen,qstart,qend,slen,sstart,send,staxids,stitle",-1);
	my $fields_space			= join(' ',@fields);
	my $fields_comma			= join(',',@fields);
	my @fields_lazy 			= @ANNOT_TABLE_FIELDS;
	my $fields_lazy_tab 		= join("\t",@fields_lazy);
	my $fields_lazy_comma	= join(',',@fields_lazy);
	my $search_name 			= "diamond.blastx";
	my $searchdb_name		= basename($db);
	my $dbtype				= "prot";
	# input params:
	my $append				= defined($args{append}) ? $args{append} : 0;	
	my $numth				= defined($args{numth}) ? $args{numth} : 8;
	my $filter_tophits		= defined($args{filter_tophits}) ? $args{filter_tophits} : 1;
	my $retain_ties			= defined($args{retain_ties}) ? $args{retain_ties} : 1;
		# diamond:
	my $evalue				= defined($args{evalue}) ? $args{evalue} : "10";
	my $min_score			= defined($args{min_score}) ? $args{min_score} : undef;
	my $min_orf				= defined($args{min_orf}) ? $args{min_orf} : 150;
	my $query_gencode		= defined($args{query_gencode}) ? $args{query_gencode} : 1;
	my $taxonlist			= defined($args{taxonlist}) ? $args{taxonlist} : undef;
	my $max_target_seqs		= defined($args{max_target_seqs}) ? $args{max_target_seqs} : 5;
	my $sensitivity			= defined($args{sensitivity}) ? $args{sensitivity} : "very-sensitive";
	my $tmpdir				= defined($args{tmpdir})? $args{tmpdir} : undef;
	my $diamond_params		= "--evalue $evalue ".
							(defined($min_score)? "--min-score $min_score ": "").
							(defined($taxonlist)? "--taxonlist $taxonlist ": "").
							(defined($tmpdir)? "--tmpdir $tmpdir ": "").
							"--min-orf $min_orf --query-gencode $query_gencode ".
							"--max-target-seqs $max_target_seqs ".
							"--$sensitivity ".
							"--outfmt 6 $fields_space ";
	# START WORKING
	# start new annotation
	if(!$append){
		system("echo \"$fields_lazy_tab\" 1> $annot");
		system("rm -f $log");
		system("touch $log");
	}
	# append requested but $annot is empty
	if($append && nlines($annot)<1){
		system("echo \"$fields_lazy_tab\" 1> $annot");
	}	
	
	if( nlines($seqs)<1){
		print STDERR "\tWARNING: $subid: empty query sequences: $seqs\n";
		return();
	}	

	# run diamond blastx
	system_call("diamond blastx --threads $numth $diamond_params --db $db --query $seqs --out $dbhits 2>> $log" );
	if( nlines($dbhits) < 1){
		print STDERR "\n\tWARNING: $subid: no hits for $seqs in $db\n\n";
		return();
	}
	
		# add headers, sort + filter bitscore
	system_call("cat $dbhits | ".
				"csvtk add-header -tn $fields_comma | ".
				"csvtk sort -j $numth -t -k qseqid:N -k bitscore:nr | ".
				"csvtk filter2 -j $numth -tf '\$bitscore>=$args{min_score}' 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");
		
		# sort and filter tophits with or without ties
	if($filter_tophits){
		filter_tophits(dbhits=>$dbhits,dbhits_flt=>"$dbhits.tmp", qcol=>'qseqid', bitscol=>'bitscore', retain_ties=>$retain_ties);
		system("mv $dbhits.tmp $dbhits");
		if( nlines($dbhits) < 2){
			print STDERR "\n\tWARNING: $subid: no hits for $seqs in $db\n\n";
			return();
		}
	}

		# add		: search,db,dbtype
		# calc		: qcov,scov
		# rename		: length>alen, staxids>staxid, stitle>sname
		# copy		: qlen > qseqlen
		# add NA		: orf, bphage
	system_call("cat $dbhits | ".
				"csvtk mutate2 -te '\"$search_name\"' -n 'search' | ".
				"csvtk mutate2 -te '\"$searchdb_name\"' -n 'db' | ".
				"csvtk mutate2 -te '\"$dbtype\"' -n 'dbtype' | ".
				"csvtk mutate2 -w5 -te '(\$qend>\$qstart) ? (\$qend-\$qstart)/\$qlen : (\$qstart-\$qend)/\$qlen ' -n qcov | ".
				"csvtk mutate2 -w5 -te '(\$send>\$sstart) ? (\$send-\$sstart)/\$slen : (\$sstart-\$send)/\$slen ' -n scov | ".
				"csvtk rename -tf length,staxids,stitle -n alen,staxid,sname | ".
				"csvtk mutate2 -te '\$qlen' -n qseqlen | ".
				"csvtk mutate2 -te '\"NA\"' -n orf | ".
				"csvtk mutate2 -te '\"NA\"' -n bphage | ".
				"csvtk mutate2 -te '\"NA\"' -n division 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");
					

	system_call("cat $dbhits | ".
				"csvtk cut -tlf $fields_lazy_comma | ".
				"csvtk del-header -t 1>> $annot 2>> $log" );
}

# Annotate with diamond blastp
# USAGE:
# annotate_diamondp(seqs=>$fasta, db=>$db, dbhits=>$dbhits, ..)
#
# Input params [optional]:
#	IN:
# seqs				: fasta sequences
# seqinfo			: tsv-file with sequence info, MUST include fields 'seqid','length'
# db					: diamond protein database
#	OUT:
# annot				: tsv-file for formatted annotations
# [append]			: append to existing $annot file, [false]
# [dbhits]			: tsv-file for raw annotations, [$seqs.diamondp.dbhits.tsv]
# [log]				: log-file, [$seqs.diamondp.log]
#	PARAMS:
# [filter_tophits]	: only report the top scoring hit for each query [true]
# [retain_ties]		: include ties when filtering the top scoring hit for each query [true]
# [numth]			: threads [8]
# [orf_finder]		: software for orf detection, [mga]
#	diamond blastp:
# [evalue]			: see diamond --evalue [10]
# [min_score]		: see diamond --min-score [undef]
# [query_gencode]	: see diamond --query_gencode [1]
# [min_orf]			: see diamond --min-orf [300]
# [taxonlist]		: see diamond --taxonlist [undef]
# [max_target_seqs]	: see diamond --max-target-seqs [5]
# [sensitivity]		: see diamond --fast/mid-sensitive/sensitive/more-sensitive/very-sensitive/ultra-sensitive [very-sensitive]
# [tmpdir]			: see diamond --tmpdir [undef]
#
sub annotate_diamondp{
	print STDERR "\n# ANNOTATE WITH DIAMOND BLASTP\n\n";

 	my (%args)			= @_;
 	my $subid			= "annotate_diamondp()";
 	
 	# Check input:
 	if( !defined($args{seqs})){
 		die "ERROR: $subid: missing input: seqs";}
	if( !(-e $args{seqs})){ 
		die "ERROR: $subid: no query sequences: $args{seqs}";}
 	if( !defined($args{seqinfo})){
 		die "ERROR: $subid: missing input: seqinfo";}
 	if( !defined($args{db})){
 		die "ERROR: $subid: missing input: db";}
 	if( !defined($args{annot})){
 		die "ERROR: $subid: missing input: annot";}
 		
	# in:
	my $seqs 				= $args{seqs};
	my $seqsinfo				= $args{seqinfo};
	my $db					= $args{db};
	# out:
	my $orfs_aa				= "$seqs.orfs.aa.fa";
	my $orfs_nt				= "$seqs.orfs.nt.fa";
	my $annot				= $args{annot};
	my $dbhits				= defined($args{dbhits}) ? $args{dbhits} : "$seqs.diamondp.dbhits.tsv";
	my $log 					= defined($args{log}) ? $args{log} : "$seqs.dimondp.log";
	# constants
	#my @fields 				= split(/,/,"qseqid,sseqid,bitscore,length,pident,qlen,qstart,qend,slen,sstart,send,qcovs,staxid,stitle",-1);
	my @fields				= split(/,/,"qseqid,sseqid,bitscore,length,pident,qlen,qstart,qend,slen,sstart,send,staxids,stitle",-1);
	my $fields_space			= join(' ',@fields);
	my $fields_comma			= join(',',@fields);
	my @fields_lazy 			= @ANNOT_TABLE_FIELDS;
	my $fields_lazy_tab 		= join("\t",@fields_lazy);
	my $fields_lazy_comma	= join(',',@fields_lazy);
	my $search_name 			= "diamond.blastp";
	my $searchdb_name		= basename($db);
	my $dbtype				= "prot";
	# input params:
	my $append				= defined($args{append}) ? $args{append} : 0;	
	my $numth				= defined($args{numth}) ? $args{numth} : 8;
	my $orf_finder			= defined($args{orf_finder}) ? $args{orf_finder}: "mga";
	my $filter_tophits		= defined($args{filter_tophits}) ? $args{filter_tophits} : 1;
	my $retain_ties			= defined($args{retain_ties}) ? $args{retain_ties} : 1;
	my $query_gencode		= defined($args{query_gencode}) ? $args{query_gencode} : 1;
	my $min_orf				= defined($args{min_orf}) ? $args{min_orf} : 300;
		# diamond:
	my $evalue				= defined($args{evalue}) ? $args{evalue} : "10";
	my $min_score			= defined($args{min_score}) ? $args{min_score} : undef;
	my $max_target_seqs		= defined($args{max_target_seqs}) ? $args{max_target_seqs} : 5;
	my $taxonlist			= defined($args{taxonlist}) ? $args{taxonlist} : undef;
	my $sensitivity			= defined($args{sensitivity}) ? $args{sensitivity} : "very-sensitive";
	my $tmpdir				= defined($args{tmpdir})? $args{tmpdir} : undef;
	my $diamond_params		= "--evalue $evalue ".
							(defined($min_score)? "--min-score $min_score ": "").
							(defined($taxonlist)? "--taxonlist $taxonlist ": "").
							(defined($tmpdir)? "--tmpdir $tmpdir ": "").
							"--max-target-seqs $max_target_seqs ".
							"--$sensitivity ".
							"--outfmt 6 $fields_space ";
	
	# START WORKING
	# start new annotation
	if(!$append){
		system("echo \"$fields_lazy_tab\" 1> $annot");
		system("rm -f $log");
		system("touch $log");
	}
	# append requested but $annot is empty
	if($append && nlines($annot)<1){
		system("echo \"$fields_lazy_tab\" 1> $annot");
	}	
	if( nlines($seqs)<1){
		print STDERR "\tWARNING: $subid: empty query sequences: $seqs\n";
		return();
	}	

	system("rm -f $orfs_aa $orfs_nt");
	detect_orfs(seqs=>$seqs, orfs_nt=>$orfs_nt, orfs_aa=>$orfs_aa, orf_finder=>$orf_finder, min_orf_length=>$min_orf, numth=>$numth);
	if( nlines($orfs_aa)<2 ){
		print STDERR "WARNING: $subid: no orfs found in $seqs\n";
		return();
	}

	# run diamond blastp
	system_call("diamond blastp --threads $numth $diamond_params --db $db --query $orfs_aa 1> $dbhits 2>> $log" );
	
	if( nlines($dbhits) < 1){
		print STDERR "\n\tWARNING: $subid: no hits for $seqs in $db\n\n";
		return();
	}
	# add headers, sort and filter by bitscore
	system_call("cat $dbhits | ".
				"csvtk add-header -tn $fields_comma | ".
				"csvtk sort -t -k qseqid:N -k bitscore:nr | ".
				"csvtk filter -tf 'bitscore>=$min_score' 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");

	# sort and filter tophits with or without ties
	if($filter_tophits){
		filter_tophits(dbhits=>$dbhits,dbhits_flt=>"$dbhits.tmp", qcol=>'qseqid', bitscol=>'bitscore', retain_ties=>$retain_ties);
		system("mv $dbhits.tmp $dbhits");
	}
	
	if( nlines($dbhits) < 2){
		print STDERR "\n\tWARNING: $subid: no hits for $seqs in $db\n\n";
		return();
	}
	# add		: search,db,dbtype
	# calc		: qcov,scov
	# rename		: length>alen, staxids>staxid, stitle>sname
	# copy		: qlen > qseqlen
	# add NA		: orf, bphage
	system_call("cat $dbhits | ".
				"csvtk mutate2 -te '\"$search_name\"' -n 'search' | ".
				"csvtk mutate2 -te '\"$searchdb_name\"' -n 'db' | ".
				"csvtk mutate2 -te '\"$dbtype\"' -n 'dbtype' | ".
				"csvtk rename -tf qseqid -n qseqid_tmp | ".
				"csvtk mutate -tf qseqid_tmp -p '^([^_]+)' -n qseqid | ".
				"csvtk mutate -tf qseqid_tmp -p 'ORF=([^_]+)\$' -n orf | ".
				"csvtk mutate2 -w5 -te '(\$qend>\$qstart) ? (\$qend-\$qstart)/\$qlen : (\$qstart-\$qend)/\$qlen ' -n qcov | ".
				"csvtk mutate2 -w5 -te '(\$send>\$sstart) ? (\$send-\$sstart)/\$slen : (\$sstart-\$send)/\$slen ' -n scov | ".
				"csvtk rename -tf length,staxids,stitle -n alen,staxid,sname | ".
				"csvtk mutate2 -te '\"NA\"' -n bphage | ".
				"csvtk mutate2 -te '\"NA\"' -n division 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");			
						
	# add fields: seqlen (length of the nucleotide query)
	system_call("csvtk join -tj $numth -f 'qseqid;seqid' -L --na 'NA' $dbhits $seqsinfo | ".
				"csvtk rename -tf length -n qseqlen 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");
								
	system_call("cat $dbhits | ".
				"csvtk cut -tlf $fields_lazy_comma | ".
				"csvtk del-header -t 1>> $annot 2>> $log" );
}


#
# USAGE:
# annotate_minimap(seqs=>$fasta, db=>$db, dbhits=>$paf, annot=>$annot, append=>!1, retain_ties=>1, log=>$log,
#					min_bits=>$min_bits,secondary=>'yes',taxonomy=>$taxonomy,$numth)
#
# Input params [optional]:
#
# seqs		: fasta sequences
# db			: minimap index
# annot		: tsv-file for formatted annotations
# [append]	: append to existing $annot file [false]
# [dbhits]	: paf-file for raw annotation [$seqs.minimap.dbhits.paf]
# [x]       : minimap2 -x option [asm20]
# [secondary] : minimap2 --secondary option [yes]
# [retain_ties] : retain ties in annotation [true]
# [tophits] : retain only tophits for each query [true]
# [log]		: log-file, [$seqs.minimap.log]
# min_bits	: minimum bitscore to keep search result
# numth		: threads
# taxonomy	: path to NCBI taxonomy dump
#
sub annotate_minimap{
	print STDERR "\n# ANNOTATE WITH MINIMAP2\n\n";

 	my (%args)			= @_;
 	my $subid			= "annotate_minimap()";
 	
 	# Check input:
 	if( !defined($args{seqs})){
 		die "ERROR: $subid: missing input: seqs";}
 	if( !defined($args{db})){
 		die "ERROR: $subid: missing input: db";}
 	if( !defined($args{annot})){
 		die "ERROR: $subid: missing input: annot";}
  	if( !defined($args{min_bits})){
 		die "ERROR: $subid: missing input: min_bits";}
 	if( defined($args{secondary})
 		&& !($args{secondary} eq 'yes' || $args{secondary} eq 'no')){
 		die "ERROR: $subid: invalid option: secondary=>$args{secondary}";
 		}
  	if( !defined($args{numth})){
 		die "ERROR: $subid: missing input: numth";}
 	if( !defined($args{taxonomy})){
 		die "ERROR: $subid: missing input: taxonomy";}

	# in:
	my $seqs 			= $args{seqs};
	my $db				= $args{db};
	my $acc2taxid		= "";
	# out:
	my $annot			= $args{annot};
	my $dbhits_paf		= defined($args{dbhits}) ? $args{dbhits} : "$seqs.minimap.dbhits.paf";
	my $log 				= defined($args{log}) ? $args{log} : "$seqs.minimap.log";
	# pars:
	my $append			= defined($args{append}) ? $args{append} : 0;	
	my $min_bits			= $args{min_bits};
	my $x				= defined($args{'x'}) ? $args{'x'} : 'asm20';
	my $secondary		= defined($args{secondary}) ? $args{secondary} : 'yes';
	my $retain_ties		= defined($args{retain_ties}) ? $args{retain_ties} : 1;
	my $filter_tophits	= defined($args{filter_tophits}) ? $args{filter_tophits} : 1;
	my $numth			= $args{numth};	
	my $taxonomy			= $args{taxonomy};
	
	my $search_name 		= "Minimap2";
	my $searchdb_name	= basename($db);
	my $dbtype			= "nucl";
	my @fields_lazy 		= @ANNOT_TABLE_FIELDS;
	my $fields_lazy_tab = join("\t",@fields_lazy);
	my $fields_lazy_comma= join(',',@fields_lazy);	
	
	# CHECK INPUT DB & SEQS
	if($db =~ m/\.mmi$/){
		$acc2taxid		= $db;
		$acc2taxid		=~ s/\.mmi$/\.acc2taxid/;}
	else{
		$acc2taxid		= "$db.acc2taxid";
		if(-e "$db.mmi"){$db = "$db.mmi"; }
	}
	if( !defined($db) || !(-e $db)){	
		die "ERROR: $subid: no minimap2 db: $db";
	}
	if( !(-e $acc2taxid)){
		die "ERROR: $subid: no acc2taxid map: $acc2taxid";
	}	
	if( !(-e $seqs)){ 
		die "ERROR: $subid: no input sequences: $seqs";
	}
	
	# start new annotation
	if(!$append){
		system("echo \"$fields_lazy_tab\" 1> $annot");
		system("rm -f $log");
		system("touch $log");
	}
	# append requested but $annot is empty
	if($append && nlines($annot)<1){
		system("echo \"$fields_lazy_tab\" 1> $annot");
	}
	
	if( nlines($seqs)<1){
		print STDERR "WARNING: $subid: empty query sequences: $seqs";
		return();
	}	
		
	## RUN MINIMAP
	system_call("minimap2 -I 6g -t $numth -x $x --secondary $secondary --cs -s $min_bits $db $seqs 1> $dbhits_paf 2>> $log" );
	if( nlines("$dbhits_paf") < 1){
		print STDERR "\n\tNO HITS FOR $seqs in $db\n\n";
		system("rm -f $dbhits_paf.*.tmp");
		return();
	}
	
	## PARSE
	my $bitscol = mcolind($dbhits_paf,"AS:i:[\-0-9]+");
	if( $bitscol<0 ){
		die "ERROR: $subid: no AS:i:num field returned by minimap2\n";
	}
	system_call("cut -f 1-11,$bitscol $dbhits_paf | ".
				"csvtk add-header -tn 'qseqid,qlen,qstart,qend,strand,sseqid,slen,sstart,send,matches,alen,bitscore.tmp' | ".
				"csvtk mutate -tf bitscore.tmp -p '^AS:i:([+-]?[0-9]+)' -n bitscore | ".
				"csvtk filter2 -tf '\$bitscore>=$min_bits' | ".
				"csvtk cut -tf '-bitscore.tmp' 1> $dbhits_paf.1.tmp 2>> $log" );
	if( nlines("$dbhits_paf.1.tmp") < 2){
		print STDERR "\n\tNO HITS FOR $seqs in $db\n\n";
		system("rm -f $dbhits_paf.*.tmp");
		return();
	}			
				
		# get top-hits including ties
	system_call("cat $dbhits_paf.1.tmp | ".
				"csvtk sort -t -j $numth -k qseqid:N -k bitscore:nr 1> $dbhits_paf.2.tmp 2>> $log");
	if($filter_tophits){
		filter_tophits(dbhits=>"$dbhits_paf.2.tmp",dbhits_flt=>"$dbhits_paf.1.tmp", qcol=>'qseqid', bitscol=>'bitscore', retain_ties=>$retain_ties);
	}
	else{
		system_call("mv $dbhits_paf.2.tmp $dbhits_paf.1.tmp");
	}
	
				
		# add fields: clen,qcov,scov,pident,search,db
	system_call("cat $dbhits_paf.1.tmp | ".
				"csvtk mutate2 -te '\"NA\"' -n 'orf' | ".
				"csvtk mutate  -tf qlen -n qseqlen | ".
				"csvtk mutate2 -w5 -te '(\$qend-\$qstart)/\$qlen' -n qcov | ".
				"csvtk mutate2 -w5 -te '(\$send-\$sstart)/\$slen' -n scov | ".
				"csvtk mutate2 -w3 -te '\$matches/\$alen*100' -n pident | ".
				"csvtk mutate2 -te '\"$search_name\"' -n 'search' | ".
				"csvtk mutate2 -te '\"$searchdb_name\"' -n 'db' | ".
				"csvtk mutate2 -te '\"$dbtype\"' -n 'dbtype' | ".
				"csvtk mutate2 -te '\"NA\"' -n bphage | ".
				"csvtk mutate2 -te '\"NA\"' -n division 1> $dbhits_paf.2.tmp 2>> $log" );
	
		# add staxid
	my $lastcol = ncol("$dbhits_paf.2.tmp") + 1;
	system_call("csvtk join -tj $numth -f 'sseqid;1'  $dbhits_paf.2.tmp $acc2taxid | ".
				"csvtk rename -tf$lastcol -n staxid 1> $dbhits_paf.1.tmp 2>> $log" );
		# add sname
	system_call("head -n1 $dbhits_paf.1.tmp | tr '\\n' '\\t' 1> $dbhits_paf.2.tmp; echo 'sname' 1>> $dbhits_paf.2.tmp 2>> $log" );
	my $staxid_col = colind("$dbhits_paf.1.tmp","staxid"); 
	system_call("csvtk del-header -t $dbhits_paf.1.tmp | ".
				"taxonkit lineage -i$staxid_col --no-lineage -n --data-dir $taxonomy -j  $numth 1>> $dbhits_paf.2.tmp 2>> $log");
		# write annotation
	system_call("cat $dbhits_paf.2.tmp | ".
				"csvtk cut -tlf $fields_lazy_comma |  ".
				"csvtk del-header -t  1>> $annot 2>> $log" );

	system("rm -f $dbhits_paf.*.tmp");
}

# USAGE:
# annotate_sans(seqs=>$seqs, dbhits=>$dbhits, annot=>$annot, append=>1, log=>$log_sans,
#				min_bits=>$min_bits,numth=>$numth,taxonomy=>$taxonomy,orf_finder=$orf_finder,min_orf_length=$length)
#
# Input params [optional]:
#
# IN:
# seqs			: fasta sequences
# seqinfo		: tsv-file with sequence info, MUST include fields 'seqid','length'
# OUT:
# annot			: tsv-file for formatted annotations
# [append]		: append to existing $annot file, [false]
# [dbhits]		: tsv-file for raw annotations, [$seqs.sans.dbhits.tsv]
# [log]			: log-file, [$seqs.sans.log]
# PARAMS:
# min_bits			: minimum bitscore to keep search result
# numth				: threads
# [orf_finder]		: software for orf detection, [mga]
# [min_orf_length]	: minimum orf length, [72]
# [retain_ties]		: include ties when filtering the top scoring hit for each query [true]
#
sub annotate_sans{
	print STDERR "\n# ANNOTATE WITH SANS\n\n";
	
 	my (%args)			= @_;
 	my $subid			= "annotate_sans()";

	# in:
	my $seqs		 		= $args{seqs}  	|| die "ERROR: $subid: missing input: seqs";
	my $seqsinfo			= $args{seqinfo} || die "ERROR: $subid: missing input: seqinfo";
	# out:
	my $orfs_nt			= "$seqs.orfs.nt.fa";
	my $orfs_aa			= "$seqs.orfs.aa.fa";
	my $annot			= $args{annot}	|| die "ERROR: $subid: missing input: annot";
	my $append			= defined($args{append}) ? $args{append} : 0;	
	my $dbhits			= defined($args{dbhits}) ? $args{dbhits} : "$seqs.sans.dbhits.tsv";
	my $log 				= defined($args{log}) ? $args{log} : "$seqs.sans.log";
	# params
	my $min_bits			= $args{min_bits} || die "ERROR: $subid: missing input: min_bits";
	my $numth			= $args{numth} || die "ERROR: $subid: missing input: numth";
	my $orf_finder		= defined($args{orf_finder}) ? $args{orf_finder}: "mga";
	my $min_orf_length	= defined($args{min_orf_length}) ? $args{min_orf_length} : 72;
	my $retain_ties		= defined($args{retain_ties})? $args{retain_ties} : 1;
	
	my $sans_params		= "-m SANStopHtaxid --SANS_H 5 -R ";
	my $search_name		= "SANSparallel";
	my $searchdb_name	= "UniProtKB";
	my $dbtype 			= "prot";
	my @fields_lazy 		= @ANNOT_TABLE_FIELDS;
	my $fields_lazy_tab	= join("\t",@fields_lazy);
	my $fields_lazy_comma= join(',',@fields_lazy);
	
	# CHECK INPUT
	if( !(-e $seqs)){ 
		die "ERROR: $subid: no input sequences: $seqs";
	}

	# start new annotation
	if(!$append){
		system("echo \"$fields_lazy_tab\" 1> $annot");
		system("rm -f $log");
		system("touch $log");
	}
	# append requested but $annot is empty
	if($append && nlines($annot)<1){
		system("echo \"$fields_lazy_tab\" 1> $annot");
	}
	
	if( nlines($seqs)<1){
		print STDERR "WARNING: $subid: empty query sequences: $seqs";
		return();
	}
	
	system("rm -f $orfs_aa $orfs_nt");
	detect_orfs(seqs=>$seqs, orfs_nt=>$orfs_nt, orfs_aa=>$orfs_aa, orf_finder=>$orf_finder, min_orf_length=>$min_orf_length, numth=>$numth);
	
	# run sans
	system_call("runsanspanz.py $sans_params -i $orfs_aa -o $dbhits &> $log" );
	if( nlines("$dbhits") < 2){
		print STDERR "\n\tSANSparallel: NO HITS FOR $seqs\n\n";
		return();
	}	
	
		# filter tophits including ties
	system_call("cat $dbhits | ".
				"csvtk filter2 -tf '\$bits>$min_bits' | ".
				"csvtk filter2 -tf '\$taxid!=\"n.d.\"' | ".
				"csvtk rename -tf bits,taxid,qpid,spid,lali,desc -n bitscore,staxid,qseqid,sseqid,alen,sname 1> $dbhits.tmp 2>> $log");
	system_call("mv $dbhits.tmp $dbhits");
	filter_tophits(dbhits=>"$dbhits", dbhits_flt=>"$dbhits.tmp", qcol=>"qseqid", bitscol=>"bitscore", retain_ties=>$retain_ties);
	system("mv $dbhits.tmp $dbhits");
	if( nlines("$dbhits") < 2){
		print STDERR "\n\tSANSparallel: NO HITS FOR $seqs\n\n";
		return();
	}		
				
		# add fields: qlen,slen,contig,orf,pident,search,db,dbtype
	system_call("cat $dbhits | ".
				"csvtk mutate2 -te '\"NA\"' -n 'qlen' | ".
				"csvtk mutate2 -te '\"NA\"' -n 'slen' | ".
				"csvtk rename -tf qseqid -n qseqid_tmp | ".
				"csvtk mutate -tf qseqid_tmp -p '^([^_]+)' -n qseqid | ".
				"csvtk mutate -tf qseqid_tmp -p 'ORF=([^_]+)\$' -n orf | ".
				"csvtk mutate2 -w3 -te '\$pide*100' -n pident | ".
				"csvtk mutate2 -te '\"$search_name\"' -n 'search' | ".
				"csvtk mutate2 -te '\"$searchdb_name\"' -n 'db' | ".
				"csvtk mutate2 -te '\"$dbtype\"' -n 'dbtype' | ".
				"csvtk mutate2 -te '\"NA\"' -n bphage | ".
				"csvtk mutate2 -te '\"NA\"' -n division 1> $dbhits.tmp 2>> $log");
	system_call("mv $dbhits.tmp $dbhits");
		
		# add fields: qseqlen (contig-length)
	system_call("csvtk join -tj $numth -f 'qseqid;seqid' -L --na 'NA' $dbhits $seqsinfo | ".
				"csvtk rename -tf length -n qseqlen 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");	
				
	system_call("cat $dbhits | ".
				"csvtk cut -tlf $fields_lazy_comma | ".
				"csvtk del-header -t 1>> $annot 2>> $log");
}


# USAGE:
# annotate_hmmscan(
#		annot=>$annot, db=>$db, dbtaxid=>$dbacc2taxid, seqs=>$seqs, seqinfo=>$seqinfo, taxonomy=>$taxonomy,
#		append=>1, dbhits=>$dbhits, log=>$log_hmmscan, 
#		min_bits=>10, max_eval=>0.001, numth=>8,
#		min_orf_length=$length, orf_finder=$orf_finder)
#
# NOTE: This function is a PROTOTYPE
#
# Input params [optional]:
#
# annot			: tsv-file for formatted annotations
# db				: indexed hmmer database (e.g. output by hmmpress)
# dbtaxid		: map all database accessions to this taxid
# dbacc2taxid	: tsv-file mapping hmmer database accessions to taxids, MUST include fields 'acc','taxid', overrides dbtaxid
# seqs			: fasta sequences
# seqinfo		: tsv-file with sequence info, MUST include fields 'seqid','length'
#
# [append]		: append to existing $annot file, [false]
# [dbhits]		: tsv-file for raw annotations, [$seqs.sans.dbhits.tsv]
# [log]			: log-file, [$seqs.sans.log]
# [min_bits]		: min bitscore [0]
# [max_eval]		: max e-value [0.01]
# [numth]		: threads, [8]
# [min_orf_length]	: minimum orf length, [0]
# [orf_finder]		: software for orf detection, [mga]
#
sub annotate_hmmscan{
	print STDERR "\n# ANNOTATE WITH HMMSCAN\n\n";
	
 	my (%args)			= @_;
 	my $subid			= "annotate_hmmscan()";
 	
 	# Check input:
 	if( !defined($args{seqs})){
 		die "ERROR: $subid: missing arguments: seqs";}
 	if( !defined($args{seqinfo})){
 		die "ERROR: $subid: missing arguments: seqinfo";}
 	if( !defined($args{db})){
 		die "ERROR: $subid: missing arguments: db";}
 	if( !defined($args{annot})){
 		die "ERROR: $subid: missing arguments: annot";}
 	if( !defined($args{dbacc2taxid}) && !(defined($args{dbtaxid})) ){
 		die "ERROR: $subid: missing arguments: dbacc2taxid OR dbtaxid";
 	}

	# in:
	my $seqs		 		= $args{seqs};
	my $seqsinfo			= $args{seqinfo};
	my $db				= $args{db};
	my $dbacc2taxid		= $args{dbacc2taxid} || undef;
	my $dbtaxid			= $args{dbtaxid}	|| undef;
	
	# out:
	my $orfs_nt			= "$seqs.orfs.nt.fa";
	my $orfs_aa			= "$seqs.orfs.aa.fa";
	my $annot			= $args{annot};
	my $dbhits			= $args{dbhits} || "$seqs.hmmscan.dbhits.tsv";
	my $log 				= $args{log} || "$seqs.hmmscan.log";
	
	# tmp/pars:
	my $append			= $args{append} || 0;
	my $min_bits			= $args{min_bits} || 10;
	my $max_eval			= $args{max_eval} || 0.01;
	my $numth			= $args{numth} || 8;
	my $orf_finder		= $args{orf_finder} || "mga";
	my $min_orf_length	= $args{min_orf_length} || 0;
	
	
	my $search_name		= "hmmscan";
	my $searchdb_name	= basename($db);
	my $dbtype			= "HMM";
	my @fields_lazy 		= @ANNOT_TABLE_FIELDS;
	my $fields_lazy_tab	= join("\t",@fields_lazy);
	my $fields_lazy_comma= join(',',@fields_lazy);
	
	# CHECK INPUT
	if( !(-e $seqs)){ 
		die "ERROR: $subid: file does not exist: $seqs";}
	if( $dbacc2taxid && !(-e $dbacc2taxid)){
		die "ERROR: $subid: file does not exist: $dbacc2taxid";}

	# start new annotation
	if(!$append){
		system("echo \"$fields_lazy_tab\" 1> $annot");
		system("rm -f $log");
		system("touch $log");
	}
	# append requested but $annot is empty
	if($append && nlines($annot)<1){
		system("echo \"$fields_lazy_tab\" 1> $annot");
	}
	
	if( nlines($seqs)<1){
		print STDERR "WARNING: $subid: empty query sequences: $seqs\n";
		return();
	}
	
	system("rm -f $orfs_aa $orfs_nt");
	detect_orfs(seqs=>$seqs, orfs_nt=>$orfs_nt, orfs_aa=>$orfs_aa, orf_finder=>$orf_finder, min_orf_length=>$min_orf_length, numth=>$numth);
	if( nlines($orfs_aa)<2 ){
		print STDERR "WARNING: $subid: no orfs found in $seqs\n";
		return();
	}
	
	# run hmmscan
	system_call("hmmscan -E $max_eval --cpu $numth --tblout $dbhits $db $orfs_aa 1> $dbhits.tmp 2>> $log");
	parse_tblout2tsv( tblout=>$dbhits, tsv => "$dbhits.tmp");
	system("mv $dbhits.tmp $dbhits");
		# filter on score>=$min_bits
	system_call("cat $dbhits | csvtk filter -tf 'score>=$min_bits' 1> $dbhits.tmp 2>> $log");
	system_call("mv $dbhits.tmp $dbhits");
		# filter tophits: filter top for each contig-ORF
	filter_tophits(dbhits=>$dbhits,dbhits_flt=>"$dbhits.tmp", qcol=>'qseqid', bitscol=>'score', retain_ties=>1);
	system("mv $dbhits.tmp $dbhits");
	

	if( nlines($dbhits) < 2){
		print STDERR "\n\tHMMSCAN: NO HITS FOR $seqs\n\n";
		return();
	}		
				
	# Match hmmer fields to Lazypipe format
	# - Lazypipe:
	#	"search,db,dbtype,qseqid,orf,qseqlen,sseqid,bitscore,alen,pident,qlen,qcov,slen,scov,staxid,sname,bphage"
	# <Lazypipe>		<Hmmscan>	
	# 	search		"hmmscan"
	#	db			basename($db)		
	#	dbtype		"HMM"
	# 	qseqid		qseqid (contig part)
	#	orf			qseqid (orf part)
	#	qseqlen		<from $seqsinfo map>
	#	sseqid		sacc
	#	bitscore		score
	#	alen			- 
	#	pident		- 
	#	qlen			-
	#	qcov			-
	#	slen			- 
	#	scov			-
	#	staxid		<from $dbacc2taxid OR $dbtaxid>
	#	sname		desc
	#	bphage		"NA"
			
	system_call("cat $dbhits | csvtk cut -tlf qseqid,sacc,score,desc | ".
				"csvtk mutate2 -te '\"$search_name\"' -n 'search' | ".
				"csvtk mutate2 -te '\"$searchdb_name\"' -n 'db' | ".
				"csvtk mutate2 -te '\"$dbtype\"' -n 'dbtype' | ".
				"csvtk rename -tf qseqid -n qseqid_tmp | ".
				"csvtk mutate -tf qseqid_tmp -p '^([^_]+)' -n qseqid | ".
				"csvtk mutate -tf qseqid_tmp -p 'ORF=([^_]+)\$' -n orf | ".
				"csvtk rename -tf sacc -n sseqid | ".
				"csvtk rename -tf score -n bitscore | ".
				
				"csvtk mutate2 -te '\"NA\"' -n 'alen' | ".
				"csvtk mutate2 -te '\"NA\"' -n 'pident' | ".
				"csvtk mutate2 -te '\"NA\"' -n 'qlen' | ".
				"csvtk mutate2 -te '\"NA\"' -n 'qcov' | ".
				"csvtk mutate2 -te '\"NA\"' -n 'slen' | ".
				"csvtk mutate2 -te '\"NA\"' -n 'scov' | ".
				"csvtk rename  -tf desc -n sname | ".
				
				"csvtk mutate2 -te '\"NA\"' -n bphage | ".
				"csvtk mutate2 -te '\"NA\"' -n division 1> $dbhits.tmp 2>> $log");
	system_call("mv $dbhits.tmp $dbhits");
	
	# ADD FIELDS: qseqlen (contig-length) + staxid
	system_call("csvtk join -tj $numth -f 'qseqid;seqid' -L --na 'NA' $dbhits $seqsinfo | ".
				"csvtk rename -tf length -n qseqlen 1> $dbhits.tmp 2>> $log" );
	system_call("mv $dbhits.tmp $dbhits");			
				
	if($dbacc2taxid){	# use $dbacc2taxid tsv-map
		system_call("csvtk join -tj $numth -f 'sseqid;acc' -L --na 'NA' $dbhits $dbacc2taxid | ".
					"csvtk rename -tf taxid -n staxid 1> $dbhits.tmp 2>> $log" );
		system_call("mv $dbhits.tmp $dbhits");
	}
	else{	# use $dbtaxid
		system_call("csvtk mutate2 -te '$dbtaxid' -n 'staxid' $dbhits 1> $dbhits.tmp 2>> $log");
		system_call("mv $dbhits.tmp $dbhits");
	}
	
	# PRINT RESULTS TO $ANNOT
	system_call("cat $dbhits | ".
				"csvtk cut -tlf $fields_lazy_comma | ".
				"csvtk del-header -t 1>> $annot 2>> $log");

}




# USAGE:
# detect_orfs(seqs=>$fasta, orfs_nt=>$orfs_nt, orfs_aa=>$orfs_aa, orf_finder=>'mga', min_orf_length=>$length);
# 
# seqs			: input seqs in fasta
# orfs_nt		: output nt orfs in fasta
# orfs_aa		: output aa orfs in fasta
# [orf_finder]	: orf finder, mga/prod/orfipy, default = mga
# min_orf_length: minimum orf length
# par_prodigal  : params for prodigal ["-p meta "]
# par_orfipy		: params for orfipy ["--table 1 --partial-3"]
#
sub detect_orfs{
	print STDERR "\n# DETECT ORFS\n\n";
	my $subid		= "detect_orfs()";

	# in:
	my (%args)		= @_;
	my $seqs 		= $args{seqs};
	# out:
	my $orfs_aa		= $args{orfs_aa};
	my $orfs_nt		= $args{orfs_nt};
	my $orfs_info	= "$orfs_nt.tsv";
	# tmp:
	my $orfs_raw		= "$orfs_nt.raw.tmp";
	my $orfs_gtf		= "$orfs_nt.gtf.tmp";
	my $orfs_gff    = "$orfs_nt.gff.tmp";
	my $seqs_fai 	= "$seqs.seqkit.fai";
	my $orfs_orfipy_dir	= "$orfs_nt.orfipy.tmp";
	# params	:
	my $orf_finder		= $args{orf_finder} || "mga";
	my $min_orf_length	= defined($args{min_orf_length}) ? $args{min_orf_length} : 10;
	my $par_prodigal		= $args{par_prodigal} || "-p meta ";
	my $par_orfipy		= $args{par_orfipy} || "--table 1 --partial-3";
	my $numth			= $args{numth} || 8;

	system("rm -f $seqs_fai");	# to ensure fai corresponds to the latest seqs
	
	if(nlines($seqs)<2){
		# empty fasta input
		print STDERR "WARNING: $subid: empty input file: exporting empy ORF fasta\n";
		system_call("touch $orfs_aa $orfs_nt $orfs_info");
		return;
	}
	
	if($orf_finder eq 'mga'){
		system_call("mga $seqs -m 1> $orfs_raw" );
		mga2gtf(mga=>$orfs_raw,gtf=>$orfs_gtf);
		system_call("seqkit subseq --gtf $orfs_gtf $seqs | ".
					"seqkit seq -w0 --min-len $min_orf_length | ".
					"seqkit replace -p '_([\\w\\:+\\-]+)\\s*\$' -r '_ORF=\$1' | ".
					"seqkit sort --natural-order --quiet -w0 1> $orfs_nt" );
		system_call("seqkit translate -w0 -f 1 $orfs_nt 1> $orfs_aa" );
		system_call("seqkit seq -n $orfs_nt | sed 's/\\([^_]\\+\\)_\\([^_]\\+\\)/\\1,\\1_\\2/' 1> $orfs_info");
	}
	elsif($orf_finder =~ m/^prod/gi){
		system_call("prodigal $par_prodigal -f sco -q -i $seqs -o $orfs_raw -d $orfs_nt.tmp" );
		system_call("grep -v '#' $orfs_raw | ".
					"sed -r 's/^>([0-9]+)_([0-9]+)_([0-9]+)_([+-])/ORF=\\2-\\3:\\4/' | nl -w1  1> $orfs_raw.map", 1); 
					# outputs "orf_num \t ORF=from-to:strand" map
		system_call("cat $orfs_nt.tmp | ".
					"seqkit seq -i  | ".
					"seqkit replace -p '_(\\w+)\$' -r '_{nr}' | ".
					"seqkit replace -p '_(\\w+)\$' -r '_{kv}' -k $orfs_raw.map -m 'NA'  | ".
					"seqkit seq -w 90 --min-len $min_orf_length | ".
					"seqkit sort --natural-order --quiet -w 90 1> $orfs_nt" );
					
		system_call("seqkit translate -w 90 -f 1 $orfs_nt 1> $orfs_aa" );
		system_call("seqkit seq -n $orfs_nt | sed 's/\\([^_]\\+\\)_\\([^_]\\+\\)/\\1,\\1_\\2/' 1> $orfs_info " );
		system("rm -f $orfs_nt.tmp $orfs_raw.map");
	}
	elsif(lc($orf_finder) eq 'orfipy'){
		system_call("orfipy $par_orfipy --procs $numth --min $min_orf_length --outdir $orfs_orfipy_dir --bed orfs.bed $seqs");
		system("rm -f $seqs.fai");
		system_call("seqkit subseq --bed $orfs_orfipy_dir/orfs.bed $seqs | ".
					"seqkit seq -i | ".
					"seqkit replace -p '_([\\w\\:+\\-]+)\\s*' -r '_ORF=\$1 ' | ".
					"seqkit sort --natural-order --quiet -w0 1> $orfs_nt" );
		system_call("seqkit subseq --bed $orfs_orfipy_dir/orfs.bed $seqs | ".
					"seqkit seq -i | ".
					"seqkit replace -p '_([\\w\\:+\\-]+)\\s*' -r '_ORF=\$1 ' | ".
					"seqkit sort --natural-order --quiet | ".
					"seqkit translate -w0 -f 1 --transl-table 1 1> $orfs_aa ");
	}
	else{
		die "ERROR: $subid: unknown orf_finder $orf_finder";
	}

	system("rm -f $seqs_fai");
	system_call("rm -fr $orfs_nt.*.tmp" );
}

# Filter host reads using bwa mem + samtools + csvtk
#
# USAGE:	
# 	filter_host_reads(r1 =>$r1, r2=>$r2, hostdb => $hostdb, res => $resdir, log=>$log,
#						numth=>16, bitscore => 50, mapq => 0, tmpdir => '.', gzip=>'pigz')
# r1	            : forward reads
# r2	            : reverce reads (can be undef for SE data)
# r1_pass       : forward reads with hostgen reads removed
# r2_pass       : reverce reads with hostgen reads removed (can be undef for SE data)
# r1_flt        : forward reads removed as host reads
# r2_flt        : reverce reads removed as host reads (can be undef for SE data)
# hostdb		    : host database in fasta[.gz] or bwa-index format
# res			: result dir for outputing results
# log			: log
#
# PARAMS:
# numth			: number of threads [16]
# bitscore		: alignment bitscore, (AS score printed by bwa mem) [0]
# mapq			: mapping quality [0]
# tmpdir			: temp dir [$res]
# gzip			: archiving utility [gzip]
#
sub filter_host_reads{
	my $SIGNATURE	= "filter_host_reads()";
	
	# in:
	my (%args)		= @_;
	# args{r1}
	# args{r2}
	# args{r1_pass}
	# args{r2_pass}
	# args{r1_flt}
	# args{r2_flt}
	# args{hostdb}
	# args{res}
	# args{log}
	if( !defined($args{r1})){
 		die "ERROR: $SIGNATURE: missing input: r1";
	}
 	if( !defined($args{r1_pass}) ){
 		die "ERROR: $SIGNATURE: missing input: r1_pass";
 	}
 	if( defined($args{r2}) && !defined($args{r2_pass}) ){
 		die "ERROR: $SIGNATURE: missing input: r2_pass";
 	}
 	if( !defined($args{r1_flt}) ){
 		die "ERROR: $SIGNATURE: missing input: r1_flt";
 	}
 	if( defined($args{r2}) && !defined($args{r2_flt}) ){
 		die "ERROR: $SIGNATURE: missing input: r2_flt";
 	}
	if( !defined($args{hostdb}) || !$args{hostdb}){
		print STDERR "\t$SIGNATURE: no hostdb specified: no filtering\n";
		return;
	}
	if( !defined($args{res}) ){
		die "ERROR: $SIGNATURE: missing input: res";
	}
	if( !defined($args{log}) ){
		die "ERROR: $SIGNATURE: missing input: log";
	}
		# create names for unpacked read files
	my $reads_gz				= 0;
	if($args{r1} =~ m/\.gz/gi){
		$reads_gz			= 1;
	}
	$args{r1_pass}			=~ s/\.gz$//i;
	$args{r1_flt}			=~ s/\.gz$//i;
	if( defined($args{r2_pass}) ){
		$args{r2_pass}		=~ s/\.gz$//i;}
	if( defined($args{r2_flt}) ){
		$args{r2_flt}		=~ s/\.gz$//i;	
	}
		
	# params
	if( !defined($args{numth}) ){
		$args{numth}		= 16;
	}
	if( !defined($args{bitscore}) ){
		$args{bitscore}	= 0;
	}
	if( !defined($args{mapq}) ){
		$args{mapq}		= 0;
	}
	if( !defined($args{tmpdir}) ){
		$args{tmpdir}	= $args{res};
	}
	if( !defined($args{gzip}) ){
		$args{gzip}		= 'gzip';
	}
	
	# tmp
	my $sam		= "$args{res}/hostgen.sam";
	my $readids	= "$args{res}/hostgen.readids";
	
	
	# START WORKING
	if($VERBAL){
		print STDERR "\n\t$SIGNATURE\n";
	}

	system("rm -f $args{log}");
	system("touch $args{log}");
	
	
	if( !((-e "$args{hostdb}.amb") && (-e "$args{hostdb}.sa"))){
		system_call("bwa index -a bwtsw $args{hostdb}");
	}

	if( !defined($args{r2}) ){ # SE-reads
		system_call("bwa mem -t $args{numth} -T $args{bitscore} $args{hostdb} $args{r1} 1> $sam 2>> $args{log}");
		system_call("sambamba view -t $args{numth} -S -F \"not(unmapped) and mapping_quality>=$args{mapq} and [AS]>=$args{bitscore}\" $sam 1> $sam.tmp 2>> $args{log}");
		system_call("mv $sam.tmp $sam");
		
		if(nlines($sam)>0){
			system_call("cut -f1 $sam | csvtk sort -Hj $args{numth} -k1 | uniq 1> $readids");	
		}
		else{
			system("touch $readids");
		}
		
		if(nlines($readids)>0){
			system_call("seqkit grep -j $args{numth} -vf $readids  --id-regexp \"^([^/\\\\s]+)\\\\s?\" $args{r1} | seqkit seq -u 1> $args{r1_pass} 2>> $args{log}");
			system_call("seqkit grep -j $args{numth} -f $readids  --id-regexp \"^([^/\\\\s]+)\\\\s?\" $args{r1} | seqkit seq -u 1> $args{r1_flt} 2>> $args{log}");		
			system_call("$args{gzip} -f $args{r1_pass} 2>> $args{log}");
			system_call("$args{gzip} -f $args{r1_flt} 2>> $args{log}");
		}
		else{
			print STDERR "\t$SIGNATURE: no host reads identified";
			if($reads_gz){
				system_call("cp $args{r1} $args{r1_pass}.gz");
			}
			else{
				system_call("cp $args{r1} $args{r1_pass}");
				system_call("$args{gzip} -f $args{r1_pass} 2>> $args{log}");
			}
		}
	}
	else{ # PE-reads
		system_call("bwa mem -t $args{numth} -T $args{bitscore} $args{hostdb} $args{r1} $args{r2} 1> $sam 2>> $args{log}");
		system_call("sambamba view -t $args{numth} -S -F \"not(unmapped) and mapping_quality>=$args{mapq} and [AS]>=$args{bitscore}\" $sam 1> $sam.tmp 2>> $args{log}");
		system_call("mv $sam.tmp $sam");
		
		if(nlines($sam)>0){
			system_call("cut -f1 $sam | csvtk sort -Hj $args{numth} -k1 | uniq 1> $readids");	
		}
		else{
			system("touch $readids");
		}
		
		if(nlines($readids)>0){
			system_call("seqkit grep -j $args{numth} -vf $readids  --id-regexp \"^([^/\\\\s]+)\\\\s?\" $args{r1} | seqkit seq -u 1> $args{r1_pass} 2>> $args{log}");
			system_call("seqkit grep -j $args{numth} -vf $readids  --id-regexp \"^([^/\\\\s]+)\\\\s?\" $args{r2} | seqkit seq -u 1> $args{r2_pass} 2>> $args{log}");
			system_call("seqkit grep -j $args{numth} -f $readids  --id-regexp \"^([^/\\\\s]+)\\\\s?\" $args{r1} | seqkit seq -u 1> $args{r1_flt} 2>> $args{log}");
			system_call("seqkit grep -j $args{numth} -f $readids  --id-regexp \"^([^/\\\\s]+)\\\\s?\" $args{r2} | seqkit seq -u 1> $args{r2_flt} 2>> $args{log}");
			system_call("$args{gzip} -f $args{r1_pass} 2>> $args{log}");
			system_call("$args{gzip} -f $args{r2_pass} 2>> $args{log}");			
			system_call("$args{gzip} -f $args{r1_flt} 2>> $args{log}");
			system_call("$args{gzip} -f $args{r2_flt} 2>> $args{log}");
		}
		else{
			print STDERR "\t$SIGNATURE: no host reads identified\n";
			if($reads_gz){
				system_call("cp $args{r1} $args{r1_pass}.gz");
				system_call("cp $args{r2} $args{r2_pass}.gz");
			}
			else{
				system_call("cp $args{r1} $args{r1_pass}");
				system_call("cp $args{r2} $args{r2_pass}");
				system_call("$args{gzip} -f $args{r1_pass} 2>> $args{log}");
				system_call("$args{gzip} -f $args{r2_pass} 2>> $args{log}");
			}
		}
	}
	system("rm -f $sam $readids");
}

# Filter host contigs using bwa mem + samtools + csvtk
#
# USAGE:	
# 	filter_host_contigs(contigs => $contigs, contigs_pass=> $contigs_pass, contigs_flt => $contigs_flt, hostdb => $hostdb, res => $resdir, log=>$log,
#						numth=>16, bitscore => 400, tmpdir => '.')
# contigs		: contig fasta
# contigs_pass	: contigs that passed filtering
# contigs_flt   : contigs filtered
# hostdb		    : host database in fasta[.gz] or bwa-index format
# res			: result dir for outputing results
# log			: log
#
# PARAMS:
# numth			: number of threads [16]
# bitscore		: alignment bitscore, (AS score printed by bwa mem) [0]
# tmpdir			: temp dir [$res]
#
# OUTPUT:
# $contigs			: input contig fasta will be replaced with non-host contigs
# $contigs.host.fa	: filtered host contigs
#
sub filter_host_contigs{
	my $SIGNATURE	= "filter_host_contigs()";
	
	# in:
	my (%args)		= @_;
	# args{contigs}
	# args{contigs_flt}
	# args{hostgen}
	# args{res}
	# args{log}
	if( !defined($args{contigs})){
 		die "ERROR: $SIGNATURE: missing input: contigs";	
	}
	if( !defined($args{contigs_pass})){
 		die "ERROR: $SIGNATURE: missing input: contigs_pass";	
	}
	if( !defined($args{contigs_flt})){
 		die "ERROR: $SIGNATURE: missing input: contigs_flt";	
	}
	if( !defined($args{hostdb}) || !$args{hostdb}){
		print STDERR "\t$SIGNATURE: no hostdb specified: no filtering\n";
		return;
	}
	if( !defined($args{res}) ){
		die "ERROR: $SIGNATURE: missing input: res"; }
	if( !defined($args{log}) ){
		die "ERROR: $SIGNATURE: missing input: log"; }
	# params
	if( !defined($args{numth}) ){
		$args{numth}		= 16; }
	if( !defined($args{bitscore}) ){
		$args{bitscore}	= 0; }
	if( !defined($args{tmpdir}) ){
		$args{tmpdir}	= $args{res}; }	
	
	# tmp:
	my $sam				= "$args{contigs}.sam.tmp";
	my $sam2				= "$args{contigs}.sam2.tmp";
	my $contigs_hostids	= "$args{contigs}.hostids.tmp";
	
	# START WORKING
	if($VERBAL){
		print STDERR "\n\t$SIGNATURE\n";
	}
	system("rm -f $args{log}");
	system("touch $args{log}");
		
	system_call("bwa mem -t $args{numth} -T $args{bitscore} $args{hostdb} $args{contigs} 1> $sam 2>> $args{log}");
	system_call("sambamba view -t $args{numth} -S -F \"not(unmapped) and [AS]>=$args{bitscore}\" $sam 1> $sam2 2>> $args{log}");
	system_call("mv $sam2 $sam");
	
	if(nlines($sam)>0){
		system_call("cut -f1 $sam | sort | uniq 1> $contigs_hostids 2>> $args{log}");
	}
	else{
		system("touch $contigs_hostids");
	}	
	
	if( nlines($contigs_hostids)>0 ){
		print STDERR "\t$SIGNATURE: found host contigs\n";
		system_call("seqkit grep -w0 -j $args{numth} -f $contigs_hostids $args{contigs}  1> $args{contigs_flt} 2>> $args{log}");
		system_call("seqkit grep -w0 -j $args{numth} -vf $contigs_hostids $args{contigs} 1> $args{contigs_pass} 2>> $args{log}");
	}
	else{
		print STDERR "\t$SIGNATURE: no host contigs identified\n";
		system("cp $args{contigs} $args{contigs_pass}");
	}	
	
	# clean
	system_call("rm -f $sam $sam2 $contigs_hostids");
}


# Filter sequences based on a supplied Annotation-table and a Condition-string
#
# USAGE:
# filter_seqs(seqs=>$fasta, seqs_flt=>$fasta_flt, annot=>$annot, seqidh=>$seqidh, condition=>$condition_string, log=>$log) 
#
# PARAMETERS: 
# 	INPUT:
# 	seqs			: Fasta-file with input sequences
#	annot	 	: TSV-file with sequence annotations. MUST INCLUDE headers used in the condition-string.
#	seqidh		: seqid header in the the annot-file
#	condition   : csvtk filter2 condition, eg \''\$division==\"Viruses\"\'. Use \-escape for extrapolated chars
# 	OUTPUT:
# 	seqs_flt		: Fasta-file with filtered sequences
# 	log			: Log
#
sub filter_seqs{
	my $SIGNATURE	= "filter_seqs()";
	
	# in:
	my (%args)			= @_;
	my $seqs 			= $args{seqs};
	my $annot			= $args{annot};
	my $seqidh			= $args{seqidh};
	my $condition		= $args{condition};
	# out:
	my $seqs_flt			= $args{seqs_flt};
	my $log 				= $args{log};
	# tmp:
	my $filter			= "$seqs_flt.flt.tmp";
	
	system("rm -f $log");
	system("touch $log");
	
	if($VERBAL){
		print STDERR "\n\t$SIGNATURE: ",basename($seqs_flt),"\n";
	}
		
	system_call("csvtk filter2 -tf $condition $annot | ".
				"csvtk cut -tlf $seqidh | tail -n+2 | uniq 1> $filter 2>> $log" );

	if(nlines($filter) > 0){
		system_call("seqkit grep -w0 -f $filter $seqs 1> $seqs_flt 2>> $log" );
	}
	else{
		print STDERR "\tWARNING: $SIGNATURE: no sequences matching $condition: exporting empty fasta\n";
		system("rm -f $seqs_flt");
		system_call("touch $seqs_flt" );
	}
	system("rm -f $filter" );
}

# Filter unannotated sequences
#
# USAGE:
# filter_unseqs(seqs=>$fasta, seqs_un=>$unfasta, annot=>$annot, seqidh=>$seqidh, log=>$log) 
#
# PARAMETERS: 
# 	INPUT:
# 	seqs			: Fasta-file with input sequences
#	annot	 	: TSV-file with sequence annotations.
#	seqidh		: seqid header in the the annot-file
# 	OUTPUT:
# 	seqs_un		: Fasta-file with sequences that have no annotation in the $annot-file
# 	log			: Log
#
sub filter_unseqs{
	my $SIGNATURE	= "filter_unseqs()";
	if($VERBAL){
		print STDERR "\n\t$SIGNATURE\n";
	}
	
	# in:
	my (%args)			= @_;
	my $seqs 			= $args{seqs};
	my $annot			= $args{annot};
	my $seqidh			= $args{seqidh};
	# out:
	my $seqs_un			= $args{seqs_un};
	my $log 				= $args{log};
	# tmp:
	my $filter			= "$seqs.flt.tmp";	
	
	system("rm -f $log");
	system("touch $log");
	
	system_call("csvtk cut -tlf $seqidh $annot | tail -n+2 | uniq 1> $filter 2>> $log");
	if( nlines($filter)>0 ){
		system_call("seqkit grep -w0 -vf $filter $seqs 1> $seqs_un 2>> $log");
	}
	else{
		# annotation is empty
		print STDERR "\tWARNING: $SIGNATURE: empty annotation file\n";
		system_call("cp $seqs $seqs_un");
	}

	system("rm -f $filter");
}


# Prints only the top scoring subject for each query in a given SAM file
# Entries are ranked based on BWA score in column that is at position >=12 and has format AS:i:num
# Enties with no AS:i:score are ignored
#
# USAGE:
# filter_tophits_SAM(in=>$samfile, out=>$samfile_flt)
#
# in		: input SAM file
# out	: filtered SAM file
#
sub filter_tophits_SAM{
	my (%args)		= @_;
	
	# Check input
	if( !defined($args{in})){
		die "ERROR: filter_tophits_SAM: missing argument: in";}
	if( !defined($args{out})){
		die "ERROR: filter_tophits_SAM: missing argument: out";}
	
		
	my $sam_in 	= $args{in};
	my $sam_out	= $args{out};
	
	open(IN,"<$sam_in") or die "Can\'t open $sam_in: $!\n";
	open(OUT,">$sam_out") or die "Can\'t open $sam_out: $!\n";
	my $l;
	my $li 			= 0;
	my $score		= -1;
	my $best_score	= -1;
	my $best_line	= !1;
	my $qname		= "";
	my $first_ali	= 1;
	my @sp;

	# The column to which bwa mem prints AS:i:score varies 
	# and can be basically any column starting from SAM OPT-field, i.e. any column starting from 12.
	# Thus to make this bulletproof we search columns 12>last for each line
	while($l=<IN>){
		$li++;
		if( $l =~ m/^@/){	# SAM HEADER SECTION
			print OUT "$l";
			next;
		}
		# SAM ALIGNMENT SECTION
		chomp($l);
		@sp 				= split(/\t/,$l,-1);
		$score			= -1;
		for(my $col=12-1; $col<scalar(@sp); $col++){
			if( $sp[$col] =~ m/AS\:i\:([+-]?[\d]+)/ ){
				$score  = $1; 
				last;
			}
		}
		if($score< 0){
			next;
		}
		#print STDERR "# found AS:i:score on line $li\n";
	
		if($first_ali){
			$qname		= $sp[0];
			$best_line 	= $l;
			$best_score	= $score;
			$first_ali	= !1;
		}
		elsif($sp[0] ne $qname){
			print OUT "$best_line\n";
			$qname 		= $sp[0];
			$best_line 	= $l;
			$best_score	= $score;
		}
		else{
			if($score > $best_score){
				$best_score 	= $score;
				$best_line 	= $l;
			}
		}
	}
	if($li > 0){
		print OUT "$best_line\n";
	}
	close(IN);
	close(OUT);
}

###
# USAGE:
# filter_tophits(dbhits=>$dbhits, dbhits_flt=>$dbhits_flt, qcol=>$qcol, bitscol=>$bitscol, retain_ties=>0)
#
# dbhits			: input tabular search result table, MUST BE SORTED by qseq, MUST INCLUDE HEADERS
#				  This can be BLAST/BLASTP/MINIMAP/SANS tabular search result table.
# dhhits_flt		: filtered dbhits
# qcol			: qseq column name
# bitscol		: bitscore column name
# [retain_ties]	: retain all top-scoring ties [false]
# 
# NOTE: #@-Comment lines are retained
#
sub filter_tophits{
	my (%args)		= @_;
	# Check input
	if( !defined($args{dbhits})){
		die "ERROR: filter_tophits: missing argument: dbhits";}
	if( !defined($args{dbhits_flt})){
		die "ERROR: filter_tophits: missing argument: dbhits_flt";}
	if( !defined($args{qcol})){
		die "ERROR: filter_tophits_SAM: missing argument: qcol";}
	if( !defined($args{bitscol})){
		die "ERROR: filter_tophits_SAM: missing argument: bitscol";}
	
	# in:
	my $dbhits		= $args{dbhits};
	my $dbhits_flt	= $args{dbhits_flt};
	my $qcol			= $args{qcol};
	my $bitscol		= $args{bitscol};
	my $retain_ties	= defined($args{retain_ties}) ? $args{retain_ties} : 0;

	# READ HEADERS 
	my $qcoli		= colind($dbhits,$qcol) -1;
	my $bitscoli		= colind($dbhits,$bitscol) -1;
	if($qcoli < 0){
		die "ERROR: filter_tophits: missing header: $qcoli";
	}
	if($bitscoli < 0){
		die "ERROR: filter_tophits: missing header: $bitscoli";
	}

	open(IN,"<$dbhits") or die "Can\'t open $dbhits: $!\n";
	open(OUT,">$dbhits_flt") or die "Can\'t open $dbhits_flt: $!\n";
	my $l=<IN>;
	print OUT "$l";

	# READ DATA
	my $score			= -1;
	my $best_score		= -1;
	my $best_line		= !1;
	my $qname			= "";
	my $first_ali		= 1;
	while($l=<IN>){
		if( $l =~ m/^[#@]/){	# comment lines
			print OUT "$l";
			next;
		}
		chomp($l);
		my @sp 			= split(/\t/,$l,-1);
		$score			= $sp[$bitscoli];
		if($first_ali){
			$qname		= $sp[$qcoli];
			$best_line 	= $l;
			$best_score	= $score;
			$first_ali	= !1;
		}
		elsif($sp[$qcoli] ne $qname){
			print OUT "$best_line\n";
			$qname  		= $sp[$qcoli];
			$best_line  = $l;
			$best_score = $score;
		}
		else{	
			if($score > $best_score){
				$best_score = $score;
				$best_line 	= $l;
			}
			elsif(($score == $best_score) && $retain_ties){
				print OUT "$best_line\n";
				$best_line 	= $l;
			}
		}
	}
	print OUT "$best_line\n";
	close(OUT);
	close(IN);
}


# CONVERT MGA ORF PREDICTIONS TO GTF2.2
# 
# USAGE:
# mga2gtf(mga=>mga_file, gtf=>gtf_file)
# 
# mga_file			: input mga prediction file
# gtf_file			: output gtf file
# [set_frame_zero]	: use 0-based frame [true]
# 
# mga_file MUST be in MGA format (ref: http://metagene.nig.ac.jp/metagene/download_mga.html):
# # [sequence name]
# # gc = [gc%], rbs = [rbs%]
# # self: [(b)acteria/(a)rchaea/(p)hage/unused(-)]
# [0:gene ID] [1:start pos.] [2:end pos.] [3:strand] [4:frame] [5:complete/partial] [6:gene score] [7:used model] [8:rbs start] [9:rbs end] [10:rbs score]
#
sub mga2gtf{
	my (%args)		= @_;
	# Check input
	if( !defined($args{mga})){
		die "ERROR: mga2gtf: missing argument: mga";}
	if( !defined($args{gtf})){
		die "ERROR: mga2gtf: missing argument: gtf";}
	
	# in:
	my $mga		= $args{mga};
	# out:
	my $gtf		= $args{gtf};
	# params:
	my $set_frame_zero	= defined($args{set_frame_zero}) ? $args{set_frame_zero} : 1;	

	my $seqid;
	my @gene_pred_list;
	my %gene_pred_hash;

	open(IN,"<$mga") or die "Can\'t open $mga: $!\n";
	while(my $l=<IN>){
		chomp($l);
        if( $l =~ m/^[\!#]/){
        		if(scalar(@gene_pred_list) > 0){	# genes predictions read for prev seq
	    			my @copy 				= @gene_pred_list;
				$gene_pred_hash{$seqid} = \@copy;
	   	     	@gene_pred_list			= ();
	    		} 
			$l 		=~ s/^[#\s]+|\s+$//g; # remove leading/trailing spaces
			my @sp	= split(/\s+/,$l);
			$seqid	= $sp[0];
			# reading two more comment lines
			$l=<IN>;
			$l=<IN>;
			next;
		}
        
		if( $l =~ m/^gene/i ){
			push(@gene_pred_list,$l);
		}
	}
	if(scalar(@gene_pred_list)>0){
		my @copy 				= @gene_pred_list;
		$gene_pred_hash{$seqid} = \@copy;
	}
	close(IN);

	# WRITE GTF2.2: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
	open(OUT,">$gtf") or die "Can\'t open $gtf: $!\n";
	
	foreach my $seq( sort keys %gene_pred_hash){
		
		my @genepred_list =  @{$gene_pred_hash{$seq}};
	
		foreach my $genepred( @genepred_list){
	
			my ($geneid,$start,$end,$strand,$frame,$complete,$score,$model,$rbs_start,$rbs_end,$rbs_score) = split(/\t/,$genepred,-1);
			if($set_frame_zero && ($strand eq '+')){
				$start  = $start + $frame;
				$frame	= 0;
			}
			if($set_frame_zero && ($strand eq '-')){
				$end	= $end - $frame;
				$frame	= 0;
			}
		
			print OUT "$seq\tbwa\tORF\t$start\t$end\t$score\t$strand\t$frame\tgene_id \"\"; transcript_id \"\"; complete $complete; model $model;\n";
		}
	}
	close(OUT);
}


# Returns query length from CIGAR-string
# CIGAR Code: M,I,D,N,S,H,P,=,X
# Codes consuming qseq: M/I/S/H/=/X
sub cigar2qlen{
	my $cigar	= shift;
	my $qlen		= 0;
	while($cigar =~ m/([0-9]+)([MISH=X]{1})/g ){
		$qlen += $1;
	}
	return $qlen;
}
# Returns query coverage as a fraction
sub cigar2qcov{
	my $cigar	= shift;
	my $qlen		= cigar2qlen($cigar);
	my $qalen	= 0;	# aligned part
	while($cigar =~ m/([0-9]+)([MI=X]{1})/g ){
		$qalen += $1;
	}
	return ($qalen/$qlen);
}
# Returns referece length from CIGAR-string
# CIGAR Code: M,I,D,N,S,H,P,=,X
# Codes consuming reference seq: M/D/N/=/X
sub cigar2rlen{
	my $cigar	= shift;
	my $rlen	= 0;
	while($cigar =~ m/([0-9]+)([MDN=X]{1})/g ){
		$rlen += $1;
	}
	return $rlen;
}
# Returns alignment length from CIGAR-string
# CIGAR Code: M,I,D,N,S,H,P,=,X
# Codes increasing alignment: MIDN=X (exclude clipping and padding)
sub cigar2alen{
	my $cigar	= shift;
	my $alen	= 0;
	while($cigar =~ m/([0-9]+)([MIDN=X]{1})/g ){
		$alen += $1;
	}
	return $alen;
}
# Returns alignment percent-identity from CIGAR-string
# CIGAR Code: M,I,D,N,S,H,P,=,X
# Assumes CIGAR is SAMv1 with '=' character denoting matches
sub cigar2pide{
	my $cigar	= shift;
	my $alen		= cigar2alen($cigar);
	my $Mlen		= 0;
	while($cigar =~ m/([0-9]+)([=]{1})/g ){
		$Mlen += $1;
	}

	return $Mlen/$alen;
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
# my $score_col = mcolind("dbhits.tsv","AS:i:[0-9]+")
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

# Parses hmmer --tblout to tsv with headers
# 
# USAGE:
# parse_tblout2tsv( tblout=> hmmer.tblout, tsv => hmmer.tsv)
# 
sub parse_tblout2tsv{
	my $subid			= "parse_tblout2tsv()";
 	my (%args)			= @_;
 	
 	# Check input:
 	$args{tblout} || die "ERROR: $subid: missing argument: tblout";
 	$args{tsv}	|| die "ERROR: $subid: missing argument: tsv";
 	(-e "$args{tblout}") || die "ERROR: $subid: file does not exist: $args{tblout}";
 		
	my @COLNAMES		= qw/sname sacc qseqid qacc eval score bias best.eval best.score best.bias dom.exp dom.reg dom.clu dom.ov dom.env dom.dom dom.rep dom.inc desc/;
	my $COLNAMES_STR	= join("\t",@COLNAMES);
	
	my $linen	= 0;
	open(IN,"<$args{tblout}") or die "Can\'t open $args{tblout}: $!\n";
	open(OUT,">$args{tsv}") or die "Can\'t open $args{tsv}: $!\n";

	print OUT "$COLNAMES_STR\n";

	while( my $l=<IN> ){
		$linen++;
		if($l =~ m/^#/){	# comment
			next;
		}
		chomp($l);
		my @sp	= split(/\s+/,$l, 19);
		print OUT join("\t",@sp),"\n";	
	}
	close(IN);
	close(OUT);	
}


# Expectation Maximization main loop
# USAGE:
# my $res = EM_loop(annot=>$annot_tsv,qseqid=>'qseqid',staxid=>'staxid',bitscore=>'bitscore', logLdiff=>0.01)
#
# PARAMS:
# annot		: TSV table with read annotations.
#			  MUST inlcude named columns $qseqid, $staxid and $bitscore
# qseqid		: column header for query-seqids
# staxid		: column header for subject-taxid
# bitscore	: column header for bitscores
# logLdiff	: minimal logL difference to continue EM
#
# OUTPUT:
# $res->{logL}	: EM log likelihood
# $res->{F}		: Hash with F-values by taxid
#
sub EM_loop{
	my $subid = "EM_loop()";
	
	# IN
	my (%args)			= @_;	
	my $annot 			= $args{annot} || die "ERROR: $subid: missing input: annot";
	my $q_prob			= $args{q_prob} || undef;
	# TMP
	my $annot_nospflt	= "$annot.nospflt.tmp";
	# OUT
	my %F				= ();
	my %Fiter			= ();	# F-values for each iteration, e.g. $Fiter{0}->{$taxid}
	
	# PARAMS
	my $par_logL_diff 	= $args{logLdiff} || 0.01;
	my $par_maxiter		= $args{maxiter} || 1000;
	my $qseqid_col		= $args{qseqid} || 'qseqid';
	my $staxid_col		= $args{staxid}	|| 'staxid';
	my $bitscore_col		= $args{bitscore} || 'bitscore';
	my $nospnorep		= $args{nospnorep} || 0;
	my $numth			= $args{numth} || 8;
	my $score2prob_temp	= $args{score2prob_temp} || 1.0;
	
	# START
	print STDERR "# STARTING $subid: logLdiff=$par_logL_diff, maxiter=$par_maxiter, temp=$score2prob_temp\n" if $VERBAL;
		# filter database hits with no species-level taxonomy
	if($nospnorep){
		printf STDERR "\t$subid: exclude DB hits with no species assigned\n" if ($VERBAL);
		system_call("csvtk grep -j $numth -tvf sname -rp '\\w+ sp\\.' $annot 1> $annot_nospflt");
		$annot	= $annot_nospflt;
	}
	
		# get log P(q|t)
	my %q_t_logprob = EM_q_t_logprob(annot=>$annot, qseqid=>$qseqid_col,staxid=>$staxid_col, bitscore=>$bitscore_col,temp=>$score2prob_temp);
										
		# init F(t) to 1/|T|
	foreach my $q(keys %q_t_logprob){
		foreach my $t(keys %{$q_t_logprob{$q}}){
			$F{$t}	= 1;
		}
	}
	my $numtaxa	= (scalar keys %F) + 0.0;
	%F			= map {$_ => 1.0 / $numtaxa } keys %F;
	$Fiter{0}	= { %F };
	my @taxids	= sort keys %F;
		
		# if P(q) is undefined, init to 1/|Q|
	my %q_prob		= ();
	if(defined($q_prob)){
		%q_prob		= %$q_prob;
	}
	else{
		my $numqseq	= (scalar keys %q_t_logprob) + 0.0;
		%q_prob		= map{$_ => 1.0 / $numqseq} keys %q_t_logprob;
	}
	my %q_logprob	= map{ $_ => log($q_prob{$_}) } keys %q_prob;

	# init logL
	my $logL			= EM_logL(q_t_logprob=>\%q_t_logprob, F=>\%F);
	my $iteration 	= 0;
	printf STDERR "\t$subid: iteration = %u, logL = %f\n",$iteration,$logL if($VERBAL);
	#printf STDERR "\tF = %s\n",join(",",@F{@taxids});
	
	for(my $i=0; $i<$par_maxiter; $i++){
		# update L(Q) = total likelihood of observing query sequences q in Q, given observed taxa
		$iteration++;
		my %F_post 		= EM_step(q_t_logprob=>\%q_t_logprob, q_logprob=>\%q_logprob, F=>\%F);
		my $logL_post	= EM_logL(q_t_logprob=>\%q_t_logprob, F=>\%F_post);
		my $logL_diff	= $logL_post - $logL;
		
		printf STDERR "\t$subid: iteration = %u, logL = %f, logLdiff = %f\n",$iteration,$logL_post,$logL_diff  if($VERBAL);
		
		$logL 				= $logL_post;
		%F					= %F_post;
		$Fiter{$iteration}	= {%F_post};
		
		if($logL_diff < $par_logL_diff && $i>1){
			# exit EM
			print STDERR "\t$subid: exiting EM\n" if($VERBAL);
			last;
		}
	}
	
	system("rm -f $annot_nospflt");

	return {logL=>$logL,F=>\%F, Fiter=>\%Fiter};
}

# Expectation Maximization step
# USAGE:
# my %F_post = EM_step(q_t_logprob=>\%q_t_logprob, q_logprob=>\%q_logprob, F=>\%F);
# 
sub EM_step{
	my $subid = "expectation_maximization_step()";
	# IN
	my (%args)			= @_;
	my $q_t_logprob		= $args{q_t_logprob} || die "ERROR: $subid: missing input: q_t_logprob";
	my %q_t_logprob		= %$q_t_logprob;
	my $q_logprob		= $args{q_logprob} || die "ERROR: $subid: missing input: q_logprob";
	my %q_logprob		= %$q_logprob;
	my $F				= $args{F} || die "ERROR: $subid: missing input: F";
	
	# OUT
	my %F				= %$F;	# updated F
	# CONST
	my $INF    			= 9**9**9;
	
	# START:
	
	# Using Baye's theorem:
	#
	# log(P(t|q)) 	= log{ P(q|t)*F(t)/SUM_t{P(q|t)*F(t)} } 
	#				= logP(q|t) + logF(t) - log SUM_t{P(q|t)*F(t)}
	# 
	# Denominator: factoring out by the largest term, P(q|t_max)*F(t_max), to avoid overflow
	# log SUM_t{P(q|t)*F(t)}
	#			= logP(tmax)+ logF(tmax)+ log SUM_t{ exp[logP(t)*F(t)]/exp[logP(tmax)*F(tmax)] }
	#			= logP(tmax)+ logF(tmax)+ log SUM_t{ exp[logP(t)-logP(tmax)] * F(t)/F(tmax)  }
	#
	my %t_q_logprob	= ();
	
	foreach my $q(keys %q_t_logprob){
		
		my @taxids		= sort keys %{$q_t_logprob{$q}};
		my @logP			= map{ $q_t_logprob{$q}->{$_} } @taxids;
		my @F			= map{ $F{$_} } @taxids;
		my @logF			= map{ ($F[$_]>0)? log($F[$_]) : -$INF } (0..$#F);
		# max_t{ logP * logF }
		my @logPlogF		= map{ $logP[$_]+$logF[$_] } (0 .. $#logP);
		my $maxi			= maxind(\@logPlogF);
		my $logP_max		= $logP[$maxi];
		my $logF_max		= $logF[$maxi];
		my $F_max		= $F[$maxi];
		
		my @sum_terms	= map{ exp($logP[$_]-$logP_max) * $F[$_]/$F_max } (0 .. $#logP);
		my $sum			= 0; map{ $sum += $_ } @sum_terms;
		my $bayes_denom	= $logP_max + $logF_max + log($sum);
		
		# calculate log(P(t|q))
		#printf STDERR "\tq: $q\n";
		
		for(my $i=0; $i<scalar(@taxids); $i++){
			my $t	= $taxids[$i];
			$t_q_logprob{$t}			= {} if(!defined($t_q_logprob{$t})); #init
			$t_q_logprob{$t}->{$q}	= $logP[$i] + $logF[$i] - $bayes_denom;
			#printf STDERR "\tlogP($t|q) = %f\n",$t_q_logprob{$t}->{$q};
		}
	}
	
	# update F(t) = SUM_q P(t|q) * P(q)
	
	#my $Q_num	 = (scalar keys %q_t_logprob) + 0.0;
	foreach my $t(keys %t_q_logprob){
		my $F_update = 0;
		foreach my $q(keys %{$t_q_logprob{$t}}){
			my $logP_t_q		= $t_q_logprob{$t}->{$q};
			my $logP_q		= $q_logprob{$q};
			$F_update 		+= exp($logP_t_q+$logP_q);
		}
		$F{$t}		= $F_update;
	}
	
	# return updated $F
	return %F;
}

# Get log-likelihood given logP(q|t) and F(t)
# USAGE: 
# my $logL	= EM_logL(q_t_logprob=>\%q_t_logprob, F=>\%F);
# 
sub EM_logL{
	my $subid = "EM_logL()";
	# IN
	my (%args)		= @_;	
	my $q_t_logprob 	= $args{q_t_logprob} || die "ERROR: $subid: missing input: q_t_logprob";
	my %q_t_logprob	= %$q_t_logprob;
	my $F			= $args{F} || die "ERROR: $subid: missing input: F";
	my %F			= %$F;
	# OUT
	my $logL			= 0;
	# CONST
	my $INF    		= 9**9**9;
	
	foreach my $q(keys %q_t_logprob){
		
		my @taxids		= sort keys %{$q_t_logprob{$q}};
		my @logP			= map{ $q_t_logprob{$q}->{$_} } @taxids;
		my @F			= map{ $F{$_} } @taxids;
		my @logF			= map{ ($F[$_]>0)? log($F[$_]) : -$INF } (0..$#F);
		# max_t{ logP * logF }
		my @logPlogF		= map{ $logP[$_]+$logF[$_] } (0 .. $#logP);
		my $maxi			= maxind(\@logPlogF);
		my $logP_max		= $logP[$maxi];
		my $logF_max		= $logF[$maxi];
		my $F_max		= $F[$maxi];
		
		my @sum_terms	= map{ exp($logP[$_]-$logP_max) * $F[$_]/$F_max } (0 .. $#logP);
		my $sum			= 0; map{ $sum += $_ } @sum_terms;
		
		$logL			+= $logP_max + $logF_max + log($sum);	
	}
	return $logL;
}

# my %q_t_logprob	= EM_q_t_logprob( annot=>$annot_tsv)
# 
# Calculates alignment probabilities, P(q|t), for each pair of query sequence, q, and reference database taxid, t.
# Returns these in a hash-hash structure.
#
# INPUT
# annot			: Path to annotation TSV. MUST have the following fields: $qseqid,$taxid,$bitscore
# [qseqid]		: Header for query-sequence-ids [qseqid]
# [staxid]		: Header for subject taxids [staxid]
# [bitscore]		: Header for bitscore [bitscore]
# OUTPUT
# %log_prob_q_t	: Hash hash of max alignment probabilities for each q-t pair.
#				  Where P(q|t) = $log_prob_q_t{$q}->{$t}
#
sub EM_q_t_logprob{
	my $subid 	= "EM_q_t_logprob()";
	
	# IN
	my (%args)			= @_;	
	my $annot 			= $args{annot} || die "ERROR: $subid: missing input: annot";
	# OUT
	my %q_t_logprob		= (); # log P(q|t), where q is query sequence, t is reference taxid
	# PARAM
	my $qseqid_col		= $args{qseqid} || 'qseqid';
	my $staxid_col		= $args{staxid}	|| 'staxid';
	my $bitscore_col		= $args{bitscore} || 'bitscore';
	my @headers_musthave	= ($qseqid_col,$staxid_col,$bitscore_col);
	my $temp				= $args{temp} || 1.0;
	print STDERR "\t$subid: temp=$temp\n" if $VERBAL;
	
	
	# START WORKING
	
	# (1) read Annotation TSV to qseqid->staxid->maxbits hash, for each staxid select alignment with max bitscore
	my %q_t_maxbits	= ();
		# read header
	open(IN,"<$annot") or die "Can\'t open $annot: $!\n";
	my $ln=0; my $l=<IN>; chomp($l);
	my @headers		= split(/\t/,$l,-1);
	my %headeri		= map {$headers[$_] => $_} 0..$#headers;
		# check that musthave headers are there
	foreach my $h(@headers_musthave){
		die "ERROR: $subid: expected header $h is missing from $annot" if(!defined($headeri{$h}));
	}
		# read data
	while($l=<IN>){
        	$ln++;
		chomp($l);
		my @sp		= split(/\t/,$l,-1);
		
		if($headeri{$qseqid_col} > $#sp){
			print STDERR "WARNING: $subid: missing qseqid-field on line $ln: skipping\n";
			next;
		}
		if($headeri{$staxid_col} > $#sp){
			print STDERR "WARNING: $subid: missing staxid-field on line $ln: skipping\n";
			next;
		}
		if($headeri{$bitscore_col} > $#sp){
			print STDERR "WARNING: $subid: missing bitscore-field on line $ln: skipping\n";
			next;
		}
		
		my $q	= $sp[$headeri{$qseqid_col}];
		my $t	= $sp[$headeri{$staxid_col}];
		my $bits	= $sp[$headeri{$bitscore_col}];
		
		if( !defined($q_t_maxbits{$q}) ){
			$q_t_maxbits{$q} 	= {};
		}
		if( !defined($q_t_maxbits{$q}->{$t})){
			$q_t_maxbits{$q}->{$t}	= $bits;
		}
		elsif( $q_t_maxbits{$q}->{$t} < $bits){
			$q_t_maxbits{$q}->{$t} 	= $bits;
		}
	}
	close(IN);
	
	# (2) Convert alignment bitscores to alignment probabilities
	#	  Use P(q|s) = exp(bitscore/t)/SUM_alignments( exp(bitscore/t) )
	#		ref: Martin Frith, 2020 (https://doi.org/10.1093/bioinformatics/btz576)
	
	# adjust for temperature
	foreach my $q(keys %q_t_maxbits){
		foreach my $t(keys %{$q_t_maxbits{$q}}){
			$q_t_maxbits{$q}->{$t}	= $q_t_maxbits{$q}->{$t}/$temp;
		}
	}
	
	%q_t_logprob		= ();
	foreach my $q(keys %q_t_maxbits){
		$q_t_logprob{$q}	= {};
		
		my $exp_sum_t	= 0;
		
		# Following are for scores for q-alignemnts
		#
		# max_t( S(t)/temp )
		my @scores		= values %{$q_t_maxbits{$q}};
		my $S_max		= max( \@scores );
		
		# exp max_t( S(t)/temp )
	
		# 	log( 1 + exp[S(t_2)]/exp[S(t_max)] +..+ exp[S(t_n)]/exp[S(t_max)] )
		#	= log( 1 + exp[S(t_2)-S(t_max)] + .. + exp[S(t_n)-S(t_max)])
		#	= log( factor_out_term_sum )
		my $factor_out_term_sum		= 0;		
		foreach my $t(keys %{$q_t_maxbits{$q}}){
			$factor_out_term_sum 	+= exp( ($q_t_maxbits{$q}->{$t}) - $S_max);
		}
		my $log_factor_out_term_sum	= log($factor_out_term_sum);
		
		foreach my $t(keys %{$q_t_maxbits{$q}}){
			
			my $S			= $q_t_maxbits{$q}->{$t}; # S(q|t)
			
			# log P(q|t)= log[  exp(S(q|t)) / SUM_t exp( S(q|t)) ]
			# 			= S(q|t) - log[ SUM_t exp( S(q|t))]
			#			= S(q|t) - S(q|t_max) - log( 1 + exp[S(t_2)-S(t_max)] + .. + exp[S(t_n)-S(t_max)])
			
			$q_t_logprob{$q}->{$t}	= $S - $S_max - $log_factor_out_term_sum;
		}
	}
	
	return %q_t_logprob;
}

# USAGE: maxind(\@values)
sub maxind{
	my @values 	= @{shift(@_)};
	my $max		= (scalar(@values)>0) ? $values[0] : !1;
	my $maxi		= (scalar(@values)>0) ? 0 : -1;
	for(my $i= 0; $i<=$#values; $i++){
		if($values[$i] > $max){
			 $max = $values[$i];
			 $maxi= $i;
		};
	}
	return $maxi;
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

