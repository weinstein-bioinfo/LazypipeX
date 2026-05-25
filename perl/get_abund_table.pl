#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $USAGE = "USAGE: $0 [-h|headers --hgabund int --hgtaxid taxid -t tail -w wmodel] annot.tsv contigs.idxstats> readn+contign+taxid\n\n".
			"-h|headers        : annot.tsv has headers\n".
			"--hgabund int     : host-genome read number\n".
			"--hgtaxid int     : host-genome taxid\n".
			"-t|conttail       : contig score tail\n".
			"-w|weights        : weighting model: taxacount|bitscore|bitscore2 [bitscore2]\n".
			"annot.tsv         : annotation tsv-file with column headers. MUST include columns in this order: contig,staxid,bitscore,contig.len\n".
			"contigs.idxstats  : SAM idxstats file output by read realignment\n\n";
			   
my $headers 		= !1;
my $hgabund 		= !1;
my $hgtaxid 		= !1;
my $wmodel  		= 'taxacount'; # options: taxacount|bitscore|bitscore2
my $cont_score_tail	= 0;

GetOptions(	'headers|h'		=> \$headers,
			'hgabund=i'		=> \$hgabund,
			'hgtaxid=i'		=> \$hgtaxid,
			'weights|w=s'	=> \$wmodel,
			'conttail|t=i'	=> \$cont_score_tail);	
$wmodel 	= lc($wmodel);


if($wmodel =~ /taxa|taxon/gi){ $wmodel = 'taxacount';}
if( !(($wmodel =~ 'taxacount') || ($wmodel eq 'bitscore') || ($wmodel eq 'bitscore2')) ){
	print STDERR "ERROR: invalid argument: --weights\n\n";
	print STDERR $USAGE; exit(1);
}
if($cont_score_tail > 10 || $cont_score_tail < 0){
	print STDERR "ERROR: invalid argument: --cont_score_tail $cont_score_tail. Use values in [0,10]\n\n";
	print STDERR $USAGE; exit(1);
}

if( scalar(@ARGV)<2){
	print STDERR "ERROR: missing input files\n\n";
	print STDERR $USAGE; exit(1);
}
my $annot	 	= shift(@ARGV);
my $idx_file 	= shift(@ARGV);


my %tax2cont;		# $taxid->$contid->0/1		: taxid mapped to a unique set of contig-ids (%hash->%hash->0/1)
my %cont2tax;		# $contid->$taxid->0/1		: contig-id mapped to a unique set of tax-ids (%hash->%hash->0/1)
my %cont2score;		# $contid->$score			: sum of bitscores for a given contig
my %cont2tax2score; # $contid->$taxid->$score	: sum of bitscores for a given contig and taxid (%hash->%hash->double)
my %cont2clen;		# $contid->int				: contig-id mapped to contig length
my %tax2score;		# $taxid->$score			: sum of bitscores for a given taxon

print STDERR "# reading $annot\n";
open(IN, "<$annot") or die "failed to open $annot\n";
my $ln  = 0;
while(my $l=<IN>) {
	$ln++;
	chomp($l);
	if($ln==1 && $headers){
		next;
	}
	my ($cont,$tax,$score,$clen)	= (0,0,0,0);
	($cont,$tax,$score,$clen) 	= split(/\t/,$l,-1);
	
	$cont2clen{$cont}	= $clen;
	
	if($tax !~ m/[0-9]+/){
		print STDERR "\tWARNING: skipping invalid taxid:$tax\n";
		next; # not valid taxid
	}
	
	if( defined($tax2cont{$tax}) ){
		$tax2cont{$tax}->{$cont} = 1;
	}
	else{
		my %tmp =($cont => 1);
		$tax2cont{$tax} = \%tmp;
	}
	
	if( defined($cont2tax{$cont}) ){
		$cont2tax{$cont}->{$tax} = 1;
	}
	else{
		my %tmp = ($tax => 1);
		$cont2tax{$cont} = \%tmp;
	}
	
	if( defined($cont2score{$cont}) ){
		$cont2score{$cont} += $score;
	}
	else{
		$cont2score{$cont} = $score;
	}
	
	if( defined($cont2tax2score{$cont})  ){
		if( defined($cont2tax2score{$cont}->{$tax}) ){
			$cont2tax2score{$cont}->{$tax} += $score;
		}
		else{
			$cont2tax2score{$cont}->{$tax} = $score;
		}
	}
	else{
		my %tmp = ($tax => $score);
		$cont2tax2score{$cont} = \%tmp;
	}
	
	if( defined($tax2score{$tax}) ){
		$tax2score{$tax} += $score;
	}
	else{
		$tax2score{$tax} = $score;
	}
}
close(IN);


# DELETING TAXA THAT ARE IN THE "TAIL" OF SCORE DIST FOR EACH CONTIG
if($cont_score_tail > 0){
  foreach my $cont(keys %cont2tax2score ){

	my $csum 	= 0;
	my %tmp_hash 	= %{$cont2tax2score{$cont}};
	
	foreach my $tax(sort { $tmp_hash{$a} <=> $tmp_hash{$b} } keys %tmp_hash ){	
		$csum 	+= $tmp_hash{$tax};
		if($csum < $cont2score{$cont}*($cont_score_tail/100.0)){
			#print STDERR "deleting\t$cont\t$tax\n";
			delete($cont2tax{$cont}->{$tax});
			delete($tax2cont{$tax}->{$cont});
			delete($cont2tax2score{$cont}->{$tax});
		}
		else{
			last;
		}
	}	
  }
}


print STDERR "# reading $idx_file\n";
my %cont2rn;
open(IN,"<$idx_file") or die "Can\'t open $idx_file: $!\n";
while(my $l=<IN>){
	chomp($l);
	my @sp = split(/\t/,$l,-1);
	if($sp[0] eq "*"){
		next; 	# unaligned reads
	}
	elsif( $sp[0] =~ m/^contig=([A-Za-z0-9\.]+)[_\s]?/ ){
		$cont2rn{$1} = $sp[2];
	}
	elsif( $sp[0] =~ m/^([A-Za-z0-9\.]+)[_\s]?/ ){
		$cont2rn{$1} = $sp[2];
	}
	else{
		$cont2rn{$sp[0]} = $sp[2];
	}	
}
close(IN);


# ESTIMATE ABUNDACIES:
my %tax2readn;
my %tax2contn;
my %tax2clen;
foreach my $tax(sort keys %tax2cont){
	my @contids		= keys %{$tax2cont{$tax}};
	my $readn		= 0;
	my $clen			= 0;
	
	foreach my $cont(@contids){
		my $contig_taxa_score 	= 0;
		foreach my $taxon(keys %{$cont2tax{$cont}}){	
			$contig_taxa_score += $tax2score{$taxon};
		}
	
		if(defined($cont2rn{$cont})){
			if($wmodel eq 'taxacount'){
				$readn += $cont2rn{$cont} / scalar(keys %{$cont2tax{$cont}} ); # assign to this taxon readn assigned to linked contig, divided by number of taxa linked to that contig
			}
			if($wmodel eq 'bitscore'){
				if($cont2score{$cont} == 0){
					print STDERR "\tWARNING: zero alignment score for contid=$cont, taxon=$tax\n";
					next;
					#print STDERR "contid\t:$cont\n";
					#print STDERR "taxon\t:$tax\n";
				}
				# for each contig
				#   weight each taxon by sum-bitscore for the contig-taxon relative to sum-bitscore for the contig
				$readn += $cont2rn{$cont} * ($cont2tax2score{$cont}->{$tax})/($cont2score{$cont});
			}
			if($wmodel eq 'bitscore2'){
				# for each contig
				#	weight each taxon by sum-bitscore for the taxon relateive to sum-bitscore for all taxa associate with the contig
				$readn += $cont2rn{$cont} * ($tax2score{$tax}/$contig_taxa_score);
			}
		}
		$clen	+= $cont2clen{$cont} || 0;
	}
	$tax2readn{$tax} 	= $readn;
	$tax2contn{$tax} 	= scalar(@contids);
	$tax2clen{$tax}		= $clen;
	#print STDERR "TAXON=$tax\tREADN=$readn\n";
}

if($hgabund && $hgtaxid){
	if(!defined($tax2readn{$hgtaxid})){
		$tax2readn{$hgtaxid}		= 0;
		$tax2contn{$hgtaxid}		= 0;
		$tax2clen{$hgtaxid}		= 0;	
	}
	$tax2readn{$hgtaxid} += $hgabund;
}

print "readn\tcontign\tassembly.size\ttaxid\n";
foreach my $tax(sort {$tax2readn{$b}<=>$tax2readn{$a}} keys %tax2readn){
	if($tax2readn{$tax} == 0){
		next;
	}
	my $readn		= $tax2readn{$tax};
	my $contign		= $tax2contn{$tax} || 0;
	my $contiglen	= $tax2clen{$tax} || 0;
	print "",join("\t",$readn,$contign,$contiglen,$tax),"\n";
}




