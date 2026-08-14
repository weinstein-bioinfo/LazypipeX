#!/usr/bin/env perl
#
# UNIT-25 of docs/testing_roadmap.md §4.2 — CIGAR helpers in Lazypipe::SeqAn.
#
# These feed query coverage, which gates min_qcov_abund / min_qcov_annot and so
# every abundance number in the reports.  Expected values live in
# tests/fixtures/cigar.cases.tsv and are derived from the SAM v1 spec by hand:
#
#     qlen  consumes query:      M I S H = X
#     rlen  consumes reference:  M D N = X
#     alen  alignment columns:   M I D N = X   (clipping and padding excluded)
#     qcov  (M I = X) / qlen
#     pide  (=) / alen           (requires a SAMv1 CIGAR that spells out '=')
#
#     prove -I perl tests/perl/cigar.t
#
use strict;
use warnings;
use Test::More;
use File::Spec;
use File::Basename qw(dirname);

my $REPO;
BEGIN { $REPO = dirname( dirname( dirname( File::Spec->rel2abs(__FILE__) ) ) ) }
use lib "$REPO/perl";

use Lazypipe::SeqAn;

my $CASES = "$REPO/tests/fixtures/cigar.cases.tsv";
open( my $fh, '<', $CASES ) or die "cannot open $CASES: $!";
my $header = <$fh>;
chomp($header);
is( $header, join( "\t", qw(cigar qlen qcov rlen alen pide) ),
	'UNIT-25 fixture header is the expected column order' );

my $EPS = 1e-9;
my $n   = 0;

while ( my $line = <$fh> ) {
	chomp($line);
	next unless length $line;
	my ( $cigar, $qlen, $qcov, $rlen, $alen, $pide ) = split( /\t/, $line, -1 );
	$n++;

	is( cigar2qlen($cigar), $qlen, "UNIT-25 cigar2qlen($cigar) == $qlen" );
	is( cigar2rlen($cigar), $rlen, "UNIT-25 cigar2rlen($cigar) == $rlen" );
	is( cigar2alen($cigar), $alen, "UNIT-25 cigar2alen($cigar) == $alen" );

	ok( abs( cigar2qcov($cigar) - $qcov ) < $EPS,
		sprintf( 'UNIT-25 cigar2qcov(%s) == %s (got %s)', $cigar, $qcov, cigar2qcov($cigar) ) );
	ok( abs( cigar2pide($cigar) - $pide ) < $EPS,
		sprintf( 'UNIT-25 cigar2pide(%s) == %s (got %s)', $cigar, $pide, cigar2pide($cigar) ) );
}
close($fh);
cmp_ok( $n, '>=', 8, 'UNIT-25 fixture exercises at least 8 CIGAR shapes' );

# Soft clipping must reduce coverage but not alignment length: this is the
# distinction that min_qcov_annot depends on.
cmp_ok( cigar2qcov('10S80M10S'), '<', cigar2qcov('100M'),
	'UNIT-25 soft clipping lowers query coverage' );
is( cigar2alen('10S80M10S'), 80, 'UNIT-25 soft clipping is excluded from alignment length' );
is( cigar2qlen('10S80M10S'), 100, 'UNIT-25 soft clipping is included in query length' );

# Hard clipping consumes query length in this implementation, exactly as soft
# clipping does; pinning it so a future change is a deliberate one.
is( cigar2qlen('5H20=5H'), 30, 'UNIT-25 hard clipping counts toward query length' );

# An unpadded 'M' CIGAR carries no match/mismatch detail, so percent identity is
# 0 by construction — callers must not read that as "0% identity".
is( cigar2pide('100M'), 0, 'UNIT-25 pide is 0 for an M-only CIGAR (no = operators)' );
cmp_ok( cigar2pide('50=2X48='), '>', 0.9, 'UNIT-25 pide is meaningful for a SAMv1 CIGAR' );

done_testing();
