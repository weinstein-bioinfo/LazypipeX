#!/usr/bin/env perl
#
# UNIT-20 .. UNIT-24 of docs/testing_roadmap.md §4.2 — Lazypipe::Utils.
#
# Pure functions, no databases, no external tools.  Run directly or via prove:
#     prove -I perl tests/perl/utils.t
#
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Basename qw(dirname);
use MIME::Base64 qw(decode_base64);

my $REPO;
BEGIN {
	$REPO = dirname(dirname(dirname(File::Spec->rel2abs(__FILE__))));
}
use File::Spec;
use lib "$REPO/perl";

use Lazypipe::Utils;
use Lazypipe::Utils qw(format_int filebin2uri);

my $FIX  = "$REPO/tests/fixtures";
my $TSV  = "$FIX/tsv/basic.tsv";
my $WORK = tempdir( CLEANUP => 1 );

# ---------------------------------------------------- UNIT-20 tsv readers ---

{
	my %h = read_tsv2hash( $TSV, 'contig', 'taxid' );
	is( scalar keys %h, 5,        'UNIT-20 read_tsv2hash: one entry per row' );
	is( $h{c1}, 11676,            'UNIT-20 read_tsv2hash: c1 -> 11676' );
	is( $h{c3}, 3048202,          'UNIT-20 read_tsv2hash: c3 -> 3048202' );

	# Documented behaviour: for a repeated key the LAST value wins.
	my %last = read_tsv2hash( $TSV, 'taxid', 'contig' );
	is( $last{1239574}, 'c4',     'UNIT-20 read_tsv2hash: repeated key keeps last value' );
}

{
	# read_tsv2kvahash keeps every value, in file order, duplicates included.
	my %h = read_tsv2kvahash( $TSV, 'taxid', 'contig' );
	is_deeply( $h{1239574}, [ 'c2', 'c4' ], 'UNIT-20 read_tsv2kvahash: all values, file order' );
	is_deeply( $h{10298},   [ 'c5' ],       'UNIT-20 read_tsv2kvahash: single-value key' );
}

{
	# read_tsv2kuvahash keeps unique values, sorted.
	my %h = read_tsv2kuvahash( $TSV, 'taxid', 'species' );
	is_deeply( $h{1239574}, [ 'Mamastrovirus 10' ],
		'UNIT-20 read_tsv2kuvahash: duplicate values collapse to one' );
}

{
	my %t = read_tsv2hashtable( $TSV, 'contig' );
	is( scalar keys %t, 5,                       'UNIT-20 read_tsv2hashtable: one row per key' );
	is( $t{c3}->{species},  'Circovirus mink',   'UNIT-20 read_tsv2hashtable: cell by column name' );
	is( $t{c5}->{bitscore}, 90,                  'UNIT-20 read_tsv2hashtable: numeric cell' );
}

{
	my @t = read_tsv2arraytable( $TSV );
	is( scalar @t, 6,                     'UNIT-20 read_tsv2arraytable: header + 5 rows' );
	is_deeply( $t[0], [qw(contig taxid species bitscore qcov)],
		'UNIT-20 read_tsv2arraytable: row 0 is the header' );
	is( $t[3]->[2], 'Circovirus mink',    'UNIT-20 read_tsv2arraytable: [3][2]' );

	# A short row is skipped with a warning rather than silently padded.
	my @warn;
	local $SIG{__WARN__} = sub { push @warn, @_ };
	my $ragged = "$FIX/tsv/ragged.tsv";
	my @r = do {
		open( my $olderr, '>&', \*STDERR ) or die;
		open( STDERR, '>', "$WORK/ragged.err" ) or die;
		my @x = read_tsv2arraytable( $ragged );
		open( STDERR, '>&', $olderr ) or die;
		@x;
	};
	is( scalar @r, 3, 'UNIT-20 read_tsv2arraytable: ragged row dropped' );
	my $err = do { local ( @ARGV, $/ ) = ("$WORK/ragged.err"); <> };
	like( $err, qr/invalid column number/, 'UNIT-20 read_tsv2arraytable: warns about the ragged row' );
}

{
	# A missing column must die with the documented ERROR: message, not return empty.
	my $err = '';
	eval { my %h = read_tsv2hash( $TSV, 'nosuchcol', 'taxid' ); 1 } or $err = $@;
	like( $err, qr/^ERROR: key col=nosuchcol/, 'UNIT-20 read_tsv2hash: missing key column dies' );

	$err = '';
	eval { my %h = read_tsv2hash( $TSV, 'contig', 'nosuchcol' ); 1 } or $err = $@;
	like( $err, qr/^ERROR: value col=nosuchcol/, 'UNIT-20 read_tsv2hash: missing value column dies' );
}

# ------------------------------------------- UNIT-21 column/line indexing ---

is( colind( $TSV, 'contig' ),   1, 'UNIT-21 colind: first column is 1-based' );
is( colind( $TSV, 'species' ),  3, 'UNIT-21 colind: middle column' );
is( colind( $TSV, 'qcov' ),     5, 'UNIT-21 colind: last column' );
is( colind( $TSV, 'nosuch' ),  -1, 'UNIT-21 colind: absent column returns -1' );

is( mcolind( $TSV, '^bit' ),    4, 'UNIT-21 mcolind: regex match' );
is( mcolind( $TSV, 'zzz' ),    -1, 'UNIT-21 mcolind: no match returns -1' );

is( ncol( $TSV ),   5, 'UNIT-21 ncol: column count' );
is( nlines( $TSV ), 6, 'UNIT-21 nlines: header + 5 rows' );

# --------------------------------------------------- UNIT-22 arithmetic ---

is( median( [ 1, 2, 3 ] ),       2,   'UNIT-22 median: odd length' );
is( median( [ 1, 2, 3, 4 ] ),    2.5, 'UNIT-22 median: even length averages the middle pair' );
is( median( [ 5, 1, 3 ] ),       3,   'UNIT-22 median: unsorted input is sorted first' );
is( median( [ 7 ] ),             7,   'UNIT-22 median: single element' );

{
	# Known rough edge: median([]) indexes an empty list, so it emits two
	# "uninitialized value" warnings and returns 0 — indistinguishable from a
	# genuine median of 0.  Marked TODO rather than failing the tier: no caller
	# in the tree passes an empty array today.  See docs/testing_roadmap.md §12.
	local $TODO = 'median() has no empty-array guard';
	my @warn;
	local $SIG{__WARN__} = sub { push @warn, @_ };
	my $m = median( [] );
	is( scalar @warn, 0, 'UNIT-22 median: empty array does not warn' );
	ok( !defined $m, 'UNIT-22 median: empty array returns undef, not 0' );
}

is( max( [ 3, 9, 2 ] ),     9,  'UNIT-22 max: positive values' );
is( max( [ -5, -2, -9 ] ), -2,  'UNIT-22 max: all negative' );
ok( !max( [] ),                 'UNIT-22 max: empty array is false' );

is( sum( [ 1, 2, 3 ] ), 6, 'UNIT-22 sum: positive values' );
is( sum( [] ),          0, 'UNIT-22 sum: empty array is 0' );
is( sum( [ -2, 2 ] ),   0, 'UNIT-22 sum: cancelling values' );

is( format_int( 1234,    ' ' ), '1 234',     'UNIT-22 format_int: four digits' );
is( format_int( 1234567, ',' ), '1,234,567', 'UNIT-22 format_int: seven digits' );
is( format_int( 123,     ',' ), '123',       'UNIT-22 format_int: no separator needed' );
is( format_int( 1000,    ',' ), '1,000',     'UNIT-22 format_int: trailing zeros' );

# ------------------------------------------------------- UNIT-23 writers ---

{
	my @a = ( 'alpha', 'beta', 'gamma' );
	my $f = "$WORK/array.txt";
	write_array2file( $f, \@a );
	ok( -s $f, 'UNIT-23 write_array2file: file written and non-empty' );
	my @back = read_file2array( $f );
	is_deeply( \@back, \@a, 'UNIT-23 write_array2file: round-trips through read_file2array' );
}

{
	my %h = ( b => 2, a => 1, c => 3 );
	my $f = "$WORK/hash.tsv";
	write_hash2file( $f, \%h );
	my @back = read_file2array( $f );
	is_deeply( \@back, [ "a\t1", "b\t2", "c\t3" ],
		'UNIT-23 write_hash2file: keys sorted, tab separated' );

	my $f2 = "$WORK/hash.csv";
	write_hash2file( $f2, \%h, ',' );
	my @back2 = read_file2array( $f2 );
	is_deeply( \@back2, [ 'a,1', 'b,2', 'c,3' ], 'UNIT-23 write_hash2file: custom separator' );
}

{
	my @table = (
		[ 'contig', 'taxid' ],
		[ 'c1',     11676   ],
		[ 'c2',     1239574 ],
	);
	my $html = '';
	open( my $fh, '>', \$html ) or die "cannot open in-memory handle: $!";
	write_table_html(
		table       => \@table,
		table_attrs => 'class="sortable"',
		th_attrs    => [ 'class="c1"', 'class="c2"' ],
		fh          => $fh,
	);
	close($fh);

	like( $html, qr/<table [^>]*class="sortable"/, 'UNIT-23 write_table_html: table attrs emitted' );
	like( $html, qr{<thead>.*contig.*</thead>}s,   'UNIT-23 write_table_html: header row in thead' );
	is( scalar( () = $html =~ /<th /g ),  2,       'UNIT-23 write_table_html: one th per column' );
	is( scalar( () = $html =~ /<tr>/g ),  3,       'UNIT-23 write_table_html: one tr per row' );
	is( scalar( () = $html =~ /<td>/g ),  4,       'UNIT-23 write_table_html: one td per body cell' );
	like( $html, qr{</table>\s*$},                 'UNIT-23 write_table_html: table closed' );
}

# ------------------------------------------------- UNIT-24 base64 data URI ---

{
	# 100 bytes covering the full 0..255 range, including NULs and high bytes,
	# so a text-mode read or an encoding slip would corrupt the round trip.
	my $bytes = join( '', map { chr( ( $_ * 7 ) % 256 ) } 0 .. 99 );
	is( length($bytes), 100, 'UNIT-24 fixture is 100 bytes' );

	my $f = "$WORK/blob.bin";
	open( my $out, '>', $f ) or die;
	binmode($out);
	print $out $bytes;
	close($out);

	my $uri = filebin2uri($f);
	like( $uri, qr{^data:application/gzip;base64,},
		'UNIT-24 filebin2uri: IGV-compatible data URI prefix' );
	unlike( $uri, qr/\n/, 'UNIT-24 filebin2uri: no embedded newlines' );

	( my $b64 = $uri ) =~ s{^data:application/gzip;base64,}{};
	is( decode_base64($b64), $bytes, 'UNIT-24 filebin2uri: base64 round-trips to the original bytes' );
}

done_testing();
