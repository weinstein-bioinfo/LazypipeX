#!/usr/bin/perl
use strict;
use warnings;
use YAML::Tiny;

#
# Config linter for the LazypipeX Tier 2 database tests.
# Implements the config-consistency half of DB-03 .. DB-05, docs/testing_roadmap.md §5.1.
#
#     config_lint.pl paths|fields|strategies [config.yaml]
#
# Prints one finding per line as "<severity>\t<message>" and exits 0 whether or
# not anything was found: the tier decides how a finding is reported, so a lint
# problem must never be confused with a linter that crashed.  Severities:
#
#   FAIL   this installation is broken and the user has to fix it
#   DRIFT  a known defect in LazypipeX itself (docs/testing_roadmap.md §12),
#          which the tier reports as a TAP TODO rather than as a failure
#
# Exit 2 is reserved for the linter failing to do its job — bad arguments, or a
# config that cannot be read or parsed.
#

my %ALLOWED_SEARCH = map { $_ => 1 }
	qw( minimap blastn blastp blastx diamondp diamondx hmmscan sans );

# Sections whose entries name a database on disk, with the fields each entry
# must carry (DB-04).
my %REQUIRED = (
	'ann.databases'  => [ qw( db name search ) ],
	'host.databases' => [ qw( db accession url ) ],
);

my $mode   = shift(@ARGV);
my $config = shift(@ARGV) // 'config.yaml';

if( !defined($mode) || $mode !~ /^(paths|fields|strategies)$/ ){
	print STDERR "USAGE: $0 paths|fields|strategies [config.yaml]\n";
	exit 2;
}
if( !(-e $config) ){
	print STDERR "ERROR: no such config file: $config\n";
	exit 2;
}

my $yaml = eval { YAML::Tiny->read($config) };
if( $@ || !defined($yaml) || !defined($yaml->[0]) ){
	print STDERR "ERROR: cannot parse $config: ".($@ || 'empty document')."\n";
	exit 2;
}
my %opt = %{ $yaml->[0] };

sub finding { printf "%s\t%s\n", $_[0], $_[1]; }

# Entries of a section, or an empty hash if the section is absent: a missing
# section is reported once by the caller rather than as a fault of every check.
sub section {
	my ($name) = @_;
	return ( ref($opt{$name}) eq 'HASH' ) ? %{ $opt{$name} } : ();
}

# ------------------------------------------------------------- DB-03 paths ---
#
# Every db path must start with an environment variable, that variable must be
# set, and a section must not mix variables.  The case this was written for:
# one host genome sat under $data/hostgen/ while the other 22 used
# $hostgenomes/, which is how a database ends up in a directory nobody backs up.
# That entry has since been moved (§12 item 4); the check guards the invariant.
#
if( $mode eq 'paths' ){
	for my $sec ( sort keys %REQUIRED ){
		my %entries = section($sec);
		next if !%entries;

		my %var_of;	# entry key -> env var its path starts with
		my %count;	# env var   -> how many entries use it

		for my $k ( sort keys %entries ){
			my $db = $entries{$k}->{db};
			next if !defined($db);			# reported by 'fields'

			# A value with no path separator is not a location on disk but the
			# name of a remote service (ann.databases:sans).  Nothing to lint.
			next if $db !~ m{/};

			if( $db !~ /^\$(\w+)/ ){
				finding( 'DRIFT', "$sec.$k.db does not start with an environment variable: $db" );
				next;
			}
			my $var = $1;
			$var_of{$k} = $var;
			$count{$var}++;

			if( !defined($ENV{$var}) ){
				finding( 'FAIL', "$sec.$k.db uses \$$var, which is not set in the environment" );
			}
		}

		# The prefix the section agrees on.  Ties break on the variable name so
		# the report is the same on every run.
		next if !%count;
		my ($majority) = sort { $count{$b} <=> $count{$a} || $a cmp $b } keys %count;
		for my $k ( sort keys %var_of ){
			next if $var_of{$k} eq $majority;
			finding( 'DRIFT', "$sec.$k.db uses \$$var_of{$k} while $count{$majority}"
				." of $sec use \$$majority: $entries{$k}->{db}" );
		}
	}
}

# ------------------------------------------------------------ DB-04 fields ---
#
# A missing field is not cosmetic: install_db.pl reads {db} and {url}, and the
# annotation steps read {name} and {search}, so an incomplete entry fails at the
# point of use, deep into a run.
#
elsif( $mode eq 'fields' ){
	for my $sec ( sort keys %REQUIRED ){
		my %entries = section($sec);
		if( !%entries ){
			finding( 'FAIL', "$sec is missing or empty" );
			next;
		}
		for my $k ( sort keys %entries ){
			my $e = $entries{$k};
			if( ref($e) ne 'HASH' ){
				finding( 'FAIL', "$sec.$k is not a mapping" );
				next;
			}
			for my $f ( @{ $REQUIRED{$sec} } ){
				my $v = $e->{$f};
				if( !defined($v) || $v eq '' ){
					finding( 'FAIL', "$sec.$k is missing the required field '$f'" );
				}
			}
			my $search = $e->{search};
			if( $sec eq 'ann.databases' && defined($search) && !$ALLOWED_SEARCH{$search} ){
				finding( 'FAIL', "ann.databases.$k.search is '$search', not one of: "
					. join( ', ', sort keys %ALLOWED_SEARCH ) );
			}
		}
	}
}

# -------------------------------------------------------- DB-05 strategies ---
#
# Strategy values are the --ann1/--ann2 command lines the pipeline substitutes
# for --anns.  A key that no longer exists in ann.databases costs a whole
# annotation round at run time, and lazypipe.pl reports it only once it gets
# there, so the dangling reference is worth catching before the run starts.
#
elsif( $mode eq 'strategies' ){
	my %ann = section('ann.databases');
	my %strat = section('ann.strategies');
	if( !%strat ){
		finding( 'FAIL', 'ann.strategies is missing or empty' );
	}
	if( !%ann ){
		finding( 'FAIL', 'ann.databases is missing or empty' );
	}
	for my $s ( sort keys %strat ){
		my $val = $strat{$s};
		next if !defined($val);
		while( $val =~ /--ann[12]\s+(\S+)/g ){
			for my $ref ( split /,/, $1 ){
				# --ann2 entries are tagged by taxon group, as in "vi:blastn.nt.abv";
				# the database key is what follows the tag.
				$ref =~ s/^\w+://;
				next if $ref eq '';
				if( !exists($ann{$ref}) ){
					finding( 'FAIL', "ann.strategies.$s references '$ref', which is not a key in ann.databases" );
				}
			}
		}
	}
}

exit 0;
