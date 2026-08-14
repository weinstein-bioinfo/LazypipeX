#!/usr/bin/perl
use strict;
use warnings;
use YAML::Tiny;

#
# Database inventory for the LazypipeX Tier 2 tests (DB-10 … DB-27, §5.2).
#
#     db_inventory.pl [config.yaml]
#
# Prints one tab-separated row per database entry:
#
#     section  key  engine  resolved_path  installed
#
# where section is ann|host, engine is the entry's `search` value (or "hostgen"
# for background filters), resolved_path has $ENV vars substituted, and
# installed is 1 when the path globs to at least one file.
#
# "Installed" uses exactly lazypipe.pl's own rule — glob("$dbpath*") is
# non-empty — so this inventory agrees with what --databases and --filters
# report, and the per-engine checks run on the same set the user sees listed.
#
# Exit 2 if the config cannot be read or parsed; the tier reports that rather
# than treating an empty inventory as "nothing installed".
#

my $config = shift(@ARGV) // 'config.yaml';

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

# $VAR and ${VAR} both appear in config.yaml paths.  An unset variable is left
# verbatim: the path then fails the -e test below and is reported as not
# installed, which is what DB-03 already flagged as a config error.
sub resolve {
	my ($path) = @_;
	return '' if !defined($path);
	$path =~ s/\$\{(\w+)\}/ defined($ENV{$1}) ? $ENV{$1} : "\${$1}" /ge;
	$path =~ s/\$(\w+)/ defined($ENV{$1}) ? $ENV{$1} : "\$$1" /ge;
	return $path;
}

sub emit {
	my ($section, $key, $engine, $raw) = @_;
	my $path = resolve($raw);
	# A value with no path separator names a remote service, not a file on
	# disk (ann.databases:sans); it can never be "installed" locally.
	my $installed = 0;
	if( $path ne '' && $path =~ m{/} ){
		my @f = glob( "$path*" );
		$installed = ( scalar(@f) > 0 ) ? 1 : 0;
	}
	printf "%s\t%s\t%s\t%s\t%d\n", $section, $key, $engine, $path, $installed;
}

my $ann = ( ref($opt{'ann.databases'}) eq 'HASH' ) ? $opt{'ann.databases'} : {};
for my $k ( sort keys %$ann ){
	emit( 'ann', $k, ( $ann->{$k}{search} // 'unknown' ), $ann->{$k}{db} );
}

my $host = ( ref($opt{'host.databases'}) eq 'HASH' ) ? $opt{'host.databases'} : {};
for my $k ( sort keys %$host ){
	emit( 'host', $k, 'hostgen', $host->{$k}{db} );
}

# The taxonomy database is not a member of either section but every taxid-aware
# step depends on it, so DB-10 … DB-13 need its location from the same place.
if( ref($opt{taxonomy}) eq 'HASH' ){
	emit( 'taxonomy', 'taxonomy', 'taxonomy', $opt{taxonomy}{db} );
}

exit 0;
