#! /usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);
use YAML::Tiny;
use Cwd qw(getcwd);

#
# DOWNLOAD AND INSTALL SEQ DATABASES FOR LAZYPIPE
#
# credit: Ilya Plyusnin, University of Helsinki, Ilja.Pljusnin@helsinki.fi
#
my $install_dir	= defined($ENV{'LAZYPIPE_INSTALL_DIR'}) ? $ENV{'LAZYPIPE_INSTALL_DIR'} : dirname(__FILE__);
my $config	= "$install_dir/config.yaml";
my $dbname	= undef;
my $url		= undef;
my $path		= undef;
my $FORCE	= 0;
my $VERBAL	= 0;
my $HELP		= 0;

my $usage= 	"\n".
			"USAGE: $0 --db dbname [-u url -p path -f -v]\n\n".
			"Download and install local databases defined in $config\n\n".
			"--config yaml   : Lazypipe config.yaml [$config]\n".
			"--db dbname     : Key or name of the annotation database to install\n".
			"                  Key, latin name or common name of the host database to install\n".
			"                  All keys/names are refences to 'ann.databases' and 'host.databases' specified in config.yaml\n".
			"Examples   :\n".
			"   --db blastn.vi.nt   : install reference database listed under this key in ann.databases\n".
			"   --db mouse          : install host database with commonName 'mouse'\n".
			"   --db hostdbs	        : install all host databases specified in config.yaml\n".
			"   --db taxonomy       : install taxonomy database\n".
			"\n".
			"--url|u url     : use this database url instead of url listed in config.yaml\n".
			"--path|p path   : use this database path instead of path listed in config.yaml\n".
			"--force|f       : Force overwrite [$FORCE]\n".
			"-v              : Verbal mode [$VERBAL]\n".
			"-h              : Print this help\n".
			"\n";

GetOptions(
	'config=s'	=> \$config,
	'db=s' 		=> \$dbname,
	'url|u=s'	=> \$url,
	'path|p=s'	=> \$path,
	'force|f'	=> \$FORCE,
	'v'			=> \$VERBAL,
	'help|h'		=> \$HELP) or die $usage;
if($HELP){
	die $usage;
}
if(!defined($dbname)){
	print STDERR "\nmissing arguments --db dbname\n";
	exit 1;
}
if(!(-e $config)){
	print STDERR "\nmissing config file: $config\n";
	exit 1;
}

my $yaml 	= YAML::Tiny->read( $config );
my %opt 		= %{$yaml->[0]};

# START WORKING

# install from supplied url/path
if( defined($url) && defined($path)){
	my $db	= {name => $dbname, db => $path, url=> $url};
	install_db($db);
	exit(0);
}

# install from config.yaml:ann.databases/host.databases
if(	defined($opt{'ann.databases'}->{"$dbname"})){
	my $db		= $opt{'ann.databases'}->{"$dbname"};
	install_db($db);
	exit(0);
}
if( defined($opt{'host.databases'}->{"$dbname"})){
	my $db		= $opt{'host.databases'}->{"$dbname"};
	install_db($db);
	exit(0);
}
foreach my $k(keys %{$opt{'ann.databases'}}){
	my $db	= $opt{'ann.databases'}->{$k};
	if(	$dbname eq $k 
		|| $dbname eq $db->{name}){
		install_db($db);
		exit(0);
	}
}
foreach my $k(keys %{$opt{'host.databases'}}){
	my $db	= $opt{'host.databases'}->{$k};
	if(	$dbname eq $k 
		|| $dbname eq $db->{latinName}
		|| $dbname eq $db->{commonName}){
		install_db($db);
		exit(0);
	}
}
if($dbname eq 'hostdbs'){
	foreach my $k(sort keys %{$opt{'host.databases'}}){
		my $db		= $opt{'host.databases'}->{$k};
		eval {
        		install_db($db);
        		1;
		} or do{
			my $e = $@;
        		print STDERR "\tINSTALLATION FAILED: $e\n";
		};
		
	}
	exit(0);
}
if($dbname eq 'taxonomy'){
	my $taxonomy	= $opt{taxonomy};
	install_taxonomy($taxonomy);
	exit(0);
}


print STDERR "\n--db $dbname did not match any database in $config\n";





# Usage install_db(\%opt)
#
# INSTALL reference/host-genome database
sub install_db{
	my $db 		= shift(@_);
	
	if($VERBAL){
		print STDERR "\n\tname      : $db->{name}\n";
		print STDERR "\tdatabase  : $db->{db}\n";
		print STDERR "\turl       : $db->{url}\n";
	}
	my $dbdir	= dirname($db->{db});
	my $dbfile	= basename($db->{db});
	my $webfile	= basename($db->{url});
		# expand env vars
	my $dbdir_exp= $dbdir;
	if($dbdir_exp =~ /\$(\w+)/g){
		my $envar	= $1;
		if(defined($ENV{$envar})){
			$dbdir_exp =~ s/\$$envar/$ENV{$envar}/g;
		}
	}
	
	my @dbfiles_ondisk	= glob("$dbdir_exp/$dbfile*");
	if(!$FORCE && scalar(@dbfiles_ondisk)>0){
		print STDERR "\tskipping: found files on disk (use --force to overwrite):\n";
		foreach my $file(@dbfiles_ondisk){
			print STDERR "\t  $file\n";
		}
		return;
	}

	system("mkdir -p $dbdir");
	system_call("wget -q $db->{url} -O $dbdir/$webfile" );
	if($webfile =~ m/\.tar\.gz$/i){
		system_call("tar -C $dbdir -xvzf $dbdir/$webfile" );
	}
}

sub install_taxonomy{
	my $taxonomy 		= shift(@_);
	
	# Check if taxonomy already up-to-date
	my $taxonomy_nodes	= "$taxonomy->{db}/nodes.dmp";
	if( (-e "$taxonomy_nodes") 
			&& ((-M "$taxonomy_nodes") <= $taxonomy->{update_time}) 
			&&  !$FORCE){
		if($VERBAL){
			print STDERR "\ttaxonomy up to date. Use --force to overwrite\n";
		}
		return;
	}
	
	system_call("wget -q $taxonomy->{url} -O $taxonomy->{db}/taxdump.tar.gz");
	system_call("tar -xzf $taxonomy->{db}/taxdump.tar.gz -C $taxonomy->{db}");
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
