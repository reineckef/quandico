use 5.008;
use strict;
use warnings;
use ExtUtils::MakeMaker;
WriteMakefile(
	NAME          => 'QUANDICO',
	AUTHOR        => q{Frank Reinecke <frank.reinecke@qiagen.com>},
	VERSION_FROM  => 'lib/QUANDICO.pm',
	ABSTRACT_FROM => 'lib/QUANDICO.pm',
	LICENSE       => 'gpl',
	EXE_FILES     => ['scripts/quandico', 'scripts/qcluster', 'scripts/qgetcounts'],

	#	PL_FILES  => {},
	PREREQ_PM => {
		'DBD::SQLite'               => 0,
		'DBI'                       => 1.622000,
		'Data::Dumper'              => 0,
		'DateTime'                  => 1.030000,
		'Env::Path'                 => 0,
		'File::Basename'            => 2.770000,
		'File::Copy'                => 2.140000,
		'File::Spec::Functions'     => 3.400000,
		'File::Temp'                => 0.230400,
		'File::chdir'               => 0.100400,
		'Getopt::Long::Descriptive' => 0.091000,
		'IPC::Run'                  => 0,
		'IPC::Cmd'                  => 0,
		'List::Util'                => 1.380000,
		'List::MoreUtils'           => 0,
		'Sort::Naturally'           => 1.020000,
		'version'                   => 0.770000,
		'Term::ANSIColor'           => 2.000000,
		'Test::More'                => 0,
		'Test::Script'              => 1
	},
	dist  => {COMPRESS => 'gzip -9f', SUFFIX => 'gz'},
	clean => {FILES    => 'QUANDICO-*'},
);
