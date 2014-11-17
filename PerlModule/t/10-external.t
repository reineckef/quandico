#!perl 
use Test::More tests => 3;
use Data::Dumper;
use Env::Path qw<PATH>;
use File::Spec;
use List::MoreUtils qw<uniq>;
use IPC::Cmd qw<can_run>;
use IPC::Run qw<run timeout>;
use Term::ANSIColor qw<color>;
my ( $in, $out, $err ) = ('','','');

BEGIN {
	use_ok('QUANDICO') || print "Bail out!\n";
}
diag("Testing external tools for QUANDICO $QUANDICO::VERSION, Perl $], $^X");

# check for a useable R installation
my @R;
my @SAMTOOLS;
my $PATH = Env::Path->PATH;
foreach my $path ( uniq $PATH->List ) {
	my $rexe = catfile( $path, 'R' );
	my $samt = catfile( $path, 'samtools' );

	# R
	if ( can_run( $rexe, "R executable" ) ) {
		run [$rexe, " --version"], \$in, \$out, \$err, timeout(3);    # or die "cat: $?"

		# diag(Dumper $in, $out, $err);
		my @V = grep { /version/i } split( /[\n\r]/, join( $in, $out, $err ) );

		# diag (Dumper \@V);
		push @R, [$rexe, $V[0]];

		# diag("Found R executable $rexe # $V[0]");
	} ## end if ( can_run( $rexe, "R executable"...))

	# SAMTOOLS
	if ( can_run( $samt, "Samtools executable" ) ) {
		run [$samt], \$in, \$out, \$err, timeout(3);    # or die "cat: $?"
		my @V = grep { /version/i } split( /[\n\r]/, join( $in, $err ) );
		@V = grep { /[Vv]ersion/ } @V;
		push @SAMTOOLS, [$samt, $V[0]];

		# diag("Found samtools: $samt # $V[0]");
	} ## end if ( can_run( $samt, "Samtools executable"...))
} ## end foreach my $path ( uniq $PATH->List)

#ok( @R,        'Testing for available R executable' );
#ok( @SAMTOOLS, 'Testing for available samtools executable' );
SKIP: {
	skip "R", 1 unless @R;
	ok(@R);
}
if (@R) {
	diag("Found at least one executables(s) of 'R'");    # in path, here is a list:");
	print STDERR @R > 1 ? color('cyan') : color('green');
	map { diag( sprintf "%s = %s\n", $_->[0], $_->[1] ) } @R;
	print STDERR color('reset');
}
else {
	print STDERR color('red');
	diag(
		<<DIAG

The 'R' programming language is required to *run* 'quandico' - but not required to just *install* 
this package. To be able to use the CNV calling algorithm implemented in the R package 'quandico', 
please download and install (or add to your path) 'R' from: http://www.r-project.org/.

Note: You will *not* be able to run 'quandico' without R!

DIAG
	);
	print STDERR color('reset');
} ## end else [ if (@R) ]
SKIP: {
	skip 'samtools', 1 unless @SAMTOOLS;
	ok(@SAMTOOLS);
}
if (@SAMTOOLS) {
	diag("Found at least one executable(s) of 'samtools'");    # in path, here is a list:");
	print STDERR @SAMTOOLS > 1 ? color('cyan') : color('green');
	map { diag( sprintf "%s = %s\n", $_->[0], $_->[1] ) } @SAMTOOLS;
	print STDERR color('reset');
}
else {
	print STDERR color('yellow');
	diag(
		<<DIAG

The executable 'samtools' is required by 'qgetcounts' - but not required to 'install' this package.
To be able to use 'qgetcounts' to extract read counts from SAM/BAM files, please download and 
install (or add to your path) the most recent version of 'samtools' from: http://www.htslib.org/.

Note: You will be able to use 'quandico' and 'qcluster' without 'samtools', but *not* 'qgetcounts'!

DIAG
	);
	print STDERR color('reset');
} ## end else [ if (@SAMTOOLS) ]
