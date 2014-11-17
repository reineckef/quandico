#!perl 

use Test::More tests => 4;
use Test::Script;
BEGIN {
    use_ok( 'QUANDICO' ) || print "Bail out!\n";
}

diag( "Testing helper scripts syntax of QUANDICO $QUANDICO::VERSION, Perl $], $^X" );

script_compiles('scripts/qgetcounts', "Helper 'getcounts' compiles OK");
script_compiles('scripts/qcluster', "Helper 'qcluster' compiles OK");
script_compiles('scripts/quandico', "Main script 'quandico' compiles OK");

