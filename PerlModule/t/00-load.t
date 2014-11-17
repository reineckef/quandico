#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'QUANDICO' ) || print "Bail out!\n";
}

diag( "Testing QUANDICO $QUANDICO::VERSION, Perl $], $^X" );
