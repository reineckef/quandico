package QUANDICO;
use 5.008;
use strict;
use warnings;
#
use Data::Dumper;
use DateTime;
use DBI;
use DBD::SQLite;
use Env::Path qw<PATH>;
use File::Basename qw<fileparse>;
use File::chdir;
use File::Copy qw<copy move>;
use File::Spec::Functions qw<rel2abs catfile>;
use File::Temp;
use Getopt::Long::Descriptive;
use List::Util qw<min max>;
use Sort::Naturally qw<nsort>;
use version 0.77;
our $VERSION = version->declare(1.12);    ## GIT_VERSION ## our $VERSION = version->declare(<%>);
our $BRANCH  = 'master';                   ## GIT_BRANCH ## our $BRANCH = '<%>';
#
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw<describe_options fileparse copy move rel2abs catfile min max nsort Dumper $CWD show_version abs_path PATH>; # symbols to export always

=head1 NAME

QUANDICO - Package with helper scripts for the R application 'quandico'

=head1 VERSION

Version 1.12

=cut

=head1 SYNOPSIS

This module holds documentation and some commonly used routines for the script 
distributed along with this package. This modules is not meant to be of any 
general use instead of providing a container for documentation and provide 
the means to package scripts for easy installation (with dependencies).

=over 3

=item * quandico - driver script to prepare and run the R application

=item * qcluster - script to group neighboring amplicon count data into groups

=item * qgetcounts - script to extract read depth for amplicon coordinates from SAM/BAM

=back

Please check the documentation of these script for detailed help. Generally, running 
these with "--help" will present usage information. Generally, quandico is the main 
driver which is also able to run the other two helper script for you.

=head1 EXPORT

Some functions of commonly used external modules are exported to the helper scripts. 
Please check the code if you need to know the details.

=head1 SUBROUTINES/METHODS

Subroutines contained in here are also only used by the scripts mentioned above, and 
these are not meant to be of any general use. 


=head2 show_version()

This prints the version number and exits.

=cut

sub show_version {
	my $name = $0;
	$name =~ s!^\.\/!!;
	printf "This is %s %s %s\nUse flags -h or --help for help.\n", $name, $VERSION->stringify, $BRANCH ? qq~($BRANCH)~ : '';
	exit;
}


=head2 abs_path() 

This will create an absolute path via rel2abs(), but always uses slashes and never 
backslash, even on Windows. The path is meant to be used inside R.

=cut

sub abs_path {
	my $abs = rel2abs( shift );
	$abs =~ tr/\\/\//;
	return $abs;
}

=head1 AUTHOR

Frank Reinecke, C<< <frank.reinecke at qiagen.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<frank.reinecke at qiage.com>, 
or through the website at L<http://code.google.com/p/quandico>. 


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc QUANDICO


You can also look for information at:

=over 4

=item * Goodle Code (report bugs there)

L<http://code.google.com/p/quandico>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Frank Reinecke.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; version 3 dated June, 2007 or at your option
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available in the source tree;
if not, write to the Free Software Foundation, Inc.,
59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.


=cut

1; # End of QUANDICO
