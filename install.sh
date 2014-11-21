#!/bin/sh
#
###############################################################################
#
# This script is written to download and install the R package 'quandico'
# and the Perl module 'QUANDICO' with required dependencies and test data.
# 
# Obviously, this does only work on Unix/Linux systems and should be seen
# as a shortcut for all manual steps described in the 'README.md' file.
#
# *** This is alpha software and not guaranteed to work! ***
# 
# Adjust these to match your needs. You may want to skip the BAM file download
# if you are limited in disk space (this needs ~500 MB of space) or your
# download speed is low.
#
# *** CONFIGURATION *** // START
#
REXE=/usr/bin/R;        # point to the R installation to use
PERL=/usr/bin/perl;     # point to the Perl installation to use
BIGDATA=0;              # set to 1 if you would like to download BAM files (~474 MB)
REFGENE=0;              # set to 0 if you would like to skip refGene.txt download (~15 MB)
TESTRUN=1;              # set to 0 if you would like to skip the test of 'quandico'
CPANM=0;                # set to 1 if you prefer App::cpanminus over classic 'make'
SAMTOOLS=~/bin/samtools # specify which samtools to use (must be >= version 1.1)
#
# *** CONFIGURATION *** // END
###############################################################################
#
# Some Shortcuts for Output Formatting
#
function action() { echo -e "\n\033[1;7m$1\033[0m\n" ; }
function skip()   { echo -e "\n\033[2;7m$1\033[0m\n" ; }
function good()   { echo -e "\033[32m$1\033[0m" ; }
function bad()    { echo -e "\033[31m$1\033[0m" ; }
function bail()   { echo -e "\033[31;7m$1\033[0m" ; exit; }

#
###############################################################################
#
# STEP 0: Basic prereqs
#
# Need R
action "[ PREREQ ] Testing for R, Perl and Samtools ... "
if [ -x $REXE ]; then 
  good "$REXE is executable and will be used"
else
  bail "$REXE is not executable"
fi
# Need Perl
if [ -x $PERL ]; then 
  good "$PERL is executable and will be used"
else
  bail "$PERL is not executable"
fi

# Need samtools for BAM
if [ -x $SAMTOOLS ]; then 
  STV=$($SAMTOOLS --version | perl -ne "if (/samtools ([\d\.]*)/gsm) { \$v = \$1 >= 1.1 ? 'ok: '. \$1 : 'old: ' . \$1 }; END { print \$v }")
  if [[ $STV == ok* ]]; then
    good "$SAMTOOLS is executable, the version is $STV (>= 1.1)."
  else
    bad "
Cannot run on BAM files without samtools version >= 1.1!
You have samtools $STV. Please consider upgrading to the most recent version 
of 'samtools' from http://www.htslib.org to use BAM files as input.
"

    if (( $BIGDATA + $TESTRUN == 2 )); then
      bail "Cannot run on BAM (BIGDATA=$BIGDATA, TESTRUN=$TESTRUN) files without samtools!"
    fi
  fi

else

  bad "
Cannot use BAM files without samtools (version >= 1.1)!
Please consider downloading the most recent version of 'samtools'
from http://www.htslib.org to use BAM files as input.
"

  if (( $BIGDATA + $TESTRUN == 2 )); then
    bail "Cannot run on BAM (BIGDATA=$BIGDATA, TESTRUN=$TESTRUN) files without samtools!"
  fi
fi

if [ -e latest ]; then 
  action "[  CLEAN  ] Removing previous 'latest' file ... "
  rm latest
fi

###############################################################################
#
# STEP 1: Find the links for the latest release
#

action "[ STEP  1  ] Looking for the latest version ... "
wget -nv https://github.com/reineckef/quandico/releases/latest -O latest

###############################################################################
#
# STEP 2: Downloads
#
RPACK=$(grep -Po "(releases/download.*.quandico.*\.tar\.gz)" < latest)
action "[ STEP  2a ] Downloading the R package: $RPACK ... "
wget -nc -nv --no-check-certificate https://github.com/reineckef/quandico/$RPACK

PERLM=$(grep -Po "(releases/download.*.QUANDICO.*\.tar\.gz)" < latest)
action "[ STEP  2b ] Downloading the Perl module: $PERLM ... "
wget -nc -nv --no-check-certificate https://github.com/reineckef/quandico/$PERLM

action "[ STEP  2c ] Downloading example data ... "
wget -nc -nv --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITQjBUZUFVcy1BNFk -O CNA902Y.bed
wget -nc -nv --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITUDZ0aXFBd2pDR0k -O M62_NA13019.tsv
wget -nc -nv --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITWU9xTndtZE5iOEE -O M62_NA12878.tsv
wget -nc -nv --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITYTduNDY4azJNZXM -O M62_NA13019.clustered
wget -nc -nv --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITWm1FS0duczVlejQ -O M62_NA12878.clustered

if [ $BIGDATA == 1 ]; then
  action "[ STEP  2d ] Downloading BAM files (474 MB), please wait ... "
  wget -nc -nv --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITMzNyakhveTh3UVE -O M62_NA13019.bam
  wget -nc -nv --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITSnU1TlVRSjRXRHM -O M62_NA12878.bam
else
  skip "[ STEP  2d ] Skipping BAM files "
fi


if [ $REFGENE == 1 ]; then
  action "[ STEP  2e ] Downloading refGene ... "
  wget -nc -nv -O- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz | gunzip > refGene.txt
  WITHREF="--cp names=refGene.txt"
else
  skip "[ STEP  2e ] Skipping refGene.txt "
  WITHREF=""
fi

###############################################################################
#
# STEP 3: Install R dependencies
#
action "[ STEP  3  ] Installing R dependencies..."
$REXE --slave --vanilla <<END_OF_R

# options(warn=-1)
# update.packages(repos="http://cran.rstudio.com/", ask=FALSE) # this will update ALL packages

for (pkg in c("MASS", "grid", "ggplot2", "gridExtra", "naturalsort", "scales")) {
  check <- require( pkg , character.only = TRUE )
  if (check) {
    message(cat(pkg, 'is already installed', sep=" "))
    update.packages(pkg, repos="http://cran.rstudio.com/", ask=FALSE)
  } else {
    message(cat("installing", pkg, "...", sep=" "))
    install.packages(pkg, repos="http://cran.rstudio.com/")
  }
}
q()

END_OF_R

echo "OK."

###############################################################################
#
# STEP 4: Install R package
#
RPACK=$(perl -e "@L = split(/\//, '$RPACK'); print \$L[-1];")
action "[ STEP  4  ] Installing 'quandico' R package '$RPACK' ... "
$REXE CMD INSTALL $RPACK
echo "OK."

###############################################################################
#
# STEP 5: Install Perl Module
#
PERLM=$(perl -e "@L = split(/\//, '$PERLM'); print \$L[-1];")
PURE=$(perl -e "print \$1 if '$PERLM' =~ /^(.*?)\.tar\.gz/;")

action "[ STEP  5  ] Installing Perl Module '$PERLM' ... "

if [ $CPANM == 1 ]; then
  echo "Installing Perl Module '$PURE' using cpanm (App::cpanminus) ..."
  cpan App::cpanminus && cpanm -v $PERLM
  echo "OK."
  echo ""
else
  echo "Installing Perl Module '$PURE' the classic way ..."
  tar xvfz $PERLM
  cd $PURE
  perl Makefile.PL
  make
  make test
  make install
  cd ..
  echo "OK."
  echo ""
fi

###############################################################################
#
# STEP 6: Run tests
#

if [ $TESTRUN != 1 ]; then
  skip "[ STEP  6  ] Test run not requested, skipping. "

else
  action "[ STEP  6  ] Testing the installation on demo data $WITHREF ... "

  if [ $BIGDATA == 1 ]; then
    echo "Running 'quandico' on BIG test data, this can take up to ten minutes to complete ..."
    echo "quandico -A CNA902Y.bed -X samtools=$SAMTOOLS --rexe $REXE -s map=M62_NA13019.bam -r map=M62_NA12878.bam -s x=2 -s y=0 -r x=2 -r y=0 $WITHREF"
    time quandico -A CNA902Y.bed -X samtools=$SAMTOOLS --rexe $REXE -s map=M62_NA13019.bam -r map=M62_NA12878.bam -s x=2 -s y=0 -r x=2 -r y=0 $WITHREF
  else
    echo "Running 'quandico' on extracted test data, this should take less than one minute ..."
    echo "quandico --rexe $REXE -s data=M62_NA13019.tsv -r data=M62_NA12878.tsv -s x=2 -s y=0 -r x=2 -r y=0 $WITHREF"
    time quandico --rexe $REXE -s data=M62_NA13019.tsv -r data=M62_NA12878.tsv -s x=2 -s y=0 -r x=2 -r y=0 $WITHREF
  fi

fi

action "[  RESULTS  ] You should see some files named 'sample_vs_reference.[ext]' "
ls -hltr sample_vs_*

good "

The R package and the Perl Module (both with dependencies)
should now be installed. Please check output for errors.

"

