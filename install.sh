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
REXE=/usr/bin/R; # point to the R installation to use
BIGDATA="1";     # set to 1 if you would like to download BAM files (~474 MB)
REFGENE="1";     # set to 1 if you would like to download refGene.txt (~15 MB)
TESTRUN="1";     # set to 1 if you would like to run 'quandico' on test data

# optional (for BIG data)
SAMTOOLS=~/bin/samtools # specify which samtools to use (must be >= version 1.1)

###############################################################################
#
# STEP 1: Find the links for the latest release
#
if [ -e latest ]; then
	echo "Cleaning previous 'latest' file..."
	# rm latest
fi

echo "Looking for the latest version..."
wget -q https://github.com/reineckef/quandico/releases/latest -O latest

###############################################################################
#
# STEP 2: Downloads
#
RPACK=$(grep -Po "(releases/download.*.quandico.*\.tar\.gz)" < latest)
echo "Downloading the R package: $RPACK"
wget -q -N --no-check-certificate https://github.com/reineckef/quandico/$RPACK

PERLM=$(grep -Po "(releases/download.*.QUANDICO.*\.tar\.gz)" < latest)
echo "Downloading the Perl module: $PERLM"
wget -q -N --no-check-certificate https://github.com/reineckef/quandico/$PERLM

echo "Downloading example data ..."
wget -q -N --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITQjBUZUFVcy1BNFk -O CNA902Y.bed
wget -q -N --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITUDZ0aXFBd2pDR0k -O M62_NA13019.tsv
wget -q -N --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITWU9xTndtZE5iOEE -O M62_NA12878.tsv
wget -q -N --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITYTduNDY4azJNZXM -O M62_NA13019.clustered
wget -q -N --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITWm1FS0duczVlejQ -O M62_NA12878.clustered

if [ $BIGDATA == 1 ]; then
  echo "Downloading BAM files (574 MB), please wait..."
#  wget -N --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITMzNyakhveTh3UVE -O M62_NA13019.bam
#  wget -N --no-check-certificate https://googledrive.com/host/0BzLnl09R3GITSnU1TlVRSjRXRHM -O M62_NA12878.bam
fi

WITHREF=""

if [ $REFGENE == 1 ]; then
  echo "Downloading refGene (from http://hgdownload.cse.ucsc.edu/goldenPath/hg19)..."
  wget -q -O- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz | gunzip > refGene.txt
  WITHREF="--cp names=refGene.txt"
fi

###############################################################################
#
# STEP 3: Install R dependencies
#
echo "Installing R dependencies..."
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
echo "Installing 'quandico' R package '$RPACK' ..."
$REXE CMD INSTALL $RPACK
echo "OK."

###############################################################################
#
# STEP 5: Install Perl dependencies
#
echo "Installing Perl dependencies (MakeMaker, cpanminus)..."
# cpan App::cpanminus ExtUtils::MakeMaker
echo "OK."

###############################################################################
#
# STEP 5: Install Perl Module
#
PERLM=$(perl -e "@L = split(/\//, '$PERLM'); print \$L[-1];")
echo "Installing Perl Module '$PERLM' ..."
cpanm -v $PERLM 
quandico -V
echo "OK."

echo ""
echo "The R package (with dependencies) and the Perl Module (with dependencies)"
echo "should now be installed (if there were no errors)."
echo ""
echo "Please consider installation of 'samtools' from: http://www.htslib.org"
echo ""

###############################################################################
#
# STEP 6: Run a test
#

if [ $TESTRUN != "1" ]; then
  echo "No test run requested, we are done."
  exit 0
fi

echo "Testing the installation on demo data $WITHREF ..."

if [ $BIGDATA == "1" ]; then
  echo "Running 'quandico' on BIG test data, this can take up to ten minutes to complete ..."
    time quandico -A CNA902Y.bed -X samtools=$SAMTOOLS --rexe $REXE -s map=M62_NA13019.bam -r map=M62_NA12878.bam -s x=2 -s y=0 -r x=2 -r y=0 $WITHREF
else
  echo "Running 'quandico' on extracted test data, this should take less than one minute ..."
  time quandico --rexe $REXE --no-cluster -s data=M62_NA13019.clustered -r data=M62_NA12878.clustered -s x=2 -s y=0 -r x=2 -r y=0 $WITHREF
fi

