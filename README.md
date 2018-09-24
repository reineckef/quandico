# quandico

[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs12859--014--0428--5-lightgrey.svg?style=flat-square)](https://doi.org/10.1186/s12859-014-0428-5)

### Abstract

The tool `quandico` applies statistical methods to detect [copy number 
variations][cnv] (CNV) in a sample by comparing the 
read counts from next generation sequencing (NGS) performed after PCR-enrichment of regions of interest, typically a 
set of genes with known or expected relevance for the sample, e.g. genes that play a role in cancer. Counts from a normal 
control (ideally matched normal from the same individual, e.g. healthy tissue) are required. 

Data are expected to be from gene-panels with primer-counts in the range of hundreds, and roughly double-digit number genes. It will not work well on whole-exome-enrichment data.

A detailed description can be found in the publication of this tool: **Quantitative analysis of differences in copy numbers using read depth obtained from PCR-enriched samples and controls**: [BMC Bioinformatics 2015, 16:17](http://www.biomedcentral.com/1471-2105/16/17 "BMC Bioinformatics"). 

## Implementation
The main script is written in `R` and there are two helper scripts and one main driver script written in Perl.

## Workflow
1. PCR-enrichment of target regions, sequencing and read-mapping (SAM/BAM file)
2. Extract counts for every PCR amplicon (`qgetcounts`)
3. Cluster read-counts into regions (`qcluster`)
4. Compare counts per cluster and perform CNV calling (`quandico`)

Steps 2, 3 and 4 can be performed with one single command line call of the Perl 
driver script `quandico`. Step 4 alone can be performed inside `R`.

## Requirements
Running `quandico` requires `R` with some commonly available packages from [CRAN][cran]. 
Please visit [The `R`-project homepage][r] for advice how to install `R` on your system.
Dependencies for `quandico` should be installed automatically when installing the package.

Accessory applications such as a command-line driver script (`quandico`), the script to extract counts of 
mapped reads (`qgetcounts`) and the region clustering script (`qcluster`) require Perl. Perl is already installed on 
most Linux systems, please check [Perl.org][perl] for details. If you are using Windows we 
recommend [Strawberry Perl][strawberry]. All dependencies should automatically be installed 
from [CPAN][cpan].

## Installation

### Fast

Download the [`install.sh`][installer] shell script for Unix/Linux. 
Check the configuration in the first section of that script, adjust 
parameters and then start it in a clean directory, e.g.:

```bash
$ mkdir quandico
$ cd quandico
$ wget https://github.com/reineckef/quandico/raw/master/install.sh
$ sh install.sh
```

This is not well tested and may fail on various systems for various 
reasons.

#### Hint for Ubuntu (maybe others)

If you are running Ubuntu linux, it is recommended to install Perl dependencies 
from repositories like so:

```bash
sudo apt-get install libenv-path-perl libgetopt-long-descriptive-perl libdatetime-perl libtest-script-perl
```
Thanks to @biocyberman #3 for this fix.

### Detailed

Install these `R` packages, either by using `R> install.packages('name')` inside an `R` session 
or by running the install command from a command line: `$ R CMD INSTALL 'name'`:

 * MASS
 * ggplot2
 * grid
 * gridExtra
 * naturalsort
 * scales

The string *n.m* will be used for the version number (major.minor) of `quandico`. This was **1.14** by the time of writing.

 * Download [the latest release][latest] and install the packaged `R` code of `quandico` (file `quandico_n.m.tar.gz`):

`$ R CMD INSTALL quandico_n.m.tar.gz`

**_optional_ but _recommended_**

 * Download [the latest release][latest] and install the Perl module `QUANDICO` (file `QUANDICO-vn.m.tar.gz`) 
 using `make` on Linux, and `dmake` on Windows:

```bash
$ tar xvfz QUANDICO-vn.m.tar.gz
$ cd QUANDICO-vn.m
$ perl Makefile.PL
$ [d]make
$ [d]make test
$ [d]make install # run this as sudo/root for system install
```

Alternatively, you can use `cpanm` from `App::cpanminus` to install to module directly from the downloaded archive. 
We recommend to use `--verbose` mode:
	
`$ cpanm --verbose QUANDICO-vn.m.tar.gz` # run this as sudo/root for system install

To start from mapped reads in SAM/BAM format, the Perl script `qgetcounts` will call `samtools` (version 1.1 or later) to 
extract the counts. Therefore, `samtools` needs to be installed, please visit [Samtools][samtools] for advice.

## Example Data

Example data files to test the tools are available on Google drive. To start from mapped reads (in SAM or BAM files, please 
note that this step additionally requires **samtools** version >= 1.1).

* [CNA902Y.bed](https://drive.google.com/open?id=0BzLnl09R3GITQjBUZUFVcy1BNFk&authuser=0) amplicon coordinates
* [M62_NA13019.bam](https://drive.google.com/open?id=0BzLnl09R3GITMzNyakhveTh3UVE&authuser=0) example case (sample) **:warning: 252 MB !**
* [M62_NA12878.bam](https://drive.google.com/open?id=0BzLnl09R3GITSnU1TlVRSjRXRHM&authuser=0) example control (reference) **:warning: 222 MB !**

Please note, these files can also be downloaded using wget (below). **Warning** Google fails to provide the file if no virus check can be done, which seems to be the case for the BAM files. You need to accept that on a per-click basis, so only via the web-browser, not `wget`:

```bash
$ wget --no-check-certificate "https://drive.google.com/uc?export=download&id=0BzLnl09R3GITQjBUZUFVcy1BNFk" -O CNA902Y.bed
$ # fails: wget --no-check-certificate "https://drive.google.com/uc?export=download&id=0BzLnl09R3GITMzNyakhveTh3UVE" -O M62_NA13019.bam
$ # fails: wget --no-check-certificate "https://drive.google.com/uc?export=download&id=0BzLnl09R3GITSnU1TlVRSjRXRHM" -O M62_NA12878.bam
```

To skip count extraction and start with clustering, the extracted counts of these samples are also available:

* [M62_NA13019.tsv](https://drive.google.com/open?id=0BzLnl09R3GITUDZ0aXFBd2pDR0k&authuser=0) control counts (sample)
* [M62_NA12878.tsv](https://drive.google.com/open?id=0BzLnl09R3GITWU9xTndtZE5iOEE&authuser=0) control counts (reference)

Please note, these files can also be downloaded using wget (sometimes these links do not work, for unclear reason):

```bash
$ wget --no-check-certificate "https://drive.google.com/uc?export=download&id=0BzLnl09R3GITUDZ0aXFBd2pDR0k" -O M62_NA13019.tsv
$ wget --no-check-certificate "https://drive.google.com/uc?export=download&id=0BzLnl09R3GITWU9xTndtZE5iOEE" -O M62_NA12878.tsv
```

A file with gene names and coordinates is required if clusters should be named using gene names. For human assemblies GRCh37 
(hg19) and GRCh38 (hg38), these file can be downloaded from [UCSC](http://hgdownload.soe.ucsc.edu/downloads.html):

* [refGene.txt.gz](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz) assembly hg19 **(5 MB)**
* [refGene.txt.gz](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz) assembly hg38 **(5 MB)**

Clustered count files ready for processing with the final step are these:

* [M62_NA13019.clustered](https://drive.google.com/open?id=0BzLnl09R3GITYTduNDY4azJNZXM&authuser=0) 
* [M62_NA12878.clustered](https://drive.google.com/open?id=0BzLnl09R3GITWm1FS0duczVlejQ&authuser=0)

Please note, these files can also be downloaded using wget (sometimes these links do not work, for unclear reason):

```bash
$ wget --no-check-certificate "https://drive.google.com/uc?export=download&id=0BzLnl09R3GITYTduNDY4azJNZXM" -O M62_NA13019.clustered
$ wget --no-check-certificate "https://drive.google.com/uc?export=download&id=0BzLnl09R3GITWm1FS0duczVlejQ" -O M62_NA12878.clustered
```

## Running

You can run all steps separately:

```bash
# extract counts
$ qgetcounts -i M62_NA13019.bam -a CNA902Y.bed > M62_NA13019.tsv
$ qgetcounts -i M62_NA12878.bam -a CNA902Y.bed > M62_NA12878.tsv

# cluster the counts
$ qcluster -i M62_NA13019.tsv [--names refGene.txt] > M62_NA13019.clustered
$ qcluster -i M62_NA12878.tsv [--names refGene.txt] > M62_NA12878.clustered

# call copy numbers
$ quandico --no-cluster \
   -s data=M62_NA13019.clustered \   # file with clustered counts
   -r data=M62_NA12878.clustered \   # file with clustered counts
   -s x=2 -s y=0 -r x=2 -r y=0   \   # sample (-s) and reference (-r) are female
   [--cp names=refGene.txt]          # optional for naming of clusters
```

Alternatively, all this can be done using one single command:

```bash
$ quandico -s map=M62_NA13019.bam -s x=2 -s y=0 \ # sample
           -r map=M62_NA12878.bam -r x=2 -r y=0 \ # reference
           -A CNA902Y.bed                       \ # amplicons
           -d results -b 13019_vs_12878         \ # output location and name
           [--cp names=refGene.txt]               # optional cluster names
```

## Example Output

In the folder 'ExampleOutput' the result files obtained by running the above mentioned steps with sample NA13019 and 
reference NA12878 from sequencing run M62 are presented. These are:

* [NA13019_overview.pdf](https://github.com/reineckef/quandico/raw/master/ExampleOutput/NA13019_overview.pdf) an overview for all clusters on one page
* [NA13019.pdf](https://github.com/reineckef/quandico/raw/master/ExampleOutput/NA13019.pdf) a detailed plot for every cluster
* [NA13019.vcf](https://github.com/reineckef/quandico/raw/master/ExampleOutput/NA13019.vcf) [VCF][vcf] of called variants
* [NA13019.csv](https://github.com/reineckef/quandico/raw/master/ExampleOutput/NA13019.csv) the called variants in tabular format

An image from one page (number 19) from the detailed PDF shows a loss (top) and a gain (bottom) event in that sample on 
chromsome X:

![Example page](/ExampleOutput/NA13019_p19.png)


## Help
The Perl command-line tools present a help and usage information when run without arguments (or by explicitly using `-h` 
or `--help`).

The `R` function provides help (after installation) by doing:

```R
R> library(quandico)
R> ?quandico
```

## Issues
Please report bugs by opening a new issued here [on GitHub](http://github.com/reineckef/quandico/issues).

## Author Contact
Fot other inquiries (no bug reports, please) you can email [me].

<!-- REFERENCES -->

[cnv]: http://en.wikipedia.org/wiki/Copy-number_variation "Copy Number Variation (Wikipedia)"
[cpan]: http://www.cpan.org "Comprehensive Perl Archive Network"
[cran]: http://cran.r-project.org "CRAN: Comprehensive R Archive Network"
[latest]: https://github.com/reineckef/quandico/releases/latest "the latest release"
[me]: mailto:frank.reinecke@qiagen.com
[perl]: http://www.perl.org "Perl Homepage"
[r]: http://www.r-project.org "R Project"
[samtools]: http://www.htslib.org "Samtools Homepage"
[strawberry]: http://www.strawberryperl.com "Strawberry Perl"
[vcf]: https://github.com/samtools/hts-specs "VCF Specifications"
[installer]: https://github.com/reineckef/quandico/raw/master/install.sh "install script"
