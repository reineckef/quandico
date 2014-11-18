# quandico

### Abstract

The tool `quandico` applies statistical methods to detect [copy number 
variations][cnv] (CNV) in a sample by comparing the 
read counts from next generation sequencing (NGS) performed after PCR-enrichment of regions of interest, typically a 
set of genes with known or expected relevance for the sample, e.g. genes that play a role in cancer. Counts from a normal 
control (ideally matched normal from the same individual, e.g. healthy tissue) are required. 

A detailed description of the method is submitted for publication.

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

The string *n.m* will be used for the version number (major.minor) of `quandico`. This was **1.12** by the time of writing.

 * Download [the latest release][latest] and install the packaged `R` code of `quandico` (file `quandico_n.m.tar.gz`):

`$ R CMD INSTALL path/to/quandico_n.m.tar.gz`

**_optional_ but _recommended_**

 * Download [the latest release][latest] and install the Perl module `QUANDICO` (file `QUANDICO-vn.m.tar.gz`) 
 using `make` on Linux, and `dmake` on Windows:

```bash
$ tar xvfz QUANDICO-vn.m.tar.gz
$ cd QUANDICO-vn.m
$ [d]make
$ [d]make test
$ [d]make install
```

Alternatively, you can use `cpanm` from `App::cpanminus` to install to module directly from the downloaded archive. 
We recommend to use `--verbose` mode:
	
`$ cpanm --verbose path/to/QUANDICO-vn.m.tar.gz`

To start from mapped reads in SAM/BAM format, the Perl script `qgetcounts` will call `samtools` (version 1.1 or later) to 
extract the counts. Therefore, `samtools` needs to be installed, please visit [Samtools][samtools] for advice.

## Example Data

Example data files to test the tools are available on Google drive. To start from mapped reads (in SAM or BAM files, please 
note that this step additionally requires **samtools** version >= 1.1), you need to download:

* [CNA902Y.bed](https://drive.google.com/open?id=0BzLnl09R3GITQjBUZUFVcy1BNFk&authuser=0) amplicon coordinates
* [M62_NA13019.bam](https://drive.google.com/open?id=0BzLnl09R3GITMzNyakhveTh3UVE&authuser=0) example case (sample) **:warning: 252 MB !**
* [M62_NA12878.bam](https://drive.google.com/open?id=0BzLnl09R3GITSnU1TlVRSjRXRHM&authuser=0) example control (reference) **:warning: 222 MB !**

To skip count extraction and start with clustering, the extracted counts of these samples are also available:

* [M62_NA13019.tsv](https://drive.google.com/open?id=0BzLnl09R3GITUDZ0aXFBd2pDR0k&authuser=0) control counts (sample)
* [M62_NA12878.tsv](https://drive.google.com/open?id=0BzLnl09R3GITWU9xTndtZE5iOEE&authuser=0) control counts (reference)

A file with gene names and coordinates is required if clusters should be named using gene names. For human assemblies GRCh37 
(hg19) and GRCh38 (hg38), these file can be downloaded from [UCSC](http://hgdownload.soe.ucsc.edu/downloads.html):

* [refGene.txt.gz](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz) assembly hg19 **(5 MB)**
* [refGene.txt.gz](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz) assembly hg38 **(5 MB)**

Clustered count files ready for processing with the final step are these:

* [M62_NA13019.clustered](https://drive.google.com/open?id=0BzLnl09R3GITYTduNDY4azJNZXM&authuser=0) 
* [M62_NA12878.clustered](https://drive.google.com/open?id=0BzLnl09R3GITWm1FS0duczVlejQ&authuser=0)

## Running

You can run all steps separately:

```bash
# extract counts
$ qgetcounts -i path/to/M62_NA13019.bam -a path/to/CNA902Y.bed -o path/to/M62_NA13019.tsv
$ qgetcounts -i path/to/M62_NA12878.bam -a path/to/CNA902Y.bed -o path/to/M62_NA12878.tsv

# cluster the counts
$ qcluster -i path/to/M62_NA13019.tsv [--names path/to/refGene.txt] > path/to/M62_NA13019.clustered
$ qcluster -i path/to/M62_NA12878.tsv [--names path/to/refGene.txt] > path/to/M62_NA12878.clustered

# call copy numbers
$ quandico --no-cluster \
   -s data=path/to/M62_NA13019.clustered \   # file with clustered counts
   -r data=path/to/M62_NA12878.clustered \   # file with clustered counts
   -s x=2 -s y=0 -r x=2 -r y=0           \   # sample (-s) and reference (-r) are female
   [--cp names=path/to/refGene.txt]          # optional for naming of clusters
```

Alternatively, all this can be done using one single command:

```bash
$ quandico -s map=path/to/M62_NA13019.bam -s x=2 -s y=0 \ # sample
           -r map=path/to/M62_NA12878.bam -r x=2 -r y=0 \ # reference
           -a /path/to/CNA902Y.bed                      \ # amplicons
           --dir results --basename 13019_vs_12878      \ # output location and name
           [--cp names=path/to/refGene.txt]               # optional cluster names
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
library(quandico)
?quandico
```

## Issues
Please report bugs to the author [Frank Reinecke][me] or file 
them on quandico's home [on GitHub](http://github.com/reineckef/quandico).


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

