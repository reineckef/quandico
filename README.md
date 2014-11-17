### quandico

**Abstract**

The tool `quandico` applies statistical methods to detect copy number variations (CNV) in a sample by comparing the 
read counts from next generation sequencing (NGS) performed after PCR-enrichment of regions of interest, typically a 
set of genes with known or expected relevance for the sample, e.g. genes that play a role in cancer. Counts from a normal 
control (ideally matched normal from the same individual, e.g. healthy tissue) are required. 

### Implementation
The main script is written in `R` and there are two helper scripts and one main driver script written in Perl.

### Workflow
0. PCR-enrichment of target regions, sequencing and read-mapping (SAM/BAM file)
1. Extract counts for every PCR amplicon (`qgetcounts`)
2. Cluster read-counts into regions (`qcluster`)
3. Compare counts per cluster and perform CNV calling (`quandico`)

Steps 1, 2 and 3 can be performed with one single command line call of the Perl 
driver script `quandico`. Step 3 alone can be performed inside R.

### Requirements
Running `quandico` requires `R` with some commonly available packages from [CRAN](http://cran.r-project.org). 
Please visit [The `R`-project homepage](http://www.r-project.org) for advice how to install `R` on your system.
Dependencies for `quandico` should be installed automatically when installing the package.

Accessory applications such as a command-line driver script (`quandico`), the script to extract counts of 
mapped reads (_qgetcounts_) and the region clustering script (`qcluster`) require Perl. Perl is already installed on 
most Linux systems, please check [Perl.org](http://www.perl.org) for details. If you are using Windows we 
recommend [Strawberry Perl](http://www.strawberryperl.com). All dependencies should automatically be installed 
from [CPAN](http://www.cpan.org).

### Installation

The string *n.m* will be used for the version number (major.minor) of `quandico`. This was **1.12** by the time of writing.

 * Download and install the packaged `R` code of `quandico` (file `quandico_n.m.tar.gz`):

`$ R CMD INSTALL path/to/quandico_n.m.tar.gz`

**_optional_ but _recommended_**

 * Download and install the Perl module `QUANDICO` (file `QUANDICO-vn.m.tar.gz`) using `make` on Linux, and `dmake` on Windows:

```
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
extract the counts. Therefore, `samtools` needs to be installed, please visit [Samtools](http://www.htslib.org) for advice.

### Running
The directory 'Data' contains two BAM files and one BED file as demo data. 

 * Download demo data

You can run all steps separately:

```
# extract counts
$ qgetcounts -i path/to/M62_NA13019.bam -a path/to/CNA902Y.bed -o path/to/M62_NA13019.counts
$ qgetcounts -i path/to/M62_NA12878.bam -a path/to/CNA902Y.bed -o path/to/M62_NA12878.counts

# cluster the counts
$ qcluster -i path/to/M62_NA13019.counts [--names path/to/refGene.txt] > path/to/M62_NA13019.clustered
$ qcluster -i path/to/M62_NA12878.counts [--names path/to/refGene.txt] > path/to/M62_NA12878.clustered

# call copy numbers
$ quandico --no-cluster \
   -s data=path/to/M62_NA13019.clustered \   # file with clustered counts
   -r data=path/to/M62_NA12878.clustered \   # file with clustered counts
   -s x=2 -s y=0 -r x=2 -r y=0           \   # sample (-s) and reference (-r) are female
   [--cp names=path/to/refGene.txt]          # optional for naming of clusters
```

Alternatively, all this can be done using one single command:

```
$ quandico -s map=path/to/M62_NA13019.bam -s x=2 -s y=0 \
           -r map=path/to/M62_NA12878.bam -r x=2 -r y=0 \
           --dir results --basename 13019_vs_12878 \
           [--cp names=path/to/refGene.txt]
```

### Help
The Perl command-line tools present a help and usage information when run without arguments (or by explicitly using `-h` 
or `--help`).

The R function provides help (after installation) by doing:

```R
library(quandico)
?quandico
```

### Issues
Please report bugs to the author [Frank Reinecke](mailto:frank.reinecke@qiagen.com) or file 
them on quandico's home [on GitHub](http://github.com/fr02081975/quandico).

