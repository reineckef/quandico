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
1. Extract counts for every PCR amplicon (_qgetcounts_)
2. Cluster read-counts into regions (_qcluster_)
3. Compare counts per cluster and perform CNV calling (_quandico_)

Steps 1, 2 and 3 can be performed with one single command line call of the Perl 
driver script `quandico`. Step 3 alone can be performed inside R.

### Requirements
Running `quandico` requires `R` with some commonly available packages from [CRAN](http://cran.r-project.org). 
These should be installed automatically when installing the package (listed dependencies).

Accessory applications such as a command-line driver script (_quandico_), the script to extract counts of 
mapped reads (_qgetcounts_) and the region clustering script (_qcluster_) require Perl and some common modules 
from [CPAN](http://www.cpan.org). All dependencies should automatically be installed. 

### Installation
 1. Install R. Please visit [The `R`-project hompage](http://www.r-project.org) for advice.
 2. Download and install the packaged `R` code of `quandico` (version n.m)

`$ R CMD INSTALL quandico_[n.m].tar.gz`

**Optional**

 3. Install Perl (only if not already installed, most Linux systems ship with Perl, check [Perl.org](http://www.perl.org) on Windows we recommend [Strawberry Perl](http://www.strawberryperl.com).
 4. Download and install the Perl module `QUANDICO` that also contains the two helper script and the main driver (version n.m):
  4.1 Download the package
  4.2 Unzip and untar (will create its own folder)
  4.3 Install using make (Linux, Unix) or dmake (Strawberry Perl on Windows)

```
$ cd QUANDICO-[n.m]
$ perl Makefile.PL
$ [d]make
$ [d]make test
$ [d]make install
```

Alternatively, you can use `cpanm` from `App::cpanminus` to install to module
directly from the downloaded archive. We recommend to use `--verbose` mode:
	
`$ cpanm --verbose QUANDICO-[n.m].tar.gz`

 5. To start from mapped reads in SAM/BAM format, the Perl script qgetcounts will call `samtools` (version 1.1 or later) 
 to extract the counts. Therefore, `samtools` needs to be installed, please visit [Samtools](http://www.htslib.org) for advice.

### Running
The Perl command-line tools present a help and usage information when run without arguments (or by explicitly using `-h` 
or `--help`). The R function provides help (after installation) by doing:

```R
library(quandico)
?quandico
```

The directory 'Data' contains two BAM files and one BED file as demo data. You can 
run all steps separately:

```
# extract counts
$ qgetcounts -i M62_NA13019.bam -a CNA902Y.bed -o M62_NA13019.counts
$ qgetcounts -i M62_NA12878.bam -a CNA902Y.bed -o M62_NA12878.counts

# cluster the counts
$ qcluster -i M62_NA13019.counts [--names refGene.txt] > M62_NA13019.clustered
$ qcluster -i M62_NA12878.counts [--names refGene.txt] > M62_NA12878.clustered

# call copy numbers
$ quandico --no-cluster \
   --sample    data=M62_NA13019.clustered \   # path to file with clustered counts
   --reference data=M62_NA12878.clustered \   # path to file with clustered counts
   -s x=2 -s y=0 -r x=2 -r y=0 \              # sample (s) and reference (r) are female
```

Alternatively, all this can be done using one single command:

```
$ quandico -s map=M62_NA13019.bam -s x=2 -s y=0 \
           -r map=M62_NA12878.bam -r x=2 -r y=0 \
           --dir results --basename 13019_vs_12878 \
           [--cp names=refGene.txt]
```

### Help
All three scripts have several detailed options, please check the output of `<script> --help` for details.

### Issues
Please report bugs to the author [Frank Reinecke](mailto:frank.reinecke@qiagen.com) or file 
them on quandico's home [on GitHub](http://github.com/fr02081975/quandico).

