\name{hg19}
\docType{data}
\alias{hg19}
\title{Lengths of human chromosomes in assembly version hg19 = GRCh37}
\description{
  This data set was created from a fasta index file 'hg19.fa.fai'.
  It contains these columns: chromosome, total_length, offset, line, bytes, 
  assembly. This kind of data is required for vcf output.
}
\usage{data(hg19)}
\format{A data frame with 25 rows}
\source{samtools -index hg19.fas}
\keyword{datasets}