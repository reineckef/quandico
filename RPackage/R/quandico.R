## // GLOBAL LIBRARY IMPORT //
##
library(grid)
library(ggplot2)
library(gridExtra)
library(MASS)
library(naturalsort)
library(scales)

## // GLOBAL CONFIGURATION SETTINGS //
##
config.default.ploidity    <- 2 # default to diploid
config.default.sample.x    <- 2 # default sample to be female
config.default.sample.y    <- 0 # default sample to be female
config.default.reference.x <- 2 # default reference to be female
config.default.reference.y <- 0 # default reference to be female

## // DATA READ FROM EXTERNAL FILES AND VERSIONING STUFF
##
## internal
config.software.name       <- 'quandico'
config.software.version    <- '1.13'  ##  GIT_VERSION ## config.software.version    <- '<%>'

## reference assembly
config.assembly.filename   <- NA
config.assembly.version    <- 'hg19'

## // SETTINGS AFFECTING NORMALIZATION
##
## absolute coverage required in either ref or sample, that means,
## count data where both reference and sample are below, will not be used
config.min_cov_for_scale     <- 20 ## absolute coverage for normalization
config.min_cov_for_call      <- 20 ## absolute coverage for calls

## // SETTINGS AFFECTING TREATMENT OF OUTLIERS AND NOISY DATA
##
## require this many amplicons per region (gene) to justify statistic testing,
## that means, at least this number of count-values must be present per 
## region (gene); this will be checked *after* removal of outliers

config.conf_level            <-  0.95
config.log2.max              <-  6.0 # max log2 (set Inf to this)
config.log2.min              <- -6.0 # min log2 (set NaN to this)

config.minimal_observations  <- 10         # (10)
config.use.weights           <- 'sum(raw)' # recommended: 'sum(raw)'
config.max.noise             <-  1

# outlier removal strategies
config.narrow.normal         <- TRUE # recommended: TRUE
config.narrow.test           <- 'shapiro' # recommended: 'shapiro'
config.narrow.pval           <- 0.05 # remove outlier untill p >= pval              (0.05)
config.narrow.lownoise       <- 1/3  # do not remove value if sd <= lownoise,        (1/3)
config.narrow.min            <- 2/3  # keep this fraction of all values in any case, (2/3)

config.flag.outliers         <- FALSE # recommended: FALSE!
config.recurse.outliers      <- FALSE # recommended: FALSE!
config.flag.extremes         <- FALSE # recommended: FALSE!

## // SETTINGS AFFECTING SENSITIVITY / SPECIFICITY => REPORTING
##
## report all copy numbers with score above the report threshold
config.score.report   <-  50 ## report above this score             (50)
config.score.check    <-  25 ## report uncertain above this score   (25)
config.score.max      <- 250 ## maximal score (hard clip)          (250)

## // SETTINGS AFFECTING THE SCORE CALCULATIONS
##
## post-process intial score from p.values (dispersion correction)
config.score.correction       <- 'sqrt(sd)*log2' #  recommended: 'sqrt(sd)*log2' = phi in the paper

## // SETTINGS AFFECTING AMOUNT OF REPORTING
##
## report output
config.report.layout         <- matrix(c(1,2,2,2,3,4,4,4), nrow = 2, byrow = TRUE)
config.report.pdf.paper      <- 'US'
config.report.pdf.colormodel <- 'cmyk'
##
## what to report
config.report.segmented       <- FALSE ## report segmented            (FALSE)
config.report.insufficient    <- TRUE  ## report insufficient          (TRUE)
config.report.noisy           <- TRUE  ## report noisy                 (TRUE)
config.report.all             <- TRUE  ## originally for dvlpmt phase  (TRUE)
config.report.genes           <- FALSE ## save a pdf for single genes (FALSE)
##
## dump data for debugging
config.write.dataset          <- FALSE  ## saves the complete data     (FALSE)

## // SETTINGS AFFECTING OPTIONAL REPORTS
##
## plots (on / off)
config.plot.rawdata      <- FALSE ## plot raw log2 counts              (FALSE)
config.plot.normalized   <- FALSE ## plot notrmalized log2 counts      (FALSE)
config.plot.overview     <- TRUE  ## plot summary per region/locus     (TRUE)
config.plot.chrom_boxes  <- FALSE ## generate boxplot per chromosome   (FALSE)
config.plot.chromosomes  <- FALSE ## generate scaled chromosome-grpahs (FALSE)

## // SETTINGS AFFECTING THE LAYOUTS OF PLOTS (also see file: themes.R)
##
## scales
config.plots.max.cn          <- 12
config.plots.ylimits         <- c(0, config.plots.max.cn + 0.9)
config.plots.ybreaks         <- c(0, 1, 2, 3, 4, 5, 6, 8, 10, 12)
config.plots.ylabels         <- c("0", "1", "2", "3", "4", "5", "6", "8", "10", "12")
config.plots.line1           <- 12.25
config.plots.line2           <- 11.50
config.plots.seg.lcolor      <- "orange"
config.plots.seg.ltype       <- "solid"
config.plots.alpha.violin    <-  0.50
config.plots.smooth.violin   <-  1.75
##
## single gene output
config.gene.layout           <- matrix(c(1,2,2,2), nrow = 1, byrow = TRUE)
config.gene.pdf.paper        <- 'USr'
##
## colors
config.color.noisy           <- "brown"
config.color.even            <- "green"
config.color.unchanged       <- "green"
config.color.insufficient    <- "gray50"
config.color.check           <- "gray50"
config.color.untested        <- "gray50"


###############################################################################
# pull_in_range
#
# input:      table, min, max, mincov
# output:     altered table 
#             - additional column outlier that has NA for out of range
#             - values in column log2 will be > min and < max
#
pull_in_range <- function(
    table, 
    min    = config.log2.min, 
    max    = config.log2.max,
    mincov = config.min_cov_for_call
  ) {
    
  x <- table$log2
  y <- x
  y[is.na(x)] <- NA
  y[is.infinite(x)] <- NA
  y[x > max] <- NA
  y[x < min] <- NA
    
  x[is.na(x)] <- min
  x[is.infinite(x)] <- max
  x[x > max] <- max
  x[x < min] <- min
  
  table$log2    <- x
  table$outlier <- y
  table
}

###############################################################################
# flag_extremes
#
# input:      x         : raw data
# output:     array of values with NA instead of outliers
#
flag_extremes <- function (x, against = 'median', target = 1, min = 10) {
  sd     <- sd(x)
  center <- ifelse(against == 'median', median(x), mean(x))
  lower <- min(x)
  upper <- max(x)
  
  while (!is.na(sd) & sd > target ) {
    # check current low and high values (non-outliers)
    low  <- min(x[x > lower])
    high <- max(x[x < upper])
    if (abs(center - low) > abs(center - high)) {
      # lowest is farther away from center than highest
      lower <- low
      
    } else {
      upper <- high
    }
    y <- x[x > lower & x < upper]
    center <- ifelse(against == 'median', median(y), mean(y))
    sd     <- sd(y)
  }
  outliers  <- rep(0, length(x))
  outliers[x < lower | x > upper] <- 1
  outliers
}

###############################################################################
# flag_outliers
#
# input:      x         : raw data
# output:     array of values with NA instead of outliers
#
flag_outliers <- function(x, na.rm = TRUE, maxnoise = 1, range = 1.5, loop = FALSE, outliers = NA) {
  qnt <- quantile(x, probs=c(.25, 0.5, .75), na.rm = na.rm)
  H <- min( range * IQR(x, na.rm = na.rm)) ## hack to force removal of extreme in noisy regions

  if (is.na(outliers)[1]) {
    outliers  <- rep(0, length(x))
  }
  seen <- sum(outliers)
  
  outliers[x < (qnt[1] - H)] <- 1         # lower than range x lower quartile
  outliers[x > (qnt[3] + H)] <- 1         # higher than range x upper quartile
  outliers[abs(x - qnt[2])>maxnoise] <- 1 # bigger distance from median than maxnoise
  
  if (loop) {
    found <- sum(outliers)
    
    if (found == 0 | found == seen) return (outliers)
    # found (additional) outliers
    x[outliers == 1] <- NA ## remove values for all flagged outliers
    return(flag_outliers(x, na.rm = TRUE, range = range, loop = TRUE, outliers = outliers))
  } else {
    outliers
  }
}

###############################################################################
# subset_cluster
#
# input:      table      : full table with amplicon data
#             gene       : name of the gene to subset (column gene)
#
# output:     subset of table
#
subset_cluster <- function(table, cluster) {
  # only use data from one cluster 
  genedata <- subset(table, gene == cluster )
  
  if (config.narrow.normal) {
    genedata$outlier <- narrow_until_normal(
        genedata$log2, 
        w         = genedata$weight, 
        test      = config.narrow.test,
        min.pval  = config.narrow.pval,
        min.left  = config.narrow.min,
        low.noise = config.narrow.lownoise
    )
  }
  else if (config.flag.outliers) {
    genedata$outlier <- flag_outliers(genedata$log2, loop = config.recurse.outliers, maxnoise = config.max.noise)
  }
  else if (config.flag.extremes) {
    genedata$outlier <- flag_extremes(genedata$log2, target = config.max.noise)
  }
  
  # additionally, flag those that fail minmal coverage in the reference
  genedata$outlier[genedata$reference < config.min_cov_for_call] <- 1
  genedata
}

###############################################################################
# test_for_change
#
# input:      table       : full table with amplicon data
#             alternative : test against this
#             call        : data frame to insert results
#
# output:     call with changed values
#
test_for_change <- function(
      table, call, 
      mu          =   0, 
      conf.level  =   0.95,
      checkscore  =   5,
      checkcolor  = 'gray50',
      minscore    =  10,
      maxscore    = 250,
      correction  =  NA
  ) {
  
  ## requirement: sufficient amount of data
  if (!call$sufficient) {
    ## defaults are already in place for this, just set report (?) and return
    if (config.report.insufficient) call$report <- TRUE
    return(call)
  }

  ## requirements are met: run the test
  wtest <- weighted_t_test( 
    table$log2,
    w          = table$weight,
    mu         = mu,
    conf.level = conf.level
  )
  
  ## insert the results into the call data frame
  call$p.val   <- wtest$p.value
  call$cn.pval <- scientific_10( wtest$p.value )
  call$cn      <- call$refcopy * 2 ^ wtest$estimate
  call$min     <- call$refcopy * 2 ^ wtest$conf.int[1]
  call$max     <- call$refcopy * 2 ^ wtest$conf.int[2]
  call$se      <- wtest$se
  call$sd      <- wtest$sd
  call$qp      <- -10 * log( wtest$se ) # * (call$usable / call$amplicons) 
  call$log2    <- wtest$estimate
  
  ## check better alternatives for this multiplier, e.g. min confident distance
  multiplier <- 1
   if (!is.na(correction)) {
     multiplier <- sqrt(call$sd) * (1 + abs( wtest$estimate) ) # \phi in the paper
  }
  # apply scoring
  call$score   <- min(c(maxscore, -10 * log10(wtest$p.value) * multiplier)) #
  
  if (is.nan(call$score)) {
    call$score  <- 0
    call$usable <- 0
  }
  
  if (call$usable < config.minimal_observations) {
    call$filter <-  paste('o', config.minimal_observations, sep="")
  }
  ## score *must* be relevant
  else if ( call$score >= config.score.report ) {
      call$cn.label <- sprintf("%.1f ~ %.1f", call$min, call$max) 
      call$cn.color <- calc_color(call$cn, call$refcopy)
      call$report   <- TRUE
      call$filter   <- 'PASS'
  } else if ( call$score >= checkscore   & 
      ( call$min   > call$refcopy | call$max < call$refcopy ) ) {
      call$cn.label <- sprintf("%.1f ~ %.1f", call$min, call$max)
      call$cn.color <- checkcolor
      call$filter   <- paste('q', config.score.report, sep="")
      call$report   <- TRUE
  }
  return(call)
}

###############################################################################
# set_filter
#
# Test if provided score is higher than pass or check and return a tag
# to put into the filter column of vcf
# 
# input:      score     : the score to check
#             pass      : threshold to set PASS
#             check     : threshold to set q<pass> (check < score < pass)
#             fail      : what to return for score < check (NA = q<check>)
#
# output:     filter    : example: PASS; q40, q20, FAIL
#
set_filter <- function (
    score, 
    obs,
    pass  = config.score.report, 
    check = config.score.check, 
    mino  = config.minimal_observations,
    fail  = NA
    ) {
  
    if (obs < mino ) {
      return(paste('o', mino, sep=""))
    }
    if (score >= pass) {
      return('PASS')
    } else if (score >= check) {
      return(paste('q', pass, sep=""))
    } else {
      if (is.na(fail)) {
        return(paste('q', check, sep=""))
      }
      else {
        return(fail)
      }
    }
}

###############################################################################
# calc_color
#
# This function will provide a color for a combination of seen and expected
# copy numbers. Unchanges (similar) values will be green, and lower seen
# values will red, higher values blue to purple.
#
# input:      copy      : number of copies found
#             refcopy   : number of copies in the reference
#             colours   : array of colour to use
#             values    : array of values for colours
#
# output:     color

calc_color <- function(
    copy, 
    refcopy = 2,
    colours = c("brown", "red", "green","blue", "purple", "black"),
    values  = c(-2,-1,0,log2(3/2),1,2)
) {
  minimal <- min(values)
  maximal <- max(values)
  ratio <- log2(copy/refcopy)
  ratio[ratio < minimal] <- minimal
  ratio[ratio > maximal] <- maximal
  
  gradient  <- gradient_n_pal(colours=colours, values=values)
  return( gradient(ratio))  
}


###############################################################################
# scientific_10
#
# input:      number              # example: 0.00021
# output:     converted format    # example: 2.1 x 10^-4
#
scientific_10 <- function(x) {
	return(parse(text=gsub("e", " %*% 10^", paste(sprintf("%3.2g", x)))))
}


###############################################################################
# vcf_header
#
# input:      assembly          : assmebly name (id), examples: hg19, GRCh38
#             filename          : name (path) of the file to read data from
#             date              : date to put in VCF in YYYYMMDD (20140331)
#             version           : version of VCF to put in header (4.1)
#
# Function to calculate the expected log2 for a given combination of inputs:
#
vcf_header <- function(
    assembly = config.assembly.version,
    filename = config.assembly.filename,
    date     = format(Sys.Date(), format="%Y%m%d"),
    version  = '4.1'
) {
  
  if (!is.na(filename)) {
    current <- read_chromosomes(assembly, filename)
  }
  else {
    data(hg19)
    current <- hg19 
    assembly <- 'hg19'
  }

  lines   <- c(
      paste("##fileformat=VCFv", version, sep=""), 
      paste("##fileDate=", date, sep=""),
      paste("##source=", config.software.name, 'V', config.software.version, sep=""),
      paste("##reference=", assembly, ".fas", sep="")  # this is a fake name!  
  )
  
  shift.by <- length(lines)
  for (row in 1:nrow(current)) {
    lines[row + shift.by] <- sprintf("##contig=<ID=%s,assembly=%s,length=%s>", 
        current$chromosome[row], 
        current$assembly[row], 
        current$total_length[row]
    )
  }  
  return(paste(lines, sep="", collapse="\n"))
}

###############################################################################
# read_chromsomes
#
# Function to read chromosome data from a file
#
# input:      assembly          : assmebly name (id), examples: hg19, GRCh38
#             filename          : name (path) of the file to read data from
#
#             * the input file can be either a custom tabular format 
#               (there is an example shipped with this code) or a .fai
#               index file generated by samtools. The type is choosen based 
#               on the presence of '.tsv' or '.fai' in the filename!
#
# output:     Data table with headers defined by the first line in file.
#
read_chromosomes <- function(
    assembly = config.assembly.version, 
    filename = config.assembly.filename
) {  
  if (is.na(filename)) {
    if (assembly == "hg19") {
      data(hg19)
      return(hg19)
    }
    else {
      print(warning("No genome index specified, and assembly not included"))
      return()
    }
  }
  
  if(length(grep(".tsv", filename))) {
  assemblies <- read.table(filename, header=TRUE)
  
  if (length(grep("^hg", assembly))) {
    names(assemblies)[4] <- 'assembly'
  } else {
    names(assemblies)[5] <- 'assembly'
  }
  
  assemblies <- assemblies[assemblies$assembly == assembly, ]
  assemblies <- assemblies
  return(assemblies)
  }  else if(length(grep(".fai", filename))) {
    assemblies <- read.table(
          filename, 
          header=FALSE, 
          sep="\t", 
          col.names = c("chromosome", "total_length", "offset", "line", "bytes")
    )
    assemblies$assembly <- rep(assembly, nrow(assemblies))
    assemblies <- assemblies
    return(assemblies)
  }
  # croak if we get here!
  warning("No assembly table generated - input data missing or not recognized!")
}


###############################################################################
# plot_overview
#
# Function to plot calls from all regions sorted by location
#
# input:      data             : the table that is also exported to csv
#             filename         : name (path) of the file to wrtie to
# output:     This will write the PDF plot to the filename

plot_overview <- function(csv, filename=NA, title=NA) {
  
  data <- read.csv(csv)
  data <- subset(data, !is.na(p.val)) # subset to keep only valid regions
  # replace 0 scores with check score for plotting (gray, not green)
  data$score[is.na(data$min)] <- config.score.check 
  
  yrange <- range(data$min, data$max, 0, 4, na.rm=TRUE)
  
  if (nrow(data)> 0) {
    
    # collect brakes, labels and lines for chromosome axis
    breaks<-c()
    labels<-c()
    lines<-c(0) # start with left line
    
    last.pos <- 0
    temp <- data.frame()
    last.row <- 0
    
    # get all chromosomes in natural order
    chromosomes <- naturalsort(unique(data$chromosome), decreasing=FALSE)
    
    # loop through all chromosomes
    for(chr in chromosomes) {    
      sub <- subset(data, chromosome==chr) # subset to keep only current chr
      
      last.row <- nrow(sub)
      if(last.row==0) next
      
      for (target in sub$locus[order(sub$start)]) {
        region <- subset(sub, locus==target)
        region <- region
        last.pos <- last.pos +1 # move on axis, first wil be 1
        region$position <- last.pos
        temp<-rbind(temp, region)
      }
      # add one more for the ending line on the right
      breaks <- c(breaks, last.pos - ((last.row-1)/2))
      last.pos <- last.pos +1
      labels <- c(labels, chr)
      lines <- c(lines, last.pos)
    }
    
    p <- ggplot(temp, aes(x=position, y=copies, colour=score, ymin=min, ymax=max)) + theme_overview()
    p <- p + geom_vline(xintercept = lines, color="white", size=0.5)
    p <- p + geom_pointrange(size=0.5)
    p <- p + scale_y_continuous('copies', limits=yrange)
    p <- p + scale_x_continuous('chromosome', breaks=breaks, labels=labels, expand=c(0.01,0.01), limits=c(0,last.pos))
    p <- p + scale_colour_gradientn(
      "score", 
      values = c(
        0,                                           # min
        0.5 * config.score.check / config.score.max, # 0.5 x check score
        config.score.check       / config.score.max, # 1.0 x check score
        config.score.report      / config.score.max, # 1.0 x call score
        config.score.report * 2  / config.score.max, # 2.0 x call score
        1                                            # max
      ), 
      colours = c(
        "green",  # min
        "green",  # 0.5 x check
        "gray",   # 1.0 x check
        "orange", # 1.0 x call
        "red",    # max
        "brown"
      ),
      limits = c(0,config.score.max),
      guide  = guide_colorbar(label.theme = element_text(size = 7, angle = 90))
    )
    
    if(!is.na(title)) p <- p + labs(title = title)
    
    if (!is.na(filename)) {
      pdf(filename, paper="US", colormodel = config.report.pdf.colormodel)
      print(p)
      dev.off()
    }
    else {
      print(p)
    }
  }
  else {
    if(!is.na(filename)) {
    pdf(filename, paper="US", colormodel = config.report.pdf.colormodel )
    plot(1, type="n", axes=F, xlab="", ylab="")
    mtext("no data for overview plot")
    dev.off()    
    }
    else {
      print("no data for overview plot")
    }
  }
  #
  ## finished
}
###############################################################################
# import_counts
#
# This function reads two tables with locations and counts, merged these into 
# a data frame and also sets the reference copy number based on the provided 
# copy numbers for X and Y in reference and sample
#
# input:      sample      = filename with sample data
#             sample.x    = expected counts for X chromosome in sample
#             sample.y    = expected counts for Y chromosome in sample
#             reference   = filename with reference data
#             reference.x = real counts (may be float) of X in reference
#             reference.y = real counts (may be float) of Y in reference
#
# output:     data.frame
#
import_counts <- function(
      sample      = NA, 
      reference   = NA, 
      sample.x    = config.default.sample.x, # default to female sample
      sample.y    = config.default.sample.y, # default to female sample
      reference.x = config.default.reference.x, # default to female reference
      reference.y = config.default.reference.y,  # default to female reference
      sample.columns    = c("chr_full", "position", "direction", "base", "gene", "sample", "orgname"),
      reference.columns = c("chr_full", "position", "direction", "base", "gene", "reference", "orgname")
  ) {
  # check if filenames are passed  
  if(is.na(sample) | is.na(reference)) {
    print(warning("Sample or reference input file not provided, using example data"))
    data(M62_NA13019)
    data(M62_NA12878)
    
    smp <- M62_NA13019 # included example data
    ref <- M62_NA12878 # included example data
  }
  else {
    # read both files
    # data for sample:
    smp <- read.table(sample, header = FALSE, col.names = sample.columns)    
    # data for reference:
    ref <- read.table(reference, header = FALSE, col.names = reference.columns) 
  }

  # merge
  depth <- merge(smp, ref, by = intersect(names(smp), names(ref)))
  
  # replace chr to make a short version of the chromosome
  depth$chr <- gsub("^chr", "", depth$chr_full)
  
  # preset refcopy number to ploidity default
  depth$refcopy <- config.default.ploidity
  
  # treat the sex chromosomes differently, depending on input. Here we 
  # will use the counts of the *sample* (although is says reference copies)
  # this is confusing, but this will be used as "expected" copy number for
  # reporting, and we expect X=1 and Y=1 for a male, even if the reference
  # is from a female!
  depth$refcopy[depth$chr == 'X'] <- sample.x # expected X copies in sample
  depth$refcopy[depth$chr == 'Y'] <- sample.y # expected Y copies in sample
  
  # calculate correction factors for x and y (for non-matched or ref-pools)
  correct.x <- ifelse(reference.x > 0, sample.x / reference.x, 0)
  correct.y <- ifelse(reference.y > 0, sample.y / reference.y, 0)
  
  # correct counts by reducing the counts that are higher than needed
  # (compared to a a truely matched reference). this is to avoid scaling up, 
  # which may add artificial precision (higher count) than it actually has
  if (correct.x < 1) {
    # X in reference is higher than in sample, reduce reference to match
    depth$reference[depth$chr == 'X'] <- depth$reference[depth$chr == 'X'] * correct.x
  }  else if (correct.x > 1) {
    # X in sample is higher than in reference, reduce sample to match
    depth$sample[depth$chr == 'X'] <- depth$sample[depth$chr == 'X'] / correct.x
  }
  if (correct.y < 1) {
    # Y in reference is higher than in sample, reduce reference to match
    depth$reference[depth$chr == 'Y'] <- depth$reference[depth$chr == 'Y'] * correct.y
  } else if (correct.y > 1) {
    # Y in sample is higher than in reference, reduce sample to match
    depth$sample[depth$chr == 'Y'] <- depth$sample[depth$chr == 'Y'] / correct.x
  }  
  # return the data.frame
  depth <- depth
  return(depth)
}

###############################################################################
# Multiple plot function                                               // START
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., 
    plotlist	=	NULL, 
    file			= NULL, 
    width     = 800, 
    height    = 600, # rendered images output
    paper     = 'special', # default for pdf output
    cols      = 1, 
    layout    = NULL
) {
	
	# Make a list from the ... arguments and plotlist
	plots    <- c(list(...), plotlist)	
	numPlots <- length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
  if (!is.null(file)) {
    pdf(file, paper=paper)
  }
  
	if (numPlots==1) {
		print(plots[[1]])		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
		}
	}
  if (!is.null(file)) {
    dev.off()
  }
}
# Multiple plot function                                                 // END
###############################################################################

###############################################################################
##
## theme_gene_left
theme_gene_left <- function(base_size = 12, base_family = "Helvetica") {
  theme(
      # Elements in this first block aren't used directly, but are inherited
      # by others
      line = element_line(colour = "black", size = 0.5, linetype = 1, lineend = "butt"),
      rect = element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
      text = element_text(
          family  = base_family, face       = "plain", 
          colour  = "black",     size       = base_size, 
          hjust   = 0.5,         vjust      = 0.5, 
          angle   = 0,           lineheight = 0.9
      ),
      axis.text  = element_text(size = rel(0.8), colour = "grey50"),
      strip.text = element_text(size = rel(0.8)),
      
      axis.line = element_blank(),
      axis.text.x = element_text(vjust = 1, size = rel(0.8)),
      axis.text.y = element_text(hjust = 1),
      axis.ticks = element_line(colour = "grey50"),
      axis.title.x = element_text(face="italic", size = rel(0.8)),
      axis.title.y = element_text(angle = 90),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks.margin = unit(0.1, "cm"),
      
      legend.background = element_rect(colour = NA),
      legend.margin = unit(1.0, "mm"),
      legend.key = element_rect(fill = "grey95", colour = "white"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = rel(0.8)),
      legend.text.align = NULL,
      legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0),
      legend.title.align = NULL,
      legend.position = "right",
      legend.direction = NULL,
      legend.justification = "center",
      legend.box = NULL,
      
      panel.background = element_rect(fill = "grey90", colour = NA),
      panel.border = element_blank(),
      panel.grid.major = element_line(colour = "white"),
      panel.grid.minor = element_line(colour = "grey95", size = 0.25),
      panel.margin = unit(0.25, "lines"),
      
      strip.background = element_rect(fill = "grey80", colour = NA),
      strip.text.x = element_text(),
      strip.text.y = element_text(angle = -90),
      
      plot.background = element_rect(colour = "white"),
      plot.title = element_text(size = rel(1.2)),
      #  	margin around entire plot (unit with the sizes of the top, right, bottom, and left margins) 
      plot.margin = unit(c(0, 2, 2, 2), "mm"),
      
      complete = TRUE
  )
}
theme_gene_left <- theme_gene_left

###############################################################################
## 
## theme_gene_right - same as left without axis.y ticks and labels

theme_gene_right <- function(base_size=12, base_family="Helvetica") {
  theme_gene_left(base_size = base_size, base_family = base_family) %+replace%
    theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(face="plain", size = rel(0.8)),
        #  	margin around entire plot (unit with the sizes of the top, right, bottom, and left margins) 
        plot.margin   = unit(c(0, 25, 2, 0), "mm"),
        ## legend.margin = unit(1.0, "mm")
        legend.position = c(1.1, 0.5),
        legend.box.just = "left"
      )
}
theme_gene_right <- theme_gene_right

###############################################################################
##
## theme_chromosome

theme_chromosome <- function(base_size = 12, base_family = "Helvetica") {
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
      theme(
          axis.text        = element_text(size = rel(0.8)),
          axis.ticks       = element_line(colour = "black"),
          legend.position  = "none",
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border     = element_rect(fill = NA, colour = "grey50"),
          panel.grid.major = element_line(colour = "grey90", size = 0.2),
          panel.grid.minor = element_line(colour = "grey98", size = 0.5),
          strip.background = element_rect(fill = "grey80", colour = "grey50"),
          strip.background = element_rect(fill = "grey80", colour = "grey50")
      )
}
theme_chromosome <- theme_chromosome

###############################################################################
##
## theme_overview

theme_overview <- function(base_size = 12,base_family = "Helvetica") {
  theme_gene_left(base_size = base_size, base_family = base_family) %+replace%
      theme(
        axis.text.x        = element_text(angle=90, hjust = 0.5, vjust = 0.5),
        axis.title.x       = element_text(face="bold", angle=0, size = rel(0.8)),
        axis.title.y       = element_text(face="bold", angle=90, size = rel(0.8)),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), # line(colour="white", size=0.5), 
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.y = element_blank(), 
        legend.position    = "bottom"
  )
}
theme_overview <- theme_overview

###############################################################################
# narrow_until_normal
#
# input:      x         : raw data
#             weights   : list of weights (one for each x)
#
# output:     list of flags with 1 (out of range) or 0 (in range)
#

#' narrow_until_normal
#' 
#' This function uses a test for standard-distribution (defaults to shapiro)
#' and 'masks' outliers (largest deviation from weighted mean) until the 
#' remaining data passes the specified test (follows normal distribution).
#' An array of 0 (usable) and 1 (flagged) will be returned.
#' 
#' @param x array of input values or measurements
#' @param w array of weights
#' @param min.pval threshold for test's p-value
#' @param min.left fraction of x to keep (between 0 and 1; defaults to 0)
#' @param low.noise stop flagging outliers if sd falls below this
#' @param test name of the test to perform ('shapiro' or 'ks')
#' @export 

narrow_until_normal <- function (x, 
    w         = NA, 
    min.pval  = 0.05, 
    min.left  = 0,
    low.noise = 1/3,
    test      = 'shapiro'
) {
  
  # preset to fail the criterion (want at least 3 values)
  exit <- FALSE
  if (length(x[!is.na(x)]) < 3) {
    exit <- TRUE
  }
  if (max(x, na.rm=TRUE) == min(x, na.rm=TRUE)) {
    exit <- TRUE
  }
  # use standard weights if missing
  if (is.na(w[1])) {
    w <- rep(1/length(x), length(x))
  }
  
  # limits
  x.lower <- min(x, na.rm = TRUE)
  x.upper <- max(x, na.rm = TRUE)
  
  # loop until we have a normal subset - or need to exit due to other reasons
  while (!exit) {
    
    # determine the center of the distribution
    x.center <- weighted.mean(x, w, na.rm=TRUE)
    
    if (test == 'shapiro') {
      sw <- shapiro.test(x)
      if (sw$p.value > min.pval) exit <- TRUE
    } else if (test == 'ks') {
      ks.result <- ks.test(x, "pnorm", mean=x.center, sd=x.sd)
      if (ks.result$p.value > min.pval) exit <- TRUE
    } else {
      print(warning("Test is not supported (none of 'shapiro' or 'ks'), will return input data as is"))
      exit <- TRUE
    }
    # check for stop criterion
    if (!exit) {
      low  <- min(x[x >= x.lower], na.rm = TRUE)
      high <- max(x[x <= x.upper], na.rm = TRUE)
      
      if ( abs(x.center - low) > abs(x.center - high) ) {
        x[x <= low] <- NA
        x.lower <- low
      } else {
        x[x >= high] <- NA
        x.upper <- high
      } 
      if (length(x[!is.na(x)]) <= length(x) * min.left) {
        exit <- TRUE
      }
      if (length(x[!is.na(x)]) < 3) {
        exit <- TRUE
      }
      if (!exit) {        
        xsd<-sd(x, na.rm=TRUE)
        if ( is.na(xsd) | xsd < low.noise) {
          exit <- TRUE
        }
      }
    }
  } # end while !exit
  
  outliers <- rep(0, length(x))
  outliers[is.na(x)] <- 1
  return(outliers)
}

#######################################################################################################################
## var_wt() - weighted variance
##
## As all of the x_i are drawn from the same distribution and the integer weights w_i indicate the number of occurrences
## ("repeat") of the observations (ratios of sample/control), the unbiased estimator of the weighted population variance
## is given by:
##
##  s^2 = \frac {\sum_{i=1}^N w_i \left(x_i - \mu^*\right)^2 } {\sum_{i=1}^n w_i - 1}
##

#' var_wt
#' 
#' This function will return the weighted variance of an array of values (x)
#' with associated weights (w). If the weighted mean (mu) is already known, 
#' it can be passed to avoid calculating it again via weighted.mean().
#' Per default, undefined values and weights are silently removed.
#' 
#' @param x array of input values or measurements
#' @param w array of weights
#' @param mu optional: weighted mean of values (x) with weights (w) 
#' @param na.rm remove undefined values and weights (TRUE)
#' @export

var_wt <- function(x, w, mu=NULL, na.rm = TRUE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  if (is.null(mu)) {
    mu <- weighted.mean(x, w, na.rm=na.rm)
  }  
  return( (sum(w*(x-mu)^2)) / (sum(w)-1) )
}

#######################################################################################################################
## weighted_t_test(x, w, mu, conf.level, alternative, na.rm) - weighted t-test 
##
## This will use weighted mean (x.w) and weighted variance (var.w) to 
## calculate a regular t-test.
##

#' weighted_t_test
#' 
#' This function performes a weighted t-test of provided values (x)
#' and weights (w) against the reference (mu). Per default, all 
#' undefined values in x are removed and those with undefined weight are 
#' also ignored. The effective sample size is used instead of the length
#' of the vector (x) to compensate for largely different weights.
#' 
#' @param x array of input values or measurements
#' @param w array of weights
#' @param mu mean to test the distribution of x against (0)
#' @param conf.level the confidence-level (0.95)
#' @param alternative alternative to use for testing (two.sided)
#' @param na.rm remove undefined values and weights (TRUE)
#' @param effective use the effective sample size (TRUE)
#' @export

weighted_t_test <- function(
      x,   w,    mu, 
      conf.level  = 0.95, 
      alternative = "two.sided", 
      na.rm       = TRUE,
      effective   = TRUE
  )  {
  
  if(!missing(conf.level) &
      ( length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  
# to achieve consistent behavior in loops, return NA-structure in case of complete missings
  if (sum(is.na(x)) == length(x)) return(list(estimate=NA, se=NA, sd=NA, conf.int=NA, statistic=NA, df=NA, p.value=NA))
  
# if only one value is present: this is the best estimate, no significance test provided
  if (sum(!is.na(x)) == 1) {
    warning("Warning weighted_t_test: only one value provided; this value is returned without test of significance!", call.=FALSE)
    return(list(estimate=x[which(!is.na(x))], se=NA, sd=NA, conf.int=NA, statistic=NA, df=NA, p.value=NA))
  }
  
  x.w   <- weighted.mean(x,w, na.rm=na.rm)
  var.w <- var_wt(x,w, x.w, na.rm=na.rm) # new version of var.w requires x.w
  n     <- length(x)

  ess   <- ifelse(effective, (sum(w)^2) / sum(w^2), n)  # effective sample size
  t     <- sqrt(ess)*((x.w-mu)/sqrt(var.w))
  se    <- sqrt(var.w)/sqrt(ess)
  
  if (alternative == "less") {
    pval <- pt(t, ess)
    cint <- c(-Inf, x.w + se*qt(conf.level, ess) )
  }
  else if (alternative == "greater") {
    pval <- pt(t, ess, lower.tail = FALSE)
    cint <- c(x.w - se * qt(conf.level, ess), Inf)
  }
  else {
    pval  <- 2 * pt(-abs(t), ess)
    alpha <- 1 - conf.level
    cint  <- x.w + se * qt(1 - alpha/2, ess)*c(-1,1)
  }
  
  names(t) <- "t"
  return(
    list(
      estimate  = x.w, 
      se        = se, 
      conf.int  = cint, 
      statistic = t, 
      df        = ess, 
      p.value   = pval, 
      sd        = sqrt(var.w)
    )
  )
}

#' quandico
#'
#' This is the main function of the package which will process read-count data 
#' and perform copy-number estimation. Results will be written to four output 
#' files:
#' \itemize{
#' \item{1. <basename>.csv}{comma separated values (CSV) and header}
#' \item{2. <basename>.vcf}{VCF file with header}
#' \item{3. <basename>.pdf}{PDF file with a graphical overview for every cluster}
#' \item{4. <basename>_overview.pdf}{summary plot with one point-range per cluster}
#' }
#' @param sample the full path to the count file containing sample data (NA)
#' @param sample.x the number of X chromosomes in the sample (2)
#' @param sample.y the number of Y chromosomes in the sample (0)
#' @param sample.columns column definition for the sample input file (NA)
#' ("chr_full", "position", "direction", "base", "gene", "sample", "orgname")
#' @param reference the full path to the count file containing reference data
#' @param reference.x the number of X chromosomes in the reference (2)
#' @param reference.y the number of Y chromosomes in the reference (0)
#' @param reference.columns column definition for the reference input file
#' ("chr_full", "position", "direction", "base", "gene", "reference", "orgname")
#' @param output.dir the full path where to write results to
#' @param output.basename prefix for output files
#' @param assembly.filename file containing chromosome details (NA = hg19.fa.fai)
#' @param assembly.version name of the underlying genome build (NA = hg19)
#' @export
#' @import MASS grid ggplot2 gridExtra naturalsort scales
#' @examples
#' \dontrun{
#' quandico(
#' 	sample            = "/path/to/sample_counts_file.tab", 
#' 	reference         = "/path/to/reference_counts_file.tab",
#' 	assembly.filename = "/path/to/fasta_index/hg19.fa.fai",
#' 	assembly.version  = "hg19",
#' 	output.dir        = "..",
#' 	output.basename   = "sample_vs_ref"
#' )
#' }

quandico <- function(
    sample            = NA, # file with sample counts to import (NA)
    sample.x          = 2,  # number of X chromosomes in sample (female = 2, male = 1)
    sample.y          = 0,  # number of Y chromosomes in sample (female = 0, male = 1)
    reference         = NA, # file with reference counts to import (NA)
    reference.x       = 2,  # number of X chromosomes in reference (female = 2, male = 1)
    reference.y       = 0,  # number of Y chromosomes in reference (female = 0, male = 1)
    sample.columns    = c("chr_full", "position", "direction", "base", "gene", "sample", "orgname"),
    reference.columns = c("chr_full", "position", "direction", "base", "gene", "reference", "orgname"),
    output.dir        = "..", # save output files to this directory
    output.basename   = 'quandico_run', # prefix for output files
    assembly.filename = config.assembly.filename, # fasta index to read chromosome lengths from
    assembly.version  = config.assembly.version   # name of the assembly/reference genome (for vcf)
  ) {

  # read the reference and sample coverage files,
  # and provide the expected counts for X and Y
  depth <- import_counts(
      sample            = sample,	# NA will fall back to example data (included)
      sample.x          = as.numeric(sample.x),
      sample.y          = as.numeric(sample.y),
      reference         = reference, # NA will fall back to example data (included)
      reference.x       = as.numeric(reference.x),
      reference.y       = as.numeric(reference.y)
  )
  
  ## get the list of chromosomes with data from the table
  chromosomes <- unique(naturalsort(depth$chr))
  
  ## // START OUTPUT - This is for visual tracking of what is being done
  ##
  options(error = quote(
      {
        sink(file=paste(output.dir , '/', output.basename, '.err', sep=""));
        dump.frames();
        print(attr(last.dump,"error.message"));
        traceback();
        sink();
        q()
      }
    )
  )
  
  ## chromosome handling
  chrlen <- read_chromosomes(assembly.version, assembly.filename) 
  
  ## pdf
  pdf( paste(output.dir , '/', output.basename, '.pdf', sep=""),
       paper       = config.report.pdf.paper,
       colormodel  = config.report.pdf.colormodel,
       title = paste(
           "quandico ",
           config.software.name, ' v', config.software.version,
           sep=""
      )
  )
  
  ###############################################################################
  ##
  ## plot the data as is - if requested in config
  ##
  if (config.plot.rawdata) {
    p <- ggplot(depth, aes(reference,sample))
    p <- p + geom_point(aes(color=chr), alpha=1/3, size=1, na.rm=TRUE) + geom_abline()
    p <- p + scale_x_log10("raw reference") + scale_y_log10("raw sample") 
    p <- p + labs(title="Raw counts scatterplot")
    p
  }
  
  ###############################################################################
  ##
  ##  normalize counts
  ##
  ## This will create the log2 from the ratios and shift by the median 
  ## to normalize the whole set (median =defined=> 0).
  # calculate the raw ratios
  depth$log2  <- log2(depth$sample / depth$reference)
  depth$mincov  <- (depth$sample + depth$reference)
  
  # filter for minimal depth in either ref or sample
  usable <- subset(depth, 
      !is.na(log2) & is.finite(log2) & mincov >= config.min_cov_for_scale
  )
  depth$log2  <- depth$log2 - median(usable$log2, na.rm = TRUE)
  
  # the next step will replace log2 values < min and > max by min and max,
  depth <- pull_in_range( 
        depth, 
        min    = config.log2.min,         # < min = min
        max    = config.log2.max
  )
  
  ## store the weights that will be used later to perform:
  ## 
  ## initialize weights to be all 1
  depth$weight <- rep(1, nrow(depth))
  
  ## use real weights (if requested)
  if (!is.na(config.use.weights)) {
      
    if  (config.use.weights == 'squared') {
      # this will divide the greater of the two normalized values by the mininmal
      # required count, then square the result (bias towards use of mainly big values)
      depth$weight  <- apply(depth, 1, function(x) as.numeric(max(x[13], x[14])))
      depth$weight <- (depth$weight / config.min_cov_for_call)^2
    } else if (config.use.weights == 'raw_max') {
      # the greater of the two raw counts
      depth$weight  <- apply(depth, 1, function(x) as.numeric(max(x[6], x[7]))) # raw
    } else if (config.use.weights == 'norm_max') {
      # the greater of the two normalized counts
      depth$weight  <- apply(depth, 1, function(x) as.numeric(max(x[13], x[14]))) # normalized
    } else if (config.use.weights == 'log2_norm_max') {
      # the greater of the two normalized counts, but log2 converted
      depth$weight  <- apply(depth, 1, function(x) as.numeric(max(x[13], x[14]))) # normalized
      depth$weight <- log2(depth$weight)
    } else if (config.use.weights == 'sum(raw)') { ## THIS IS RECOMMENDED!
      # raw counts summed
      depth$weight <- (depth$sample + depth$reference)     # 
    } else if (config.use.weights == 'log2(sum(raw))') {
      # log2 of raw counts summed
      depth$weight <- log2(depth$sample + depth$reference) # 
    } else if (config.use.weights == 'norm_sum') {
      # normalized counts summed
      depth$weight <- (depth$sample_norm + depth$reference_norm) / 2 
    }
  
  }
  
  depth$cn <- depth$refcopy * (2 ^ depth$log2)
  maxcov   <- max(depth$sample, depth$reference, na.rm = TRUE)
  plots    <- list()
  counter  <- 0
  
  ## // create lists to store results for export after the loop
  ## regions
  genes.all <- data.frame()
  vcf       <- data.frame()
  
  # sort the vcf table by chromosome / position
  depth    <- depth[naturalorder(paste(depth$chr, depth$position)),] 
  allgenes <- depth$gene
  lastgene <- '_NOT_A_GENE_'
  
  ## loop all regions
  for (cluster in unique(depth$gene)) {
  
    # subset gene data
    genedata <- subset_cluster(depth, cluster)
    n.total  <- nrow(genedata)
    
    # smoothed data (no outliers)
    majority  <- subset(genedata, outlier < 1)
    majority  <- majority
    n.usable  <- nrow(majority)
    #majority.dip <- dip.test(sort(majority$log2))  
    
    # use this data frame to centrally store values to report
    call <- data.frame(
        
        ## gene specific data, tunneled through for reporting
        genename    = cluster,
        chr         = genedata$chr[1], 
        from        = min(genedata$position),
        to          = max(genedata$position),
        refcopy     = genedata$refcopy[1],
        
        ## summary of data, to make decisions based on these
        amplicons   = n.total,
        outliers    = n.total - n.usable,
        usable      = n.usable,
        sufficient  = ifelse( n.usable > config.minimal_observations, TRUE, FALSE),
  
        ## in case of copy number changes, collect statistics for that
        alternative = NA,
        p.val       = NA,
        conf.level  = config.conf_level,
        score       = 0,
        qp          = 0,
        filter      = set_filter(0, n.usable, config.score.report, config.score.check, config.minimal_observations),
  
        ## provide an alternative copy number call from weighted_t_test
        cn          = genedata$refcopy[1],
        min         = NA,
        max         = NA,
        se          = NA,
        log2        = NA,
        sd          = NA,
        
        ## test for segmentation of the gene region
        segmented   = FALSE,
        seq.p.value = NA,#1/majority.sd,
        segments    = NA,
        
        ## collect strings and colors for cn output, defaults to "unchanged"
        cn.color    = config.color.unchanged,
        cn.label    = "unchanged",                ## first  line of left panel
        cn.pval     = 1,                         ## second line of left panel
        
        ## collect strings and colors for segmentation, defaults to "unimodal"
        sg.color    = config.color.even,
        sg.label    = "unimodal",           ## first  line of right panel
        sg.pval     = 1, #majority.dip$p.value,  ## second line of right panel
        
        ## flag to set for PDF reporting
        report      = config.report.all
    )
    
    # check for general copy number change (different from reference)
    call <- test_for_change(
          majority, call, 
          conf.level = config.conf_level, 
          minscore   = config.score.report,
          maxscore   = config.score.max,
          correction = config.score.correction,
          checkscore = config.score.check,
          checkcolor = config.color.check
      )  
      
    if (!call$sufficient) {
      call$cn.color <- config.color.insufficient
      call$cn.label <- "insufficient"
      call$cn.pval  <- ''
      call$filter   <- paste('o', config.minimal_observations, sep="")
      call$score    <- 0
    }
    
    xlabel<-sprintf("%s ~ %.1f Mbp", call$chr, mean(genedata$position)/1000000)
    
    if (call$report & call$usable > 0 ) {
      
      ## // adjust the plots y-axis limits for extremely high copy numbers
      these.ybreaks <- config.plots.ybreaks
      these.ylimits <- config.plots.ylimits
      these.ylabels <- config.plots.ylabels
      
      if (!is.na(call$max)) {
        ## increase limits if necessary
        current <- these.ybreaks[length(these.ybreaks)]
        # current
        while (current < (call$max + 2)) {
          current <- current + 4
          these.ybreaks <- c(these.ybreaks, current)
          these.ylabels  <- c(these.ylabels, current )
          these.ylimits[2] <- current + 0.9
        }
      }
      ## // the left panel showing the violin plot for all amplicons of the gene
      panel.left <- ggplot(majority, aes(factor(gene), cn, weight=weight/sum(weight))) + theme_gene_left() # , label=call$cn.label
      ## add layers
      panel.left <- panel.left + geom_abline(intercept=call$refcopy, slope=0, colour="green", size=4, alpha=0.25)
      panel.left <- panel.left + geom_violin(fill=call$cn.color, adjust=config.plots.smooth.violin, alpha=config.plots.alpha.violin, adjust=1.5)
      panel.left <- panel.left + scale_y_continuous(limits=these.ylimits, breaks=these.ybreaks, labels=these.ylabels)  
      panel.left <- panel.left + labs(list(y = "copies", x = sprintf("n = %.0f (%.0f)", call$usable, call$amplicons)), title=xlabel)
      panel.left <- panel.left + geom_text(x=1, y=config.plots.line1, label=paste(call$cn.label), vjust=0, size=3.5, color=call$cn.color, parse=FALSE)
      panel.left <- panel.left + geom_text(x=1, y=config.plots.line2, label=paste(round(call$score,0), " / ", round(call$qp, 0)), vjust=0, size=3.5, color=call$cn.color, parse=FALSE)
      
      ## // the right panel showing the amplicons in chromosomal context (zoom)
      panel.right <- ggplot(genedata, aes(position / 1000000, cn)) + theme_gene_right()
      ## add layers
      panel.right <- panel.right + geom_abline(intercept=call$refcopy, slope=0, colour="green", size=4, alpha=0.25)
      panel.right <- panel.right + scale_size("coverage", trans="log", range=c(1, 3), limits=c(config.min_cov_for_call,maxcov), breaks=c(config.min_cov_for_call,config.min_cov_for_call*10,config.min_cov_for_call*100))
      
      # no segmentation checked, so outliers may be removed if requested
      if(config.flag.outliers | config.narrow.normal | config.flag.extremes) {
        # plot outliers and indicate by color change
        panel.right <- panel.right + geom_point(aes(size=mincov, colour=factor(outlier)), alpha=0.25, na.rm=TRUE) 
        panel.right <- panel.right + scale_colour_manual("outlier", values = c("black", "red"), labels=c("no", "yes"))
      } else {
        # plot all values
        panel.right <- panel.right + geom_point(aes(size=mincov), alpha=0.25, na.rm=TRUE) 
      }
      # label x axis and scale y
      panel.right <- panel.right + labs(list(x = sprintf("chr%s [Mbp]", call$chr)))
      panel.right <- panel.right + scale_y_continuous(limits=these.ylimits, breaks=these.ybreaks, labels=these.ylabels) 
          
      if (config.report.genes) {
        # save the individual gene to a single file
        multiplot(
            plotlist = list(panel.left, panel.right), 
            file     = paste(output.dir , '/', output.basename, '_', cluster, '.pdf', sep=""),
            paper    = config.gene.pdf.paper,
            layout   = config.gene.layout
        )
      }
    } # end if call$report and usable > 0
    else {
      ## this is purely for regions without any data
      panel.left <- ggplot(genedata, aes(factor(gene), cn, weight=weight/sum(weight))) + theme_gene_left() # , label=call$cn.label
      ## add layers
      panel.left <- panel.left + geom_abline(intercept=call$refcopy, slope=0, colour="green", size=4, alpha=0.25)
      panel.left <- panel.left + geom_blank()
      panel.left <- panel.left + scale_y_continuous(limits=config.plots.ylimits, breaks=config.plots.ybreaks, labels=config.plots.ylabels)  
      panel.left <- panel.left + labs(list(y = "copies", x = sprintf("n = %.0f (%.0f)", call$usable, call$amplicons)), title=xlabel)
      panel.left <- panel.left + geom_text(x=1, y=config.plots.line1, label=paste("no data"), vjust=0, size=3.5, color="gray50", parse=FALSE)
      panel.left <- panel.left + geom_text(x=1, y=config.plots.line2, label="- / -", vjust=0, size=3.5, color="gray50", parse=FALSE)
      
      panel.right <- ggplot(genedata, aes(position / 1000000, cn)) + theme_gene_right() + geom_blank()
      ## add layers
      panel.right <- panel.right + geom_abline(intercept=call$refcopy, slope=0, colour="green", size=4, alpha=0.25)
      panel.right <- panel.right + scale_size("coverage", trans="log", range=c(1, 3), limits=c(config.min_cov_for_call,maxcov), breaks=c(config.min_cov_for_call,config.min_cov_for_call*10,config.min_cov_for_call*100))
      panel.right <- panel.right + labs(list(x = sprintf("chr%s [Mbp]", call$chr)))
      panel.right <- panel.right + scale_y_continuous(limits=config.plots.ylimits, breaks=config.plots.ybreaks, labels=config.plots.ylabels)
    }
  
    # collect this plot	
    counter <- counter + 1
    plots[[counter]] <- panel.left
    counter <- counter + 1
    plots[[counter]] <- panel.right
    
    # do the multiplotting
    if (counter == 4) {
      multiplot(plotlist = plots, layout=config.report.layout)
      plots<-list()
      counter = 0
    }
    
    # assemble results for final output
    call$cn.pval <- NULL
    genes.all    <- rbind(genes.all, call)
  
    # dispatch results to lists for export
    if (call$report) {
      bp      <- call$to -  call$from + 1
      refline <- subset(depth, chr == call$chr & position == call$from  )
  
      if (call$filter == 'PASS') {
        vcf <- rbind(
            vcf, data.frame(
                CHROM    = call$chr, #paste('chr', call$chr, sep=""),
                POS      = call$from, 
                ID       = '.',
                REF      = refline$base[1],
                ALT      = ifelse( call$cn < call$refcopy , '<DEL>', '<DUP>'),
                QUAL     = round(call$score),
                FILTER   = call$filter,
                INFO     = paste(
                  'IMPRECISE', ';',
                  'SVTYPE=', ifelse( call$cn < call$refcopy , 'DEL', 'DUP'), ';',
                  'END=', call$to, ';',
                  'SVLEN=', ifelse( call$cn < call$refcopy , -bp, bp), ';',
                  'LOCUS=', call$genename,
                   sep=""
                ),
                FORMAT    = 'CN:CNMIN:CNMAX:CNQ:PVAL',
                SAMPLE    = paste(
                  round(call$cn,  2), 
                  round(call$min, 2), 
                  round(call$max, 2),
                  round(call$qp,  0),
                  format(call$p.val, digits=2),
                  sep=":"
                )
            )
        )
      } else if (call$score >= config.score.check) {
        vcf <- rbind(
            vcf, data.frame(
                CHROM   = call$chr, #paste('chr', call$chr, sep=""),
                POS    = call$from, 
                ID      = '.',
                REF     = refline$base[1],
                ALT      = ifelse( call$cn < call$refcopy , '<DEL>', '<DUP>'),
                QUAL     = round(call$score),
                FILTER   = call$filter,
                INFO     = paste(
                    'IMPRECISE', ';',
                    'SVTYPE=', ifelse( call$cn < call$refcopy , 'DEL', 'DUP'), ';',
                    'END=', call$to, ';',
                    'SVLEN=', ifelse( call$cn < call$refcopy , -bp, bp), ';',
                    'LOCUS=', call$genename,
                    sep=""
                ),
                FORMAT    = 'CN:CNMIN:CNMAX:CNQ:PVAL',
                SAMPLE    = paste(
                    round(call$cn,  2), 
                    round(call$min, 2), 
                    round(call$max, 2),
                    round(call$qp,  0),
                    format(call$p.val, digits=2),
                    sep=":"
                )
            )
        )
      } else {
        
        info <- paste('UNCHANGED;END=', call$to,';LOCUS=', call$genename, sep="")
        if (!call$sufficient) {
          info <- paste('NOCALL;END=', call$to,';LOCUS=', call$genename, sep="")
        }
        vcf <- rbind(
            vcf, data.frame(
                CHROM   = call$chr, #paste('chr', call$chr, sep=""),
                POS     = call$from, 
                ID      = '.',
                REF     = refline$base[1],
                ALT     = '.',
                QUAL    = round(call$score),
                FILTER  = call$filter,
                INFO    = info,
                FORMAT  = 'CN:CNMIN:CNMAX:CNQ:PVAL',
                SAMPLE  = paste(
                    round(call$cn,  2), 
                    round(call$min, 2), 
                    round(call$max, 2),
                    round(call$qp,  0),
                    format(call$p.val, digits=2),
                    sep=":"
                )
            )
        )
        }
    }
  }
  
  if (counter == 2) {
    # add the last plot to the pdf (if any)
    multiplot(plotlist = plots, layout=config.report.layout)
  }
  
  if (nrow(vcf) > 0) {
    
  # export the vcf
  writeLines(
    paste(vcf_header(assembly.version, assembly.filename), '
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=0,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=UNCHANGED,Number=0,Type=Flag,Description="No copy number variation in this region">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SEGMENTED,Number=0,Type=Flag,Description="Structural variation segments named locus group (e.g. gene)">
##INFO=<ID=NOCALL,Number=0,Type=Flag,Description="Copy number analysis was not possible (see FILTER for reason)">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=LOCUS,Number=1,Type=String,Description="Locus description (e.g. gene name in this region)">
##ALT=<ID=CNV,Description="Copy number variable region">
##ALT=<ID=DEL,Description="Deletion relative to reference (decreased copy number, loss)">
##ALT=<ID=DUP,Description="Duplication relative to reference (increased copy number, gain)">
##FILTER=<ID=q', config.score.check ,',Description="Region not changed in copy number: Score below ', config.score.check,'">
##FILTER=<ID=q', config.score.report ,',Description="Region potentially changed in copy number: Score below ', config.score.report,'">
##FILTER=<ID=o', config.minimal_observations , ',Description="Less than ' , config.minimal_observations ,' amplicons have suffcient read depth (<' , config.min_cov_for_call , ')">
##FORMAT=<ID=CN,Number=1,Type=Float,Description="Estimated copy number for imprecise events">
##FORMAT=<ID=CNMIN,Number=1,Type=Float,Description="Minimal copy number supported by statistical test (conf=', config.conf_level ,')">
##FORMAT=<ID=CNMAX,Number=1,Type=Float,Description="Maximal copy number supported by statistical test (conf=', config.conf_level ,')">
##FORMAT=<ID=CNQ,Number=1,Type=Integer,Description="Copy number quality score">
##FORMAT=<ID=PVAL,Number=1,Type=Float,Description="Probability for normal copy number">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE',
    sep=""),
    paste(output.dir , '/', output.basename, '.vcf', sep="")
  )
  
  # sort the vcf table by chromosome / position
  vcf <- vcf[naturalorder(paste(vcf$CHROM, vcf$POS)),]
  vcf$CHROM <- paste('chr', vcf$CHROM, sep="")
  
  # export the vcf to a file
  write.table(
      vcf,
      paste(output.dir , '/', output.basename, '.vcf', sep=""),
      sep       = "\t", 
      quote     = FALSE,
      append    = TRUE,
      row.names = FALSE,
      col.names = FALSE
  )
  }
  
  ## export excel file (tab separated for now)
  ##
  ## first select columns to include in export:
  ## available:
  ## genename	chr	from	to	refcopy	amplicons	outliers	usable	noise	sufficient	noisy	alternativep.val	conf.level	score	filter	cn	min	max	se	ml.cn	ml.prob	segmented	seq.p.value	segments	cn.color	cn.label	sg.color	sg.label	sg.pval	report
  if (nrow(genes.all)> 0) {
  
    export <- genes.all[c("chr", "from", "to", "genename",  "amplicons", "outliers", "usable", "refcopy", "log2", "cn", "min", "max", "qp", "sd", "score", "p.val", "filter")]
    names(export) <- c("chromosome", "start", "end", "locus", "amplicons", "outliers", "usable", "expected", "log2", "copies", "min", "max", "qp", "sd", "score", "p.val", "filter")
    
    # round to integer
    export$score <- round(export$score, 0)
    export$qp    <- round(export$qp, 0)
    
    # round to two decimal
    export$log2   <- round(export$log2,   2)
    export$copies <- round(export$copies, 2)
    export$min    <- round(export$min, 2)
    export$max    <- round(export$max, 2)
    export$sd     <- round(export$sd, 2)
    
    # format p.values
    export$p.val <- format(export$p.val, digits=2)
    
    # remove NA
    export$p.val[is.na(export$copies)]  <- ''
    export$score[is.na(export$copies)]  <- ''
    export$min[is.na(export$copies)]    <- ''
    export$max[is.na(export$copies)]    <- ''
    export$log2[is.na(export$copies)]   <- ''
    export$sd[is.na(export$copies)]     <- ''
    export$copies[is.na(export$copies)] <- ''
    
    # sort by chromosome and position
    export            <- export[naturalorder(paste(export$chromosome, export$start)),] 
    export$chromosome <- paste('chr', export$chromosome, sep="")
    
    # save the subset
    write.table(
        export,
        paste(output.dir , '/', output.basename, '.csv', sep=""), 
        sep=",", 
        quote = FALSE,
        row.names = FALSE
    )
  
  }
  
  # create the overview with final reported ranges (separate page)
  plot_overview(
      csv      = paste(output.dir , '/', output.basename, '.csv', sep=""), 
      filename = paste(output.dir , '/', output.basename, '_overview.pdf', sep="")
  )
  
  # save the full dataset (for debugging)
  if (config.write.dataset) {
    write.table(
        depth,
        paste(output.dir , '/', output.basename, '_full_data.tsv', sep=""), 
        sep="\t", 
        quote = FALSE,
        row.names = FALSE
    )
  }
  
  # all finished
  devoff <- dev.off()
}
