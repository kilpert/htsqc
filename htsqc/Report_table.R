library(gdata)
library(Rsamtools)
library(R.utils)

.Library
.Library.site
.libPaths()

sessionInfo()

print("Generating project report...")

args = commandArgs(TRUE)

## For debugging only!!! #######################################################
# args = c('/data/processing/kilpert/test/htsqc/data/A237_Trompouki_Trompouki',
#        '/data/processing/kilpert/test/htsqc/output/A237_Trompouki_Trompouki',
#        '1000000',
#        'False')
# setwd(args[2])
################################################################################
main_indir = args[1]
main_outdir = args[2]
total_reads = as.numeric(args[3])
paired = as.logical(args[4]=="True")

## BAM #########################################################################

indir = file.path(main_outdir,"HISAT")
if ( file.exists( indir ) ){
  files <- list.files(indir, pattern="*.bam$", full.names=T)
  files
  
  names = c()
  total = c()
  mapped = c()
  lowqual = c()
  discon = c()
  spliced = c()
  tlen_median = c()
  tlen_1k = c()
  for (file in files) {
    # file = files[[15]]
    name = gsub(".bam$","",basename(file))
    name
    names = c(names, name)
    
    ##     scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA, 
    ##                 hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
    ##                 isFirstMateRead = NA, isSecondMateRead = NA, isNotPrimaryRead = NA,
    ##                 isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
    ##                 isDuplicate = NA)
    
    ## usual (concordant mapping) ##############################################
    if (paired==T) {
      flag = scanBamFlag( isPaired = T,
                          isProperPair = T,
                          isFirstMateRead=T,
                          isUnmappedQuery=F,
                          isSecondaryAlignment=F
      )
    } else {
      flag = scanBamFlag( isUnmappedQuery=F,
                          isSecondaryAlignment=F
      )
    }
    ##p = ScanBamParam(what=scanBamWhat(), flag=flag)
    p = ScanBamParam(what=c("qname", "mapq", "cigar", "isize"), flag=flag)
    
    ## mapped ##################################################################
    filt <- list(function(x) x$mapq >= 5)
    bam_mapped <- filterBam(file, tempfile(), param=p, filter=FilterRules(filt))
    num_mapped = countBam(bam_mapped)$records
    mapped = c(mapped, num_mapped)
    
    ## lowqual #################################################################
    filt_lowqual <- list(function(x) x$mapq < 5)
    bam_lowqual <- filterBam(file, tempfile(), param=p, filter=FilterRules(filt_lowqual))
    num_lowqual = countBam(bam_lowqual)$records
    lowqual = c(lowqual, num_lowqual)

    ## TLEN median
    if (paired==T) {
      b = scanBam(bam_mapped)[[1]]
      num_tlen_median = summary(abs(b$isize))[["Median"]]
    } else {
      num_tlen_median = NA
    } 
    tlen_median = c(tlen_median, num_tlen_median)
    
    ## TLEN 1k #################################################################
    filt <- list(function(x) x$isize >= 1000)
    if (paired==T) {
      bam_mapped_tlen_1k <- filterBam(bam_mapped, tempfile(), filter=FilterRules(filt))
      num_tlen_1k = countBam(bam_mapped_tlen_1k)$records
    } else {
      num_tlen_1k = NA
    }
    tlen_1k = c(tlen_1k, num_tlen_1k)
    
    ## spliced #################################################################
    filt <- list(function(x) grepl("N", x$cigar))
    bam_mapped_spliced <- filterBam(bam_mapped, tempfile(), filter=FilterRules(filt))
    num_spliced = countBam(bam_mapped_spliced)$records
    spliced = c(spliced, num_spliced)
    
    ## disconcordant ###########################################################
    if (paired==T) {
      flag = scanBamFlag( isPaired = T,
                          isProperPair = F,
                          isFirstMateRead=T,
                          isUnmappedQuery=F,
                          isSecondaryAlignment=F
      )
      p = ScanBamParam(what=c("qname", "mapq", "cigar", "isize"), flag=flag)
      filt <- list(function(x) x$mapq >= 5)
      discon_bam_mapped <- filterBam(file, tempfile(), param=p, filter=FilterRules(filt))
      num_discon = countBam(discon_bam_mapped)$records
    } else {
      num_discon = NA
    }
    discon = c(discon, num_discon)
    
    ############################################################################
    
    cat(name, num_mapped, num_lowqual, num_spliced, num_discon, num_tlen_1k, num_tlen_median, "\n")
  }
}  

report = data.frame(row.names=names)
report$project = basename(main_indir)
report$sample = sapply(rownames(report), function (x) unlist(strsplit(x, "___"))[[2]] )
report$genome = sapply(rownames(report), function (x) unlist(strsplit(x, "___"))[[1]] )
if (paired==T) { 
  report$lib_type = "PE"
} else { 
  report$lib_type = "SE"
}
report$total = total_reads
report$total_perc = 100

report$mapped = mapped
report$mapped_perc = report$mapped / report$total * 100

report$lowqual = lowqual
report$lowqual_perc = report$lowqual / report$total * 100

report$spliced = spliced
report$spliced_perc = report$spliced / report$total * 100

report$discon = discon
report$discon_perc = report$discon / report$total * 100

report$tlen_1k = tlen_1k
report$tlen_1k_perc = report$tlen_1k / report$total * 100

report$tlen_median = tlen_median
rownames(report) = NULL

report$project = sub("Project_", "", report$project)

report

write.table(report,"Report.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA) # save to file

