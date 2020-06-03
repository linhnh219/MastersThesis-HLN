source("https://bioconductor.org/biocLite.R")
#install.packages("RCurl", repos='http://cran.us.r-project.org')
#library(RCurl)
#biocLite("DNAcopy")
library(DNAcopy)
#biocLite("genoset")
library(genoset)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
	stop("Argument missing", call.=FALSE)
}

data = read.table(args[1], header=FALSE)
#dim(data)
names(data) = c('CHROM','POS', 'REF', 'ALT', 'DEPTH', 'GQUAL', 'BAF')
var.names = paste("v", 1:nrow(data), sep="")
sample.names='Sample'
num.samples = 1
num.vars = length(var.names)

locs = GRanges(ranges= IRanges(start=data$POS,end=data$POS),seqnames=factor(data$CHROM, levels=unique(data$CHROM)))

pData=data.frame(matrix(sample.names,nrow=1,ncol=1,dimnames=list(sample.names,'a')))

baf = matrix(data$BAF,nrow=nrow(data),ncol=1,dimnames=list(var.names,sample.names))
#mbaf = matrix(abs(data$BAF-0.5),nrow=nrow(data),ncol=1,dimnames=list(var.names,sample.names))

baf.ds = GenoSet(rowRanges=locs, assays=list(baf=baf), colData=pData )
#baf.ds[, , "mbaf"] = baf2mbaf(baf.ds[, , "baf"], hom.cutoff = 0.95)

baf.ds[, ,"mbaf"] = matrix(baf2mbaf(baf.ds[, , "baf"], hom.cutoff = 1))
baf.ds[,,"baf.segs"] = runCBS(baf.ds[, ,"mbaf"], rowRanges(baf.ds))

segments = as.data.frame(runCBS(baf.ds[, ,"mbaf"], rowRanges(baf.ds), return.segs=TRUE))

for (i in 1:nrow(segments)) {
	if (segments[i,6] > 0.9) segments[i,7] = 'LOH'
	else if (segments[i,6] > 0.65) segments[i,7] = 'AI'
	else segments[i,7] = 'N'
}

write.table(segments[,2:7],args[2],sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
