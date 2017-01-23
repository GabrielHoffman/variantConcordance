#!/bin/env Rscript

# Gabriel Hoffman
# August 19, 2016

# DATA: VCF of genotypes

# OBJECTIVE: Show similarity of variants based on VCF

# METHOD: bcftools gtcheck and hierarchical clustering

#######
# RUN #
#######
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(readr))

spec = matrix(c(
	'vcf', 		'v', 1, "character",
	'regions', 	'r', 1, "character",
	'main', 	'm', 1, "character",
	'height', 	'h', 1, "integer",
	'width', 	'w', 1, "integer",
	'rmargin', 	'i', 1, "integer",
	'plotOnly', 'p', 0, "logical",
	'out', 		'o', 1, "character"
), byrow=TRUE, ncol=4);

opt = getopt(spec);

# Read in sites file
regions = ''
if( ! is.null(opt$regions) ){
	regions = paste("-R", opt$regions)
}
cmd = ""

opt$discordance_pairs = gsub(".pdf$", "", opt$out)

cmd = paste0(cmd, "bcftools view ", regions, " ", opt$vcf, " | bcftools gtcheck -G1 -p", opt$discordance_pairs)

if( ! is.null(opt$plotOnly) ){
	cat("Computing concordance...\n")
	# run command
	system(cmd)
}

# grep sites
system(paste0('echo -e "CN\tDiscordance\\tnsites\\tAMD\\tsample_i\\tsample_j" >', opt$discordance_pairs, "_2.tab"))
system(paste0('grep \"^CN\" ', opt$discordance_pairs, ".tab | sed 's/pholder/_/g' >>", opt$discordance_pairs, "_2.tab"))

cat("Reading concordance results...\n")
# read discordance scores
data = read_tsv(paste0(opt$discordance_pairs, "_2.tab"))

# create matrix storing all pairs of discordance
# use sample names as row/col names
sampleIDs = sort(unique(c(data$sample_i, data$sample_j)))
data$sample_i = factor(data$sample_i, sampleIDs)
data$sample_j = factor(data$sample_j, sampleIDs)

discorMat = matrix(NA, nlevels(data$sample_j), nlevels(data$sample_j))
rownames(discorMat) = levels(data$sample_j)
colnames(discorMat) = levels(data$sample_j)
diag(discorMat) = 1

# loop thru all pairs
for(k in 1:nrow(data)){
	discorMat[data$sample_i[k], data$sample_j[k]] = data$Discordance[k] / data$nsites[k]
	discorMat[data$sample_j[k], data$sample_i[k]] = data$Discordance[k] / data$nsites[k]
}
n_sites = paste(range(data$nsites), collapse=' - ')


# image and subset of matrix
# image(discorMat)
# discorMat[20:25, 20:25]

# hierarchical clustering
hcl = hclust( as.dist( discorMat ) )

if( is.null(opt$height)){
	opt$height = 10
}

if( is.null(opt$width)){
	opt$width = 7
}

if( is.null(opt$rmargin)){
	opt$rmargin = 13
}

graphics.off()
pdf( file=opt$out, height=opt$height, width=opt$width)
par(mar=c(5,1,1,opt$rmargin))
plot(as.dendrogram(hcl), horiz=TRUE, xlab="Discordance", main=opt$main)
abline(v=0.01, lty=2, col="grey", lwd=2)
abline(v=0.05, lty=3, col="grey", lwd=2)
legend("topleft", legend=c("1%", "5%"), lty=2:3, col="grey", lwd=2, bty='n')
legend("bottomleft", legend=paste0("# sites:", n_sites), bty='n')
dev.off()


