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
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(readr))

spec = matrix(c(
	'vcf', 		'v', 1, "character",
	'regions', 	'r', 1, "character",
	'main', 	'm', 1, "character",
	'class', 	'l', 1, "character",
	'height', 	'h', 1, "integer",
	'width', 	'w', 1, "integer",
	'rmargin', 	'i', 1, "integer",
	'plotOnly', 'p', 0, "logical",
	'cex', 		'c', 1, "double",
	'out', 		'o', 1, "character"
), byrow=TRUE, ncol=4);

opt = getopt(spec);

if( is.null(opt$cex) ){
	opt$cex = 1
}

# Read in sites file
regions = ''
if( ! is.null(opt$regions) ){
	regions = paste("-R", opt$regions)
}
cmd = ""

opt$discordance_pairs = gsub(".pdf$", "", opt$out)

cmd = paste0(cmd, "bcftools view ", regions, " ", opt$vcf, " | bcftools gtcheck -G1 -p", opt$discordance_pairs)

if( is.null(opt$plotOnly) ){
	cat("Computing concordance...\n")
	# run command
	system(cmd)

	# grep sites
	system(paste0('echo -e "CN\tDiscordance\\tnsites\\tAMD\\tsample_i\\tsample_j" >', opt$discordance_pairs, "_2.tab"))
	system(paste0('grep \"^CN\" ', opt$discordance_pairs, ".tab | sed 's/pholder/_/g' >>", opt$discordance_pairs, "_2.tab"))
}


cat("Reading concordance results...\n")

getDistanceMatrix = cmpfun(function( file ){

	# read discordance scores
	data = read_tsv( file )

	# create matrix storing all pairs of discordance
	# use sample names as row/col names
	sampleIDs = sort(unique(c(data$sample_i, data$sample_j)))
	data$sample_i = factor(data$sample_i, sampleIDs)
	data$sample_j = factor(data$sample_j, sampleIDs)

	discorMat = matrix(.5, nlevels(data$sample_j), nlevels(data$sample_j))
	rownames(discorMat) = levels(data$sample_j)
	colnames(discorMat) = levels(data$sample_j)
	diag(discorMat) = 1

	idx = which(data$nsites != 0)
	n_rows = length(idx)

	# loop thru all pairs
	for(k in idx ){

		if( k %% 3000 == 0 ){
			cat("\r", match(k, idx), "/", n_rows, "  ", round(match(k, idx) / n_rows * 100, 1), "%")
		}

		value = data$Discordance[k] / data$nsites[k]

		value = ifelse( is.nan(value), .5, value)

		discorMat[data$sample_i[k], data$sample_j[k]] = value
		discorMat[data$sample_j[k], data$sample_i[k]] = value
	}
	cat("\r", match(k, idx), "/", n_rows, "  ", round(n_rows / n_rows * 100, 1), "%", "\n\n")

	# discorMat[data$sample_i, data$sample_j] = with(data, Discordance / nsites)

	n_sites = paste(range(data$nsites), collapse=' - ')

	return( list(discorMat = discorMat, n_sites = n_sites, data=data) )
})

res = getDistanceMatrix( paste0(opt$discordance_pairs, "_2.tab") )

discorMat = res$discorMat
n_sites = res$n_sites
data = res$data

cat("# NAN: ", sum(is.nan(discorMat)), "\n")
cat("# 0.5: ", sum(discorMat == 0.5), "\n\n")


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

# opt$cex = .6
# opt$height = 150

graphics.off()
pdf( file=opt$out, height=opt$height, width=opt$width)
par(mar=c(5,1,1,opt$rmargin), cex=opt$cex)
plot(as.dendrogram(hcl), horiz=TRUE, xlab="Discordance", main=opt$main)
abline(v=0.01, lty=2, col="grey", lwd=2)
abline(v=0.05, lty=3, col="grey", lwd=2)
legend("topleft", legend=c("1%", "5%"), lty=2:3, col="grey", lwd=2, bty='n')
legend("bottomleft", legend=paste0("# sites:", n_sites), bty='n')
dev.off()

if( ! is.null(opt$class) & file.exists(opt$class) ){

	cat("Writing summary file...\n")
	info = read.table( opt$class, header=TRUE)
	info$Sample = as.character(info$Sample)

	result = c()

	for( i in 1:ncol(discorMat)){
		if( i %% 50 ==0 ){
			cat("\r", i, "/", ncol(discorMat), "   ",  paste(format(i / ncol(discorMat) * 100, digits=3), "%"), "\t")
		}

		for( class in levels(info$Class) ){

			idx = match(with(info, Sample[Class==class]), colnames(discorMat))

			if( sum(!is.na(idx)) == 0){
				next
			}

			k = which.min( discorMat[i,idx])

			id = colnames(discorMat)[i]
			id2 = colnames(discorMat)[idx[k]]

			nSites = with(data, nsites[(sample_i %in% c(id, id2)) &  (sample_j %in% c(id, id2))])

			value = c(Sample = colnames(discorMat)[i], CompareToClass = class, BestMatchInClass=id2, 
				Discordance = format(discorMat[i,idx[k]], digits=2), nSites = nSites)

			result = rbind(result, value)
		}
	}
	cat("\r", i, "/", ncol(discorMat), "   ",  paste(format(i / ncol(discorMat) * 100, digits=3), "%"), "\n\n")

	write.table( result, paste0(opt$discordance_pairs, "_results.tsv"), quote=FALSE, row.names=FALSE, sep="\t")
}
