# variantConcordance
Compute hierarchical clustering from RNA-seq/WES/genotyping data


# Running analysis
~~~~
module load R
module load bcftools tabix

VCF={Your favorite VCF}

# Create hierarchical clustering using all sites
./variant_similarity_tree.R \
--vcf $VCF \
--main "Example" \
--out example_tree.pdf


# Create hierarchical clustering using only sites in example.regions
./variant_similarity_tree.R \
--vcf $VCF \
--regions example.regions \
--main "Example: regions" \
--out example_regions_tree.pdf
~~~~

If the plot is not readable, you can change the dimensions with the --height and --width arguments.  Defaults are 10 and 7 inches, respectively

# Merging VCF
~~~~
# Merge VCF's and extract SNP's
bcftools merge -m none *.vcf.gz | bcftools view -v snps | bgzip > combined.vcf.gz
~~~~
