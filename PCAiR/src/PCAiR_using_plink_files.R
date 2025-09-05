# conda activate R4.4
library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(Matrix)

setwd('/data100t1/home/wanying/BioVU/20250711_BP_in_hispanic_AGD/data/PCAiR')

# ##### Step 1. Conver plink file to GDS file, and read in GDS data #####
cat('# Step 1. Conver plink file to GDS file, and read in GDS data\n')
plink_bed = "/data100t1/home/wanying/BioVU/20250711_BP_in_hispanic_AGD/data/plink_files_for_pc_and_null_model/merged_all_chr.bed"
plink_bim = "/data100t1/home/wanying/BioVU/20250711_BP_in_hispanic_AGD/data/plink_files_for_pc_and_null_model/merged_all_chr.bim"
plink_fam = "/data100t1/home/wanying/BioVU/20250711_BP_in_hispanic_AGD/data/plink_files_for_pc_and_null_model/merged_all_chr.fam"
gds_fn = 'genotype.gds'
SNPRelate::snpgdsBED2GDS(bed.fn = plink_bed, 
                         bim.fn = plink_bim, 
                         fam.fn = plink_fam, 
                         out.gdsfn = gds_fn)
gds <- SNPRelate::snpgdsOpen(gds_fn)

# ##### Step 2. LD pruning #####
cat('\n# Step 2. LD pruning\n')
snpset <- SNPRelate::snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, ld.threshold=sqrt(0.1), verbose=T, num.thread=20)
pruned <- unlist(snpset, use.names=FALSE)
write.table(pruned, 'pruned_SNPs_for_PCAiR.txt', col.names=F, row.names=F, sep='\t', quote=F)
cat(length(pruned))

# ##### Step 3. Pairwise Measures of Ancestry Divergence #####
cat('\n# Step 3. Pairwise Measures of Ancestry Divergence\n')
# Calculate the KING-robust estimates directly from a GDS file
king <- SNPRelate::snpgdsIBDKING(gds, type='KING-robust', num.thread=20, verbose =T)
KINGmat <- GENESIS::kingToMatrix(king)
SNPRelate::snpgdsClose(gds)
print(KINGmat[1:5,1:5])

# ##### Step 4. Run PC-AiR #####
cat('\n\n# Step 4. Run PC-AiR\n')
HapMap_geno <- GWASTools::GdsGenotypeReader(filename = gds_fn)

# Create a GenotypeData class object
HapMap_genoData <- GWASTools::GenotypeData(HapMap_geno)

# run PC-AiR on pruned SNPs
mypcair <- pcair(HapMap_genoData,
                 kinobj = KINGmat,
                 divobj = KINGmat,
                 snp.include = pruned)

summary(mypcair)


cat('\n\n# Save results\n')
# Name column headers
colnames(mypcair$vectors) = sapply(1:ncol(mypcair$vectors), function(x) paste0('PCAiR',x))
# Save PCs and other values
# Make the matrix a df so row names can have column header (Use "SampleID" as row name header)
df = data.frame(SampleID = rownames(mypcair$vectors), mypcair$vectors, row.names = NULL, check.names = FALSE)
write.table(df, "PCAiR.eigenvector.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(mypcair$values, 'PCAiR.eigenvalue.txt', col.names=T, row.names=T, sep='\t', quote=F)
write.table(mypcair$varprop, 'PCAiR.varpop.txt', col.names=T, row.names=T, sep='\t', quote=F)


