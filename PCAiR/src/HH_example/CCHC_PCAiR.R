source activate genesis
library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(Matrix)

snpgdsBED2GDS(bed.fn = "/vgipiper04/CCHC/QC_preimputed_cleandata_APR22/preimputation_b38_TOPMed_aligned/CCHC_rescan_nodups_v2_nobatchsnps_b38-updated.bed",
	     bim.fn = "/vgipiper04/CCHC/QC_preimputed_cleandata_APR22/preimputation_b38_TOPMed_aligned/CCHC_rescan_nodups_v2_nobatchsnps_b38-updated.bim",
	     fam.fn = "/vgipiper04/CCHC/QC_preimputed_cleandata_APR22/preimputation_b38_TOPMed_aligned/CCHC_rescan_nodups_v2_nobatchsnps_b38-updated.fam", 
	     out.gdsfn = "cchc.gds")

gds = snpgdsOpen("cchc.gds")
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, ld.threshold=sqrt(0.1), verbose=T, num.thread=20)
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)
#[1] 184427
write.table(pruned, 'CCHC_pruned_SNPs.txt', col.names=F, row.names=F, sep='\t', quote=F)

king <- snpgdsIBDKING(gds, type='KING-robust', num.thread=20, verbose =T)
KINGmat <- kingToMatrix(king)
snpgdsClose(gds)
gds_geno <- GdsGenotypeReader(filename = 'cchc.gds')
genoData <- GenotypeData(gds_geno)
genoData
mypcair <- pcair(genoData, kinobj = KINGmat, divobj = KINGmat, snp.include = pruned, num.cores=20, verbose=T)
summary(mypcair)
#Call:
#.pcair(gdsobj = gdsobj, kinobj = kinobj, divobj = divobj, kin.thresh = kin.thresh,
#           div.thresh = div.thresh, unrel.set = unrel.set, sample.include = sample.include,
#	       snp.include = snp.include, num.cores = num.cores, verbose = verbose)
#PCA Method: PC-AiR
#Sample Size: 4841
#Unrelated Set: 3055 Samples
#Related Set: 1786 Samples
#Kinship Threshold: 0.02209709
#Divergence Threshold: -0.02209709
#Principal Components Returned: 32
#Eigenvalues: 28.286 6.665 3.596 2.146 2.032 2.028 1.905 1.881 1.866 1.857 ...
#SNPs Used: 184427

colnames(mypcair$vectors) = sapply(1:32, function(x) paste0('PCAiR',x))
write.table(mypcair$vectors, 'CCHC_PCAiR.eigenvector.txt', col.names=T, row.names=T, sep='\t', quote=F)
write.table(mypcair$values, 'CCHC_PCAiR.eigenvalue.txt', col.names=T, row.names=T, sep='\t', quote=F)
write.table(mypcair$varprop, 'CCHC_PCAiR.varpop.txt', col.names=T, row.names=T, sep='\t', quote=F)
