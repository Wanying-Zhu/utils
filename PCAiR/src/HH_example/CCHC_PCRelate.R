library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(Matrix)
gds = snpgdsOpen("cchc.gds")
pruned = read.table('CCHC_pruned_SNPs.txt', header=F, stringsAsFactors=F)[,1]
king <- snpgdsIBDKING(gds, type='KING-robust', num.thread=20, verbose =T)
#IBD analysis (KING method of moment) on genotypes:
#Excluding 0 SNP on non-autosomes
#Excluding 0 SNP (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
#    # of samples: 4,841
#    # of SNPs: 714,911
#    using 20 threads
#No family is specified, and all individuals are treated as singletons.
#Relationship inference in the presence of population stratification.
#KING IBD:    the sum of all selected genotypes (0,1,2) = 1839266145
#CPU capabilities: Double-Precision SSE2
#Sat Jun 25 00:43:45 2022    (internal increment: 20480)
#[==================================================] 100%, completed, 1.2m
#Sat Jun 25 00:45:01 2022    Done.
KINGmat <- kingToMatrix(king)
#Using 4841 samples provided
#Identifying clusters of relatives...
#    4816 relatives in 3 clusters; largest cluster = 4812
#Creating block matrices for clusters...
#25 samples with no relatives included
#Putting all samples together into one block diagonal matrix
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

CCHC_genoData = GenotypeBlockIterator(genoData, snpInclude=pruned)
mypcrelate <- pcrelate(CCHC_genoData, pcs = mypcair$vectors[,1:5], training.set = mypcair$unrels, BPPARAM = BiocParallel::SerialParam())
#Using 1 CPU cores
#4841 samples to be included in the analysis...
#Betas for 5 PC(s) will be calculated using 3055 samples in training.set...
#Running PC-Relate analysis for 4841 samples using 184427 SNPs in 19 blocks...
#Performing Small Sample Correction...

write.table(mypcrelate$kinBtwn, 'CCHC_PCRelate.kinBtwn.txt', col.names=T, row.names=F, sep='\t', quote=F)
write.table(mypcrelate$kinSelf, 'CCHC_PCRelate.kinSelf.txt', col.names=T, row.names=F, sep='\t', quote=F)

primus = mypcrelate$kinBtwn[,c(1,2,4)]
primus$k1 = 1 - mypcrelate$kinBtwn$k0 - mypcrelate$kinBtwn$k2
primus$k2 =  mypcrelate$kinBtwn$k2
primus$pi = primus$k2 + 0.5*primus$k1
primus$FID1 = 0
primus$FID2 = 0
primus$RT = 'UN'
primus$EZ = NA
primus = primus[,c('FID1','ID1','FID2','ID2','RT','EZ','k0','k1','k2','pi')]
colnames(primus) = c('FID1','IID1','FID2','IID2','RT','EZ','IDB0','IBD1','IBD2','PI_HAT')
write.table(primus, 'CCHC_PCRelate.genome', col.names=T, row.names=F, sep='\t', quote=F)

cchc_grm2  = pcrelateToMatrix(mypcrelate, scaleKin = 2)
#Using 4841 samples provided
#Identifying clusters of relatives...
#    4841 relatives in 1 clusters; largest cluster = 4841
#Creating block matrices for clusters...
#0 samples with no relatives included
write.table(as.matrix(cchc_grm2), 'CCHC_PCRelate.GRM.2.csv', col.names=T, row.names=T, sep=',', quote=F)
write.table(as.matrix(cchc_grm2), 'CCHC_PCRelate.GRM.2.txt', col.names=T, row.names=F, sep='\t', quote=F)


