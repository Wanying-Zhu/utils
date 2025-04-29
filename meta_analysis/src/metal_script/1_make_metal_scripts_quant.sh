### METAL script
# WZ notes:
# Check reference: https://genome.sph.umich.edu/wiki/METAL_Documentation
# METAl script can be run as: metal some_metal_script


## Fill in directory:
mydir="/proj/epi/CVDGeneNas/xinruo/page_hl_htn/06_meta_METAL_mvp" # No it doesn't end with a /

mkdir -p ${mydir}/metal_jobs
mkdir -p ${mydir}/metal_res


for trait in sbp dbp; do

# Check if the submit file exists and remove it if it does:
  submit_file="submitMETAL_${trait}.sh"
  if [ -f "$submit_file" ]; then
    rm "$submit_file"
  fi

# SCHEME SAMPLESIZE
# default approach, uses p-value and direction of effect, weighted according to sample size
echo "
# === START ANALYSIS ===
SCHEME STDERR 			# Use effect size estimates and standard errors
GENOMICCONTROL ON 		# Automatically correct test stats to account for small amounts of population stratification or unaccounted for relatedness

# Track Allele Frequencies
AVERAGEFREQ ON
MINMAXFREQ ON
#FLIP

VERBOSE OFF
#ADDFILTER EffN > 30
#ADDFILTER INFO > 0.4

CUSTOMVARIABLE TotalSampleSize




# === DESCRIBE AND PROCESS THE 1st (MEGA) INPUT FILE ===
SEPARATOR TAB
MARKER VCF_ID
ALLELE ALT REF  # Test allele first
FREQ ALT_AF
EFFECT BETA
STDERR SE 
PVALUE PVALUE 
WEIGHT N_INFORMATIVE
LABEL TotalSampleSize as N_INFORMATIVE

# Quant trait file:
PROCESSFILE /proj/epi/CVDGeneNas/xinruo/page_hl_htn/01_page/06_studyQC/${trait}/CLEANED.${trait}_w.invn.MEGA_HA.out.metal.txt






# === THE 2nd (BioMe_OMNI) INPUT FILE HAS THE SAME FORMAT AND CAN BE PROCESSED IMMEDIATELY ===

# Quant trait file:
PROCESSFILE /proj/epi/CVDGeneNas/xinruo/page_hl_htn/01_page/06_studyQC/${trait}/CLEANED.${trait}_w.invn.BioMe_OMNI_HA.out.metal.txt





# === DESCRIBE AND PROCESS THE 3rd (CCHC) INPUT FILE ===
SEPARATOR TAB
MARKER MarkerID
ALLELE Allele2 Allele1  # Test allele first
FREQ AF_Allele2
EFFECT BETA
STDERR SE 
PVALUE p.value 
WEIGHT N
LABEL TotalSampleSize as N

# Quant trait file:
#   ***The cleaned CCHC SBP and DBP files were renamed to remove "_adj" from their filenames so the code can run properly.***
PROCESS /proj/epi/CVDGeneNas/xinruo/page_hl_htn/02_cchc/step3_easystrata_output/CLEANED.${trait}_w.invn_allchr.txt




# === DESCRIBE AND PROCESS THE 4th (SLS) INPUT FILE ===
SEPARATOR TAB
MARKER VCF_ID
ALLELE ALT REF  # Test allele first
FREQ ALT_AF
EFFECT BETA
STDERR SE 
PVALUE PVALUE 
WEIGHT N_INFORMATIVE
LABEL TotalSampleSize as N_INFORMATIVE

# Quant trait file:
PROCESS /proj/epi/CVDGeneNas/xinruo/page_hl_htn/03_SLS_22yr/step1_gwas_results/gwas_results_combined/EasyStrata/${trait}/CLEANED.${trait}_allchr.txt





# === DESCRIBE AND PROCESS THE 5th (MVP) INPUT FILE ===
SEPARATOR TAB
MARKER VCF_ID
ALLELE alt ref  # Test allele first
FREQ af
EFFECT beta
STDERR sebeta 
PVALUE pval 
WEIGHT num_samples
LABEL TotalSampleSize as num_samples

# Quant trait file:
PROCESS /proj/epi/CVDGeneNas/xinruo/page_hl_htn/05_mvp_summary_stats/MVP_R4_${trait}_HIS_VCFID.txt









# === METAL output file ===
OUTFILE ${mydir}/metal_res/${trait}.meta.out .tbl
ANALYZE HETEROGENEITY
CLEAR

QUIT



"> ${mydir}/metal_jobs/metal.${trait}.sh


## Submitting METAL
echo "module add metal

sbatch -o ${mydir}/metal_jobs/${trait}.slurmlog --mem=\"40GB\" -t 8:00:00 --wrap=\"metal ${mydir}/metal_jobs/metal.${trait}.sh\" --job-name=\"${trait}.metal\"
" >> submitMETAL_${trait}.sh


done


