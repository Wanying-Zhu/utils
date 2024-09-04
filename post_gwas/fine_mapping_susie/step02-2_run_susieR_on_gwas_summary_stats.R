# Run susieR using GWAS summary stats
# Works in R4.4.1 on Below lab server
# R3.6 had some issue, but might be solved with newer version of susie
# Check HH's code: /vgipiper04/CCHC/local_ancestry/phaser/mixQTL/finemapping/SuSie/run_susie.R
# Official documentation: https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html

# Example call
# Rscript /belowshare/vumcshare/data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/code/utils/step02-2_run_susieR_on_gwas_summary_stats.R \
# --gwas_file "/data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/output/susieR_fine_mapping/gwas_results_by_region/TG-O-54:3-_[NL-18:1].chr22.region_1.filter_by_p0.05.merged_region" \
# --cols CHR SNP POS N BETA SE P \
# --phenotype_fn "/data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/lipid_species_pheno/TG-O-54:3-_[NL-18:1].pheno" \
# --ld_matrix_fn "/data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/output/ld_matrix/snps_filtered_by_pval0.05/TG-O-54:3-_[NL-18:1].chr22.region_1.ld" \
# --zscore_col z_score \
# --min_cs_corr 0.3

library(argparse)
library(susieR)
library(data.table)
print(paste('# susieR version:', packageVersion("susieR")))

# ############### Helper function ###############

# Process argument
parse_args <- function() {
    # Create a parser object
    parser <- argparse::ArgumentParser(description = "Parse command-line arguments")

    # Add arguments
    parser$add_argument("--output_path", help="Path to the outout files", default='./') # Omit required=TRUE if default value is set
    parser$add_argument("--output_prefix", help="Prefix to the outout files", default='output')
    parser$add_argument("--gwas_file",
                        help="Path to the gwas region file (sumstat), must contain columns specified by --cols flag",
                        required=TRUE)
    parser$add_argument("--cols", help="Columns of the GWAS region file. At least need sample size N; beta and se, or z score", nargs='*',
                        default=c('SNP', 'N', 'BETA', 'SE', 'P'))
    parser$add_argument("--beta_col", default='BETA',
                        help="Column name of the beta. Only needed if z score is not provided")
    parser$add_argument("--se_col", default='SE',
                        help="Column name of the standard error. Only needed if z score is not provided")
    parser$add_argument("--phenotype_fn",
                        help="Path to the phenotype file used in GWAS",
                        required=TRUE)
    parser$add_argument("--ld_matrix_fn",
                        help="Path to the ld matrix file created by plink",
                        required=TRUE)
    parser$add_argument("--sample_size_col", default='N',
                        help="Column name of the sample size in GWAS summary stats")
    parser$add_argument("--zscore_col", default='z_score',
                        help="Column name of the z score. Only needed if BETA and SE are not provided")
    parser$add_argument("--min_cs_corr", default = 0.5, type = "double",
                        help='min_abs_corr of susie_get_cs. A "purity" threshold for the CS. Any CS that contains a pair of variables with correlation less than this threshold will be filtered out and not reported.')

    # Parse the command-line arguments
    args <- parser$parse_args()
    return(args)
}
# ############### End of helper function ###############


# Parse the command-line arguments
args <- parse_args()

# Manually add '/' if outout_path does not end with '/'
if (substr(args$output_path, nchar(args$output_path), nchar(args$output_path)) != "/") {
    args$output_path = paste0(args$output_path, '/', sep='')
}

log_fn = paste(args$output_path, args$output_prefix, '.log', sep='')
logfile <- file(log_fn, open = "w")
sink(logfile, type = "output", split = TRUE)
sink(logfile, type = "message")
print('# Arguments used:')
print(args)

# Load GWAS summary stats file
print('# Load GWAS summary stats')
sumstat = data.table::fread(args$gwas_file)
print(dim(sumstat))

# Load phenotype used in GWAS
print('# Load phenotype fn, get variance(y)')
y = data.table::fread(args$phenotype_fn, header=FALSE)
print(dim(y))
y = y[, colnames(y)[ncol(y)], with=FALSE] # last column in the phenotype file contains phenotype measures
var_y = var(y)
var_y = as.numeric(var_y) # Need a single value

# Load LD matrix
print('# Load LD matrix')
ld_matrix = data.table::fread(args$ld_matrix_fn, header=FALSE)
print(dim(ld_matrix))

# Get number of samples
print('# Get number of samples')
n = as.integer(sumstat[1, args$sample_size_col, with=FALSE])
print(n)

cat('\n# Run susie_rss using summary stats\n')
# Check if z score is in GWAS summary stats
if (args$zscore_col %in% colnames(sumstat)) {
    print('# Use Z score in susie_rss')
    z_score = sumstat[[args$zscore_col]]
    # By default SuSiE assumes at most 10 causal variables, with L = 10
    # LD matrix needs to be a data matrix
    result = susieR::susie_rss(z=z_score, R=as.matrix(ld_matrix), n=n, var_y=var_y, L=10)
} else {
    print('# Use BETA and SE in susie_rss')
    beta = sumstat[[args$beta_col]]
    se = sumstat[[args$se_col]]
    result = susieR::susie_rss(R=as.matrix(ld_matrix), n=n, var_y=var_y, bhat=beta, shat=se, L=10)
}

# Get PIP and credible sets (95% and 99% cresdibe sets)
pips = susieR::susie_get_pip(result)
# Default min_abs_corr = 0.5
# Sex Xcorr = as.matrix(ld_matrix) so that susie_get_cs can filter based on min_abs_corr provided
credible_sets_95 = susieR::susie_get_cs(result, coverage = 0.95, Xcorr = as.matrix(ld_matrix),
                                        min_abs_corr = args$min_cs_corr)

# credible_sets_99 = susieR::susie_get_cs(result, coverage = 0.99,
#                                         min_abs_corr = args$min_cs_corr) # Default min_abs_corr = 0.5

# Add fine-mapping info to GWAS summary stat
sumstat[, PIP:=susie_get_pip(result)]
# Number of credible sets
n_cs = length(names(credible_sets_95$cs))


if (is.null(result$sets$cs)) {
    # Move on if no credible set to output
    print('# No credible set to output')
} else{
    # Else add credible set and PIP to GWAS summary stats and save
    print('# Add credible set and PIP to GWAS summary stats and save')
    for (x in 1:n_cs) {
        col_name = paste('95_credible_set_', x, sep='')
        element_name = paste('L', x, sep='') # Each credible sets are stored as L1, L2,...
        # Label SNPs in a given credible set to be 1, else 0
        sumstat[, (col_name):=0] # Initilize column to be all zeros
    
        # Indices of snps in current credible set
        # Must use [[]] not $, since element name is stored in a variable
        cs_index = credible_sets_95$cs[[element_name]]
        print(paste('# Number of SNPs in credible set ', element_name, ': ', length(cs_index), sep=''))
        sumstat[cs_index, (col_name):=1] # Initilize column to be all zeros
    }
    # Save GWAS summary stats with fine-mapping results
    print('# Save GWAS summary stats with fine-mapping results')
    output_fn = paste(args$output_path, args$output_prefix, '.susie.txt', sep='')
    print(paste('# -', output_fn))
    fwrite(sumstat, output_fn, sep='\t')
    
    # Print other info to log file in case needed
    cat('\n# ########## Other info in the fine-mapping result ##########\n')
    cat('# purity\n')
    print(result$sets$purity)
    cat('# coverage\n')
    print(result$sets$coverage)
    # print(names(result$sets))
}

cat('\n# DONE\n')
sink() # Close log file