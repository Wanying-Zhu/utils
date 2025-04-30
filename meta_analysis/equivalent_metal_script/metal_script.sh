### METAL script from Alice
# WZ notes:
# Check reference: https://genome.sph.umich.edu/wiki/METAL_Documentation
# METAl script can be run as: metal metal_script.sh
# This is the equivalent METAL settings of the python meta-analysis code

# This METAL script will run sample size-weighted meta-analysis across CCHC proteomics batch 1 an 2 adult samples

# To use this METAL script:
# (1) It is assumed that the column names are the same across all input files.
# (2) You need to use 'idmap_bridge_proteomics.txt' (pinned in slack channel) to get the unique ID 'OIDHT_OID3072' for each protein.
# (3) You need to create 2 fake allele columns:
#       - One is called 'dummy_ref' that has T (or anyone in A, T, G, C) for all rows.
#       - The other called 'dummy_alt' that has G (or anyone in A, T, G, C) for all rows.
# (4) Run this code as: metal metal_script.sh


# Leave these as they are:
SEPARATOR TAB           # The delimiter of your file; default is WHITESPACE, also accepts COMMA
SCHEME SAMPLESIZE       # Since we are running N-weighted
USESTRAND OFF
VERBOSE OFF


# Provide the following according to your input files:
MARKER protein    # Name of unique protein ID
ALLELE dummy_ref dummy_alt     # Names of the fake allele columns
PVALUE pval                # Name of P value column
EFFECT beta                # Name of beta column
WEIGHT sample_size               # Name of sample size column



# List the file paths for each input:
PROCESSFILE ../example_data/result1.txt
PROCESSFILE ../example_data/result2.txt 

# Leave the space before .tbl as it is
OUTFILE ../example_output/METAL_meta_output .tbl

ANALYZE HETEROGENEITY

CLEAR
QUIT


