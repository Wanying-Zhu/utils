# Modifed from an example provided by Alice and Misa (UNC)
# plink --bfile genotype_file \
# --extract SNP.extract \
# --threads 1 \
# --update-ref-allele SNP.allele1.txt 2 1 \ # Make sure alt allele is the same in LD calculation and GWAS
# --r square \ # make a matrix!
# --out output.LD \
# --make-just-bim


# Example call:
# code_path=/data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/code/utils
# project_path=/data100t1/home/wanying/CCHC/lipidomics
# bash ${code_path}/step02_get_ld_matrix.sh \
# ${project_path}/input_docs/lipidomic_sample_plink_imputed/chr11 \
# ${project_path}/post_gwas_finemapping/output/gwas_regions_split_for_plink/regions_merged/PC-16:0_20:4-.region_2.snp_list \
# /data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/output/susieR_fine_mapping/ld_matrix/all_snps
# PC-16:0_20:4-

genotype_file=$1 # Genotype file
snplist=$2 # SNPs in the region to be extracted from genotype file
output_path=$3 # Output path
output_prefix=$4 # output file prefix

cd ${output_path}
echo $PWD

# Check if output file already exists to avoid overwriting
if [ ! -f ${output_prefix}.ld ]
then
    echo "File not found, run plink!"
    # Outout file not found, need to run plink
    # Make sure alt allele is the same in LD calculation and GWAS using:
    # --update-ref-allele or --keep-allele-order
    # Becasue some version of plink flips allele when read in vcf.
    # Should be fine if GWAS was done using the same plink binary file.
    
    plink --bfile ${genotype_file} \
    --extract ${snplist} \
    --threads 8 \
    --keep-allele-order \
    --r square \
    --out ${output_prefix} \
    --make-just-bim # Save a bim file as are a reference of header and column names
else
    echo "Output LD matrix file exists. Skip plink LD calculation"
fi

