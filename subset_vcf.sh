# Subset vcf to only include ADVENTO sampless
# Usage: subset_vcf.sh <chromosome_number> <output_path> <sample_list>

# Refer to ~/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_test_vcfs/subset_vcfs.sh
chr_num=$1
output_dir=$2
sample_list=$3 # Sample list, one sampel perline without header
# sample_list=/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_training_vcfs/genome_maximum_independent_set_3rd_degree.txt

vcf=/data100t1/share/advento/brazil-cvd.${chr_num}.vcf.gz

output=${output_dir}/ADVENTO.chr${chr_num}.vcf

# Print chromosome, position, ref allele, the first alternate allele, AF field of INFO column and genotype
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%INFO/AF\t%FORMAT\n' -S ${sample_list} -o ${output} ${vcf}
echo "# Run bcftools on chr${chr_num}"
echo "# - File saved to ${output}"

bcftools view -o ${output} -S ${sample_list} ${vcf} # Include samples in the list
# bcftools view -o ${output} -S ^${sample_list} ${vcf} # Exclude samples in the list

echo "# bgzip output"
bgzip ${output}
echo "# tabix output"
tabix ${output}.gz