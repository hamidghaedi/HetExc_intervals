#!/bin/bash
#SBATCH --job-name=gnomAD15
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2-00:00:00
#SBATCH --mem=50gb
#SBATCH --output=gnomad15.%J.out
#SBATCH --error=gnomad15.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qaedi65@gmail.com

#module load StdEnv/2020 gcc/9.3.0
module load vcflib/1.0.3 bcftools/1.16 bedtools/2.31.0

cd /home/ghaedi/projects/def-gooding-ab/ghaedi/gnomad

# Define the BED file
# the gtf2bed.ipynb shows the process of making this file
# Then used bed tools recommendation to sort the bed file: sort -k1,1 -k2,2n exome_coordinates.bed > in.sorted.bed
BED_FILE="in.sorted.bed"

# Additional parameters
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
output_dir="exome_result"
mkdir -p "$output_dir"
all_intervals_file="$output_dir/all_intervals_with_length.bed"

for chromosome in "${chromosomes[@]}"; do
    EXOME_FILE="exome_vcf/gnomad.exome.v4.0.sites.chr${chromosome}.vcf.bgz"
    
    if [ ! -f "$EXOME_FILE" ]; then
        echo "Downloading exome VCF file for chromosome $chromosome..."
        wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/$EXOME_FILE" -O "$EXOME_FILE"
    fi

    echo "Filtering exome VCF file for chromosome $chromosome based on inbreeding coefficient..."
    bcftools view -i 'INFO/inbreeding_coeff <= -0.3' "$EXOME_FILE" -Oz -o "$output_dir/chr${chromosome}_exome_filtered_by_InbreedCoef.vcf.gz"

    echo "Intersecting exome VCF file for chromosome $chromosome with exons..."
    bedtools intersect -a "$output_dir/chr${chromosome}_exome_filtered_by_InbreedCoef.vcf.gz" -b "$BED_FILE" -sorted -header > "$output_dir/chr${chromosome}_exome_intersected.vcf"

    echo "Converting exome intersected VCF file for chromosome $chromosome to TSV..."
    vcf2tsv "$output_dir/chr${chromosome}_exome_intersected.vcf" > "$output_dir/chr${chromosome}_exome_intersected.tsv"

    # Additional code for interval processing
    exome_file="$output_dir/chr${chromosome}_exome_intersected.vcf"
    bedtools merge -i "$exome_file" -d 300 -c 1 -o count > "$output_dir/merged_intervals_d300_${chromosome}.bed"
    awk '$4 != 1' "$output_dir/merged_intervals_d300_${chromosome}.bed" > "$output_dir/filtered_intervals_chr${chromosome}.bed"
    awk '{print $0, $3 - $2}' "$output_dir/filtered_intervals_chr${chromosome}.bed" > "$output_dir/intervals_with_length_chr${chromosome}.bed"
    cat "$output_dir/intervals_with_length_chr${chromosome}.bed" >> "$all_intervals_file"
    echo "Processing for exome chromosome $chromosome completed."
done

echo "All intervals with length have been appended to $all_intervals_file."
