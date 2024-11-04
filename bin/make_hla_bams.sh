#!/bin/bash -eu

hla_bam=$1
scripts_dir=$2
max_mismatch=$3
sample_id=$4
seq=$5
aligner=$6
paired_end=$7

echo "hla_bam=$hla_bam"
echo "scripts_dir=$scripts_dir"
echo "max_mismatch=$max_mismatch"
echo "sample_id=$sample_id"
echo "seq=$seq"
echo "aligner=$aligner"
echo "paired_end=$paired_end"


# Get patient alleles to loop through
alleles=$(grep ^'>' *reference.fa | sed 's/^>*//' | sort -u )

# Loop through alleles
for i in ${alleles}; do
    echo "Generating allele-specific BAM for: ${i}"

    # Step 1: Get reads that align to the allele
    echo "Step 1: Getting reads that align to ${i}"
    samtools view -b -o "${sample_id}_${seq}_${aligner}.${i}.bam" "${hla_bam}" "${i}"
    echo "samtools view -b -o ${sample_id}_${seq}_${aligner}.${i}.bam ${hla_bam} ${i}"

    # Step 2: Sort allele-specific BAM
    echo "Step 2: Sorting allele-specific BAM"
    samtools sort "${sample_id}_${seq}_${aligner}.${i}.bam" -o "${sample_id}_${seq}_${aligner}.${i}.sorted.bam"

    # Step 3: Index allele-specific BAM
    echo "Step 3: Indexing allele-specific BAM"
    samtools index "${sample_id}_${seq}_${aligner}.${i}.sorted.bam"

    # Step 4: Get size of allele BAM before filtering
    echo "Step 4: Getting size of allele BAM before filtering"
    echo "bam="${sample_id}_${seq}_${aligner}.${i}.sorted.bam
    num_reads_in_unfiltered_allele_bam=$(samtools view -c -F 2304 "${sample_id}_${seq}_${aligner}.${i}.sorted.bam")
    echo "num_reads_in_unfiltered_allele_bam="$num_reads_in_unfiltered_allele_bam
    if [[ "${num_reads_in_unfiltered_allele_bam}" -eq "0" ]]; then
        echo "${sample_id}_${seq}_${aligner}.${i}.sorted.bam is empty - skipping to next allele"
        continue
    fi

    # Step 5: Filter reads with more than max_mismatch mismatches
    echo "Step 5: Filtering reads with more than ${max_mismatch} mismatches"
    Rscript "${scripts_dir}quantify_mismatches.R" \
    --bam_path "${sample_id}_${seq}_${aligner}.${i}.sorted.bam" \
    --max_mismatch "${max_mismatch}" \
    --allele "${i}" \
    --sample_name "${sample_id}_${seq}_${aligner}" \
    --paired_end "$paired_end"


    # Step 6: Generate final allele-specific BAM if passed reads file is not empty
    if [ -s "${sample_id}_${seq}_${aligner}.${i}.passed_reads.txt" ]; then
        echo "Step 6: Generating final allele-specific BAM for ${i}"
        picard FilterSamReads -I "${sample_id}_${seq}_${aligner}.${i}.sorted.bam" \
            -FILTER includeReadList \
            -RLF "${sample_id}_${seq}_${aligner}.${i}.passed_reads.txt" \
            -OUTPUT "${sample_id}_${seq}_${aligner}.${i}.sorted.filtered.bam"

        # Index final BAM
        samtools index "${sample_id}_${seq}_${aligner}.${i}.sorted.filtered.bam"

        # Add number of reads before and after filtering to the CSV
        num_reads_in_filtered_allele_bam=$(samtools view -c -F 2304 "${sample_id}_${seq}_${aligner}.${i}.sorted.filtered.bam")
        echo "${i},${num_reads_in_unfiltered_allele_bam},${num_reads_in_filtered_allele_bam}" >> "${sample_id}_${seq}_${aligner}.hla_bam_read_count.csv"

    else
    
        num_reads_in_filtered_allele_bam=0
        echo "${i},${num_reads_in_unfiltered_allele_bam},${num_reads_in_filtered_allele_bam}" >> "${sample_id}_${seq}_${aligner}.hla_bam_read_count.csv"

        echo "No reads passed mismatch filter for allele ${i} - No allele-specific BAM can be generated!"
    fi

    rm "${sample_id}_${seq}_${aligner}.${i}.bam" "${sample_id}_${seq}_${aligner}.${i}.sorted.bam" "${sample_id}_${seq}_${aligner}.${i}.sorted.bam.bai" 

done