#!/bin/bash -eu

# Declare vars. 
fq_prefix=""
task_cpus=""
sample_name=""
splice_junction_path=""
sort_ram=""
reference_dir=""
original_genome_size=""
fasta_file=""
gtf_file=""
star_overhang=""

while getopts "f:t:s:j:b:r:l:a:g:o:" opt; do
    case "${opt}" in
        f)
            fq_prefix="${OPTARG}";
            ;;
        t)
            task_cpus="${OPTARG}";
            ;;                       
        s)
            sample_name="${OPTARG}";
            ;; 
        j)
            splice_junction_path="${OPTARG}";
            ;; 
        b)
            sort_ram="${OPTARG}";
            ;; 
        r)
            reference_dir="${OPTARG}";
            ;;          
        l)
            original_genome_size="${OPTARG}";
            ;;                         
        a)
            fasta_file="${OPTARG}";
            ;;    
        g)
            gtf_file="${OPTARG}";
            ;;     
        o)
            star_overhang="${OPTARG}";
            ;;                                                             
        \?)
            echo -e "Error: Invalid option: -${OPTARG}\n" >&2
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))


# Check if all required arguments are present 
if [[ -z "$fq_prefix" ]]; then
    echo "Error: Missing required argument: -f <f>"
    exit 1
fi

if [[ -z "$task_cpus" ]]; then
    echo "Error: Missing required argument: -t <task_cpus>"
    exit 1
fi

if [[ -z "$sample_name" ]]; then
    echo "Error: Missing required argument: -s <sample_name>"
    exit 1
fi

if [[ -z "$splice_junction_path" ]]; then
    echo "Error: Missing required argument: -j <splice_junction_path>"
    exit 1
fi

if [[ -z "$sort_ram" ]]; then
    echo "Error: Missing required argument: -b <sort_ram>"
    exit 1
fi

if [[ -z "$reference_dir" ]]; then
    echo "Error: Missing required argument: -r <reference_dir>"
    exit 1
fi

if [[ -z "$original_genome_size" ]]; then
    echo "Error: Missing required argument: -l <original_genome_size>"
    exit 1
fi

if [[ -z "$fasta_file" ]]; then
    echo "Error: Missing required argument: -l <fasta_file>"
    exit 1
fi

if [[ -z "$gtf_file" ]]; then
    echo "Error: Missing required argument: -g <gtf_file>"
    exit 1
fi

if [[ -z "$star_overhang" ]]; then
    echo "Error: Missing required argument: -o <star_overhang>"
    exit 1
fi

echo fq_prefix=$fq_prefix
echo task_cpus=$task_cpus
echo sample_name=$sample_name
echo splice_junction_path=$splice_junction_path
echo sort_ram=$sort_ram
echo reference_dir=$reference_dir
echo original_genome_size=$original_genome_size
echo fasta_file=$fasta_file
echo gtf_file=$gtf_file
echo star_overhang=$star_overhang
echo PWD=$PWD
echo ""

echo "Running STAR"
set +e
# STAR second align
    STAR \
    --runThreadN $task_cpus \
    --genomeDir $reference_dir \
    --readFilesIn $fq_prefix".1.fq.gz" $fq_prefix".2.fq.gz" \
    --readFilesCommand gunzip -c \
    --outSAMunmapped None \
    --outFilterType BySJout \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNmax 1 \
    --alignSJoverhangMin 3 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --outFileNamePrefix $sample_name".subset.sorted." \
    --outSAMattributes All \
    --sjdbFileChrStartEnd $splice_junction_path \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM $sort_ram 

star_error_code=$?
set -e


if [[ "$star_error_code" -eq 139 ]]; then
    echo "STAR had 139 error, trying with bigger genome"

    next_genome_length=$(expr $original_genome_size + 2)
    new_genome_lengths=($(seq $next_genome_length 2 15))

    # make sure max is in array
    if [[ ! " ${new_genome_lengths[*]} " =~ " 15 " ]]; then
        new_genome_lengths+=(15)
    fi

    for new_genome_length in "${new_genome_lengths[@]}";do
        echo "Trying genome length: $new_genome_length"

        echo "Making new genome"

        STAR \
        --runThreadN $task_cpus \
        --runMode genomeGenerate \
        --genomeDir "star_genome_"$new_genome_length \
        --genomeFastaFiles $fasta_file \
        --sjdbGTFfile $gtf_file \
        --sjdbOverhang $star_overhang \
        --genomeSAindexNbases $new_genome_length

        echo "Running STAR"
        set +e
        # STAR second align
        STAR \
        --runThreadN $task_cpus \
        --genomeDir "star_genome_"$new_genome_length \
        --readFilesIn $fq_prefix".1.fq.gz" $fq_prefix".2.fq.gz" \
        --readFilesCommand gunzip -c \
        --outSAMunmapped None \
        --outFilterType BySJout \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNmax 1 \
        --alignSJoverhangMin 3 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --outFileNamePrefix $sample_name".subset.sorted." \
        --outSAMattributes All \
        --sjdbFileChrStartEnd $splice_junction_path \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM $sort_ram 

    star_error_code=$?
    set -e

    if [[ ! "$star_error_code" -eq 139 ]]; then
      echo "Quitting with error "$star_error_code
      exit $star_error_code
    fi

    done
    
fi

echo "Quitting with error "$star_error_code
exit $star_error_code
