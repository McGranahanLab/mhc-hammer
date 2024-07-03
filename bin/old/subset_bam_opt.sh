#!/bin/bash -eu

# Declare vars. 
bam=""
kmer_file=""
fish_reads=""
contig_reads=""
unmapped_reads=""
hla_coords=""
out_prefix=""
task_cpus=""
sort_mem=""
contig_file=""
only_fish=""

while getopts "b:k:f:c:u:h:p:t:m:d:o:" opt; do
    case "${opt}" in
        b)
            bam="${OPTARG}";
            ;;
        k)
            kmer_file="${OPTARG}";
            ;;   
        f)
            fish_reads="${OPTARG}";
            ;;                         
        c)
            contig_reads="${OPTARG}";
            ;; 
        u)
            unmapped_reads="${OPTARG}";
            ;; 
        h)
            hla_coords=$(cat "${OPTARG}");
            ;;
        p)
            out_prefix="${OPTARG}";
            ;;
        t)
            task_cpus="${OPTARG}";
            ;;                       
        m)
            sort_mem="${OPTARG}";
            ;; 
        d)
            contig_file="${OPTARG}";
            ;; 
        o)
            only_fish="${OPTARG}";
            ;; 
        \?)
            echo -e "Error: Invalid option: -${OPTARG}\n" >&2
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))


# Check if all required arguments are present except for kmer_file
if [[ -z "$bam" ]]; then
    echo "Error: Missing required argument: -b <bam>"
    exit 1
fi

if [[ -z "$fish_reads" ]]; then
    echo "Error: Missing required argument: -f <fish_reads>"
    exit 1
fi

if [[ -z "$contig_reads" ]]; then
    echo "Error: Missing required argument: -c <contig_reads>"
    exit 1
fi

if [[ -z "$unmapped_reads" ]]; then
    echo "Error: Missing required argument: -u <unmapped_reads>"
    exit 1
fi

if [[ -z "$hla_coords" ]]; then
    echo "Error: Missing required argument: -h <hla_coords>"
    exit 1
fi

if [[ -z "$out_prefix" ]]; then
    echo "Error: Missing required argument: -p <out_prefix>"
    exit 1
fi

if [[ -z "$task_cpus" ]]; then
    echo "Error: Missing required argument: -t <task_cpus>"
    exit 1
fi

if [[ -z "$sort_mem" ]]; then
    echo "Error: Missing required argument: -m <sort_mem>"
    exit 1
fi

if [[ -z "$contig_file" ]]; then
    echo "Error: Missing required argument: -d <contig_file>"
    exit 1
fi

# Check if kmer_file is required and present when fish_reads is set to true
if [[ "$fish_reads" == "true" ]] && [[ -z "$kmer_file" ]]; then
    echo "Error: kmer_file is required when fish_reads is set to true."
    exit 1
fi

# Check if kmer_file is required and present when only_fish is set to true
if [[ "$only_fish" == "true" ]] && [[ -z "$kmer_file" ]]; then
    echo "Error: kmer_file is required when only_fish is set to true."
    exit 1
fi

# create tmp dir 
tmp_dir="${out_prefix}_tmpDir"

echo creating temporary directory ${tmp_dir} - this will be deleted after the job finishes 
mkdir -p ${tmp_dir}

# assign output filenames 
out_hla_sam_name="${tmp_dir}/${out_prefix}.subset.sam"
uniq_out_hla_sam_name="${tmp_dir}/${out_prefix}.subset.uniq.sam"
out_hla_bam_name="${tmp_dir}/${out_prefix}.subset.bam"
tmp_sam_name="${tmp_dir}/${out_prefix}.tmp.sam"
read_names_path="${tmp_dir}/${out_prefix}.reads.tmp.txt"
read_names_unique_path="${tmp_dir}/${out_prefix}.reads.uniq.tmp.txt"
bam_header="${tmp_dir}/${out_prefix}.header"
sorted_out_hla_bam_name="${out_prefix}.subset.sorted.bam"
unmapped_read_count=0
region_read_count=0
kmer_read_count=0
contig_read_count=0

# Get the reads in the desired region and store them in a temporary BAM file
echo Extracting reads from ${bam}

# create file with read counts which will be appended to by each task
echo "read_origin,read_count" > "${out_prefix}.read_counts.csv"

# first, check the chromosome name in hla_coords is present in the bam file
# if not, exit with error
chr_name=$(echo ${hla_coords} | cut -d ":" -f 1)
if ! samtools view -H ${bam} | grep -q ${chr_name}; then
    echo "Error: chromosome name ${chr_name} not found in bam file ${bam}"
    exit 1
fi

echo Extracting reads aligning to ${hla_coords}
samtools view -@ ${task_cpus} -bh -o ${tmp_dir}/region.bam ${bam} ${hla_coords}
region_read_count=$(samtools view -@ ${task_cpus} -c ${tmp_dir}/region.bam)
echo hla_region,${region_read_count} >> "${out_prefix}.read_counts.csv"

# only extract kmer reads if only_fish is set to true 
if [ ${only_fish} == true ]; then
    echo "Extracting reads that contain HLA reference kmers"
    # Extract header from input BAM file
    samtools view -@ ${task_cpus} -H ${bam} > ${tmp_dir}/header.sam
    # Filter reads containing HLA reference kmers and concatenate with header
    samtools view -@ ${task_cpus} ${bam} | grep -F -f ${kmer_file} | cat ${tmp_dir}/header.sam - | samtools view -@ ${task_cpus} -bS - > ${tmp_dir}/kmer.bam
    
    kmer_read_count=$(samtools view -@ ${task_cpus} -c ${tmp_dir}/kmer.bam)

    echo fish_kmers,${kmer_read_count} >> "${out_prefix}.read_counts.csv"
    
    samtools merge -c -p -@ ${task_cpus} ${tmp_dir}/merged_input.bam ${tmp_dir}/region.bam ${tmp_dir}/kmer.bam
    # remove kmer.bam and input.bam and header.sam
    rm ${tmp_dir}/region.bam ${tmp_dir}/kmer.bam ${tmp_dir}/header.sam
    # rename merged_input.bam to input.bam
    mv ${tmp_dir}/merged_input.bam ${tmp_dir}/input.bam

else
    # Extract unmapped reads if params == true
    if [ ${unmapped_reads} == true ]; then
        echo "Extracting unmapped reads" 
        samtools view -@ ${task_cpus} -bf 4 -o ${tmp_dir}/unmapped.bam ${bam}
        unmapped_read_count=$(samtools view -@ ${task_cpus} -c ${tmp_dir}/unmapped.bam)
        echo unmapped_reads,${unmapped_read_count} >> "${out_prefix}.read_counts.csv"
    fi

    # Merge region.bam and unmapped.bam if both exist, otherwise use the appropriate file as input for the next steps
    if [ ${unmapped_reads} == true ]; then
        samtools merge -c -p -@ ${task_cpus} ${tmp_dir}/input.bam ${tmp_dir}/region.bam ${tmp_dir}/unmapped.bam
        # remove temporary files
        rm ${tmp_dir}/region.bam ${tmp_dir}/unmapped.bam
    else
        mv ${tmp_dir}/region.bam ${tmp_dir}/input.bam
    fi

    # Extract reads aligned to contigs if params == true
    if [ ${contig_reads} == true ]; then
        echo "Extracting reads that map to contigs"
       
        # Generate a list of standard chromosome names (1-22, MT, X, Y, M) with and without "chr" prefix
        chroms=$(echo {1..22} MT X Y M | awk '{for (i=1; i<=NF; i++) print $i "\nchr" $i}')

        # Extract the header from the input BAM file
        samtools view -H ${bam} > ${bam_header}
        
        # Extract the sequence names (contig names) from the BAM header and save them to a temporary file
        grep @SQ <(samtools view -H ${bam}) | cut -f 2 | cut -d ":" -f2 > ${tmp_dir}/contigs.txt

        # Count the number of standard chromosomes present in the contig list
        count=`grep -E '^chr[1-9]{1}$|chr1[0-9]{1}$|chr2[0-2]{1}$|^[0-9]{1}$|^1[0-9]{1}$|2[0-2]{1}$|^chr[XYMT]{1,2}$|^[XYMT]{1,2}$' ${tmp_dir}/contigs.txt | wc -l`
        
        # If no standard chromosomes are found in the contig list, print a warning and exit
        if [ ${count} == 0 ]; then
            echo "No standard chromosomes in the contigs of the BAM - quitting job"
            exit
        fi

        # check if the contigs_placeholder.txt file was provided and replace it with a list of contigs from the bam header if so
        if [[ "$contig_file" == "contigs_placeholder.txt" ]]; then
            echo "No contig file provided. generating contig list from bam header."

            # overwrite contig_file with a temporary file
            contig_file="${tmp_dir}/extract.txt"

            # Filter out standard chromosomes from the contig list, keeping only non-standard contigs in the output file
            grep -v -E '^chr[1-9]{1}$|chr1[0-9]{1}$|chr2[0-2]{1}$|^[0-9]{1}$|^1[0-9]{1}$|^2[0-2]{1}$|^chr[XYMT]{1,2}$|^[XYMT]{1,2}$' ${tmp_dir}/contigs.txt > ${contig_file}

        else 
            echo "Using contig file provided by user: ${contig_file}"
        fi   

        # Iterate through the non-standard contigs in the ${contig_file} file
        cat ${contig_file} | while read line || [[ -n $line ]];
        do  
            # If the contig is a standard chromosome, skip extraction
            if [[ "${chroms[*]}\n" =~ "${line}" ]]; then
                echo ${line} will not be extracted
            else
                # Extract reads mapping to the non-standard contig and append them to the contig.bam file
                echo Extracting reads mapping to ${line}
                samtools view -@ ${task_cpus} -bh -o ${tmp_dir}/${line}_tmp_contig.bam ${bam} ${line}
                # Count the number of reads in the temporary contig.bam file
                contig_read_count=$(samtools view -@ ${task_cpus} -c ${tmp_dir}/${line}_tmp_contig.bam)
                # Append the contig name and read count to the contig_read_counts.csv file
                echo ${line},${contig_read_count} >> "${out_prefix}.read_counts.csv"
            fi
        done

        tmp_contig_bams=$(ls ${tmp_dir}/*tmp_contig.bam)

        # Merge input.bam and tmp_contig.bams
        samtools merge -c -p -@ ${task_cpus} ${tmp_dir}/contigs.bam ${tmp_contig_bams}

        # remove temporary contig.bam files
        rm ${tmp_contig_bams}
        
        # Count the number of reads in the contigs.bam file
        contig_read_count=$(samtools view -@ ${task_cpus} -c ${tmp_dir}/contigs.bam)

        # Merge input.bam and contigs.bam
        samtools merge -c -p -@ ${task_cpus} ${tmp_dir}/merged_input.bam ${tmp_dir}/input.bam ${tmp_dir}/contigs.bam

        # remove input.bam and contigs.bam
        rm ${tmp_dir}/contigs.bam ${tmp_dir}/input.bam
        
        # rename merged_input.bam to input.bam
        mv ${tmp_dir}/merged_input.bam ${tmp_dir}/input.bam

    fi

    # Extract reads that contain HLA reference kmers if params == true 
    if [ ${fish_reads} == true ]; then
        echo "Extracting reads that contain HLA reference kmers"
        # Extract header from input BAM file
        samtools view -@ ${task_cpus} -H ${bam} > ${tmp_dir}/header.sam
        # Filter reads containing HLA reference kmers and concatenate with header
        samtools view -@ ${task_cpus} ${bam} | grep -F -f ${kmer_file} | cat ${tmp_dir}/header.sam - | samtools view -@ ${task_cpus} -bS - > ${tmp_dir}/kmer.bam
        
        kmer_read_count=$(samtools view -@ ${task_cpus} -c ${tmp_dir}/kmer.bam)

        # add kmer read count to read_counts.csv
        echo fish_kmers,${kmer_read_count} >> "${out_prefix}.read_counts.csv"
        
        samtools merge -c -p -@ ${task_cpus} ${tmp_dir}/merged_input.bam ${tmp_dir}/input.bam ${tmp_dir}/kmer.bam
        # remove kmer.bam and input.bam and header.sam
        #rm ${tmp_dir}/input.bam ${tmp_dir}/kmer.bam ${tmp_dir}/header.sam
        # rename merged_input.bam to input.bam
        mv ${tmp_dir}/merged_input.bam ${tmp_dir}/input.bam
    fi

fi

echo "Extracting unique reads from ${tmp_dir}/input.bam"

# some lines in the bam file may have been added multiple times, so we need to remove duplicates
# Extract header from input BAM file
samtools view -@ ${task_cpus} -H ${tmp_dir}/input.bam > ${tmp_dir}/header.sam
# Sort and get unique reads, concatenate with header, and convert to BAM
samtools view -@ ${task_cpus} ${tmp_dir}/input.bam | sort -S ${sort_mem} -T ${tmp_dir} -u --parallel=${task_cpus} | cat ${tmp_dir}/header.sam - | samtools view -@ ${task_cpus} -bS - > ${tmp_dir}/uniq_input.bam

# remove input.bam
rm ${tmp_dir}/input.bam

# sort subset bam
echo "Sorting subset bam"
samtools sort -@ ${task_cpus} -T ${tmp_dir}/samtools ${tmp_dir}/uniq_input.bam -o ${sorted_out_hla_bam_name}

# Remove remaining temporary files
rm -rf ${tmp_dir}

echo "Indexing bam"
samtools index ${sorted_out_hla_bam_name}

final_bam_read_count=$(samtools view -@ ${task_cpus} -c ${sorted_out_hla_bam_name})
echo final_bam_read_count,${final_bam_read_count} >> "${out_prefix}.read_counts.csv"


echo "Done! Here's a summary of reads coming from each step:"
echo "---------------------------------------------"
printf " %-55s %10s\n" "Number of reads extracted from ${hla_coords}:" "${region_read_count}"
printf " %-55s %10s\n" "Number of reads extracted from unmapped reads:" "${unmapped_read_count}"
printf " %-55s %10s\n" "Number of reads extracted from contigs:" "${contig_read_count}"
printf " %-55s %10s\n" "Number of reads extracted from reads containing HLA reference kmers:" "${kmer_read_count}"
printf " %-55s %10s\n" "Final number of unique reads in ${sorted_out_hla_bam_name}:" "${final_bam_read_count}"
echo "---------------------------------------------"