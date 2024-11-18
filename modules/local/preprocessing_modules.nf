process FLAGSTAT {

    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'
    
    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), \
          path("${meta.sample_id}_${meta.seq}.library_size.txt")              , emit: library_size
    path "${meta.sample_id}_${meta.seq}.flagstat"
    path "${meta.sample_id}_${meta.seq}.library_size_with_unmapped.txt"       , emit: library_size_with_unmapped
    path "${meta.sample_id}_${meta.seq}.library_size_without_unmapped.txt"    , emit: library_size_without_unmapped
    path "versions.yml"                                                       , emit: versions

    script: 
    """
    # run flagstat
    samtools flagstat ${bam_file[0]} > ${meta.sample_id}_${meta.seq}.flagstat

    if [ "${meta.paired_end}" = "true" ]; then
        samtools view -c -f 1 -F 2304 ${bam_file[0]} > ${meta.sample_id}_${meta.seq}.library_size_with_unmapped.txt
        samtools view -c -f 1 -F 2308 ${bam_file[0]} > ${meta.sample_id}_${meta.seq}.library_size_without_unmapped.txt
    else
        samtools view -c -F 2304 ${bam_file[0]} > ${meta.sample_id}_${meta.seq}.library_size_with_unmapped.txt
        samtools view -c -F 2308 ${bam_file[0]} > ${meta.sample_id}_${meta.seq}.library_size_without_unmapped.txt
    fi
    if [ ${params.include_unmapped_reads_in_library_size} == true ]; then
        cp ${meta.sample_id}_${meta.seq}.library_size_with_unmapped.txt ${meta.sample_id}_${meta.seq}.library_size.txt
    else
        cp ${meta.sample_id}_${meta.seq}.library_size_without_unmapped.txt ${meta.sample_id}_${meta.seq}.library_size.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

process SUBSET_BAMS {

    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_low'
    
    label 'error_retry'

    input:
    tuple val(meta), path(bam_file)
    path mhc_coords
    path contigs_file_ch

    output:
    tuple val(meta), path("${meta.sample_id}_${meta.seq}.subset.sorted.bam*") , emit: subset_bam
    path "${meta.sample_id}_${meta.seq}.read_counts.csv"                      , emit: read_counts
    path "versions.yml"                                                       , emit: versions

    script: 
    def args = task.ext.args ?: ''
    sort_mem = (task.memory.giga*0.8).intValue()
    """
    echo ${contigs_file_ch}
    # subset bam 
    subset_bam_opt.sh -b ${bam_file[0]} ${args} \
    -f ${params.fish_reads} -c ${params.contig_reads} \
    -d ${contigs_file_ch} \
    -u ${params.unmapped_reads} -h ${mhc_coords} \
    -p ${meta.sample_id}_${meta.seq} -t ${task.cpus} \
    -m ${sort_mem}G -o ${params.fish_reads_only} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

process GENERATE_HLA_FQS {

    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'

    label 'error_retry'

    input:
    tuple val(meta), path(subset_bam_file)

    output:
    tuple val(meta), path("${meta.sample_id}_${meta.seq}*.fq.gz")             , emit: fqs
    path "versions.yml"                                                       , emit: versions

    script: 
    """
    # Run samtools collate and fastq to generate paired fastq files
    if [ "${meta.paired_end}" = "true" ]; then
        samtools collate -u -O ${subset_bam_file[0]} | \
        samtools fastq -1 ${meta.sample_id}_${meta.seq}.1.fq.gz \
                       -2 ${meta.sample_id}_${meta.seq}.2.fq.gz \
                       -s /dev/null \
                       -0 /dev/null \
                       -n
    else
        samtools collate -u -O ${subset_bam_file[0]} | \
        samtools fastq -0 ${meta.sample_id}_${meta.seq}.fq.gz \
                       -n
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}