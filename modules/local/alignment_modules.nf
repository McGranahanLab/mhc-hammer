process NOVOALIGN {
    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'
    label 'process_long'
    label 'error_retry'

    input:
    tuple val(meta), path(fqs), path(personalised_reference), path(novoindex)

    output:
    tuple val(meta), path(personalised_reference), path("${meta.sample_id}_${meta.seq}.hla.rehead.bam*"), emit: make_hla_bam_input
    path ("versions.yml"), emit: versions

    script:
    if (fqs instanceof List) {
        if (fqs.size() == 2) {
            // Paired-end reads
            fq1 = fqs[0]
            fq2 = fqs[1]
            fq = ""  // Not used for paired-end
        } else if (fqs.size() == 1) {
            // Single-end reads as a List with one element
            fq = fqs[0]
            fq1 = ""
            fq2 = ""
        } else {
            throw new IllegalArgumentException("Unexpected List input format for fqs: ${fqs}")
        }
    } else if (fqs instanceof String || fqs instanceof nextflow.processor.TaskPath) {
        // Single-end read as a String or TaskPath
        fq = fqs
        fq1 = ""
        fq2 = ""
    } else {
        throw new IllegalArgumentException(
            "Unexpected input format for fqs: ${fqs} " +
            "(type: ${fqs.getClass().getName()}). " +
            "Expected either a List with 2 elements (for paired-end reads) " +
            "or a String/ List with 1 element (for single-end reads)."
        )
    }
    """
    # Create named pipes to avoid writing uncompressed FASTQs to disk
    if [ "${meta.paired_end}" = "true" ]; then
        mkfifo fq1_uncompressed fq2_uncompressed
        gzip -cdf ${fq1} > fq1_uncompressed &
        gzip -cdf ${fq2} > fq2_uncompressed &
    else
        mkfifo fq_uncompressed
        gzip -cdf ${fq} > fq_uncompressed &
    fi

    # Align reads using novoalign
    if [ "${meta.paired_end}" = "true" ]; then
        novoalign -d ${novoindex} -f fq1_uncompressed fq2_uncompressed -F STDFQ -R 0 -r All 9999 -o SAM -o FullNW 1> ${meta.sample_id}.sam 2> ${meta.sample_id}.metrics
    else
        novoalign -d ${novoindex} -f fq_uncompressed -F STDFQ -R 0 -r All 9999 -o SAM -o FullNW 1> ${meta.sample_id}.sam 2> ${meta.sample_id}.metrics
    fi

    # Convert SAM to BAM
    samtools view -b -o ${meta.sample_id}.bam ${meta.sample_id}.sam

    # Sort BAM and remove duplicates
    samtools sort -o ${meta.sample_id}.sorted.bam ${meta.sample_id}.bam

    # Keep properly paired reads
    if [ "${meta.paired_end}" = "true" ]; then
        samtools view -f 2 -b -o ${meta.sample_id}.hla.bam ${meta.sample_id}.sorted.bam
    else
        cp ${meta.sample_id}.sorted.bam ${meta.sample_id}.hla.bam
    fi

    # Add sample ID to BAM header
    samtools addreplacerg -r ID:${meta.sample_id} -r SM:${meta.sample_id} -o ${meta.sample_id}_${meta.seq}.hla.rehead.bam ${meta.sample_id}.hla.bam
    samtools index ${meta.sample_id}_${meta.seq}.hla.rehead.bam  

    # Remove intermediate files
    rm ${meta.sample_id}.sam ${meta.sample_id}.sorted.bam ${meta.sample_id}.bam ${meta.sample_id}.hla.bam
  
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        novoalign: \$(novoalign --version | sed -e "s/V//g")
    END_VERSIONS
    """
}


process MAKE_HLA_ALLELE_BAMS {
    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'

    input:
    tuple val(meta), path(personalised_reference), path(hla_bam), val(aligner)

    output:
    tuple val(meta.patient_id), val(meta), \
    path("${meta.sample_id}_${meta.seq}_${aligner}.*sorted.filtered*"),                   emit: mosdepth_input
    tuple val(meta), path("${meta.sample_id}_${meta.seq}_${aligner}.*sorted.filtered*"),
    path(personalised_reference),                                                         emit: hla_allele_bams, optional: true
    path("${meta.sample_id}_${meta.seq}_${aligner}.hla_bam_read_count.csv"),                   emit: hla_bam_read_count
    tuple val(meta), path("${meta.sample_id}_passed_hla_genes.txt"),                      emit: passed_hla_genes, optional: true
    tuple val(meta), path("${meta.sample_id}_passed_heterozygous_hla_genes.txt"),         emit: passed_heterozygous_hla_genes, optional: true
    tuple val(meta), path("${meta.sample_id}_passed_hla_alleles.txt"),                    emit: passed_hla_alleles, optional: true
    tuple val(meta), path("${meta.sample_id}_passed_heterozygous_hla_alleles.txt"),       emit: passed_heterozygous_hla_alleles, optional: true
    
    path "versions.yml",                                                                  emit: versions

    script: 
    fasta = personalised_reference[0]
    """
    echo running script to produce allele-specific bams... 

    make_hla_bams.sh ${hla_bam[0]} ${projectDir}/bin/ \
        ${params.max_mismatch} ${meta.sample_id} ${meta.seq} ${aligner} ${meta.paired_end}

    touch ${meta.sample_id}_passed_heterozygous_hla_genes.txt
    touch ${meta.sample_id}_passed_heterozygous_hla_alleles.txt
    touch ${meta.sample_id}_passed_hla_genes.txt
    touch ${meta.sample_id}_passed_hla_alleles.txt
    
    # get list of hla_genes in the fasta file
    hla_genes=\$(grep '^>' ${fasta} | sed 's/^>//' | cut -d '_' -f 1,2 | sort -u)
    echo "Checking for genes that pass"
    for hla_gene in \${hla_genes}; do  
        echo \${hla_gene}
        alleles=(\$(grep '^>' ${fasta} | sed 's/^>//' | grep \${hla_gene} ))

        allele_count=\${#alleles[@]}
        echo \${allele_count}
        if [ \${allele_count} -eq "2" ]; then
            
            allele1_bam=${meta.sample_id}_${meta.seq}_${aligner}.\${alleles[0]}.sorted.filtered.bam
            allele2_bam=${meta.sample_id}_${meta.seq}_${aligner}.\${alleles[1]}.sorted.filtered.bam

            if [ -f \${allele1_bam} ] && [ -f \${allele2_bam} ] ; then
                echo \${hla_gene} >> ${meta.sample_id}_passed_heterozygous_hla_genes.txt
                echo \${hla_gene} >> ${meta.sample_id}_passed_hla_genes.txt
            fi

            if [ -f \${allele1_bam} ] ; then
                echo \${alleles[0]} >> ${meta.sample_id}_passed_heterozygous_hla_alleles.txt
            fi

            if [ -f \${allele2_bam} ] ; then
                echo \${alleles[1]} >> ${meta.sample_id}_passed_heterozygous_hla_alleles.txt
            fi


        elif [ \${allele_count} -eq "1" ]; then
            allele1_bam=${meta.sample_id}_${meta.seq}.\${alleles[0]}.sorted.filtered.bam
            if [ -f \${allele1_bam} ] ; then
                  echo \${hla_gene} >> ${meta.sample_id}_passed_hla_genes.txt
            fi
        fi
    done


    echo "Check for alleles that pass"
    # create file of hla alleles with non-empty bam files
    alleles=\$(grep ^'>' ${fasta} | sed 's/^>//' | sort -u)
    for hla_allele in \${alleles}; do
        
        bam_path=${meta.sample_id}_${meta.seq}_${aligner}.\${hla_allele}.sorted.filtered.bam
        echo \$bam_path
        if [ -f \${bam_path} ]; then
            echo "BAM exists!"
            echo \$hla_allele >> ${meta.sample_id}_passed_hla_alleles.txt
        fi
    done 

    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    RSAM_VERSION=\$(Rscript -e "cat(paste(packageVersion('Rsamtools'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        picard: \$(picard SamToFastq --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
        R: \${R_VERSION}
        Rsamtools: \${RSAM_VERSION}
    END_VERSIONS
    """
}

process STAR_ALIGN_FIRST_PASS {

    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_low'
    label 'process_long'
    label 'error_retry'

    input:
    tuple val(meta), path(bam), path(gtf), path(fasta), val(genome_size)
    path kmer_file

    output:
    path ("${meta.sample_id}_${meta.seq}.subset.sortedSJ.out.tab"), emit: splicing_table

    tuple val(meta.patient_id), \
         val(meta), \
         path ("${meta.sample_id}_${meta.seq}.subset.sortedSJ.out.tab"), emit: kmer_input

    tuple val(meta.sample_id), \
          val(meta), \
          path(bam), \
          path("${meta.sample_id}_reference"), \
          path(fasta), \
          path(gtf), \
          val(genome_size), \
          env(overhang)                                               , emit : second_pass_input

    path ("versions.yml")                                         , emit: versions

    script: 
    """
    mkdir ${meta.sample_id}_reference  

    # convert bams to fastq and align to reference

    fq_name=${meta.sample_id}_${meta.seq}.subset.sorted

    # subset bam
    small_sam=${meta.sample_id}_${meta.seq}_small.sam
    samtools view -H ${bam[0]} > \${small_sam}

    samtools view ${bam[0]} | grep -F -f ${kmer_file} >> \${small_sam}

    # generate fastqs
    if [ "${meta.paired_end}" = "true" ]; then
        samtools collate -u -O \${small_sam} | \
        samtools fastq -1 \${fq_name}.1.fq.gz \
                       -2 \${fq_name}.2.fq.gz \
                       -s /dev/null \
                       -0 /dev/null \
                       -n
    else
        samtools collate -u -O \${small_sam} | \
        samtools fastq -0 \${fq_name}.fq.gz \
                       -n
    fi

    # rm sam
    rm \${small_sam}

    # run fastqc
    if [ "${meta.paired_end}" = "true" ]; then
        fastqc \${fq_name}.1.fq.gz \${fq_name}.2.fq.gz
    else
        fastqc \${fq_name}.fq.gz
    fi

    unzip \\*zip

    # Get max read length from fastqc report
    if [ "${meta.paired_end}" = "true" ]; then
        cd \${fq_name}.1_fastqc
    else
        cd \${fq_name}_fastqc
    fi
    read_len=\$(grep length fastqc_data.txt | awk '{print \$NF}' | cut -f2 -d-)
    overhang=\$(( \${read_len} - 1 ))

    cd ../

    # Generate STAR ref
    STAR \
    --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --genomeDir ${meta.sample_id}_reference \
    --genomeFastaFiles ${fasta} \
    --sjdbGTFfile ${gtf} \
    --sjdbOverhang \${overhang} \
    --genomeSAindexNbases ${genome_size}

    if [ "${meta.paired_end}" = "true" ]; then
        read_files_in="\${fq_name}.1.fq.gz \${fq_name}.2.fq.gz"
    else
        read_files_in="\${fq_name}.fq.gz"
    fi

    # STAR first run 
    STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${meta.sample_id}_reference \
    --readFilesIn \${read_files_in}  \
    --readFilesCommand gunzip -c \
    --outSAMunmapped None \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --alignIntronMin 20  \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 3 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 2 \
    --outFileNamePrefix \${fq_name} \
    --outSAMtype None

    # remove fastq files
    if [ "${meta.paired_end}" = "true" ]; then
        rm \${fq_name}.1.fq.gz \${fq_name}.2.fq.gz
    else
        rm \${fq_name}.fq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
} 

process GENERATE_ALL_SPLICE_TABLE {

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'
    label 'error_retry'

    input:
    path(splicing_tables)

    output:
    path("all_splice_sites.tab"), emit: cohort_splice_table

    script:
    """
    cat ${splicing_tables} | awk '(\$5 > 0 && \$7 >= ${params.uniq_num_across_junc} && \$6==0)' | cut -f1-6 | sort | uniq > all_splice_sites.tab
    """
}

process ALT_SPLICE_KMER {

    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'
    label 'error_retry'

    input:
    tuple val(meta), path(sj_tab), path(gtf), path(fasta)
    path(kmer_file)

    output:
    tuple val(meta.sample_id), path ("*_kmers.txt"), emit: sample_kmer
    path ("versions.yml"), emit: versions

    script:
    """
    # Run Rscript to detect novel splice junctions
    Rscript ${projectDir}/bin/make_region_kmer_file.R \
        --region_sj_table ${sj_tab} \
        --hla_genome_fasta_path ${fasta} \
        --gtf ${gtf} \
        --sample_id ${meta.sample_id} \
        --original_kmer ${kmer_file}

    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}      
    END_VERSIONS
    """

}

process STAR_ALIGN_SECOND_PASS {

    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_low'
    label 'process_long'

    input:
    tuple val(meta), path(bam), path(star_reference), path(fasta), path(gtf), val(genome_size), val(overhang), path(kmer)
    path(cohort_splice_table)

    output:
    tuple val(meta.patient_id), val(meta), \
          path("${meta.sample_id}.subset.sorted.bam*")      , emit: star_aligned_bams
    tuple val(meta.patient_id), val(meta), \
          path("${meta.sample_id}.subset.sorted.SJ.out.tab"), emit: sj_table
    path ("versions.yml")                                   , emit: versions

    script:
    """
    # First need to subset the cohort splicing table so it only contains HLA alleles matching this patient
    awk 'NR==FNR{a[\$0]; next} \$1 in a' ${star_reference}/chrName.txt all_splice_sites.tab > ${meta.patient_id}.splice.tab

    # subset bam
    small_sam=${meta.sample_id}_${meta.seq}_small.sam
    samtools view -H ${bam[0]} > \${small_sam}

    samtools view ${bam[0]} | grep -F -f ${kmer} >> \${small_sam}

    # generate fastqs
    if [ "${meta.paired_end}" = "true" ]; then
        samtools collate -u -O \${small_sam} | \
        samtools fastq -1 ${meta.sample_id}.1.fq.gz \
                       -2 ${meta.sample_id}.2.fq.gz \
                       -s /dev/null \
                       -0 /dev/null \
                       -n
    else
        samtools collate -u -O \${small_sam} | \
        samtools fastq -0 ${meta.sample_id}.fq.gz \
                       -n
    fi

    # rm sam
    rm \${small_sam}

    # This script runs star, and if it errors with code 139 reruns with a larger genome
    
    star_second_pass.sh \
    -f ${meta.sample_id} \
    -t ${task.cpus} \
    -s ${meta.sample_id} \
    -j ${meta.patient_id}.splice.tab \
    -b ${task.memory.toBytes()} \
    -r ${meta.sample_id}_reference \
    -p ${meta.paired_end} \
    -l ${genome_size} \
    -a ${fasta} \
    -g ${gtf} \
    -o ${overhang}

    # rename
    mv ${meta.sample_id}.subset.sorted.Aligned.sortedByCoord.out.bam ${meta.sample_id}.subset.sorted.bam

    # index
    samtools index ${meta.sample_id}.subset.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

process MOSDEPTH {

    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bed), val(aligner)

    output:
    path  "${meta.sample_id}.${meta.seq}.${aligner}.mosdepth.bed", emit: mosdepth_bed
    path("versions.yml"), emit: versions

    script:
    """
    
    for file in *bam ; do
        echo "Runninng mosdepth for" \$file

        # get allele name
        allele=\${file##*_novoalign.}   
        allele=\${allele##*_star.}   
        allele=\${allele%%.sorted.filtered.bam}

        mosdepth --no-per-base --chrom \$allele --by ${bed} ${meta.sample_id} \$file
        gunzip ${meta.sample_id}".regions.bed.gz"
        mv ${meta.sample_id}.regions.bed ${meta.sample_id}.\$allele.${meta.seq}.${aligner}.mosdepth.bed
    done

    cat *.mosdepth.bed > ${meta.sample_id}.${meta.seq}.${aligner}.mosdepth.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """

}