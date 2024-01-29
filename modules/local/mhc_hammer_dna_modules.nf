process DETECT_CN_AND_AIB {
    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'

    input:
    tuple val(patient_id), val(meta), path(tumour_bams), path(personalised_reference), path(normal_bams), \
          path(normal_flagstat), val(passed_heterozygous_hla_genes), val(purity_ploidy), path(tumour_flagstat), \
          path(snp_positions), path(gtf)

    val(min_depth)
    val(snp_type) 
    val(aligner)

    output: 
    path("*.pdf"), emit: plots, optional: true
    path("*_dna_analysis.csv"), emit: dna_analysis_output, optional: true
    path("*.coverage.csv"), optional: true

    path("versions.yml"), emit: versions

    script: 
    purity = purity_ploidy[0]
    ploidy = purity_ploidy[1]

    empty_purity_check = purity.isEmpty()
    empty_ploidy_check = ploidy.isEmpty()

    patient_fasta = personalised_reference[0]
    hla_genes_to_run = passed_heterozygous_hla_genes.collect().join(" ")
    """  
    echo "rerun"  
    for gene in ${hla_genes_to_run}
    do
        echo Doing \${gene}
        # Declare inputs
        hla_gene_alleles=\$(grep ^'>' ${patient_fasta} | sed 's/^>//' | grep \${gene} | sort -u)
        allele1=\$(echo \${hla_gene_alleles} | cut -f 1 -d ' ')
        allele2=\$(echo \${hla_gene_alleles} | cut -f 2 -d ' ')
        allele1_snp_path=\${allele1}_genome.snp_pos.bed
        allele2_snp_path=\${allele2}_genome.snp_pos.bed

        if [ ${snp_type} = "exon_snps" ]; then
            echo Filtering for exon snps
            Rscript ${baseDir}/bin/filter_for_gene_exon_snps.R \
            --gtf_path ${gtf} \
            --allele1_snp_path \$allele1_snp_path \
            --allele2_snp_path \$allele2_snp_path

            # update the snp paths to be the one that contains the exon snps
            allele1_snp_path=\${allele1}_genome.exon_snp_pos.bed
            allele2_snp_path=\${allele2}_genome.exon_snp_pos.bed
        fi
        
        allele1_tumour_filtered_bam=${meta.sample_id}_${meta.seq}_${aligner}.\${allele1}.sorted.filtered.bam
        allele1_gl_filtered_bam=${meta.normal_sample_name}_${meta.seq}_${aligner}.\${allele1}.sorted.filtered.bam
        
        allele2_tumour_filtered_bam=${meta.sample_id}_${meta.seq}_${aligner}.\${allele2}.sorted.filtered.bam
        allele2_gl_filtered_bam=${meta.normal_sample_name}_${meta.seq}_${aligner}.\${allele2}.sorted.filtered.bam

        echo Step 1: run mpileup on tumour and germline bams 

        allele1_tumour_mpileup=${meta.sample_id}_${meta.seq}.\${allele1}.mpileup
        allele1_gl_mpileup=${meta.normal_sample_name}_${meta.seq}.\${allele1}.mpileup

        allele2_tumour_mpileup=${meta.sample_id}_${meta.seq}.\${allele2}.mpileup
        allele2_gl_mpileup=${meta.normal_sample_name}_${meta.seq}.\${allele2}.mpileup

        samtools mpileup \$allele1_tumour_filtered_bam -f ${patient_fasta} > \$allele1_tumour_mpileup
        samtools mpileup \$allele1_gl_filtered_bam -f ${patient_fasta} > \$allele1_gl_mpileup
        samtools mpileup \$allele2_tumour_filtered_bam -f ${patient_fasta} > \$allele2_tumour_mpileup
        samtools mpileup \$allele2_gl_filtered_bam -f ${patient_fasta} > \$allele2_gl_mpileup
        
        echo Step 2: convert mpileup to coverage csvs
        allele1_tumour_coverage=${meta.sample_id}_${meta.seq}_${aligner}.\${allele1}.coverage.csv
        allele1_gl_coverage=${meta.normal_sample_name}_${meta.seq}_${aligner}.\${allele1}.coverage.csv

        allele2_tumour_coverage=${meta.sample_id}_${meta.seq}_${aligner}.\${allele2}.coverage.csv
        allele2_gl_coverage=${meta.normal_sample_name}_${meta.seq}_${aligner}.\${allele2}.coverage.csv

        Rscript --vanilla ${baseDir}/bin/mpileup_to_csv.R \
        --mpileup_path \$allele1_tumour_mpileup \
        --hla_ref_path ${patient_fasta} \
        --allele \$allele1 \
        --out_csv_name \$allele1_tumour_coverage

        Rscript --vanilla ${baseDir}/bin/mpileup_to_csv.R \
        --mpileup_path \$allele1_gl_mpileup \
        --hla_ref_path ${patient_fasta} \
        --allele \$allele1 \
        --out_csv_name \$allele1_gl_coverage

        Rscript --vanilla ${baseDir}/bin/mpileup_to_csv.R \
        --mpileup_path \$allele2_tumour_mpileup \
        --hla_ref_path ${patient_fasta} \
        --allele \$allele2 \
        --out_csv_name \$allele2_tumour_coverage

        Rscript --vanilla ${baseDir}/bin/mpileup_to_csv.R \
        --mpileup_path \$allele2_gl_mpileup \
        --hla_ref_path ${patient_fasta} \
        --allele \$allele2 \
        --out_csv_name \$allele2_gl_coverage

        echo Step 3: getting positions in germline that are greater than the minimum depth filter of ${min_depth} for each allele

        allele1_filtered_positions=${meta.normal_sample_name}_${meta.seq}.\${allele1}.filtered_positions.bed
        allele2_filtered_positions=${meta.normal_sample_name}_${meta.seq}.\${allele2}.filtered_positions.bed

        Rscript --vanilla ${baseDir}/bin/get_filtered_pos_bed.R \
        --coverage_path \$allele1_gl_coverage \
        --min_depth ${min_depth} \
        --out_bed_name \$allele1_filtered_positions \
        --allele \$allele1

        Rscript --vanilla ${baseDir}/bin/get_filtered_pos_bed.R \
        --coverage_path \$allele2_gl_coverage \
        --min_depth ${min_depth} \
        --out_bed_name \$allele2_filtered_positions \
        --allele \$allele2

        echo Step 4: generating coverage files that only contain positions that pass minimum depth in the germline
        
        allele1_tumour_coverage_at_filtered_positions=${meta.sample_id}_${meta.seq}.\${allele1}.coverage_at_filtered_positions.csv
        allele1_gl_coverage_at_filtered_positions=${meta.normal_sample_name}_${meta.seq}.\${allele1}.coverage_at_filtered_positions.csv
        allele2_tumour_coverage_at_filtered_positions=${meta.sample_id}_${meta.seq}.\${allele2}.coverage_at_filtered_positions.csv
        allele2_gl_coverage_at_filtered_positions=${meta.normal_sample_name}_${meta.seq}.\${allele2}.coverage_at_filtered_positions.csv
        
        Rscript --vanilla ${baseDir}/bin/get_filtered_coverage.R \
        --coverage_path \$allele1_tumour_coverage \
        --filtered_positions_path \$allele1_filtered_positions \
        --out_csv \$allele1_tumour_coverage_at_filtered_positions

        Rscript --vanilla ${baseDir}/bin/get_filtered_coverage.R \
        --coverage_path \$allele2_tumour_coverage \
        --filtered_positions_path \$allele2_filtered_positions \
        --out_csv \$allele2_tumour_coverage_at_filtered_positions

        Rscript --vanilla ${baseDir}/bin/get_filtered_coverage.R \
        --coverage_path \$allele1_gl_coverage \
        --filtered_positions_path \$allele1_filtered_positions \
        --out_csv \$allele1_gl_coverage_at_filtered_positions

        Rscript --vanilla ${baseDir}/bin/get_filtered_coverage.R \
        --coverage_path \$allele2_gl_coverage \
        --filtered_positions_path \$allele2_filtered_positions \
        --out_csv \$allele2_gl_coverage_at_filtered_positions

        echo Step 5: Getting snp positions that pass minimum depth in the germline in both alleles
    
        Rscript --vanilla ${baseDir}/bin/get_gene_filtered_snp_positions.R \
        --allele1_snp_bed \$allele1_snp_path \
        --allele2_snp_bed \$allele2_snp_path \
        --allele1 \$allele1 \
        --allele2 \$allele2 \
        --allele1_filtered_pos_bed \$allele1_filtered_positions \
        --allele2_filtered_pos_bed \$allele2_filtered_positions \
        --sample_name ${meta.normal_sample_name}_${meta.seq} 

        if [[ ${empty_purity_check} == "true" || ${empty_ploidy_check} == "true" ]]; then
            echo Either purity or ploidy information is missing for ${meta.sample_id}, cannot get allele specific CN without this!           
        else
            echo Step 6: Getting allele specific copy number 
            
            cn_output="${meta.sample_id}.\${gene}.${snp_type}.cn.csv"
            cn_plots_prefix="${meta.sample_id}.\${gene}.${snp_type}.cn"

            allele1_gene_filtered_snp_positions=${meta.normal_sample_name}_${meta.seq}.\${allele1}.filtered_snp_positions.bed
            allele2_gene_filtered_snp_positions=${meta.normal_sample_name}_${meta.seq}.\${allele2}.filtered_snp_positions.bed

            Rscript --vanilla ${baseDir}/bin/get_cn.R \
            --allele1 \$allele1 \
            --allele2 \$allele2 \
            --allele1_gl_coverage_file \$allele1_gl_coverage_at_filtered_positions \
            --allele2_gl_coverage_file \$allele2_gl_coverage_at_filtered_positions \
            --allele1_tumour_coverage_file \$allele1_tumour_coverage_at_filtered_positions \
            --allele2_tumour_coverage_file \$allele2_tumour_coverage_at_filtered_positions \
            --allele1_snp_bed \$allele1_gene_filtered_snp_positions \
            --allele2_snp_bed \$allele2_gene_filtered_snp_positions \
            --purity ${purity} \
            --ploidy ${ploidy} \
            --tumour_library_size_path ${tumour_flagstat} \
            --gl_library_size_path ${normal_flagstat} \
            --gtf_path ${gtf} \
            --cn_output_path \$cn_output \
            --cn_plots_prefix \$cn_plots_prefix \
            --scripts_dir ${projectDir}/bin/

            echo Step 7: Getting allele expected depth and generating plots
            expected_depth_output="${meta.sample_id}.\${gene}.${snp_type}.expected_depth.csv"
            expected_depth_plots_prefix="${meta.sample_id}.\${gene}.${snp_type}.expected_depth"

            Rscript --vanilla ${baseDir}/bin/get_expected_depth.R \
            --allele1 \$allele1 \
            --allele2 \$allele2 \
            --allele1_gl_coverage_file \$allele1_gl_coverage_at_filtered_positions \
            --allele2_gl_coverage_file \$allele2_gl_coverage_at_filtered_positions \
            --allele1_snp_bed \$allele1_gene_filtered_snp_positions \
            --allele2_snp_bed \$allele2_gene_filtered_snp_positions \
            --purity ${purity} \
            --tumour_library_size_path ${tumour_flagstat} \
            --gl_library_size_path ${normal_flagstat} \
            --expected_depth_output_path \$expected_depth_output \
            --expected_depth_plots_prefix \$expected_depth_plots_prefix \
            --scripts_dir ${projectDir}/bin/

        fi

        echo Step 8: Get coverage where reads only count once to a snp

        allele1_tumour_snp_reads_overlap=${meta.sample_id}_${meta.seq}.\${allele1}.snp_reads_overlap.bed
        allele2_tumour_snp_reads_overlap=${meta.sample_id}_${meta.seq}.\${allele2}.snp_reads_overlap.bed
        allele1_gl_snp_reads_overlap=${meta.normal_sample_name}_${meta.seq}.\${allele1}.snp_reads_overlap.bed
        allele2_gl_snp_reads_overlap=${meta.normal_sample_name}_${meta.seq}.\${allele2}.snp_reads_overlap.bed

        bedtools intersect -loj -bed -b \$allele1_tumour_filtered_bam \
                        -a \$allele1_gene_filtered_snp_positions \
                        > \$allele1_tumour_snp_reads_overlap

        bedtools intersect -loj -bed -b \$allele2_tumour_filtered_bam \
                        -a \$allele2_gene_filtered_snp_positions \
                        > \$allele2_tumour_snp_reads_overlap

        bedtools intersect -loj -bed -b \$allele1_gl_filtered_bam \
                        -a \$allele1_gene_filtered_snp_positions \
                        > \$allele1_gl_snp_reads_overlap

        bedtools intersect -loj -bed -b \$allele2_gl_filtered_bam \
                        -a \$allele2_gene_filtered_snp_positions \
                        > \$allele2_gl_snp_reads_overlap
        
        # generally, there will be 9 columns, but if there is a NULL event,
        # i.e. no reads overlap a position, then there will be 15 columns
        # we only want the first 7 columns
        # first save first 7 columns, so we can use fread which is a lot faster when
        # the tables get big

        allele1_tumour_snp_reads_overlap_7cols=${meta.sample_id}_${meta.seq}.\${allele1}.snp_reads_overlap_7cols.bed
        allele2_tumour_snp_reads_overlap_7cols=${meta.sample_id}_${meta.seq}.\${allele2}.snp_reads_overlap_7cols.bed
        allele1_gl_snp_reads_overlap_7cols=${meta.normal_sample_name}_${meta.seq}.\${allele1}.snp_reads_overlap_7cols.bed
        allele2_gl_snp_reads_overlap_7cols=${meta.normal_sample_name}_${meta.seq}.\${allele2}.snp_reads_overlap_7cols.bed

        cut -f1-7 \$allele1_tumour_snp_reads_overlap > \$allele1_tumour_snp_reads_overlap_7cols
        cut -f1-7 \$allele2_tumour_snp_reads_overlap > \$allele2_tumour_snp_reads_overlap_7cols
        cut -f1-7 \$allele1_gl_snp_reads_overlap > \$allele1_gl_snp_reads_overlap_7cols
        cut -f1-7 \$allele2_gl_snp_reads_overlap > \$allele2_gl_snp_reads_overlap_7cols

        allele1_tumour_reads_count_once_coverage=${meta.sample_id}_${meta.seq}.\${allele1}.coverage_at_filtered_snps_reads_count_once.csv
        allele2_tumour_reads_count_once_coverage=${meta.sample_id}_${meta.seq}.\${allele2}.coverage_at_filtered_snps_reads_count_once.csv
        allele1_gl_reads_count_once_coverage=${meta.normal_sample_name}_${meta.seq}.\${allele1}.coverage_at_filtered_snps_reads_count_once.csv
        allele2_gl_reads_count_once_coverage=${meta.normal_sample_name}_${meta.seq}.\${allele2}.coverage_at_filtered_snps_reads_count_once.csv

        Rscript --vanilla ${baseDir}/bin/count_reads_once.R \
        --snp_reads_overlap_bed \$allele1_tumour_snp_reads_overlap_7cols \
        --snp_path \$allele1_snp_path \
        --out_csv \$allele1_tumour_reads_count_once_coverage

        Rscript --vanilla ${baseDir}/bin/count_reads_once.R \
        --snp_reads_overlap_bed \$allele2_tumour_snp_reads_overlap_7cols \
        --snp_path \$allele2_snp_path \
        --out_csv \$allele2_tumour_reads_count_once_coverage

        Rscript --vanilla ${baseDir}/bin/count_reads_once.R \
        --snp_reads_overlap_bed \$allele1_gl_snp_reads_overlap_7cols \
        --snp_path \$allele1_snp_path \
        --out_csv \$allele1_gl_reads_count_once_coverage

        Rscript --vanilla ${baseDir}/bin/count_reads_once.R \
        --snp_reads_overlap_bed \$allele2_gl_snp_reads_overlap_7cols \
        --snp_path \$allele2_snp_path \
        --out_csv \$allele2_gl_reads_count_once_coverage

        echo Step 11: Get AIB
        logr_aib_output="${meta.sample_id}.\${gene}.${snp_type}.logr_aib.csv"
        logr_aib_prefix="${meta.sample_id}.\${gene}.${snp_type}.logr_aib"

        Rscript --vanilla ${baseDir}/bin/get_logr_aib.R \
        --allele1 \$allele1 \
        --allele2 \$allele2 \
        --allele1_gl_reads_count_once_coverage \$allele1_gl_reads_count_once_coverage \
        --allele2_gl_reads_count_once_coverage \$allele2_gl_reads_count_once_coverage \
        --allele1_tumour_reads_count_once_coverage \$allele1_tumour_reads_count_once_coverage \
        --allele2_tumour_reads_count_once_coverage \$allele2_tumour_reads_count_once_coverage \
        --allele1_snp_bed \$allele1_gene_filtered_snp_positions \
        --allele2_snp_bed \$allele2_gene_filtered_snp_positions \
        --tumour_library_size_path ${tumour_flagstat} \
        --gl_library_size_path ${normal_flagstat} \
        --logr_aib_output_path \$logr_aib_output \
        --logr_aib_plots_prefix \$logr_aib_prefix \
        --scripts_dir ${projectDir}/bin/

    done

    Rscript --vanilla ${baseDir}/bin/concatenate_dna_analysis_tables.R \
    --genes ${hla_genes_to_run} \
    --snp_type ${snp_type} \
    --sample_name ${meta.sample_id} \
    --aligner ${aligner}

    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")
    SEQ_VERSION=\$(Rscript -e "cat(paste(packageVersion('seqinr'), collapse = ''))")
    GGPUBR_VERSION=\$(Rscript -e "cat(paste(packageVersion('ggpubr'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}
        seqinr: \${SEQ_VERSION}
        ggpubr: \${GGPUBR_VERSION}
    END_VERSIONS
    """
}


process DETECT_MUTS {
    tag "${meta.sample_id}"

    container "library://tpjones15/mhchammer/mhchammer_detect_muts:latest"

    label 'process_single'
 
    input:
    tuple val(patient_id), val(meta), path(tumour_bams), path(normal_bams), \
          val(passed_hla_alleles), path(personalised_reference), path(gtf)

    output: // for now - only saving the sample level vep tables for mutation depth script
    path("*.vcf")                           optional true
    tuple val(meta), path(normal_bams), \
                     path(tumour_bams), \
                     path("*vep.txt") ,     emit: mutation_output
    path "versions.yml"                    ,emit: versions

    script:
    patient_fasta = personalised_reference[1]
    hla_alleles_to_run = passed_hla_alleles.collect().join(" ")
    // get zipped_gtf from [zipped_gtf, zipped_gtf.tbi]
    zipped_gtf = gtf[0]
    """
    for allele in ${hla_alleles_to_run}
    do
        echo "Detecting mutations in \${allele}"
        normal_bam=${meta.normal_sample_name}_${meta.seq}.\${allele}.sorted.filtered.bam
        tumour_bam=${meta.sample_id}_${meta.seq}.\${allele}.sorted.filtered.bam

        gatk --java-options '-Xmx${task.memory.toGiga()}g -Xms1g' Mutect2 \
        -R ${patient_fasta} \
        -I \${normal_bam} \
        -I \${tumour_bam} \
        -normal ${meta.normal_sample_name} \
        --f1r2-tar-gz ${meta.sample_id}.\${allele}.f1r2.tar.gz \
        --output ${meta.sample_id}.\${allele}.vcf 
    
        if tail -n 1 ${meta.sample_id}.\${allele}.vcf | grep -q CHROM
        then 
            echo "No mutations detected for allele \${allele} skipping to next"
            continue
        else
            echo "mutations detected for allele \${allele} - Will now run FilterMutectCalls"
        fi

        gatk --java-options '-Xmx${task.memory.toGiga()}g -Xms1g' LearnReadOrientationModel \
        -I ${meta.sample_id}.\${allele}.f1r2.tar.gz \
        -O ${meta.sample_id}.\${allele}.read-orientation-model.tar.gz 


        gatk --java-options '-Xmx${task.memory.toGiga()}g -Xms1g' FilterMutectCalls \
        -V ${meta.sample_id}.\${allele}.vcf \
        -R ${patient_fasta} \
        --ob-priors ${meta.sample_id}.\${allele}.read-orientation-model.tar.gz \
        -O ${meta.sample_id}.\${allele}.filt.vcf 

        # split multi allelics first 
        bcftools norm -m-any ${meta.sample_id}.\${allele}.filt.vcf \
        --output ${meta.sample_id}.\${allele}.norm.filt.vcf

        # run VEP
        echo "running vep"
        
        vep -i ${meta.sample_id}.\${allele}.norm.filt.vcf --gtf ${zipped_gtf} --fasta ${patient_fasta} --vcf -o ${meta.sample_id}.\${allele}.norm.filt.vep.vcf \
        --fields "Allele,Consequence,IMPACT,Feature_type,Feature,EXON,INTRON,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,DISTANCE,STRAND,FLAGS"

        gatk --java-options '-Xmx${task.memory.toGiga()}g -Xms1g' VariantsToTable \
        -V ${meta.sample_id}.\${allele}.norm.filt.vep.vcf  \
        -O ${meta.sample_id}.\${allele}.norm.filt.vep.txt \
        -F CHROM -F POS -F REF -F ALT -F FILTER  -GF AD -GF DP -F CSQ   \
        --show-filtered

        # move the norm.filt.vcf to just filt.vcf 
        mv ${meta.sample_id}.\${allele}.norm.filt.vcf ${meta.sample_id}.\${allele}.filt.vcf
    done

    # check if no mutations were detected
    if [ ! -f ${meta.sample_id}.*.vep.txt ]
    then
        echo "No mutations detected for any alleles - creating empty vep file to prevent pipeline stalling"
        touch ${meta.sample_id}.empty.vep.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/Using.*\$//')
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//') 
    END_VERSIONS
    """
}

process PARSE_MUTATIONS {
    tag "${patient_id}"

    container "library://tpjones15/mhchammer/mhchammer_detect_muts:latest"

    label 'process_single'

    input:
    tuple val(patient_id), val(normal_ids), \
    path(normal_bams), path(tumour_bams), \
    path(vep_tables)

    path(inventory)

    output:
    path "${patient_id}_mutations.csv", emit: patient_mutations 

    path "versions.yml"               , emit: versions

    script: 
    // get list of normal samples 
    normal_samples = normal_ids.collect().unique().join(" ")
    // get rid of bai files from bam channel
    tumour_bam_files = tumour_bams.findAll{ it.getName().endsWith(".bam") }.collect().join(" ")
    // get list of non-empty vep tables 
    vep_tables = vep_tables.findAll{ !it.getName().endsWith("empty.vep.txt") }.collect().join(" ")
    """
    echo parsing mutation output

    Rscript ${baseDir}/bin/make_mutation_table.R \
            --vep_tables ${vep_tables} --wxs_tumour_bam_files ${tumour_bam_files} \
            --gl_sample_ids ${normal_samples} --mutation_save_path "${patient_id}_mutations.csv" \
            --scripts_dir ${projectDir}/bin/ --inventory ${inventory}

    # get R and library versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")
    SNV_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}
        deepSNV: \${SNV_VERSION}
    END_VERSIONS
    """
}