process DETECT_HLA_REPRESSION {
    
    tag "${meta.sample_id}"

    container "library://tpjones15/default/final_lohhla:latest"

    label 'process_single'

    input:
    tuple val(patient_id), val(meta), path(tumour_bams), path(personalised_reference), \
          path(normal_bams), path(normal_flagstat), val(passed_alleles), \
          path(tumour_flagstat), path(snp_positions), path(gtf)
    val(snp_type)          
    val(reference_type)
    val(aligner)

    output:
    path("*.pdf"), emit: plots, optional: true
    path("*_rna_repression.csv"), emit: hla_repression_tables, optional: true

    val(meta.patient_id), emit: patient_ids
    
    path("versions.yml"), emit: versions

    script: 
    patient_fasta = personalised_reference[0]

    hla_alleles_to_run = passed_alleles.collect().join(" ")
    """ 

    for allele in ${hla_alleles_to_run}; do
        echo Doing allele \${allele}
        
        # Declaring input files
        allele_snp_path=\${allele}_${reference_type}.snp_pos.bed
        allele_tumour_filtered_bam=${meta.sample_id}_${meta.seq}_${aligner}.\${allele}.sorted.filtered.bam
        allele_gl_filtered_bam=${meta.normal_sample_name}_${meta.seq}_${aligner}.\${allele}.sorted.filtered.bam

        if [ ${snp_type} = "exon_snps" ]; then
            echo Filtering for exon snps
            Rscript ${baseDir}/bin/filter_for_allele_exon_snps.R \
            --gtf_path ${gtf} \
            --snp_path \$allele_snp_path

            # update the allele snp path to be the one that contains the exon snps
            allele_snp_path=\${allele}_${reference_type}.exon_snp_pos.bed
        fi

        echo Step 1: Make sure reads only count onces towards a SNP

        allele_tumour_snp_reads_overlap=${meta.sample_id}_${meta.seq}.\${allele}.snp_reads_overlap.bed
        allele_gl_snp_reads_overlap=${meta.normal_sample_name}_${meta.seq}.\${allele}.snp_reads_overlap.bed

        allele_tumour_snp_reads_overlap_7cols=${meta.sample_id}_${meta.seq}.\${allele}.snp_reads_overlap_7cols.bed
        allele_gl_snp_reads_overlap_7cols=${meta.normal_sample_name}_${meta.seq}.\${allele}.snp_reads_overlap_7cols.bed

        allele_tumour_reads_count_once_coverage=${meta.sample_id}_${meta.seq}.\${allele}.coverage_at_filtered_snps_reads_count_once.csv
        allele_gl_reads_count_once_coverage=${meta.normal_sample_name}_${meta.seq}.\${allele}.coverage_at_filtered_snps_reads_count_once.csv

        bedtools intersect -loj -bed -b \$allele_tumour_filtered_bam \
                        -a \$allele_snp_path \
                        > \$allele_tumour_snp_reads_overlap

        bedtools intersect -loj -bed -b \$allele_gl_filtered_bam \
                        -a \$allele_snp_path \
                        > \$allele_gl_snp_reads_overlap
        
        # generally, there will be 9 columns, but if there is a NULL event,
        # i.e. no reads overlap a position, then there will be 15 columns
        # we only want the first 7 columns
        # first save first 7 columns, so we can use fread which is a lot faster when
        # the tables get big
        cut -f1-7 \$allele_tumour_snp_reads_overlap > \$allele_tumour_snp_reads_overlap_7cols
        cut -f1-7 \$allele_gl_snp_reads_overlap > \$allele_gl_snp_reads_overlap_7cols

        Rscript --vanilla ${baseDir}/bin/count_reads_once.R \
        --snp_reads_overlap_bed \$allele_tumour_snp_reads_overlap_7cols \
        --out_csv \$allele_tumour_reads_count_once_coverage

        Rscript --vanilla ${baseDir}/bin/count_reads_once.R \
        --snp_reads_overlap_bed \$allele_gl_snp_reads_overlap_7cols \
        --out_csv \$allele_gl_reads_count_once_coverage

        # get repression 
        echo Step 4: Compare tumour and normal samples to detect HLA allelic repression 
        Rscript --vanilla ${baseDir}/bin/get_rna_tumour_normal_comparison.R \
        --allele \$allele \
        --sample_name ${meta.sample_id} \
        --patient ${meta.patient_id} \
        --allele_snp_bed \$allele_snp_path \
        --allele_gl_coverage \$allele_gl_reads_count_once_coverage \
        --allele_tumour_coverage \$allele_tumour_reads_count_once_coverage \
        --tumour_library_size ${tumour_flagstat} \
        --gl_library_size ${normal_flagstat} \
        --output_path ${meta.sample_id}.\${allele}.${aligner}.${snp_type}.tumour_normal_comparison.csv \
        --plots_prefix ${meta.sample_id}.\${allele}.${aligner}.${snp_type}.tumour_normal_comparison \
        --scripts_dir ${projectDir}/bin/ 

    done

    Rscript --vanilla ${baseDir}/bin/concatenate_repression_tables.R \
    --alleles ${hla_alleles_to_run} \
    --snp_type ${snp_type} \
    --tumour_sample_name ${meta.sample_id} \
    --normal_sample_name ${meta.normal_sample_name} \
    --aligner ${aligner}

    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")
    GGBEE_VERSION=\$(Rscript -e "cat(paste(packageVersion('ggbeeswarm'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard FilterSamReads --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}
        ggbeeswarm: \${GGBEE_VERSION}
    END_VERSIONS
    """
}

process GET_HLA_ALLELIC_IMBALANCE {
    tag "${meta.sample_id}"

    container "library://tpjones15/default/final_lohhla:latest"

    label 'process_single'

    input:
    tuple val(meta), path(hla_bams), path(personalised_reference), \
          val(passed_heterozygous_genes), \
          path(single_allele_reads), \
          path(snp_positions), path(gtf)
    val(snp_type)          
    val(reference_type)
    val(aligner)          

    output: 
    path("*.pdf"),                         emit: plots, optional: true
    path("*rna_aib.csv"),                  emit: rna_aib_tables, optional: true
    val(meta.patient_id),                  emit: patient_ids
    path("versions.yml"),                  emit: versions

    script: 

    heterozygous_hla_genes = passed_heterozygous_genes.collect().join(" ")
    """
    for gene in ${heterozygous_hla_genes}; do

        echo Calculating RNA AIB for \${gene}

        # Declare inputs
        hla_gene_alleles=\$(grep ^'>' ${personalised_reference} | sed 's/^>//' | grep \${gene} | sort -u)
        allele1=\$(echo \${hla_gene_alleles} | cut -f 1 -d ' ')
        allele2=\$(echo \${hla_gene_alleles} | cut -f 2 -d ' ')
        allele1_snp_path=\${allele1}_${reference_type}.snp_pos.bed
        allele2_snp_path=\${allele2}_${reference_type}.snp_pos.bed

        if [ ${snp_type} = "exon_snps" ]; then
            echo Filtering for exon snps
            Rscript ${baseDir}/bin/filter_for_gene_exon_snps.R \
            --gtf_path ${gtf} \
            --allele1_snp_path \$allele1_snp_path \
            --allele2_snp_path \$allele2_snp_path

            # update the snp paths to be the one that contains the exon snps
            allele1_snp_path=\${allele1}_${reference_type}.exon_snp_pos.bed
            allele2_snp_path=\${allele2}_${reference_type}.exon_snp_pos.bed
        fi

        allele1_filtered_bam=${meta.sample_id}_${meta.seq}_${aligner}.\${allele1}.sorted.filtered.bam       
        allele2_filtered_bam=${meta.sample_id}_${meta.seq}_${aligner}.\${allele2}.sorted.filtered.bam

        # First calculate RNA AIB without considering reads that map to multiple alleles
        # Step 1 - If a read overlaps more than one SNP, make sure it is only counted once

        allele1_snp_reads_overlap=${meta.sample_id}_${meta.seq}.\${allele1}.snp_reads_overlap.bed
        allele2_snp_reads_overlap=${meta.sample_id}_${meta.seq}.\${allele2}.snp_reads_overlap.bed

        bedtools intersect -loj -bed -b \$allele1_filtered_bam \
                        -a \$allele1_snp_path \
                        > \$allele1_snp_reads_overlap

        bedtools intersect -loj -bed -b \$allele2_filtered_bam \
                        -a \$allele2_snp_path \
                        > \$allele2_snp_reads_overlap
        
        # generally, there will be 9 columns, but if there is a NULL event,
        # i.e. no reads overlap a position, then there will be 15 columns
        # we only want the first 7 columns
        # first save first 7 columns, so we can use fread which is a lot faster when
        # the tables get big
        
        allele1_snp_reads_overlap_7cols=${meta.sample_id}_${meta.seq}.\${allele1}.snp_reads_overlap_7cols.bed
        allele2_snp_reads_overlap_7cols=${meta.sample_id}_${meta.seq}.\${allele2}.snp_reads_overlap_7cols.bed

        cut -f1-7 \$allele1_snp_reads_overlap > \$allele1_snp_reads_overlap_7cols
        cut -f1-7 \$allele2_snp_reads_overlap > \$allele2_snp_reads_overlap_7cols

        allele1_reads_count_once_coverage=${meta.sample_id}_${meta.seq}.\${allele1}.coverage_at_filtered_snps_reads_count_once.csv
        allele2_reads_count_once_coverage=${meta.sample_id}_${meta.seq}.\${allele2}.coverage_at_filtered_snps_reads_count_once.csv
        
        Rscript --vanilla ${baseDir}/bin/count_reads_once.R \
        --snp_reads_overlap_bed \$allele1_snp_reads_overlap_7cols \
        --out_csv \$allele1_reads_count_once_coverage

        Rscript --vanilla ${baseDir}/bin/count_reads_once.R \
        --snp_reads_overlap_bed \$allele2_snp_reads_overlap_7cols \
        --out_csv \$allele2_reads_count_once_coverage

        # Step 2 - Calculate allelic imbalance
        
        aib_output="${meta.sample_id}.\${gene}.${aligner}.${snp_type}.aib.csv"
        aib_prefix="${meta.sample_id}.\${gene}.${aligner}.${snp_type}.aib"

        Rscript --vanilla ${baseDir}/bin/get_rna_aib.R \
        --patient ${meta.patient_id} \
        --sample_name ${meta.sample_id} \
        --allele1 \$allele1 \
        --allele2 \$allele2 \
        --allele1_reads_count_once_coverage \$allele1_reads_count_once_coverage \
        --allele2_reads_count_once_coverage \$allele2_reads_count_once_coverage \
        --allele1_snp_bed \$allele1_snp_path \
        --allele2_snp_bed \$allele2_snp_path \
        --aib_output_path \$aib_output \
        --aib_plots_prefix \$aib_prefix \
        --plot_title "RNA depth AIB" \
        --scripts_dir ${projectDir}/bin/

        # Now calculate RNA AIB where multimapping reads are removed
        
        # Step 1: making hla allele bams that only contain reads mapping to a single allele
        
        reads_mapping_single_allele_list=${meta.sample_id}_${aligner}_reads_mapping_one_allele.csv
        allele1_reads_mapping_single_allele_bam=${meta.sample_id}_${meta.seq}.\${allele1}.reads_mapping_single_allele.bam 
        allele2_reads_mapping_single_allele_bam=${meta.sample_id}_${meta.seq}.\${allele2}.reads_mapping_single_allele.bam 
        
        picard FilterSamReads \
        I=\$allele1_filtered_bam \
        O=\$allele1_reads_mapping_single_allele_bam \
        READ_LIST_FILE=\$reads_mapping_single_allele_list \
        FILTER=includeReadList 

        picard FilterSamReads \
        I=\$allele2_filtered_bam \
        O=\$allele2_reads_mapping_single_allele_bam \
        READ_LIST_FILE=\$reads_mapping_single_allele_list \
        FILTER=includeReadList 

        # Step 2: If a read overlaps more than one SNP, make sure it is only counted once
        
        allele1_snp_reads_overlap_single_allele=${meta.sample_id}_${meta.seq}.\${allele1}.snp_reads_overlap_single_allele.bed
        allele2_snp_reads_overlap_single_allele=${meta.sample_id}_${meta.seq}.\${allele2}.snp_reads_overlap_single_allele.bed

        bedtools intersect -loj -bed -b \$allele1_reads_mapping_single_allele_bam \
                        -a \$allele1_snp_path \
                        > \$allele1_snp_reads_overlap_single_allele

        bedtools intersect -loj -bed -b \$allele2_reads_mapping_single_allele_bam \
                        -a \$allele2_snp_path \
                        > \$allele2_snp_reads_overlap_single_allele
        
        # generally, there will be 9 columns, but if there is a NULL event,
        # i.e. no reads overlap a position, then there will be 15 columns
        # we only want the first 7 columns
        # first save first 7 columns, so we can use fread which is a lot faster when
        # the tables get big

        allele1_snp_reads_overlap_7cols_single_allele=${meta.sample_id}_${meta.seq}.\${allele1}.snp_reads_overlap_7cols_single_allele.bed
        allele2_snp_reads_overlap_7cols_single_allele=${meta.sample_id}_${meta.seq}.\${allele2}.snp_reads_overlap_7cols_single_allele.bed

        cut -f1-7 \$allele1_snp_reads_overlap_single_allele > \$allele1_snp_reads_overlap_7cols_single_allele
        cut -f1-7 \$allele2_snp_reads_overlap_single_allele > \$allele2_snp_reads_overlap_7cols_single_allele

        allele1_reads_count_once_coverage_single_allele=${meta.sample_id}_${meta.seq}.\${allele1}.coverage_at_filtered_snps_reads_count_once_single_allele.csv
        allele2_reads_count_once_coverage_single_allele=${meta.sample_id}_${meta.seq}.\${allele2}.coverage_at_filtered_snps_reads_count_once_single_allele.csv
        
        Rscript --vanilla ${baseDir}/bin/count_reads_once.R \
        --snp_reads_overlap_bed \$allele1_snp_reads_overlap_7cols_single_allele \
        --out_csv \$allele1_reads_count_once_coverage_single_allele

        Rscript --vanilla ${baseDir}/bin/count_reads_once.R \
        --snp_reads_overlap_bed \$allele2_snp_reads_overlap_7cols_single_allele \
        --out_csv \$allele2_reads_count_once_coverage_single_allele

        aib_output="${meta.sample_id}.\${gene}.${aligner}.${snp_type}.single_mapping.aib.csv"
        aib_prefix="${meta.sample_id}.\${gene}.${aligner}.${snp_type}.single_mapping.aib"

        Rscript --vanilla ${baseDir}/bin/get_rna_aib.R \
        --patient ${meta.patient_id} \
        --sample_name ${meta.sample_id} \
        --allele1 \$allele1 \
        --allele2 \$allele2 \
        --allele1_reads_count_once_coverage \$allele1_reads_count_once_coverage_single_allele \
        --allele2_reads_count_once_coverage \$allele2_reads_count_once_coverage_single_allele \
        --allele1_snp_bed \$allele1_snp_path \
        --allele2_snp_bed \$allele2_snp_path \
        --aib_output_path \$aib_output \
        --aib_plots_prefix \$aib_prefix \
        --plot_title "RNA depth AIB" \
        --scripts_dir ${projectDir}/bin/
    done

    Rscript --vanilla ${baseDir}/bin/concatenate_rna_aib_tables.R \
    --genes ${heterozygous_hla_genes} \
    --snp_type ${snp_type} \
    --sample_name ${meta.sample_id} \
    --aligner ${aligner}

    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")
    GGBEE=\$(Rscript -e "cat(paste(packageVersion('ggbeeswarm'), collapse = ''))")
    RSAM=\$(Rscript -e "cat(paste(packageVersion('Rsamtools'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard FilterSamReads --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}
        ggbeeswarm: \${GGBEE}
        Rsamtools: \${RSAM}
    END_VERSIONS
    """
}


process GET_HLA_ALLELE_EXPRESSION {
    tag "${meta.sample_id}"

    container "library://tpjones15/default/final_lohhla:latest"

    label 'process_single'

    input:
    tuple val(meta), path(hla_bams), path(personalised_reference), \
          path(flagstat), val(aligner)

    output: 
    tuple val(meta.sample_id), \
    path("*_reads_mapping_one_allele.csv"),             emit: tables_for_aib, optional: true

    path("*_rpkm.csv"),                                 emit: rpkm_table, optional: true

    val(meta.patient_id),                               emit: patient_ids
    
    path("versions.yml"),                               emit: versions

    """
    echo First, get reads that map to a single or multiple alleles, or multiple genes
    # Get reads that map only to a single allele, and count them


    Rscript --vanilla ${baseDir}/bin/count_multimapping_reads.R \
    --hla_bam_paths $hla_bams \
    --sample_name ${meta.sample_id} \
    --reference_fasta_path ${personalised_reference} \
    --patient ${meta.patient_id} \
    --out_file_prefix ${meta.sample_id}_${aligner}

    echo Step 2: Quantifying RPKM
    Rscript ${projectDir}/bin/get_rpkm.R \
    --library_size_path ${flagstat} \
    --allele_read_count_path ${meta.sample_id}_${aligner}_allele_read_count.csv \
    --sample_name ${meta.sample_id} \
    --reference_fasta_path ${personalised_reference} \
    --out_file_prefix ${meta.sample_id}_${aligner}

    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")
    RSAM=\$(Rscript -e "cat(paste(packageVersion('Rsamtools'), collapse = ''))")
    SEQINR=\$(Rscript -e "cat(paste(packageVersion('seqinr'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard FilterSamReads --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}
        Rsamtools: \${RSAM}
        seqinr: \${SEQINR}
    END_VERSIONS
    """
}


// process GENERATE_FINAL_RPKM_TABLE {
    
//     tag "GENERATE_FINAL_RPKM_TABLE"

//     container "library://tpjones15/default/final_lohhla:latest"

//     label 'process_single'

//     input:
//     path(rpkm_tables)

//     output:
//     path("cohort_allele_rpkm.csv"), emit: allele_rpkm_table
//     path("cohort_gene_rpkm.csv"), emit: gene_rpkm_table

//     path("versions.yml"), emit: versions

//     script:
//     // make object with rpkm tables separated by a space
//     rpkm_tables = rpkm_tables.join(" ")
//     """
//     # Run Rscript to generate RPKM table
//     Rscript --vanilla ${baseDir}/bin/create_final_rpkms.R \
//     --rpkm_tables ${rpkm_tables} 

//     # Get R version and package versions
//     R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
//     DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
//     AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         R: \${R_VERSION}
//         data.table: \${DT_VERSION}
//         argparse: \${AP_VERSION}
//     END_VERSIONS
//     """

// }

// process GENERATE_COHORT_RNA_TABLES {
    
//     tag "GENERATE_COHORT_RNA_TABLES"

//     container "library://tpjones15/default/final_lohhla:latest"

//     label 'process_single'

//     input: 
//     path  rna_csvs 
//     path  rna_csvs_file
//     path  inventory
//     val aligner
//     val snp_type

//     output:
//     path("*.csv"),  emit: rna_cohort_table

//     path("versions.yml"), emit: versions

//     """
//     # Run Rscript to generate AIB table
//     Rscript --vanilla ${baseDir}/bin/create_cohort_rna_tables.R \
//     --rna_aib_file_path ${rna_csvs_file} \
//     --inventory_path ${inventory} \
//     --min_frac_mapping_uniquely ${params.min_frac_mapping_uniquely} \
//     --max_frac_mapping_multi_gene ${params.max_frac_mapping_multi_gene} \
//     --out_file_prefix mhc_hammer_rna_${aligner}_${snp_type}

//     # Get R version and package versions
//     R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
//     DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
//     AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         R: \${R_VERSION}
//         data.table: \${DT_VERSION}
//         argparse: \${AP_VERSION}
//     END_VERSIONS
//     """

// }