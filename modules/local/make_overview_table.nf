process CREATE_MHC_HAMMER_TABLE {

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'

    input: 
    path  output_csvs 
    path  inventory
    path  input_csvs_file
    path  hlahd_germline_samples
    val(min_depth)

    output:
    path "cohort_mhc_hammer_gene_table.csv", emit: cohort_mhc_hammer_gene_table
    path "cohort_mhc_hammer_allele_table.csv", emit: cohort_mhc_hammer_allele_table
    path "versions.yml", emit: versions

    // when:
    // only run if there are tumour samples in the cohort
    // total_tumour_count > 0 

    script:
    """
    Rscript ${projectDir}/bin/make_cohort_overview_table.R \
    --inventory_path ${inventory} \
    --csv_tables_path ${input_csvs_file} \
    --hlahd_germline_samples_path ${hlahd_germline_samples} \
    --max_cn_range ${params.max_copy_number_range} \
    --min_n_snps ${params.min_number_of_snps} \
    --min_expected_depth ${params.min_expected_depth} \
    --min_frac_mapping_uniquely ${params.min_frac_mapping_uniquely} \
    --max_frac_mapping_multi_gene ${params.max_frac_mapping_multi_gene} \
    --dna_snp_min_depth ${min_depth}

    # get R version and data.table version
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

process CREATE_MOSDEPTH_COHORT_TABLE {

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'

    input: 
    path  output_csvs 
    path  inventory
    path  input_csvs_file

    output:
    path "mosdepth_*.csv", emit: overview_table
    path "versions.yml", emit: versions

    script:
    """
    Rscript ${projectDir}/bin/make_cohort_mosdepth_table.R \
    --inventory_path ${inventory} \
    --csv_tables_path ${input_csvs_file} 

    # get R version and data.table version
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

process CREATE_LIBRARY_SIZE_TABLE {

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'

    input: 
    path  output_csvs 
    path  inventory
    path  input_csvs_file

    output:
    path "cohort_library_size.csv", emit: library_size_table
    path "versions.yml", emit: versions

    // when:
    // only run if there are tumour samples in the cohort
    // total_tumour_count > 0 

    script:
    """
    echo "running..."
    Rscript ${projectDir}/bin/make_cohort_library_size_table.R \
    --inventory_path ${inventory} \
    --csv_tables_path ${input_csvs_file} 

    # get R version and data.table version
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

process CREATE_ALTERNATIVE_SPLICING_COHORT_TABLE {

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'

    input: 
    path  output_csvs 
    path  inventory
    path  input_csvs_file

    output:
    path "novel_splicing_events.csv", emit: novel_splicing_events
    path "known_splicing_events.csv", emit: known_splicing_events
    path "versions.yml", emit: versions

    script:
    """
    Rscript ${projectDir}/bin/make_cohort_alternative_splicing_table.R \
    --inventory_path ${inventory} \
    --csv_tables_path ${input_csvs_file} 

    # get R version and data.table version
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