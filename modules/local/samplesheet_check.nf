process SAMPLESHEET_CHECK {
    tag "${samplesheet}"

    executor 'local' 
    memory 100.MB

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    script: 
    """
    Rscript --vanilla ${baseDir}/bin/check_samplesheet.R \
        --sample_sheet_path ${samplesheet} \
        --validated_sample_sheet_path samplesheet.valid.csv \
        --run_hlahd ${params.run_hlahd}
    
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


process CHECK_HLA_TYPE_INPUT {
    tag "CHECK_HLA_TYPE_INPUT"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'

    input:
    path(mhc_fasta)
    path(hla_allele_files)
    path(samplesheet)

    output:
    path "versions.yml" , emit: versions
    path "samplesheet.valid.csv", emit: validated_samplesheet

    script:
    """
    # Check if HLA fasta file exists
    Rscript ${projectDir}/bin/check_hla_fasta.R \\
        --hla_fasta_path ${mhc_fasta} \\
        --hla_allele_files ${hla_allele_files} \\
        --parsed_csv ${samplesheet}
    
    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")
    SEQINR_VERSION=\$(Rscript -e "cat(paste(packageVersion('seqinr'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}
        seqinr: \${SEQINR_VERSION}
    """
}
