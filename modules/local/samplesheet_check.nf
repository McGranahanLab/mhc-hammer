process SAMPLESHEET_CHECK {
    tag "${samplesheet}"

    executor 'local' 
    memory 100.MB

    container 'library://tpjones15/default/final_lohhla:latest'

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in McGranahanLab/mhc_hammer/bin/
    """
    Rscript --vanilla ${baseDir}/bin/check_samplesheet.R \
        --sample_sheet_path ${samplesheet} \
        --validated_sample_sheet_path samplesheet.valid.csv
    
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
