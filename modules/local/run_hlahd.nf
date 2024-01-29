process HLAHD_LOCAL { 
    tag "${meta.patient_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_low'
    label 'error_retry'

    input:
    tuple val(meta), path(fqs)
    path mhc_gtf

    output:
    tuple val(meta.patient_id), \
          path("${meta.patient_id}_hla_alleles.csv"), emit: genotype
    path "${meta.patient_id}_hla_alleles.csv"       , emit: allele_table
    path "log/"                                     , emit: hlahd_logs
    path "result/"                                  , emit: hlahd_results
    path "versions.yml"                             , emit: versions

    script: 
    def args = task.ext.args ?: ''
    """
    
    # export paths to local binaries
    export PATH=\${PATH}:\${BOWTIE2_PATH}:\${HLAHD_BIN_PATH}
    bash \${HLAHD_PATH}/bin/hlahd.sh \
    ${args} \
    -t ${task.cpus} \
    -f \${HLAHD_PATH}/freq_data \
    ${fqs[0]} \
    ${fqs[1]} \
    \${HLAHD_PATH}/HLA_gene.split.txt \
    \${HLAHD_PATH}/dictionary \
    ${meta.patient_id} \
    "."

    # clean up directory - remove high mem intermediate files
    rm -rf ${meta.patient_id}/exon
    rm -rf ${meta.patient_id}/intron
    rm -rf ${meta.patient_id}/mapfile
    rm -rf ${meta.patient_id}/maplist
    rm ${meta.patient_id}/pickup.sh
    rm ${meta.patient_id}/estimation.sh

    # move the contents of the patient directory into the current directory 
    mv ${meta.patient_id}/* ./

    # remove the patient directory
    rm -rf ${meta.patient_id}/
    
    # Sometimes HLAHD will terminate due to memory requirements, but the job continues to run
    # To ensure this doesn't happen - check if HLA A, B and C estimates exist
    # if not, don't produce required outputs.
    # The patient hla allele types csv is not produced and the job will be restarted with more resources
    if [ -f result/${meta.patient_id}_A.est.txt ] \
        && [ -f result/${meta.patient_id}_B.est.txt ] \
        && [ -f result/${meta.patient_id}_C.est.txt ] \
        && [ -f result/${meta.patient_id}_E.est.txt ] \
        && [ -f result/${meta.patient_id}_F.est.txt ] \
        && [ -f result/${meta.patient_id}_G.est.txt ]; then
    
        Rscript ${projectDir}/bin/hlahd_parse_output.R \
        --hlahd_folder result \
        --gtf_path ${mhc_gtf} \
        --sample_id ${meta.patient_id} \
        --genes A B C E F G

        mv result/${meta.patient_id}_hla_alleles.csv ./

    else 
        echo "HLAHD failed to produce HLA A, B and C estimates"
        echo "This process will be restarted with more resources"
    fi
    
    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")

    # Write versions to YAML file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        HLAHD: \$(head -n 1 .command.log | awk '{print \$3}')
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}
    END_VERSIONS
    """
}

process HLAHD { 
    tag "${meta.patient_id}"

    container "${params.singularityDir}hlahd.sif"

    label 'process_low'
    label 'error_retry'

    input:
    tuple val(meta), path(fqs)
    path mhc_gtf

    output:
    tuple val(meta.patient_id), \
          path("${meta.patient_id}_hla_alleles.csv"), emit: genotype
    path "${meta.patient_id}_hla_alleles.csv"       , emit: allele_table
    path "log/"                                     , emit: hlahd_logs
    path "result/"                                  , emit: hlahd_results
    path "versions.yml"                             , emit: versions

    script: 
    def args = task.ext.args ?: ''
    """
    # get hlahd path
    echo "non local"
    hlahd_path="/bin/hlahd/"

    bash \${hlahd_path}bin/hlahd.sh \
    ${args} \
    -t ${task.cpus} \
    -f \${hlahd_path}/freq_data \
    ${fqs[0]} \
    ${fqs[1]} \
    \${hlahd_path}/HLA_gene.split.txt \
    \${hlahd_path}/dictionary \
    ${meta.patient_id} \
    "."

    # clean up directory - remove high mem intermediate files
    rm -rf ${meta.patient_id}/exon
    rm -rf ${meta.patient_id}/intron
    rm -rf ${meta.patient_id}/mapfile
    rm -rf ${meta.patient_id}/maplist
    rm ${meta.patient_id}/pickup.sh
    rm ${meta.patient_id}/estimation.sh

    # move the contents of the patient directory into the current directory 
    mv ${meta.patient_id}/* ./

    # remove the patient directory
    rm -rf ${meta.patient_id}/
    
    # Sometimes HLAHD will terminate due to memory requirements, but the job continues to run
    # To ensure this doesn't happen - check if HLA A, B and C estimates exist
    # if not, don't produce required outputs.
    # The patient hla allele types csv is not produced and the job will be restarted with more resources
    if [ -f result/${meta.patient_id}_A.est.txt ] \
        && [ -f result/${meta.patient_id}_B.est.txt ] \
        && [ -f result/${meta.patient_id}_C.est.txt ] \
        && [ -f result/${meta.patient_id}_E.est.txt ] \
        && [ -f result/${meta.patient_id}_F.est.txt ] \
        && [ -f result/${meta.patient_id}_G.est.txt ]; then
    
        Rscript ${projectDir}/bin/hlahd_parse_output.R \
        --hlahd_folder result \
        --gtf_path ${mhc_gtf} \
        --sample_id ${meta.patient_id} \
        --genes A B C E F G

        mv result/${meta.patient_id}_hla_alleles.csv ./

    else 
        echo "HLAHD failed to produce estimates for all genes"
        echo "This process will be restarted with more resources"
    fi
    
    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")

    # Write versions to YAML file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        HLAHD: \${HLAHD_VERSION}
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}
    END_VERSIONS
    """
}
