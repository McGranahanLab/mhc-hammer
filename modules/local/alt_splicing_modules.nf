
process DETECT_ALT_SPLICING {

    tag "${meta.sample_id}"

    container "library://tpjones15/default/final_lohhla:latest"

    label 'process_single'

    input:
    tuple val(meta), path(sample_splice_table), path(patient_gtf), path(fasta)
    path(codon_file_path)

    output:
    tuple val(meta), \
        path("${meta.sample_id}_novel_splice_junctions.csv"), \
        path("${meta.sample_id}_known_splice_junctions.csv"), emit: sjs_table
    path ("*_novel_splice_junctions.csv"), emit: novel_splice_junctions_tables
    path ("*_known_splice_junctions.csv"), emit: known_splice_junctions_tables
    path ("versions.yml"), emit: versions

    script: // scripts used are bundled with the pipeline, in McGranahanLab/MHChammer/bin/
    """
    echo "rerun"
    echo ${sample_splice_table}
    # Run Rscript to detect novel splice junctions
    Rscript ${projectDir}/bin/annotate_star_splice_junctions.R \
        --sj_tab_path ${sample_splice_table} \
        --hla_genome_fasta_path ${fasta} \
        --gtf_path ${patient_gtf} \
        --sample_id ${meta.sample_id} \
        --codon_table_path ${codon_file_path} \
        --scripts_dir ${projectDir}/bin/

    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    SQ_VERSION=\$(Rscript -e "cat(paste(packageVersion('seqinr'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}      
        seqinr: \${SQ_VERSION}      
    END_VERSIONS
    """

}

process TUMOUR_NORMAL_ALT_SPLICING_ENRICHMENT {

    tag "${meta.sample_id}"

    container "library://tpjones15/default/final_lohhla:latest"

    label 'process_single'

    input:
    tuple val(meta), path(tumour_novel_sjs), path(tumour_known_sjs), path(normal_novel_sjs), path(normal_known_sjs)

    output:
    path ("*_tumour_normal_splice_junctions.csv"), emit: tumour_normal_splice_junctions_tables
    path ("versions.yml"), emit: versions

    script: // scripts used are bundled with the pipeline, in McGranahanLab/MHChammer/bin/
    """

    # Run Rscript to detect novel splice junctions
    Rscript ${projectDir}/bin/splice_junction_tumour_normal_enrichment.R \
        --tumour_novel_sjs_path ${tumour_novel_sjs} \
        --tumour_known_sjs_path ${tumour_known_sjs} \
        --normal_novel_sjs_path ${normal_novel_sjs} \
        --normal_known_sjs_path ${normal_known_sjs} \
        --sample_name ${meta.sample_id} \
        --scripts_dir ${projectDir}/bin/

    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    SQ_VERSION=\$(Rscript -e "cat(paste(packageVersion('seqinr'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}      
        seqinr: \${SQ_VERSION}      
    END_VERSIONS
    """

}

