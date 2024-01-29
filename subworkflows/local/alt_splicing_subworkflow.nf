
include { DETECT_ALT_SPLICING; TUMOUR_NORMAL_ALT_SPLICING_ENRICHMENT } from '../../modules/local/alt_splicing_modules'
include { filterByMetadataField } from '../../lib/core_functions'

workflow ALT_SPLICING {
    
    take:
    splice_junction_ch
    reference_ch
    codon_file

    main:
    // make empty versions channel
    versions = Channel.empty()

    // merge star align second pass channel with patient gtf
    alt_splicing_input_ch = splice_junction_ch
                              .combine(reference_ch, by: 0)
                              .map{ patient_id, meta, sj_table, gtf, personalised_reference, genome_size ->
                                    tuple(meta, sj_table, gtf, personalised_reference)}
    
    DETECT_ALT_SPLICING ( 
      alt_splicing_input_ch,
      codon_file
    )

    // get tumour normal enrichment if the tumour has a matched normal
    tumour_sjs_ch = DETECT_ALT_SPLICING.out.sjs_table
            .filter(filterByMetadataField("sample_type", "tumour"))
            .map{ meta, novel_sjs, known_sjs -> tuple(meta.normal_sample_name, meta, novel_sjs, known_sjs) }
            
    normal_sjs_ch = DETECT_ALT_SPLICING.out.sjs_table
            .filter(filterByMetadataField("sample_type", "normal"))
            .map{ meta, novel_sjs, known_sjs -> tuple(meta.sample_id, novel_sjs, known_sjs) }


    combined_tumour_normal_sjs_ch = tumour_sjs_ch
            .combine(normal_sjs_ch, by:0) // combine normal hla bams channel
            .map{ normal_sample_name, meta, tumour_novel_sjs, tumour_known_sjs, normal_novel_sjs, normal_known_sjs ->
                tuple(meta, tumour_novel_sjs, tumour_known_sjs, normal_novel_sjs, normal_known_sjs) 
                }

    TUMOUR_NORMAL_ALT_SPLICING_ENRICHMENT ( 
          combined_tumour_normal_sjs_ch
        )

    versions = versions.mix(DETECT_ALT_SPLICING.out.versions.first())

    emit: 
    tumour_normal_sjs_ch = TUMOUR_NORMAL_ALT_SPLICING_ENRICHMENT.out.tumour_normal_splice_junctions_tables
    novel_sjs_ch = DETECT_ALT_SPLICING.out.novel_splice_junctions_tables
    known_sjs_ch = DETECT_ALT_SPLICING.out.known_splice_junctions_tables
    versions
}