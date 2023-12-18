//
// Subset bam files and produce fastqs in preparation for next steps 
//

include { GENERATE_HLA_FQS } from '../../modules/local/preprocessing_modules'
include { SUBSET_BAMS } from '../../modules/local/preprocessing_modules'
include { FLAGSTAT } from '../../modules/local/preprocessing_modules'
include { filterByMetadataField } from '../../lib/core_functions'

workflow PREPROCESSING {
    take:
    input_bams
    gl_sample_count
    mhc_coords
    contigs_file_ch

    main:
    versions = Channel.empty()

    // Subset input bam files for reads that could align to the HLA class I genes - only if params.run_bam_subsetting == true
        
    FLAGSTAT ( input_bams )
    SUBSET_BAMS ( input_bams, mhc_coords, contigs_file_ch )
    GENERATE_HLA_FQS ( SUBSET_BAMS.out.subset_bam )

    versions = versions.mix(FLAGSTAT.out.versions.first()) 
    versions = versions.mix(SUBSET_BAMS.out.versions.first())
    versions = versions.mix(GENERATE_HLA_FQS.out.versions.first()) 
    
    // Channel containing only the wxs normal sample fqs
    normal_fqs = GENERATE_HLA_FQS.out.fqs
               .filter(filterByMetadataField("seq", "wxs"))
               .filter(filterByMetadataField("sample_type", "normal"))


    // Only need to run HLAHD on one germline sample per patient - use gl_sample_count to prevent stalling
    // Combine with germline_sample_count_ch [patient_id, sample_count]
    hlahd_sample = normal_fqs
                .map{ meta, fqs -> tuple(meta.patient_id, meta, fqs) }
                .combine(gl_sample_count, by: 0)
                // Map the key, meta, fqs, and num_germlines fields into a single tuple
                .map { key, meta, fqs, num_germlines -> tuple(groupKey(key, num_germlines), meta, fqs) }
                // Group tuples by the first field (groupKey)
                .groupTuple(by: 0)
                // Take the first sample of each group
                .map { key, meta, fqs -> meta.sort()[0] }
                // Map the sample_id field into a new tuple
                .map { meta -> tuple(meta.sample_id) }
                
    emit:

    // Channel containing all fqs to be aligned to personalised references
    all_fqs = GENERATE_HLA_FQS.out.fqs

    // RNA outputs
    // Channel containing all subset rna bams
    rna_bams = SUBSET_BAMS.out.subset_bam
                .filter(filterByMetadataField("seq", "rnaseq"))

    // Channel containing only RNA fqs -> This will be piped into mhchammer RNA
    rna_fqs = GENERATE_HLA_FQS.out.fqs
                .filter(filterByMetadataField("seq", "rnaseq"))
    
    rna_flagstat = FLAGSTAT.out.flagstat
                .filter(filterByMetadataField("seq", "rnaseq"))                 
                .map{ meta, flagstat -> tuple(meta.sample_id, flagstat)}

    // DNA outputs                         
    // combine with normal_fqs so we only take one normal sample into hlahd
    hlahd_input = normal_fqs
                .map{ meta, fqs -> tuple(meta.sample_id, meta, fqs) }
                .join(hlahd_sample, by:0)
                .map{ sample_id, meta, fqs -> tuple(meta, fqs) }

    // Channel containing only WXS fqs -> This will be piped into mhchammer WXS
    wxs_fqs = GENERATE_HLA_FQS.out.fqs
                .filter(filterByMetadataField("seq", "wxs")) 

    dna_flagstat = FLAGSTAT.out.flagstat
                .filter(filterByMetadataField("seq", "wxs"))
                .map{ meta, flagstat -> tuple(meta.sample_id, flagstat)} 

    // Channels contain the mapped and unmapped rate to go into the final table
    unmapped_count_ch = FLAGSTAT.out.library_size_with_unmapped
    mapped_count_ch = FLAGSTAT.out.library_size_without_unmapped

    
    versions
    
}


workflow SUBSET_BAM_PREPROCESSING {
    take:
    input_bams
    gl_sample_count
    mhc_coords

    main:
    versions = Channel.empty()

    // perform flagstat on non-subsetted bam and generate fqs
        
    GENERATE_HLA_FQS ( input_bams )

    versions = versions.mix(GENERATE_HLA_FQS.out.versions.first()) 
    
    // Channel containing only the wxs normal sample fqs
    normal_fqs = GENERATE_HLA_FQS.out.fqs
               .filter(filterByMetadataField("seq", "wxs"))
               .filter(filterByMetadataField("sample_type", "normal"))


    // Only need to run HLAHD on one germline sample per patient - use gl_sample_count to prevent stalling
    // Combine with germline_sample_count_ch [patient_id, sample_count]
    hlahd_sample = normal_fqs
                .map{ meta, fqs -> tuple(meta.patient_id, meta, fqs) }
                .combine(gl_sample_count, by: 0)
                // Map the key, meta, fqs, and num_germlines fields into a single tuple
                .map { key, meta, fqs, num_germlines -> tuple(groupKey(key, num_germlines), meta, fqs) }
                // Group tuples by the first field (groupKey)
                .groupTuple(by: 0)
                // Take the first sample of each group
                .map { key, meta, fqs -> meta.sort()[0] }
                // Map the sample_id field into a new tuple
                .map { meta -> tuple(meta.sample_id) }
    emit:

    // Channel containing all fqs to be aligned to personalised references
    all_fqs = GENERATE_HLA_FQS.out.fqs

    // RNA outputs
    rna_bams = input_bams
                .filter(filterByMetadataField("seq", "rnaseq"))

    // Channel containing only RNA fqs -> This will be piped into mhchammer RNA
    rna_fqs = GENERATE_HLA_FQS.out.fqs
                .filter(filterByMetadataField("seq", "rnaseq"))

    // DNA outputs                         
    // combine with normal_fqs so we only take one normal sample into hlahd
    hlahd_input = normal_fqs
                .map{ meta, fqs -> tuple(meta.sample_id, meta, fqs) }
                .join(hlahd_sample, by:0)
                .map{ sample_id, meta, fqs -> tuple(meta, fqs) }

    // Channel containing only WXS fqs -> This will be piped into mhchammer WXS
    wxs_fqs = GENERATE_HLA_FQS.out.fqs
                .filter(filterByMetadataField("seq", "wxs")) 
    
    versions
    
}