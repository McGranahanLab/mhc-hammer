//
// SUBWORKFLOW 1: Detecting HLA allele specific repression in RNA tumour samples vs matched normals
//

include { DETECT_HLA_REPRESSION as DETECT_HLA_REPRESSION_ALL_SNPS;
        DETECT_HLA_REPRESSION as DETECT_HLA_REPRESSION_EXON_SNPS;
        GET_HLA_ALLELE_EXPRESSION;
        GET_HLA_ALLELIC_IMBALANCE as GET_HLA_ALLELIC_IMBALANCE_ALL_SNPS;
        GET_HLA_ALLELIC_IMBALANCE as GET_HLA_ALLELIC_IMBALANCE_EXON_SNPS} from '../../modules/local/mhc_hammer_rna_modules'
include { filterByMetadataField } from '../../lib/core_functions'

workflow RNA_ANALYSIS {

    take: 
    rna_flagstat_ch
    snp_positions
    patient_allele_tables_ch
    hla_allele_bams_ch
    passed_heterozygous_hla_genes_ch
    passed_heterozygous_hla_alleles_ch
    hla_bam_read_count_ch
    inventory
    reference_type
    mhc_reference
    aligner
    
    main:
    versions = Channel.empty()
    rna_exon_snps_output_ch = Channel.empty()
    rna_all_snps_output_ch = Channel.empty()

    // MODULE 1: Detecting HLA allele specific repression in tumour samples vs matched normals
    
    // Filter the input channels for tumour and normal samples
    tumour_hla_genes_ch = passed_heterozygous_hla_alleles_ch
            .filter(filterByMetadataField("sample_type", "tumour"))
            .map{ meta, passed_alleles -> 
                tuple(meta.normal_sample_id, meta.sample_id, passed_alleles)
            }
    
    normal_hla_genes_ch =  passed_heterozygous_hla_alleles_ch
            .filter(filterByMetadataField("sample_type", "normal"))
            .map{ meta, passed_alleles -> 
                tuple(meta.sample_id, passed_alleles)
            }
    // Join the tumour and normal channels on sample_id and their matched normal
    // This keeps tumour samples with matched normals
    tumour_normal_ch = tumour_hla_genes_ch.combine(normal_hla_genes_ch, by:[0]) 

    // Compare passed heterozygous genes for tumour and normal samples
    tumour_normal_passed_allele_ch = tumour_normal_ch.map { normal_id, sample_id, tumour_passed_het_alleles, normal_passed_het_alleles ->
        normal_alleles = normal_passed_het_alleles.text.readLines().findAll { it =~ /^hla_/ }.sort()
        tumour_alleles = tumour_passed_het_alleles.text.readLines().findAll { it =~ /^hla_/ }.sort()
        

    common_alleles = tumour_alleles.intersect(normal_alleles)
    
    // // Can only detect repression in alleles that pass in both tumour and normal
        if (common_alleles) {
            [sample_id, tuple(common_alleles)]
        } else {
            [''] // make common alleles an empty channel so it doesnt combine - wont run MHChammer 
        }
    }

    tumour_bams_ch = hla_allele_bams_ch
                .filter(filterByMetadataField("sample_type", "tumour"))
                .map{ meta, hla_bams, reference -> tuple(meta.normal_sample_id, meta, hla_bams, reference) }

    normal_bams_ch = hla_allele_bams_ch
                .filter(filterByMetadataField("sample_type", "normal"))
                .map{ meta, hla_bams, reference -> tuple(meta.sample_id, hla_bams) }

    repression_ch = tumour_bams_ch
                .combine(normal_bams_ch, by:0) // This will only keep tumour samples with a matched normal
                .combine(rna_flagstat_ch, by:0) // This will merge in the normal flagstat
                .map{ normal_sample_id, meta, tumour_hla_bams, reference, normal_hla_bams, normal_flagstat ->
                    tuple(meta.sample_id, meta, tumour_hla_bams, reference, normal_hla_bams, normal_flagstat) 
                    }
                .combine(tumour_normal_passed_allele_ch, by: 0) // This adds the alleles we can measure repression in
                .combine(rna_flagstat_ch, by:0) // This adds the tumour bam flagstat csv
                .map{ tumour_sample_id, meta, tumour_hla_bams, reference, normal_hla_bams, normal_flagstat, passed_alleles, tumour_flagstat ->
                    tuple(meta.patient_id, meta, tumour_hla_bams, reference, normal_hla_bams, normal_flagstat, passed_alleles, tumour_flagstat) 
                    }
                .combine(snp_positions, by:0) // This adds in the patient snp positions
                .combine(mhc_reference, by:0) 
                .map{ tumour_sample_id, meta, tumour_hla_bams, reference, normal_hla_bams, normal_flagstat, passed_alleles, tumour_flagstat, snp_positions, gtf, fasta, genome_size ->
                    tuple(meta.patient_id, meta, tumour_hla_bams, reference, normal_hla_bams, normal_flagstat, passed_alleles, tumour_flagstat, snp_positions, gtf) 
                    }
    
    if (reference_type == "genome") {
        
        // DETECT_HLA_REPRESSION_EXON_SNPS ( repression_ch, "exon_snps", reference_type, aligner )
        DETECT_HLA_REPRESSION_ALL_SNPS ( repression_ch, "all_snps", reference_type, aligner )

        // rna_exon_snps_output_ch = rna_exon_snps_output_ch
        //     .mix(DETECT_HLA_REPRESSION_EXON_SNPS.out.hla_repression_tables)

        rna_all_snps_output_ch = rna_all_snps_output_ch
            .mix(DETECT_HLA_REPRESSION_ALL_SNPS.out.hla_repression_tables)

    } else {
        
        DETECT_HLA_REPRESSION_ALL_SNPS ( repression_ch, "all_snps", reference_type, aligner )
        rna_all_snps_output_ch = rna_all_snps_output_ch
            .mix(DETECT_HLA_REPRESSION_ALL_SNPS.out.hla_repression_tables)
            
    }

    // MODULE 2: Detect allele specific expression
     
    allele_bam_ref_ch = hla_allele_bams_ch
            .map{ meta, hla_bams, reference -> tuple(meta.sample_id, meta, hla_bams, reference) }          
            
    detect_hla_expression_ch = allele_bam_ref_ch
            .combine(rna_flagstat_ch, by:0)
            .map{ sample_id, meta, bams, reference, flagstat -> 
                tuple(meta, bams, reference, flagstat, aligner)}

    GET_HLA_ALLELE_EXPRESSION ( detect_hla_expression_ch )

    rna_all_snps_output_ch = rna_all_snps_output_ch
            .mix(GET_HLA_ALLELE_EXPRESSION.out.rpkm_table)
    rna_exon_snps_output_ch = rna_exon_snps_output_ch
            .mix(GET_HLA_ALLELE_EXPRESSION.out.rpkm_table)        

    // MODULE 3: Detect allelic imbalance 

    // Get the passed heterozygous genes as a value instead of a path
    het_hla_genes_merge_ch = passed_heterozygous_hla_genes_ch
            .map{ meta, passed_heterozygous_hla_genes -> 
                passed_heterozygous_genes = passed_heterozygous_hla_genes.text.readLines().findAll { it =~ /^hla_/ }.sort()
                tuple(meta.sample_id, passed_heterozygous_genes)
            }

    allelic_imbalance_ch = allele_bam_ref_ch // channel with sample_name, bams and reference
            .combine(het_hla_genes_merge_ch, by:0) 
            .combine(GET_HLA_ALLELE_EXPRESSION.out.tables_for_aib, by:0) 
            .map{ sample_id, meta, bams, reference, passed_heterozygous_hla_genes, reads_mapping_one_allele ->
            tuple(meta.patient_id, meta, bams, reference, passed_heterozygous_hla_genes, reads_mapping_one_allele) }
            .combine(snp_positions, by:0)
            .combine(mhc_reference, by:0)             
            .map{ patient_id, meta, bams, reference, passed_heterozygous_hla_genes, reads_mapping_one_allele, snp_positions, gtf, fasta, genome_size ->
            tuple(meta, bams, reference, passed_heterozygous_hla_genes, reads_mapping_one_allele, snp_positions, gtf) }

    if (reference_type == "genome") {
        
        // GET_HLA_ALLELIC_IMBALANCE_EXON_SNPS ( allelic_imbalance_ch, "exon_snps", reference_type, aligner )
        GET_HLA_ALLELIC_IMBALANCE_ALL_SNPS ( allelic_imbalance_ch, "all_snps", reference_type, aligner )

        // rna_exon_snps_output_ch = rna_exon_snps_output_ch
        //     .mix(GET_HLA_ALLELIC_IMBALANCE_EXON_SNPS.out.rna_aib_tables)

        rna_all_snps_output_ch = rna_all_snps_output_ch
            .mix(GET_HLA_ALLELIC_IMBALANCE_ALL_SNPS.out.rna_aib_tables)

    } else {
        
        GET_HLA_ALLELIC_IMBALANCE_ALL_SNPS ( allelic_imbalance_ch, "all_snps", reference_type, aligner )
        rna_all_snps_output_ch = rna_all_snps_output_ch
            .mix(GET_HLA_ALLELIC_IMBALANCE_ALL_SNPS.out.rna_aib_tables)
            
    }          

    // versions = versions.mix(DETECT_HLA_REPRESSION_EXON_SNPS.out.versions.first())
    // versions = versions.mix(DETECT_HLA_REPRESSION_ALL_SNPS.out.versions.first())
    // versions = versions.mix(GET_HLA_ALLELIC_IMBALANCE.out.versions.first())
    // versions = versions.mix(GENERATE_COHORT_RNA_TABLES.out.versions)

    emit:
    // versions
    rna_all_snps_output_ch
    rna_exon_snps_output_ch
}