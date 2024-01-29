//
// Align fqs to personalised reference and generate HLA allele specific bams
//

include { 
    NOVOALIGN; 
    MAKE_HLA_ALLELE_BAMS;
    MOSDEPTH
    } from '../../modules/local/alignment_modules'

workflow GENERATE_RNA_NOVOALIGN_BAMS {
    take: 
    fqs_ch
    personalised_reference_ch
    mosdepth_bed_ch

    main:
    versions = Channel.empty()

    // Combine fqs with patient specific fasta channel
    // map to create common key to combine channels 
    fqs_tmp = fqs_ch.map { meta, fqs ->
        tuple(meta.patient_id, meta, fqs)
    }

    alignment_input_ch = fqs_tmp.combine(personalised_reference_ch, by: 0).map { patient_id, meta, fqs, reference, novoindex ->
        tuple(meta, fqs, reference, novoindex)
    }
    
    // Make rna bams from novoalign
    NOVOALIGN( alignment_input_ch) 
    
    make_hla_bams_input_ch =  NOVOALIGN.out.make_hla_bam_input
                    .map{meta, reference, bam -> tuple(meta, reference, bam, "novoalign")}

    MAKE_HLA_ALLELE_BAMS ( make_hla_bams_input_ch )

    mosdepth_input_ch = MAKE_HLA_ALLELE_BAMS.out.mosdepth_input
                                    .combine(mosdepth_bed_ch, by: 0)
                                    .map{ patient_id, meta, bam, bed -> tuple(meta, bam, bed, "novoalign") }

    MOSDEPTH( mosdepth_input_ch )

    emit:

    // Channel containing file with passed allele pairs 
    passed_heterozygous_hla_genes_ch = MAKE_HLA_ALLELE_BAMS.out.passed_heterozygous_hla_genes
    passed_heterozygous_hla_alleles_ch = MAKE_HLA_ALLELE_BAMS.out.passed_heterozygous_hla_alleles
    hla_allele_bams_ch = MAKE_HLA_ALLELE_BAMS.out.hla_allele_bams
    hla_bam_read_count_ch = MAKE_HLA_ALLELE_BAMS.out.hla_bam_read_count
    mosdepth_ch = MOSDEPTH.out.mosdepth_bed
    versions = versions.mix(NOVOALIGN.out.versions.first())              // channel: [ versions.yml ]
    versions = versions.mix(MAKE_HLA_ALLELE_BAMS.out.versions.first())    // channel: [ versions.yml ]
    versions = versions.mix(MOSDEPTH.out.versions.first())    // channel: [ versions.yml ]

}