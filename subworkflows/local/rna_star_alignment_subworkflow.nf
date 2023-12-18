//
// Align fqs to personalised reference and generate HLA allele specific bams
//

include { 
    STAR_ALIGN_FIRST_PASS; GENERATE_ALL_SPLICE_TABLE;
    STAR_ALIGN_SECOND_PASS; ALT_SPLICE_KMER;
    MAKE_HLA_ALLELE_BAMS;
    MOSDEPTH
    } from '../../modules/local/alignment_modules'

workflow GENERATE_RNA_STAR_BAMS {
    take: 
    mosdepth_bed_ch
    subset_bam_ch
    star_reference_ch
    kmer_file_ch

    main:
    versions = Channel.empty()

    // Make rna bams from star
    // merge subset_bam channel and star input channel 
    star_input_ch = subset_bam_ch.map{meta, bam -> tuple(meta.patient_id, meta, bam)}
                                .combine(star_reference_ch, by: 0)
                                 .map{ patient_id, meta, bam, gtf, personalised_reference, genome_size -> 
                                    tuple(meta, bam, gtf, personalised_reference, genome_size) 
                                     }
    
    // MODULE: subset bam for kmer reads, then SamToFastq, then fastqc to get read length, then STAR align
    STAR_ALIGN_FIRST_PASS ( star_input_ch, kmer_file_ch )

    // MODULE: make new kmer file
    // create channel containing only meta.patient_id patient_gtf and personalised_reference
    patient_gtf_ch = star_reference_ch 
                      .map{ patient_id, gtf, personalised_reference, genome_size ->
                            tuple(patient_id, gtf, personalised_reference) }

    kmer_input_ch = STAR_ALIGN_FIRST_PASS.out.kmer_input
                                         .combine(patient_gtf_ch, by: 0)
                                         .map{patient_id, meta, sj_tab, gtf, personalised_reference ->
                                              tuple(meta, sj_tab, gtf, personalised_reference)}

    
    ALT_SPLICE_KMER (
      kmer_input_ch,
      kmer_file_ch
    )

    // MODULE: simple 1 liner to generate the cohort level splice table 
    STAR_ALIGN_FIRST_PASS.out.splicing_table.collect() | GENERATE_ALL_SPLICE_TABLE 

    // MODULE: STAR align second pass - now utilising cohort level splicing table (filtered for patient specific HLAs)
    star_align_second_pass_ch = STAR_ALIGN_FIRST_PASS.out.second_pass_input
                                                     .combine(ALT_SPLICE_KMER.out.sample_kmer,by: 0)
                                                     .map{ sample_id, meta, bam, star_reference, fasta, gtf, genome_size, overhang, kmer ->
                            tuple(meta, bam, star_reference, fasta, gtf, genome_size, overhang, kmer) }

    // star_align_second_pass_ch.view()

    STAR_ALIGN_SECOND_PASS ( 
      star_align_second_pass_ch,
      GENERATE_ALL_SPLICE_TABLE.out.cohort_splice_table
    )

    make_hla_bams_input_ch = STAR_ALIGN_SECOND_PASS.out.star_aligned_bams
                              .combine(star_reference_ch, by:0)
                              .map{patient_id, meta, bam, gtf, reference, genome_size ->
                              tuple(meta, reference, bam, "star")}

    MAKE_HLA_ALLELE_BAMS(make_hla_bams_input_ch)                              

    mosdepth_input_ch = MAKE_HLA_ALLELE_BAMS.out.mosdepth_input
                                    .combine(mosdepth_bed_ch, by: 0)
                                    .map{ patient_id, meta, bam, bed -> tuple(meta, bam, bed, "star") }

    MOSDEPTH( mosdepth_input_ch )

    emit:
    hla_allele_bams_ch = MAKE_HLA_ALLELE_BAMS.out.hla_allele_bams
    splice_junction_ch = STAR_ALIGN_SECOND_PASS.out.sj_table
    passed_heterozygous_hla_genes_ch = MAKE_HLA_ALLELE_BAMS.out.passed_heterozygous_hla_genes
    passed_heterozygous_hla_alleles_ch = MAKE_HLA_ALLELE_BAMS.out.passed_heterozygous_hla_alleles
    hla_bam_read_count_ch = MAKE_HLA_ALLELE_BAMS.out.hla_bam_read_count
    mosdepth_ch = MOSDEPTH.out.mosdepth_bed
    
    // hla_allele_bams_ch = MAKE_HLA_ALLELE_BAMS.out.hla_allele_bams
    // passed_heterozygous_hla_genes_ch = MAKE_HLA_ALLELE_BAMS.out.passed_heterozygous_hla_genes
    // passed_heterozygous_hla_alleles_ch = MAKE_HLA_ALLELE_BAMS.out.passed_heterozygous_hla_alleles
    // hla_bam_read_count_ch = MAKE_HLA_ALLELE_BAMS.out.hla_bam_read_count

    // // Channel containing bam aligned to personalised HLA reference
    // hla_bam_ch = STAR_ALIGN_SECOND_PASS.out.make_bam_input

    // Channel containing HLA allele specific bams
    // hla_allele_bams_ch = MAKE_HLA_ALLELE_BAMS.out.hla_allele_bams

    // Channel containing the number of reads before and after filtering
    // aligned_read_count_ch = MAKE_HLA_ALLELE_BAMS.out.aligned_read_count

    // Channel containing file with passed allele pairs between gl and tumour sample
    // passed_hla_genes_ch = MAKE_HLA_ALLELE_BAMS.out.passed_hla_genes
    // passed_hla_alleles_ch = MAKE_HLA_ALLELE_BAMS.out.passed_hla_alleles
    // passed_heterozygous_hla_alleles_ch = MAKE_HLA_ALLELE_BAMS.out.passed_heterozygous_hla_alleles
    // hla_bam_read_count_ch = MAKE_HLA_ALLELE_BAMS.out.hla_bam_read_count

    // versions = versions.mix(NOVOALIGN.out.versions.first())              // channel: [ versions.yml ]
    // versions = versions.mix(MAKE_HLA_ALLELE_BAMS.out.versions.first())    // channel: [ versions.yml ]
    // versions = versions.mix(MOSDEPTH.out.versions.first())    // channel: [ versions.yml ]

}