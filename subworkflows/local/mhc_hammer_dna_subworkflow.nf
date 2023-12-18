//
// SUBWORKFLOW 1: Detecting HLA allele specific CN and AIB in DNA 
//

include { DETECT_MUTS; PARSE_MUTATIONS; 
            DETECT_CN_AND_AIB as DETECT_CN_AND_AIB_ALL_SNPS;
            DETECT_CN_AND_AIB as DETECT_CN_AND_AIB_EXON_SNPS } from '../../modules/local/mhc_hammer_dna_modules'
include { filterByMetadataField; mutations_detected_check } from '../../lib/core_functions'

workflow DNA_ANALYSIS {
    take: 
    dna_flagstat_ch
    genome_snp_positions_ch
    gtf_ch
    purity_ploidy_ch
    hla_allele_bams_ch
    passed_heterozygous_hla_genes_ch
    
    main:
    versions = Channel.empty()

    // Filter the input channels for tumour and normal samples
    tumour_hla_genes_ch = passed_heterozygous_hla_genes_ch
            .filter(filterByMetadataField("sample_type", "tumour"))
            .map{ meta, passed_heterozygous_genes -> 
                tuple(meta.normal_sample_name, meta.sample_id, passed_heterozygous_genes)
            }
    
    normal_hla_genes_ch =  passed_heterozygous_hla_genes_ch
            .filter(filterByMetadataField("sample_type", "normal"))
            .map{ meta, passed_heterozygous_genes -> 
                tuple(meta.sample_id, passed_heterozygous_genes)
            }
    // Join the tumour and normal channels on sample_id and their matched normal
    tumour_normal_hla_genes_ch = tumour_hla_genes_ch.combine(normal_hla_genes_ch, by:[0])

    // Compare passed heterozygous genes for tumour and normal samples
    common_passed_heterozygous_hla_genes_ch = tumour_normal_hla_genes_ch
                                .map { normal_id, sample_id, tumour_passed_heterozygous_genes, normal_passed_heterozygous_genes ->
        tumour_heterozygous_genes = tumour_passed_heterozygous_genes.text.readLines().findAll { it =~ /^hla_/ }.sort()
        normal_heterozygous_genes = normal_passed_heterozygous_genes.text.readLines().findAll { it =~ /^hla_/ }.sort()

        common_heterozygous_genes = tumour_heterozygous_genes.intersect(normal_heterozygous_genes)

        if (common_heterozygous_genes) {
            [sample_id, tuple(common_heterozygous_genes)]
        } else {
            [''] // make common heterozygous genes an empty channel so it doesnt combine - wont run MHChammer 
        }
    }

    tumour_bams_ch = hla_allele_bams_ch
            .filter(filterByMetadataField("sample_type", "tumour"))
            .map{ meta, hla_bams, reference -> tuple(meta.normal_sample_name, meta, hla_bams, reference) }

    normal_bams_ch = hla_allele_bams_ch
            .filter(filterByMetadataField("sample_type", "normal"))
            .map{ meta, hla_bams, reference -> tuple(meta.sample_id, hla_bams) }

    cn_and_aib_input_ch = tumour_bams_ch
            .combine(normal_bams_ch, by:0) // combine normal hla bams channel
            .combine(dna_flagstat_ch, by:0) // combine normal bam flagstat csv
            .map{ normal_sample_name, meta, tumour_hla_bams, reference, normal_hla_bams, normal_flagstat ->
                tuple(meta.sample_id, meta, tumour_hla_bams, reference, normal_hla_bams, normal_flagstat) }
            .combine(common_passed_heterozygous_hla_genes_ch, by:0) // combine channel with genes to run MHChammer on
            .combine(purity_ploidy_ch, by:0) // Combine purity ploidy channel
            .combine(dna_flagstat_ch, by:0) // combine tumour bam flagstat csv
            .map{ tumour_sample_id, meta, tumour_hla_bams, reference, normal_hla_bams, normal_flagstat, passed_hetero_hla_genes, purity_ploidy, tumour_flagstat ->
                tuple(meta.patient_id, meta, tumour_hla_bams, reference, normal_hla_bams, normal_flagstat, passed_hetero_hla_genes, purity_ploidy, tumour_flagstat) 
                }
            .combine(genome_snp_positions_ch, by:0) // combine patient snp positions channel
            .combine(gtf_ch, by:0) // combine patient gtf 
    
    // Detect LOH and AIB 
    DETECT_CN_AND_AIB_ALL_SNPS ( cn_and_aib_input_ch, params.min_depth, "all_snps", "novoalign")
    // DETECT_CN_AND_AIB_EXON_SNPS ( cn_and_aib_input_ch, params.min_depth, "exon_snps", "novoalign" )

    emit:
    versions = versions.mix(DETECT_CN_AND_AIB_ALL_SNPS.out.versions.first())
    dna_all_snps_output_ch = DETECT_CN_AND_AIB_ALL_SNPS.out.dna_analysis_output
    // dna_exon_snps_output_ch = DETECT_CN_AND_AIB_EXON_SNPS.out.dna_analysis_output
}

//
// SUBWORKFLOW 2: Detecting HLA allele specific mutations in DNA 
//
workflow DETECT_MUTATIONS {

    take:
    hla_allele_bams_ch
    passed_hla_alleles_ch
    personilsed_reference_and_gtf
    wxs_sample_count_ch
    germline_sample_count_ch
    gl_samples
    inventory

    main:

    versions = Channel.empty()

    // Filter the input channels for tumour and normal samples
    tumour_hla_alleles_ch = passed_hla_alleles_ch
            .filter(filterByMetadataField("sample_type", "tumour"))
            .map{ meta, passed_alleles -> 
                tuple(meta.normal_sample_name, meta.sample_id, passed_alleles)
            }
    
    normal_hla_alleles_ch =  passed_hla_alleles_ch
            .filter(filterByMetadataField("sample_type", "normal"))
            .map{ meta, passed_alleles -> 
                tuple(meta.sample_id, passed_alleles)
            }
    // Join the tumour and normal channels on sample_id and their matched normal
    joined_hla_alleles_ch = tumour_hla_alleles_ch.combine(normal_hla_alleles_ch, by:[0])

    // Compare passed alleles for tumour and normal samples
    common_passed_hla_alleles_ch = joined_hla_alleles_ch.map { normal_id, sample_id, tumour_passed_hla_alleles, normal_passed_hla_alleles ->
        tumour_alleles = tumour_passed_hla_alleles.text.readLines().findAll { it =~ /^hla_/ }.sort()
        normal_alleles = normal_passed_hla_alleles.text.readLines().findAll { it =~ /^hla_/ }.sort()

        common_alleles = tumour_alleles.intersect(normal_alleles)
    // Output warning message when their is missmatch between the passed hla alleles between tumour and normal samples
    // passed here is a non-empty hla-allele-specific bam 
    // Can only detect mutations where there is a matched normal
        if (common_alleles) {
            if (tumour_alleles != normal_alleles) {
                log.warn("Mismatch in passed hla alleles for tumour and normal samples:" +
                "Normal_id: ${normal_id}, Tumour_id: ${sample_id}\n" +
                "Tumour alleles: ${tumour_alleles.join(', ')}\n" +
                "Normal alleles: ${normal_alleles.join(', ')}\n" +
                "Mutation calling will only run on the alleles with a matched normal")
            }
            [sample_id, tuple(common_alleles)]
        } else {
            log.warn("No common passed hla alleles for tumour and normal samples: Normal_id: ${normal_id}, Tumour_id: ${sample_id}\n" +
            "Mutation calling will not run for this sample")
            [''] // make commoon alleles an empty channel so it doesnt combine - wont run MHChammer 
        }
    }

    tumour_bams_ch = hla_allele_bams_ch
            .filter(filterByMetadataField("sample_type", "tumour"))
            .map{ meta, hla_bams, reference -> tuple(meta.normal_sample_name, meta, hla_bams) }

    normal_bams_ch = hla_allele_bams_ch
            .filter(filterByMetadataField("sample_type", "normal"))
            .map{ meta, hla_bams, reference -> tuple(meta.sample_id, hla_bams) }

    combined_tumour_normal_mutation_input_ch = tumour_bams_ch
            .combine(normal_bams_ch, by:0) // combine normal hla bams channel
            .map{ normal_sample_name, meta, tumour_hla_bams, normal_hla_bams ->
                tuple(meta.sample_id, meta, tumour_hla_bams, normal_hla_bams) 
                }
            .combine(common_passed_hla_alleles_ch, by:0) // combine channel with alleles to run MHChammer on
            .map{ tumour_sample_id, meta, tumour_hla_bams, normal_hla_bams, passed_hla_alleles ->
                tuple(meta.patient_id, meta, tumour_hla_bams, normal_hla_bams, passed_hla_alleles) 
                }
            .combine(personilsed_reference_and_gtf, by:0) // combine reference and zipped_gtf

    // MODULE: run mutect2 in single sample mode
    DETECT_MUTS ( combined_tumour_normal_mutation_input_ch ) 

    // Merge sample level output for mutation calls
    patient_mutations = DETECT_MUTS.out.mutation_output
            // move patient_id outside meta and change meta to normal_sample_name
            .map{ meta, normal_bams, tumour_bams, mutations ->
                tuple(meta.patient_id, meta.normal_sample_name, normal_bams, tumour_bams, mutations)
            }
            // combine with wxs_sample_count_ch [patient_id, sample_count]
            .combine(wxs_sample_count_ch, by:0)
            // groupKey prevents stalling
            .map{ key, normal_sample_name, normal_bams, tumour_bams, mutations, sample_count -> 
                tuple ( groupKey(key, sample_count), normal_sample_name, normal_bams, tumour_bams, mutations )
            }
            .groupTuple(by:0)
            // sort normal_bams list alphabetically to ensure consistent order for unique() function
            // This is to ensure -resume works correctly when multiple normal bams are present 
            .map{ key, normal_sample_name, normal_bams, tumour_bams, mutations ->
                sorted_normal_bams = normal_bams.sort()
                unique_normal_bams = sorted_normal_bams.unique { it.toString().tokenize('/').last() }
                tuple ( key, normal_sample_name, unique_normal_bams, tumour_bams, mutations )
            }
            // flatten bam and mutation fields so multiple sample bams and mutations are merged into single fields
            .map{ key, normal_sample_name, normal_bams, tumour_bams, mutations -> 
                tuple ( key, normal_sample_name, normal_bams.flatten(), tumour_bams.flatten(), mutations.flatten() )
            }
            // remove any cases in mutations where the vep table ends in empty.vep.txt
            .map{ key, normal_sample_name, normal_bams, tumour_bams, mutations -> 
                tuple ( key, normal_sample_name, normal_bams, tumour_bams, mutations.findAll{ !it.getName().endsWith("empty.vep.txt") } )
            }
 
    def criteria = branchCriteria {
        samples_with_no_muts: !mutations_detected_check(it)
        samples_with_muts: mutations_detected_check(it)
    }

    mutation_input_by_status = patient_mutations.branch(criteria)

    // MODULE: Parse mutations
    // PARSE_MUTATIONS( mutation_input_by_status.samples_with_muts, inventory )
    
    versions = versions.mix(DETECT_MUTS.out.versions.first()) // channel: [ versions.yml ]
    // versions = versions.mix(PARSE_MUTATIONS.out.versions.first())

    emit:
    versions

}