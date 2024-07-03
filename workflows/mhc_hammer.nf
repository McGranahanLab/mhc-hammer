/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
nextflow.enable.dsl=2

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Stage placeholder files to be used as an optional input where required
contigs_placeholder = file("${projectDir}/assets/contigs_placeholder.txt", checkIfExists: true)
// transcriptome_placeholder = file("${projectDir}/assets/transcriptome_placeholder.txt", checkIfExists: true)

// Check if mhc_kmer file is provided if params.fish_reads = true
if (params.kmer_file) { 
    mhc_kmer_ch = file(params.kmer_file, checkIfExists: true)
    if (mhc_kmer_ch.isEmpty()) {exit 1, "File provided with --kmer_file is empty: ${mhc_kmer_ch.getName()}!"} 
}

// Check if contig_reads is true, if yes, check if a contigs_file is provided, if yes, check if it exists and is not empty and check that there is only a single contig per line
if (params.contig_reads) {
    if (params.contigs_file) {
        contigs_file_ch = file(params.contigs_file, checkIfExists: true)
        if (contigs_file_ch.isEmpty()) {exit 1, "File provided with --contigs_file is empty: ${contigs_file_ch.getName()}!"}

        // Read the lines of the file into a list
        def lines = contigs_file_ch.readLines()

        // Check that there is only a single element per line
        lines.each { line ->
            def elements = line.split(/[\s,]+/) //split by any whitespace or comma
            if (elements.size() != 1) {
                exit 1, "File provided with --contigs_file should have only a single element per line. Offending line: '${line}'. File: ${contigs_file_ch.getName()}!"
            }
        }
    } else {
        contigs_file_ch = contigs_placeholder
    }
} else {
    contigs_file_ch = contigs_placeholder
}

// Check mhc_coords file if params.mhc_coords is true, also check the file exists and is not empty
if (params.mhc_coords) {
    mhc_coords_ch = file(params.mhc_coords, checkIfExists: true)
    if (mhc_coords_ch.isEmpty()) {exit 1, "File provided with --mhc_coords is empty: ${mhc_coords_ch.getName()}!"}
}

// Check mhc fasta file 
if (params.mhc_fasta) { 
     mhc_fasta_ch = file(params.mhc_fasta, checkIfExists: true)
    if (mhc_fasta_ch.isEmpty()) {exit 1, "File provided with --mhc_fasta is empty: ${mhc_fasta_ch.getName()}!"} 
}

// Check mhc gtf file 
if (params.mhc_gtf) { 
     mhc_gtf_ch = file(params.mhc_gtf, checkIfExists: true)
    if (mhc_gtf_ch.isEmpty()) {exit 1, "File provided with --mhc_gtf is empty: ${mhc_gtf_ch.getName()}!"} 
}

// Check mhc transcriptome fasta file
if (params.mhc_transcriptome_fasta) { 
     mhc_transcriptome_fasta_ch = file(params.mhc_transcriptome_fasta, checkIfExists: true)
    if (mhc_transcriptome_fasta_ch.isEmpty()) {exit 1, "File provided with --mhc_fasta is empty: ${mhc_transcriptome_fasta_ch.getName()}!"} 
} 

if (params.codon_table) { 
     codon_table_ch = file(params.codon_table, checkIfExists: true)
    if (codon_table_ch.isEmpty()) {exit 1, "File provided with --mhc_gtf is empty: ${codon_table_ch.getName()}!"} 
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOWS
//

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPROCESSING; PREPROCESSING_WITHOUT_BAM_SUBSETTING } from '../subworkflows/local/preprocessing_subworkflow'
include { DNA_ANALYSIS; DETECT_MUTATIONS } from '../subworkflows/local/mhc_hammer_dna_subworkflow'
include { GENERATE_DNA_NOVOALIGN_BAMS } from '../subworkflows/local/dna_novoalign_alignment_subworkflow'
include { GENERATE_RNA_NOVOALIGN_BAMS } from '../subworkflows/local/rna_novoalign_alignment_subworkflow'
include { GENERATE_RNA_STAR_BAMS } from '../subworkflows/local/rna_star_alignment_subworkflow'
include { ALT_SPLICING } from '../subworkflows/local/alt_splicing_subworkflow'
include { RNA_ANALYSIS as NOVOALIGN_RNA_ANALYSIS;
            RNA_ANALYSIS as STAR_RNA_ANALYSIS } from '../subworkflows/local/mhc_hammer_rna_subworkflow'

//
// MODULES 
//

include { HLAHD; HLAHD_LOCAL } from '../modules/local/run_hlahd'
include { GENERATE_REFERENCES } from '../modules/local/generate_references'
include { CREATE_TUMOUR_NORMAL_SPLICING_COHORT_TABLE; CREATE_LIBRARY_SIZE_TABLE;
        CREATE_MOSDEPTH_COHORT_TABLE; CREATE_MHC_HAMMER_TABLE; 
        CREATE_KNOWN_SPLICING_COHORT_TABLE; CREATE_NOVEL_SPLICING_COHORT_TABLE } from '../modules/local/make_overview_table'

/*
========================================================================================
    IMPORT NF-CORE MODULE
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

// import custom functions
include { check_hlahd_results; check_tumour_results } from '../lib/core_functions'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/


workflow MHC_HAMMER {

    ch_versions = Channel.empty()

    // Define empty channels here so the RNA analysis can be run independently of the DNA analysis
    cohort_mosdepth_input_ch = Channel.empty()
    rna_all_snps_ch = Channel.empty()
    rna_bam_read_count_ch = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK ( ch_input, mhc_fasta_ch )
    
    // Add software used in INPUT_CHECK
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
 
    // These channels will be used to prevent pipeline stalling
    tumour_wxs_sample_count_ch = INPUT_CHECK.out.tumour_wxs_sample_count 

    rna_sample_count_ch = INPUT_CHECK.out.rna_sample_count
    
    germline_sample_count_ch = INPUT_CHECK.out.germline_sample_count

    normal_rna_count = INPUT_CHECK.out.normal_rna_sample_count

    gl_samples = INPUT_CHECK.out.gl_samples

    purity_ploidy_ch = INPUT_CHECK.out.purity_ploidy

    total_tumour_count = INPUT_CHECK.out.total_tumour_sample_count

    //
    // SUBWORKFLOW: Run Preprocessing steps - gererates bams / fqs containing putative HLA aligned reads, runs flagstat on input bams
    // 

    if ( params.run_bam_subsetting ) {

    PREPROCESSING ( 
        INPUT_CHECK.out.bams_ch, 
        germline_sample_count_ch,
        mhc_coords_ch,
        contigs_file_ch
    )

    // update downstream input channels
    preprocessing_versions = PREPROCESSING.out.versions
    hlahd_input_ch = PREPROCESSING.out.hlahd_input
    wxs_fqs_ch = PREPROCESSING.out.wxs_fqs
    dna_library_size_ch = PREPROCESSING.out.dna_library_size
    rna_fqs_ch = PREPROCESSING.out.rna_fqs
    rna_bams_ch = PREPROCESSING.out.rna_bams
    rna_library_size_ch = PREPROCESSING.out.rna_library_size

    } else {
    
    PREPROCESSING_WITHOUT_BAM_SUBSETTING ( 
        INPUT_CHECK.out.bams_ch, 
        germline_sample_count_ch,
        mhc_coords_ch
    )
    // update downstream input channels 
    preprocessing_versions = PREPROCESSING_WITHOUT_BAM_SUBSETTING.out.versions
    hlahd_input_ch = PREPROCESSING_WITHOUT_BAM_SUBSETTING.out.hlahd_input
    wxs_fqs_ch = PREPROCESSING_WITHOUT_BAM_SUBSETTING.out.wxs_fqs
    dna_library_size_ch = INPUT_CHECK.out.dna_library_size
    rna_fqs_ch = PREPROCESSING_WITHOUT_BAM_SUBSETTING.out.rna_fqs
    rna_bams_ch = PREPROCESSING_WITHOUT_BAM_SUBSETTING.out.rna_bams
    rna_library_size_ch = INPUT_CHECK.out.rna_library_size
    }

    // Add software used in PREPROCESSING
    ch_versions = ch_versions.mix(preprocessing_versions)

    //
    // MODULE: Run HLAHD on WXS normal sample
    //

    // only run the rest of the pipeline if the run_preprocessing_only flag is false
    if ( !params.run_preprocessing_only ) {

        if ( params.run_hlahd) {

        if ( params.hlahd_local_install ) {
            HLAHD_LOCAL ( hlahd_input_ch,  mhc_gtf_ch )

            // define outputs
            hlahd_genotype_ch = HLAHD_LOCAL.out.genotype
            ch_versions = ch_versions.mix(HLAHD_LOCAL.out.versions)

        } else {
            HLAHD ( hlahd_input_ch,  mhc_gtf_ch )

            // define outputs
            hlahd_genotype_ch = HLAHD.out.genotype
            ch_versions = ch_versions.mix(HLAHD.out.versions)
        }

    } else { // if not running hlahd, use user provided hla alleles

        // define outputs
        hlahd_genotype_ch = INPUT_CHECK.out.hla_alleles

    }

    // Now we check the output of HLAHD to see if it has failed for any samples
    // Failed samples will be excluded from the rest of the pipeline
    // Failed samples will have no typed HLA genes in their HLAHD output

    // Define branch criteria
    def criteria = branchCriteria {
        failed_hlahd_samples: check_hlahd_results(it)
        successful_hlahd_samples: !check_hlahd_results(it)
    }

    // Apply the branch criteria to the HLAHD genotype output
    hlahd_samples_by_status = hlahd_genotype_ch
        .branch(criteria)

    // output a warning if any samples failed HLAHD
    hlahd_samples_by_status.failed_hlahd_samples.view { 
    patient_id, file -> 
    log.warn("HLAHD failed to produce HLA allele estimates for sample ${patient_id}.\n" +
         "As a result, MHCHammer will not be executed for this sample.\n" +
         "The issue may be due to overly strict subset_bam process parameters. Consider allowing more reads to be included.")
    }
    
    // MODULE: generates sample specific reference files for the rest of the pipeline
    GENERATE_REFERENCES ( 
        hlahd_samples_by_status.successful_hlahd_samples, 
        mhc_fasta_ch, 
        mhc_gtf_ch, 
        mhc_transcriptome_fasta_ch 
    )

    ch_versions = ch_versions.mix(GENERATE_REFERENCES.out.versions.first())

    // SUBWORKFLOW: align wxs fqs to personalised reference and generate allele specific bams
    GENERATE_DNA_NOVOALIGN_BAMS (
        wxs_fqs_ch,
        GENERATE_REFERENCES.out.genome_reference,
        GENERATE_REFERENCES.out.mosdepth_reference
    )
    
    // ch_versions = ch_versions.mix(GENERATE_DNA_NOVOALIGN_BAMS.out.versions)

    // SUBWORKFLOW: detect HLA allele copy number and allelic imbalance
    DNA_ANALYSIS (
        dna_library_size_ch,
        GENERATE_REFERENCES.out.genome_snp_positions,
        GENERATE_REFERENCES.out.patient_gtf,
        purity_ploidy_ch,
        GENERATE_DNA_NOVOALIGN_BAMS.out.hla_allele_bams_ch,
        GENERATE_DNA_NOVOALIGN_BAMS.out.passed_heterozygous_hla_genes_ch
    )
    
    ch_versions = ch_versions.mix(DNA_ANALYSIS.out.versions)

    // SUBWORKFLOW: detect HLA allele mutations
    DETECT_MUTATIONS (
        GENERATE_DNA_NOVOALIGN_BAMS.out.hla_allele_bams_ch,
        GENERATE_DNA_NOVOALIGN_BAMS.out.passed_heterozygous_hla_alleles_ch,
        GENERATE_REFERENCES.out.mutation_calling_reference,
        tumour_wxs_sample_count_ch,
        germline_sample_count_ch,
        gl_samples,
        INPUT_CHECK.out.checked_inventory
    )

    ch_versions = ch_versions.mix(DETECT_MUTATIONS.out.versions)

    if ( params.run_rna_analysis ) {

        // SUBWORKFLOW: align rna fqs to personalised reference and generate allele specific bams
               
        GENERATE_RNA_NOVOALIGN_BAMS (
            rna_fqs_ch,
            GENERATE_REFERENCES.out.transcriptome_reference,
            GENERATE_REFERENCES.out.mosdepth_exons_reference
        )
        
        ch_versions = ch_versions.mix(GENERATE_RNA_NOVOALIGN_BAMS.out.versions)

        //
        // SUBWORKFLOW: Detect HLA allele expression and allelic imbalance in all samples
        //              Detect HLA allelic repression in samples with a match normal
        //

        NOVOALIGN_RNA_ANALYSIS (
            rna_library_size_ch,
            GENERATE_REFERENCES.out.transcriptome_snp_positions,
            GENERATE_REFERENCES.out.transcriptome_allele_tables,
            GENERATE_RNA_NOVOALIGN_BAMS.out.hla_allele_bams_ch,
            GENERATE_RNA_NOVOALIGN_BAMS.out.passed_heterozygous_hla_genes_ch,
            GENERATE_RNA_NOVOALIGN_BAMS.out.passed_heterozygous_hla_alleles_ch,
            GENERATE_RNA_NOVOALIGN_BAMS.out.hla_bam_read_count_ch,
            INPUT_CHECK.out.checked_inventory,
            "transcriptome",
            GENERATE_REFERENCES.out.star_input,
            "novoalign"
        )

        ch_versions = ch_versions.mix(NOVOALIGN_RNA_ANALYSIS.out.versions)

        // overwrite rna_all_snps_ch with the output from NOVOALIGN_RNA_ANALYSIS when the rnaseq analysis is run
        rna_all_snps_ch = NOVOALIGN_RNA_ANALYSIS.out.rna_all_snps_output_ch

        GENERATE_RNA_STAR_BAMS (
            GENERATE_REFERENCES.out.mosdepth_reference,
            rna_bams_ch,
            GENERATE_REFERENCES.out.star_input,
            mhc_kmer_ch
        )

        STAR_RNA_ANALYSIS (
            rna_library_size_ch,
            GENERATE_REFERENCES.out.genome_snp_positions,
            GENERATE_REFERENCES.out.genome_allele_tables,
            GENERATE_RNA_STAR_BAMS.out.hla_allele_bams_ch,
            GENERATE_RNA_STAR_BAMS.out.passed_heterozygous_hla_genes_ch,
            GENERATE_RNA_STAR_BAMS.out.passed_heterozygous_hla_alleles_ch,
            GENERATE_RNA_STAR_BAMS.out.hla_bam_read_count_ch,
            INPUT_CHECK.out.checked_inventory,
            "genome",
            GENERATE_REFERENCES.out.star_input,
            "star"
        )

        // updates the mosdepth cohort channel only if RNA samples are being run
        cohort_mosdepth_input_ch = cohort_mosdepth_input_ch
            .mix(GENERATE_RNA_STAR_BAMS.out.mosdepth_ch)
            .mix(GENERATE_RNA_NOVOALIGN_BAMS.out.mosdepth_ch)

        // overwrite rna_bam_read_count_ch with the output from GENERATE_RNA_NOVOALIGN_BAMS when the rnaseq analysis is run
        rna_bam_read_count_ch = GENERATE_RNA_NOVOALIGN_BAMS.out.hla_bam_read_count_ch

        if ( params.run_alt_splicing_analysis ) {

            // SUBWORKFLOW: Generate STAR aligned bams to detect exon skipping and intron retention

            ALT_SPLICING ( 
                GENERATE_RNA_STAR_BAMS.out.splice_junction_ch,
                GENERATE_REFERENCES.out.star_input,
                codon_table_ch
            )

            // Collate known splicing cohort table
            cohort_known_splicing_input_ch = ALT_SPLICING.out.known_sjs_ch.collect()
            cohort_known_splicing_input_file = cohort_known_splicing_input_ch
                    .flatten()
                    .map { it.getName() }
                    .collectFile(name: 'input_csvs.txt', newLine: true)
            
            CREATE_KNOWN_SPLICING_COHORT_TABLE(
                cohort_known_splicing_input_ch,
                INPUT_CHECK.out.checked_inventory,
                cohort_known_splicing_input_file
            )    

            // Collate novel splicing cohort table
            cohort_novel_splicing_input_ch = ALT_SPLICING.out.novel_sjs_ch.collect()
            cohort_novel_splicing_input_file = cohort_novel_splicing_input_ch
                    .flatten()
                    .map { it.getName() }
                    .collectFile(name: 'input_csvs.txt', newLine: true)
            CREATE_NOVEL_SPLICING_COHORT_TABLE(
                cohort_novel_splicing_input_ch,
                INPUT_CHECK.out.checked_inventory,
                cohort_novel_splicing_input_file
            )           

            // Collate tumour/normal enrichment splicing cohort table
            // ALT_SPLICING.out.tumour_normal_sjs_ch.first().view()
            cohort_tumour_normal_splicing_input_ch = ALT_SPLICING.out.tumour_normal_sjs_ch.collect()
            cohort_tumour_normal_splicing_input_file = cohort_tumour_normal_splicing_input_ch
                    .flatten()
                    .map { it.getName() }
                    .collectFile(name: 'input_csvs.txt', newLine: true)
            CREATE_TUMOUR_NORMAL_SPLICING_COHORT_TABLE(
                cohort_tumour_normal_splicing_input_ch,
                INPUT_CHECK.out.checked_inventory,
                cohort_tumour_normal_splicing_input_file
            )           

            ch_versions = ch_versions.mix(ALT_SPLICING.out.versions)
        }
    }
    
    // Collate mosdepth output
   cohort_mosdepth_input_ch = cohort_mosdepth_input_ch
            .mix(GENERATE_DNA_NOVOALIGN_BAMS.out.mosdepth_ch)
            .flatten()
            .collect()
            .map{ csvs ->
                sorted_csvs = csvs.sort()
                unique_csvs = sorted_csvs.unique { it.toString().tokenize('/').last() }
                tuple(unique_csvs)
                }

    cohort_mosdepth_input_text_file = cohort_mosdepth_input_ch
            .flatten()
            .map { it.getName() }
            .collectFile(name: 'input_csvs.txt', newLine: true)

    CREATE_MOSDEPTH_COHORT_TABLE(
        cohort_mosdepth_input_ch,
        INPUT_CHECK.out.checked_inventory,
        cohort_mosdepth_input_text_file
    )

    // collate library sizes table - only if flagstat is run
    if ( params.run_bam_subsetting ) {

        library_size_input_ch = PREPROCESSING.out.unmapped_count_ch
                    .mix(PREPROCESSING.out.mapped_count_ch)
                    .flatten()
                    .collect()
                    .map{ csvs ->
                        sorted_csvs = csvs.sort()
                        unique_csvs = sorted_csvs.unique { it.toString().tokenize('/').last() }
                        tuple(unique_csvs)
                        }

            // Create text file with input csvs for final table
            library_size_text_file = library_size_input_ch
                    .flatten()
                    .map { it.getName() }
                    .collectFile(name: 'input_csvs.txt', newLine: true)

            CREATE_LIBRARY_SIZE_TABLE(
                library_size_input_ch,
                INPUT_CHECK.out.checked_inventory,
                library_size_text_file
            )

    }
    
    // create cohort mutations table
        

    // collate overview table
    genome_allele_tables_ch = GENERATE_REFERENCES.out.genome_allele_tables
            .map{patient_id, genome_allele_tables -> tuple(genome_allele_tables)}

    transcriptome_allele_tables_ch = GENERATE_REFERENCES.out.transcriptome_allele_tables
            .map{patient_id, transcriptome_allele_tables -> tuple(transcriptome_allele_tables)}
    
    dna_library_size_for_output = dna_library_size_ch
            .map{sample_name, library_size_path -> tuple(library_size_path)}
    rna_library_size_for_output = rna_library_size_ch
            .map{sample_name, library_size_path -> tuple(library_size_path)}            

    allele_table_input_ch = rna_bam_read_count_ch
            .mix(GENERATE_DNA_NOVOALIGN_BAMS.out.hla_bam_read_count_ch)
            .mix(DNA_ANALYSIS.out.dna_all_snps_output_ch)
            .mix(rna_all_snps_ch)
            .mix(genome_allele_tables_ch)
            .mix(transcriptome_allele_tables_ch)
            .mix(rna_library_size_for_output)
            .mix(dna_library_size_for_output)
            .flatten()
            .collect()
            .map{ csvs ->
                sorted_csvs = csvs.sort()
                unique_csvs = sorted_csvs.unique { it.toString().tokenize('/').last() }
                tuple(unique_csvs)
                }
    
    // Create text file with input csvs for final table
    allele_input_text_file = allele_table_input_ch
            .flatten()
            .map { it.getName() }
            .collectFile(name: 'input_csvs.txt', newLine: true)
    
    // We also want a text file that has the germline samples that were used for hlahd
    hlahd_gl_samples_ch = hlahd_input_ch.map{meta, fastqs -> tuple(meta.sample_id)}
                .flatten()
                .collect()
                .map{
                    gl_sample_names -> tuple(gl_sample_names)
                }
    
    hlahd_germline_samples_text_file = hlahd_gl_samples_ch
            .flatten()
            .collectFile(name: 'hlahd_germline_samples.txt', newLine: true)

    CREATE_MHC_HAMMER_TABLE(
        allele_table_input_ch,
        INPUT_CHECK.out.checked_inventory,
        allele_input_text_file,
        hlahd_germline_samples_text_file,
        params.min_depth
    )

    }
    

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

}

/*
========================================================================================
    COMPLETION SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/