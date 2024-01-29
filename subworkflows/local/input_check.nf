//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK; CHECK_HLA_TYPE_INPUT } from '../../modules/local/samplesheet_check'
include { filterByMetadataField } from '../../lib/core_functions'

workflow INPUT_CHECK {

    take:
    samplesheet // file: /path/to/samplesheet.csv
    mhc_fasta_ch

    main:
    // These channels are specific for when the input is subset bam files  
    dna_library_size = Channel.empty()
    rna_library_size = Channel.empty()
    hla_alleles = Channel.empty()
    versions = Channel.empty()
    
    // Run process
    SAMPLESHEET_CHECK ( samplesheet )

    versions = versions.mix(SAMPLESHEET_CHECK.out.versions)

    if ( !params.run_bam_subsetting ) {
    // flagstat channel w matched norm
    dna_library_size = SAMPLESHEET_CHECK.out.csv
                            .splitCsv ( header:true, sep:',' )
                            .filter { row -> row.sequencing_type == 'wxs' }
                            .map { create_library_size_channels(it) }   
                            .unique()

    // flagstat channel without matched norm
    rna_library_size = SAMPLESHEET_CHECK.out.csv
                            .splitCsv ( header:true, sep:',' )
                            .filter { row -> row.sequencing_type == 'rnaseq' }
                            .map { create_library_size_channels(it) }   
                            .unique()
    }

    if ( !params.run_hlahd ) {

        // get hla alleles for each sample
        hla_alleles_to_check = SAMPLESHEET_CHECK.out.csv
                .splitCsv ( header:true, sep:',' )
                .map { create_hlahd_channels(it) }   
                .unique()

        hla_allele_files = hla_alleles_to_check
                .map { patient_id, hla_alleles -> tuple(hla_alleles) }
                .collect()
                .flatten()

        CHECK_HLA_TYPE_INPUT ( mhc_fasta_ch, hla_allele_files, SAMPLESHEET_CHECK.out.csv )

        hla_alleles = CHECK_HLA_TYPE_INPUT.out.validated_samplesheet
                .splitCsv(header: true, sep: ',')
                .map { create_hlahd_channels(it) }
                .unique()

        versions = versions.mix(CHECK_HLA_TYPE_INPUT.out.versions)
    }

    emit:

    // Emit empty channel or flagstat channel [sample_id, flagstat] if using subset bams as input
    dna_library_size
    rna_library_size

    // Emit empty channel or hla alleles channel [patient_id, hla_alleles.csv] if user hla alleles are provided
    hla_alleles

    // parsed input.csv
    checked_inventory = SAMPLESHEET_CHECK.out.csv

    // get purity and ploidy information for each sample
    purity_ploidy = SAMPLESHEET_CHECK.out.csv
            .splitCsv ( header:true, sep:',' )
            .filter { row -> row.sequencing_type == 'wxs' }
            .map{ row -> tuple(row.sample_name, tuple(row.purity, row.ploidy)) }

    // generate bam file channel
    bams_ch = SAMPLESHEET_CHECK.out.csv
            .splitCsv ( header:true, sep:',' )
            .map { create_bam_channels(it) }
            .unique()

    // generate channel with list of wxs gl samples for each patient [patient_id, [gl_id_1, gl_id_2 ...]]
    gl_samples = bams_ch
            .filter(filterByMetadataField("seq", "wxs"))
            .filter(filterByMetadataField("sample_type", "normal"))
            .map{ meta, bam -> tuple(meta.patient_id, meta.sample_id) }
            .groupTuple()

    // get number of rna samples per patient
    rna_sample_count = SAMPLESHEET_CHECK.out.csv
            .splitCsv ( header:true, sep:',' )
            .filter({row -> row.sequencing_type == 'rnaseq'})
            .map{row -> tuple(row.patient, file(row.sample_name))}
            .groupTuple(by:0)
            .map{ patient, list_F -> tuple( patient, list_F.size() ) }

    // get number of tumour wxs samples per patient
    tumour_wxs_sample_count = bams_ch
            .filter(filterByMetadataField("seq", "wxs"))
            .filter(filterByMetadataField("sample_type", "tumour"))
            .map {meta, bam -> tuple(meta.patient_id, file(bam[0])) }
            .groupTuple()
            .map{ patient, list_F -> tuple( patient, list_F.size() ) }

    // get number of gl wxs samples per patient
    germline_sample_count = bams_ch
            .filter(filterByMetadataField("seq", "wxs"))
            .filter(filterByMetadataField("sample_type", "normal"))
            .map {meta, bam -> tuple(meta.patient_id, file(bam[0])) }
            .groupTuple()
            .map{ patient, list_F -> tuple( patient, list_F.size() ) }

    normal_rna_sample_count = bams_ch
            .filter(filterByMetadataField("seq", "rnaseq"))
            .filter(filterByMetadataField("sample_type", "normal"))
            .map {meta, bam -> tuple(meta.patient_id, file(bam[0])) }
            .groupTuple()
            .map{ patient, list_F -> tuple( patient, list_F.size() ) }

        // get number of tumour samples 
        total_tumour_sample_count = SAMPLESHEET_CHECK.out.csv
                .splitCsv ( header:true, sep:',' )
                .filter({row -> row.sample_type == 'tumour'})
                // count the number of samples
                .count()        
    versions
}

def create_bam_channels(LinkedHashMap row) {
    def meta = [:]
    meta.patient_id     = row.patient
    meta.sample_id      = row.sample_name
    meta.sample_type    = row.sample_type
    meta.seq            = row.sequencing_type
    meta.normal_sample_name  = row.normal_sample_name

    def bam_array = []
    if (!file(row.bam_path).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> BAM file does not exist!\n${row.bam_path}"
    }

    // Check bam is indexed
    def dir_name =  file(row.bam_path).getParent()
    def simple_name = file(row.bam_path).getBaseName()
    def bai_files = [file("${dir_name}/${simple_name}.bai"), file("${dir_name}/${simple_name}.bam.bai")]

    if (!bai_files.any { it.exists() }) {
        exit 1, "ERROR: Please check if BAM file is indexed!\n${bai_files.join(' / ')} not found!"
    }

    def bai_file = bai_files.findResult { it.exists() ? it : null }
    bam_array = [ meta, [ file(row.bam_path), file("${bai_file}") ] ]

    return bam_array

}

def create_library_size_channels(LinkedHashMap row) {
    def library_size_array = []
    if (!file(row.library_size_path).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> You've run the pipeline using --run_bam_subsetting false, but the following library size file does not exist!\n${row.library_size_path}"
    }
    library_size_array = [ row.sample_name, file(row.library_size_path) ]

    return library_size_array

}

def create_hlahd_channels(LinkedHashMap row) {
        def hla_alleles_array = []
        if (!file(row.hla_alleles_path).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> You've run the pipeline using --run_hlahd false, but the following hla alleles file does not exist!\n${row.hla_alleles_path}"
        }
        hla_alleles_array = [ row.patient, file(row.hla_alleles_path) ]
        
        return hla_alleles_array
}
