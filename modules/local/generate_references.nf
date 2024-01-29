process GENERATE_REFERENCES {
    tag "${patient_id}"

    container "library://tpjones15/mhchammer/mhchammer_core:latest"

    label 'process_single'

    input:
    tuple val(patient_id), path(genotype)
    path mhc_fasta
    path mhc_gtf
    path mhc_transcriptome_fasta

    output:
    // mhchammer dna references
    tuple val(patient_id), path("*genome.snp_pos.bed")                    , optional:true, emit: genome_snp_positions
    tuple val(patient_id), path("${patient_id}.gtf")                      , emit: patient_gtf
    tuple val(patient_id), path("${patient_id}_mhc_genome_reference.fa*"),  
                           path("${patient_id}.mhc_genome_fasta.nix")     , emit: genome_reference
    tuple val(patient_id), path("${patient_id}_mhc_genome_reference*"),
                           path("${patient_id}_sorted.gtf.gz*")           , emit: mutation_calling_reference
    tuple val(patient_id), path("${patient_id}_genome_allele_table.csv")  , emit: genome_allele_tables

    // mhchammer rna references
    tuple val(patient_id), path("${patient_id}.gtf"), \
                           path("${patient_id}_mhc_genome_reference.fa"), 
                           env(genome_size)                                    , emit: star_input

    tuple val(patient_id), path("${patient_id}.bed")                           , optional:true, emit: mosdepth_reference
    tuple val(patient_id), path("${patient_id}.exons.bed")                     , optional:true, emit: mosdepth_exons_reference
    tuple val(patient_id), path("*transcriptome.snp_pos.bed")        , optional:true, emit: transcriptome_snp_positions
    tuple val(patient_id), path("${patient_id}_mhc_transcriptome_reference.fa"), 
                           path("${patient_id}.mhc_transcriptome_fasta.nix")   , optional:true, emit: transcriptome_reference 
    tuple val(patient_id), path("${patient_id}_transcriptome_allele_table.csv"), optional:true, emit: transcriptome_allele_tables

    path ("versions.yml")                                                      , emit: versions
    
    run_rna_analysis = "${params.run_rna_analysis}"
    script: 
    """
    if [ "${params.run_rna_analysis}" == "true" ]; then
        # Run Rscript to generate transcriptome allele mismatches
        Rscript --vanilla ${baseDir}/bin/allele_mismatch.R \
        --genome_or_transcriptome transcriptome \
        --patient_id ${patient_id} \
        --patient_hla_alleles ${genotype} \
        --mhc_fasta ${mhc_transcriptome_fasta}

        echo "Creating personalised transcriptome reference"
        # generate personalised transcriptome fasta
        Rscript --vanilla ${baseDir}/bin/generate_personalised_transcriptome_fasta.R \
        --patient_id ${patient_id} \
        --patient_hla_alleles ${patient_id}_transcriptome_allele_table.csv \
        --mhc_fasta_path ${mhc_transcriptome_fasta}

        echo "Indexing transcriptome reference"
        # index reference
        samtools faidx ${patient_id}_mhc_transcriptome_reference.fa
        echo "Done"

        echo "Creating transcriptome sequence dictionary"
        # create reference dictionary
        picard CreateSequenceDictionary -R ${patient_id}_mhc_transcriptome_reference.fa
        echo "Done"

        echo "Creating transcriptome novoindex"
        # create novoalign index
        novoindex ${patient_id}.mhc_transcriptome_fasta.nix ${patient_id}_mhc_transcriptome_reference.fa
        echo "Done"

        echo "Creating bed file"
        # create bed file for mosdepth 
        Rscript --vanilla ${baseDir}/bin/make_bed_file.R \
        --patient_hla_alleles ${genotype} \
        --patient_id ${patient_id} \
        --mhc_gtf ${mhc_gtf}
        echo "Done"
 
    fi

    # Run Rscript to generate genome allele mismatches
    Rscript --vanilla ${baseDir}/bin/allele_mismatch.R \
    --genome_or_transcriptome genome \
    --patient_id ${patient_id} \
    --patient_hla_alleles ${genotype} \
    --mhc_fasta ${mhc_fasta}
    
    echo "Creating personalised genome reference"
    # Run Rscript to generate genome gtf and fasta for star
    Rscript --vanilla ${baseDir}/bin/create_gtf_and_ref.R \
    --gtf_out_path ${patient_id}.gtf \
    --fa_out_path ${patient_id}_mhc_genome_reference.fa \
    --genome_out_path ${patient_id}_genome_size.txt \
    --hla_path ${genotype} \
    --mhc_fasta ${mhc_fasta} \
    --mhc_gtf ${mhc_gtf}
    echo "Done"

    genome_size=`cat ${patient_id}_genome_size.txt` 

    echo "Sorting, zipping and indexing personalised GTF"
    # sort gtf 
    sort -k1,1 -k4,4n -k5,5n -t\$'\\t' ${patient_id}.gtf > ${patient_id}_sorted.gtf

    # zip gtf
    bgzip ${patient_id}_sorted.gtf

    # index gtf
    tabix -p gff ${patient_id}_sorted.gtf.gz
    echo "Done"

    echo "Indexing genome reference"
    # index genome reference
    samtools faidx ${patient_id}_mhc_genome_reference.fa
    echo "Done"

    echo "Creating genome sequence dictionary"
    # create genome reference dictionary
    picard CreateSequenceDictionary -R ${patient_id}_mhc_genome_reference.fa
    echo "Done"

    echo "Creating genome novoindex"
    # create novoalign index
    novoindex ${patient_id}.mhc_genome_fasta.nix ${patient_id}_mhc_genome_reference.fa
    echo "Done"

    # Get R version and package versions
    R_VERSION=\$(Rscript -e "cat(as.character(getRversion()))")
    DT_VERSION=\$(Rscript -e "cat(paste(packageVersion('data.table'), collapse = ''))")
    AP_VERSION=\$(Rscript -e "cat(paste(packageVersion('argparse'), collapse = ''))")
    SEQINR_VERSION=\$(Rscript -e "cat(paste(packageVersion('seqinr'), collapse = ''))")
    BIOSTRING_VERSION=\$(Rscript -e "cat(paste(packageVersion('Biostrings'), collapse = ''))")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*tabix //; s/Copyright.*\$//')
        novoalign: \$(novoalign --version)
        R: \${R_VERSION}
        data.table: \${DT_VERSION}
        argparse: \${AP_VERSION}
        seqinr: \${SEQINR_VERSION}
        Biostrings: \${BIOSTRING_VERSION}
    END_VERSIONS
    """
}