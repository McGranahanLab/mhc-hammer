// This is a slightly modified version of the original process from nf-core
process CUSTOM_DUMPSOFTWAREVERSIONS { 
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    container 'https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0'

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}