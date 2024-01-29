# Creating the MHC reference files

The reference files used by MHC Hammer are created from data downloaded from the IMGT database and avaliable in the assets folder of the MHC Hammer git repository. The instructions to recreate these files are below.

### Clone this repo
The scripts to create the MHC Hammer references are included in this repo.

```bash
git clone git@github.com:McGranahanLab/mhc-hammer.git
cd mhc-hammer
project_dir=${PWD}
```

### Download the preprocessing singularity container
This singularity container contains the tools needed to create the MHC Hammer references.

```bash
mkdir -p ${project_dir}/singularity_images
cd ${project_dir}/singularity_images
singularity pull --arch amd64 library://tpjones15/default/mhc_hammer_preprocessing:latest
mhc_hammer_preprocessing_sif="${project_dir}/singularity_images/mhc_hammer_preprocessing_latest.sif"
```

### Download the IMGT data

Download the IMGT data from github (https://github.com/ANHIG/IMGTHLA):
```bash
cd ${project_dir}"/assets/"
singularity exec -B ${project_dir}:${project_dir} $mhc_hammer_preprocessing_sif git clone https://github.com/ANHIG/IMGTHLA.git
```
The IMGT data is regularity updated and the above command will automatically clone the latest release. If you want to use a specific release you can `git checkout` the release branch. For example to use v3.38:

```bash
cd IMGTHLA
singularity exec -B ${project_dir}:${project_dir} $mhc_hammer_preprocessing_sif git checkout 3380
```

### Create the HLA reference files from the IMGT data

Once you have downloaded the IMGT database you then need to run two R scripts that create the HLA allele reference files used in the MHC Hammer pipeline. Note that these scripts have been run on 

**Step 1:** Run `scripts/convert_hla_dat.R`. This script parses the `hla.dat` file downloaded with the IMGT data and converts it into two csv files.
```bash
singularity exec -B ${project_dir}:${project_dir} $mhc_hammer_preprocessing_sif Rscript ${project_dir}/scripts/convert_hla_dat.R \
    --path_to_hla_dat ${project_dir}/assets/IMGTHLA/hla.dat \
    --save_dir ${project_dir}/assets/mhc_references \
    --functions_file ${project_dir}/scripts/mhc_reference_functions.R
```
The script `convert_hla_dat.R` will convert the information in the `hla.dat` file to two csv files containing information taken from the `hla.dat` file:
- `all_allele_features.csv`. This contains the avaliable sequence information of different features (intron, exon, UTR) for each HLA allele. This includes the columns:
    - start - the start position of the feature.
    - end - the end position of the feature.
    - type - either intron, exon or UTR.
    - number - The exon/intron numbers. These are counted from the 5' end.
    - length - The length of the feature.
    - allele_name - The name of the allele.
    - gene - The name of the gene.
    - feature_seq. The sequence of the feature.
- `all_allele_info.csv`. This contains information from  on the HLA alleles. This includes the columns:
    - allele_name - name of the allele. E.g. A*01:01:01:01
    - cell_lines - Cell lines from which the sequence was obtained.
    - gene - the gene name. E.g. HLA-A
    - ethnic - Possible ethnic origin of the cells sequenced for this allele.
    - partial - Differentiates between complete and partial region.
    - cds_line - The line from hla.dat that defined the known positions of the CDS sequence. 
    - translation - The translated protein. 
    - seq_length - the length of the sequence. 
    - seq - the avaliable sequence of the allele. 

For more information on hla.dat see https://github.com/ANHIG/IMGTHLA/blob/Latest/Manual.md.

**Step 2:** Run `create_mhc_hammer_references.R`. This script creates the reference files required to run MHC Hammer.

```bash
singularity exec -B ${project_dir}:${project_dir} $mhc_hammer_preprocessing_sif Rscript ${project_dir}/scripts/create_mhc_hammer_references.R \
    --imgt_dir ${project_dir}/assets//IMGTHLA/ \
    --all_allele_features_path ${project_dir}/assets/mhc_references/all_allele_features.csv \
    --strand_info ${project_dir}/assets/strand_info.txt \
    --functions_file ${project_dir}/scripts/mhc_reference_functions.R \
    --outdir ${project_dir}/assets/mhc_references/
```
This creates the following:
* ./assets/mhc_references/gtf/mhc.gtf - GTF file with all HLA alleles in the IMGT database
* ./assets/mhc_references/genome/mhc_genome.fasta - the complete genome (DNA) sequence for all alleles in the IMGT database
* ./assets/mhc_references/genome/mhc_genome_strand.fasta - the complete genome (DNA) sequence for all alleles in the IMGT database. This version is strand specific.
* ./assets/mhc_references/transcriptome/ the complete transcriptome (RNA) sequence for all alleles in the IMGT database.

### Creating the HLA kmer files
The HLA kmer files can be used by MHC Hammer to subset the input BAM files. To create the HLA kmer files you first run an R script that reads in all IMGT `*nuc` and `*gen` fasta files in the IMGT database directory and saves them in one single fasta file:
```bash
singularity exec -B ${project_dir}:${project_dir} $mhc_hammer_preprocessing_sif Rscript ${project_dir}/scripts/make_fasta_file.R \
    --imgt_dir ${project_dir}/assets/IMGTHLA \
    --out_dir ${project_dir}/assets/kmer_files
```

Next run jellyfish on that fasta file to make the kmer file. 
```bash
singularity exec -B ${project_dir}:${project_dir} $mhc_hammer_preprocessing_sif jellyfish count \
    -m 30 \
    --size 100M \
    --output ${project_dir}/assets/kmer_files/jellyfish_counts \
    ${project_dir}/assets/kmer_files/all_fasta.fasta

singularity exec -B ${project_dir}:${project_dir} $mhc_hammer_preprocessing_sif jellyfish dump \
    ${project_dir}/assets/kmer_files/jellyfish_counts | grep -v '^>' > ${project_dir}/assets/kmer_files/imgt_30mers.fa

```
