## Introduction

MHC Hammer is a bioinformatics pipeline designed for the analysis of the class I HLA genes using WES and RNAseq data. 

MHC Hammer requires for each patient a germline WES sample. From this, the patient's HLA allele types are predicted using [HLA-HD](https://pubmed.ncbi.nlm.nih.gov/28419628/). Then, depending on what is else input to the pipeline, the following analysis can be run:

To estimate DNA HLA allelic imbalance and somatic mutations:
- A tumour WES BAM file.
To estimate DNA HLA copy number (and LOH):
- A tumour WES BAM file with corresponding purity and ploidy estimates. 
To estimate RNA HLA allelic expression, allelic imbalance and alternative splicing:
- A tumour or normal RNAseq BAM file.
To estimate RNA HLA allelic repression or tumour/normal enrichment of alternative splicing events:
- A tumour and normal tissue RNAseq BAM file. The normal sample should be from the same patient and tissue as the tumour.

## Steps before running the pipeline
 
### Installing Nextflow and singularity
1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/)

### Make an inventory file

You need to create a inventory file with the following columns:
- patient - the patient name. MHC Hammer will replace spaces in the patient name with underscores. Required.
- sample_name - the sample name. MHC Hammer will replace spaces in the sample name with underscores. Required.
- sample_type - either `tumour` or `normal`.  Required.
- bam_path - full path to the wxs or rnaseq BAM file. Required.
- sequencing_type - either `wxs` or `rnaseq`.  Required.
- purity - the purity of the tumour region. Can be left empty.
- ploidy - the ploidy of the tumour region. Can be left empty.
- normal_sample_name - when sequencing_type is `wxs` this is the matched germline wxs. When sequencing_type is `rnaseq` this is the matched rnaseq normal name. 

The inventory should be a csv file and is input to the pipeline with the `--input` parameter. 

The following is an example inventory for a single patient with:
- two tumour regions with wxs (sample_name1 and sample_name2), one of which has rnaseq (sample_name1)
- one germline wxs sample (sample_name3)
- one normal rnaseq sample (sample_name4)

|  patient |  sample_name | sample_type | bam_path                 | sequencing_type | purity | ploidy | normal_sample_name |
| :------: | :----------: | :---------: | :----------------------: | :-------------: | :----: | :----: | :--------------: |
| patient1 | sample_name1 |   tumour    | path/to/sample_name1.bam |       wxs       |   0.5  |    3   |   sample_name3   |
| patient1 | sample_name2 |   tumour    | path/to/sample_name2.bam |       wxs       |   0.3  |   2.5  |   sample_name3   |
| patient1 | sample_name3 |   normal    | path/to/sample_name3.bam |       wxs       |        |        |                  |
| patient1 | sample_name1 |   tumour    | path/to/sample_name4.bam |       rnaseq    |        |        |   sample_name4   |
| patient1 | sample_name4 |   normal    | path/to/sample_name5.bam |       rnaseq    |        |        |                  |

### Clone this repo
```bash
git clone git@github.com:McGranahanLab/mhc-hammer.git
cd mhc-hammer
project_dir=${PWD}
```
### Pull the preprocessing singularity image
This singularity image contains software that is needed to make the input files for MHC Hammer, including git lfs, R and jellyfish. 

```bash
mkdir ${project_dir}/singularity_images
cd singularity_images
singularity pull --arch amd64 library://tpjones15/default/mhc_hammer_preprocessing:latest
mhc_hammer_preprocessing_sif="${project_dir}/singularity_images/mhc_hammer_preprocessing_latest.sif"
```

### Download the IMGT data
Download the IMGT data from github (https://github.com/ANHIG/IMGTHLA):
```bash
mkdir ${project_dir}/assets/IMGT/
cd ${project_dir}/assets/IMGT/
singularity exec -B ${PWD}:${PWD} $mhc_hammer_preprocessing_sif git clone https://github.com/ANHIG/IMGTHLA.git
```
The IMGT data is regularity updated and the above command will automatically clone the latest release. If you want to use a specific release you can `git checkout` the release branch. For example to use v3.38:

```bash
cd ${project_dir}/assets/IMGTHLA
singularity exec -B ${PWD}:${PWD} $mhc_hammer_preprocessing_sif git checkout 3380
```

### Create the HLA reference files from the IMGT data

Once you have downloaded the IMGT database you then need to run two R scripts that create the HLA allele reference files used in the MHC Hammer pipeline. 

First run `scripts/convert_hla_dat.R`:
```bash
singularity exec -B ${PWD}:${PWD} $mhc_hammer_preprocessing_sif Rscript ${project_dir}/scripts/convert_hla_dat.R \
    --path_to_hla_dat ${project_dir}/assets/IMGT/IMGTHLA/hla.dat \
    --save_dir ${project_dir}/assets/mhc_references \
    --functions_file ${project_dir}/scripts/mhc_reference_functions.R
```
The script `convert_hla_dat.R` will convert the information in the `hla.dat` file to two csv files:
- `all_allele_features.csv` - this contains TODO
- `all_allele_info.csv` - this contains TODO

Then create the references using `create_mhc_hammer_references.R`

```bash
singularity exec -B ${PWD}:${PWD} $mhc_hammer_preprocessing_sif Rscript ${project_dir}/scripts/create_mhc_hammer_references.R \
    --imgt_dir ${project_dir}/assets/IMGT/IMGTHLA/ \
    --all_allele_features_path ${project_dir}/assets/mhc_references/all_allele_features.csv \
    --strand_info ${project_dir}/assets/strand_info.txt \
    --functions_file ${project_dir}/scripts/mhc_reference_functions.R \
    --outdir ${project_dir}/assets/mhc_references/
```
This creates the following:
* ./assets/mhc_references/gtf/mhc.gtf - GTF file with all HLA alleles in the IMGT database
* ./assets/mhc_references/genome/mhc_genome.fasta - the complete genome (DNA) sequence for all alleles in the IMGT database
* ./assets/mhc_references/genome/mhc_genome_strand.fasta - the complete genome (DNA) sequence for all alleles in the IMGT database. This version is strand specific which means that alleles on the reverse strand (HLA-B and HLA-C) TODO
* ./assets/mhc_references/transcriptome/ the complete transcriptome (RNA) sequence for all alleles in the IMGT database.

### Creating the HLA kmer files
The HLA kmer files can be used by MHC Hammer to subset the input BAM files. To create the HLA kmer files you first run an R script that reads in all IMGT `*nuc` and `*gen` fasta files in the IMGT database directory and saves them in one single fasta file:
```bash
singularity exec -B ${PWD}:${PWD} $mhc_hammer_preprocessing_sif Rscript scripts/make_fasta_file.R \
    --imgt_dir ${project_dir}/assets/IMGT/IMGTHLA/ \
    --out_dir ${project_dir}/assets/kmer_files
```

Then you run jellyfish on that fasta file to make the kmer file. 
<!--- TODO but probably doesnt make sense to save as fa file? --->
```bash
singularity exec -B ${PWD}:${PWD} $mhc_hammer_preprocessing_sif jellyfish count \
    -m 30 \
    --size 100M \
    --output jellyfish_counts \
    all_fasta.fasta

singularity exec -B ${PWD}:${PWD} $mhc_hammer_preprocessing_sif jellyfish dump jellyfish_counts | grep -v '^>' > imgt_30mers.txt

```

### Downloading HLA-HD and creating a singularity container
We use [HLA-HD](https://pubmed.ncbi.nlm.nih.gov/28419628/) within MHC Hammer to predict the HLA allele types of each sample. 

**We are unable to provide a singularity container for this tool.**

Instead, we provide two methods for installing HLA-HD and its dependencies locally, before running the pipeline. 

#### Option 1: Install HLA-HD and its dependencies locally (recommended)

The setps are as follows:
1. On the HLA-HD website fill in the [download request form](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/download-request/) to get a download link for HLA-HD
2. After downloading HLA-HD and cloning this repo, place your downloaded hlahd.version.tar.gz file into the project bin directory.
```bash
mv /path/to/hlahd_download ${project_dir}/bin/
```
3. Edit the [install_hlahd.sh](scripts/install_hlahd.sh) script for your desired IMGT database. If you are using the most up to date version of IMGT, comment out line 97. If you want to install a particular version, replace
```
https://media.githubusercontent.com/media/ANHIG/IMGTHLA/v3.38.0-alpha/hla.dat
```
on line 97 with the path to the specific IMGT version hla.dat file. The line above downloads the hla.dat file from IMGT version 3.38

4. Run the hlahd installation script

**This requires g++ and wget to be installed**
```bash
bash ${project_dir}/scripts/install_hlahd.sh -p ${project_dir} -h ${hlahd_download}
```
This script will install bowtie2 (2.5.1) and HLA-HD (whichever version was downloaded) and store them in ```${project_dir}/bin/``` directory. 

5. When running the pipeline ensure you run with ```--local_hlahd_install true```

#### Option 2: Create your own HLA-HD singularity container

The steps are as follows:

1. On the HLA-HD website fill in the [download request form](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/download-request/) to get a download link for HLA-HD
2. After downloading HLA-HD, edit the [hlahd_container.def](assets/hlahd_container_template.def) file in the assets directory as follows:
   1. Update the `/path/to/downloaded/hlahd.version.tar.gz` in the `%files` section 
   2. Update the `HLAHD_VERSION` variable in the `%post` section 
   3. Update the hla.dat file used by HLA-HD 
      - **Ensure the same IMGT hla.dat file is being used by HLA-HD as was used to generate the mhc_references**
      - If you are using the most up-to-date version of the IMGT database, then comment out line 41.
      - If you are using a specific version of the IMGT database, replace 
        ```
        https://media.githubusercontent.com/media/ANHIG/IMGTHLA/v3.38.0-alpha/hla.dat
        ```
        on line 41 with the path to the specific IMGT version hla.dat file. The line above downloads the hla.dat file from IMGT version 3.38
3. Build the singularity image:
   ```bash
   singularity build hlahd.sif hlahd_container_template.def
   ``` 
4. Move the image into the singularity_images directory
   ```bash
   mv hlahd.sif singularity_images
   ```
5. When running the MHC Hammer pipeline ensure you run with ```--hlahd_local_install false```.

### Update the config files

In the `.conf/hpc/config` file you can set parameters specific to your HPC system. For example:
* `workDir` - this is where intermediate files are saved. 
* `queueSize` - how many jobs are run at a given time.
* `max_memory` - the maximum memory avaliable on your HPC system.
* `max_cpus` - the maximum number of cpus avaliable on your HPC system.

In the `nextflow.config` file you can change the pipeline parameters from the default parameters. These parameters include:

**BAM subsetting options** 

The following parameters control what reads are included in the subsetted BAM file.
- `unmapped_reads`: can be `true` or `false`. If `true` unmapped reads in the input BAM file are included in the subsetted BAM file.
- `contig_reads`: can be `true` or `false`. If `true` reads that map to either:
  - any contig in the input BAM files (if the `contigs_file` parameter is `null`), or
  - specific contigs listed in the `contigs_file`
  are included in the subsetted BAM file.
- `contigs_file`: can be `null` or the path to a file. This file should contain a list of contigs present in the input BAM. Reads that map to these contigs will be included in the subsetted BAM file. 
- `fish_reads`: can be `true` or `false`. If `true` reads in the input BAM files that contain a kmer from the `imgt_30mers.txt` file are included in the subsetted BAM.
- `kmer_file`: The path the the HLA kmer file. [see the section "Creating the HLA kmer files".](###-Creating-the-HLA-kmer-files)
- `fish_reads_only`: can be `true` or `false`, If `true` only reads that contain a kmer from the `imgt_30mers.txt` file are included in the subsetted BAM.
- `mhc_coords`: path to the file containing a set of genomic coordinates (chr:star-stop). Any read in the input BAM file that lies within these coordinates is included in the subsetted BAM.
- `save_subset_bams`: Save the subsetted BAM files.

**Generate FQ job params**
- `save_hla_fqs`: Save the HLA fastq files. 

**HLA-HD parameters** 

The following are HLA-HD specific parameters. See [HLA-HD](https://pubmed.ncbi.nlm.nih.gov/28419628/).
- `hlahd_ignored_read_length`: A read whose length is shorter than this parameter is ignored.
- `hlahd_trimming_ratio`: Trimming option. If a match sequence is not found in the dictionary, trim the read until some sequence is matched to or reaches this ratio.
- `save_hlahd_predictions`: Save the HLA-HD predictions.
- `hlahd_local_install`- `true` if you have installed HLA-HD locally (option 1), `false` if you have created a singularity container.

**MHC Hammer reference files** 

These parameters provide the path to the reference files that MHC Hammer requires. They are created in the section "Create the HLA reference files from the IMGT data"
- `mhc_gtf`: Path to the GTF file.
- `mhc_fasta`: Path to the strand specific FASTA file.
- `mhc_transcriptome_fasta`: Path to the transcriptome FASTA file.
- `save_references`: if `true` then the patient specific HLA reference files will be saved

**Library size calculation**
- `include_unmapped_reads_in_library_size`: can be `true` or `false`. If true, unmapped reads are included in library size (e.g `samtools view -c -f 1 -F 2304`). If false, unmapped reads are not included in library size (`samtools view -c -f 1 -F 2308`).
- `save_flagstat`

**HLA BAM alignment parameters**
- `max_mismatch` - The maximum number of mismatches allowed in a sequencing read aligned to the patient specific HLA reference. If a read alignment has more than `max_mismatch` mismatches, then it will be removed from the filtered HLA BAM file.
- `save_unfiltered_novoalign_bam`: Save the unfiltered novoalign BAM file.
- `save_unfiltered_star_bam`: Save the unfiltered STAR bam file.
- `save_filtered_bams`: Save the filtered BAM file
- `save_filtered_read_counts`: Save the filtered read count file.

**DNA allelic AIB and copy number parameters**
- `min_depth` - SNPs with less than `min_depth` in the germline will not be used to calculate DNA AIB and copy number. 

**Mutation calling params**
- `save_vcfs`: Save the mutect VCF files.
- `save_vep_output`: Save the VEP output files.

**Alternative splicing parameters**
- `uniq_num_across_junc`: STAR parameter specifying the minimum number of unique reads supporting a junction in the first STAR run for the junction to be called.
- `codon_table` - Path to the file that contains the mapping between codons and amino acids. This should be included in the assests directory when the mhc_hammer repo is cloned. 
- `save_sample_sj_tab`: 
- `save_cohort_sj_tab`: 
- `save_novel_splice_junctions`: 
- `save_sample_kmer`: 

**Parameters for filtering the samples in the output table**
- `min_frac_mapping_uniquely` - Minimum fraction of reads mapping uniquely to each allele.
- `max_frac_mapping_multi_gene` - Maximum fraction of reads mapping to multiple genes.
- `min_number_of_snps` - Minimum number of SNPs that pass coverage.
- `max_copy_number_range` - Maximum range in the 95% confidence interval.
- `min_expected_depth` - Minimum expected depth.

## Outputs
Depending on the input files, the pipeline will output a number of different files.

Tables in the cohort_tables folder include:
* `novel_splicing_events.csv` - this contains the novel splicing events, i.e. splicing events that are not in the GTF file.
* `known_splicing_events.csv` - this contains the known splicing events, i.e. splicing events that are in the GTF file.
* `cohort_library_size.csv` - this contains the library size for every sample in the cohort.
* `mosdepth_rnaseq_star.csv` - this contains the mosdepth output for every allele specific STAR HLA BAM file.
* `mosdepth_novoalign_rnaseq_star.csv` - this contains the mosdepth output for every allele specific RNAseq NOVOALIGN HLA BAM file.
* `mosdepth_novoalign_wes_star.csv` - this contains the mosdepth output for every allele specific WES NOVOALIGN HLA BAM file.
* `cohort_mhc_hammer_gene_table.csv` - this table has a single line per gene and sample, and contains:
    * the DNA AIB, copy number and LOH calls
    * the RNA AIB and repression calls
    * Filters based on the WES and RNAseq data
* `cohort_mhc_hammer_allele_table.csv` - this table contains the same information as `cohort_mhc_hammer_gene_table.csv` but with a single line per allele and sample.

## Test dataset
A test dataset is provided. The input BAMs and inventory are in `./test_data`. To run the pipeline with the test dataset:
```bash
nextflow run main.nf -profile test,singularity --fish_reads_only true --run_alt_splicing_analysis false
```

## Citations
This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

 > The nf-core framework for community-curated bioinformatics pipelines.
 >
 > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
 >
 > Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
 >

An extensive list of references for the tools used by the pipeline can be found in the [CITATIONS.md](CITATIONS.md) file.
