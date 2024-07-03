# MHC Hammer outputs

By defult, the output is saved in the working directory in a folder called `mhc_hammer_results`. Depending on the input files, the pipeline will output the following:

The output tables in the `mhc_hammer_results/cohort_tables` folder are explained below.

***

**`cohort_mhc_hammer_gene_table.csv`**

This table combines the results of the WXS and RNAseq analysis (excluding the alternative splicing results). Each row represents a sample/gene. Every sample included in the inventory file should be included in this table. 

Sample information columns:
* `patient` - the patient name, defined in the input inventory file.
* `sample_name` - the sample name, defined in the input inventory file.
* `sample_type` - either 'tumour' or 'normal', as defined in the input inventory file.
* `gene` - either `HLA-A`, `HLA-B` or `HLA-C`.
* `allele1` - the full allele1 name. 
* `allele2` - the full allele2 name. 
* `short_allele1` - the allele1 name shortened to two digits.
* `short_allele2` - the allele2 name shortened to two digits.
* `homozygous` - `TRUE` if `allele1` is the same as `allele2`.
* `purity` - the purity of the tumour sample, as defined in the input inventory file.
* `ploidy` - the ploidy of the tumour sample, as defined in the input inventory file.
* `sample_has_wxs` - TRUE if the sample has WXS, as defined in the input inventory file.
* `sample_has_rna` - TRUE if the sample has RNAseq, as defined in the input inventory file.
* `wxs_fail` - Does the sample/gene pass the WXS quality filters? This column is `TRUE` if any of the following columns are `TRUE`: `fail_homozygous`, `fail_n_wxs_snps`, `fail_expected_depth`, `fail_cn_range`. If `sample_has_wxs` is `FALSE` this column is empty.
* `rna_fail` - Does the sample/gene pass the RNAseq quality filters? This column is `TRUE` if any of the following columns are `TRUE`: `fail_homozygous`, `fail_rna_n_snps`, `fail_multi_gene_mapping`, `fail_unique_allele_mapping`. If `sample_has_rna` is `FALSE` this column is empty.
* `rna_normal_fail` - Does the sample in `rnaseq_normal_sample_name` pass the sample/gene RNAseq quality filters? This column is `TRUE` if any of the following columns are `TRUE` for the `rnaseq_normal_sample_name` sample: `fail_homozygous`, `fail_rna_n_snps`, `fail_multi_gene_mapping`, `fail_unique_allele_mapping`. If `sample_has_rna` is `FALSE` or `rnaseq_normal_sample_name` is empty, this column is empty.
* `wxs_germline_id` - the name of the WXS germline sample, as defined in the inventory file.
* `hlahd_germline_sample_name` - the name of the WXS germline sample that is input to HLA-HD to predict HLA alleles.
* `rnaseq_normal_sample_name` - the name of the RNAseq matched normal sample, as defined by the inventory.

WXS results:
* `n_genome_snps` - the number of mismatches between the two alleles, measured using the entire genome sequence defined in `assets/mhc_references/genome/mhc_genome.fasta`.
* `wxs_library_size` - the library size used in the DNA analysis.
* `allele1_dna_n_reads_before_filtering` - The number of reads in the WXS HLA allele1 BAM file before the BAM file is filtered to remove reads with more than the input parameter `max_mismatch`.
* `allele1_dna_n_reads_after_filtering` - The number of reads in the WXS HLA allele1 BAM file after the BAM file is filtered to remove reads with more than the input parameter `max_mismatch`.
* `allele2_dna_n_reads_before_filtering` - The number of reads in the WXS HLA allele2 BAM file before the BAM file is filtered to remove reads with more than the input parameter `max_mismatch`.
* `allele2_dna_n_reads_after_filtering` - The number of reads in the WXS HLA allele2 BAM file after the BAM file is filtered to remove reads with more than the input parameter `max_mismatch`.
* `cn1_binned` - the allele1 copy number.
* `cn1_binned_lower` - the lower value of the 95% confidence interval for the cn1_binned.
* `cn1_binned_upper` - the upper value of the 95% confidence interval for the cn1_binned.
* `allele1_cn_range` - `cn1_binned_upper - cn1_binned_lower`
* `cn2_binned` - the allele2 copy number.
* `cn2_binned_lower` - the lower value of the 95% confidence interval for the cn2_binned.
* `cn2_binned_upper` - the upper value of the 95% confidence interval for the cn2_binned.
* `allele2_cn_range` - `cn2_binned_upper - cn2_binned_lower`
* `cn_n_bins` - the number of bins used to calculate the copy number.
* `cn_n_snps` - the number of filtered SNPs used to calculate the copy number.
* `allele1_expected_depth` - the expected depth of allele 1.
* `allele2_expected_depth` - the expected depth of allele 2.
* `logr_aib_paired_t_test` - the p-value of the log R t-test
* `logr_aib_paired_wilcoxon_test` - the p-value of the log R Wilcoxon test.
* `logr_aib_n_snps` - the number of filtered SNPs used to calculate the Log R AIB.
* `dna_aib` - is there DNA AIB in this sample/gene? `TRUE` if `logr_aib_paired_wilcoxon_test < 0.01`.
* `loh` - is there HLA LOH in this sample/gene? `TRUE` if `logr_aib_paired_wilcoxon_test < 0.01` and `cn1_binned < 0.5` and/or `cn2_binned < 0.5`.
* `minor_allele` - the allele with the lowest copy number.
* `allele1_loh` - is there HLA LOH of allele1 in this sample? `TRUE` if `dna_aib = TRUE` and `cn1_binned < 0.5`
* `allele2_loh` - is there HLA LOH of allele2 in this sample? `TRUE` if `dna_aib = TRUE` and `cn2_binned < 0.5`

RNAseq results
* `n_transcriptome_snps` -  - the number of mismatches between the two alleles, measured using the exome sequence defined in `assets/mhc_references/transcriptome/mhc_cds.fasta`.
* `rna_library_size` - the library size used in the RNAseq analysis.
* `allele1_rna_n_reads_before_filtering` - the number of reads in the RNAseq HLA allele1 BAM file before the BAM file is filtered to remove reads with more than the input parameter `max_mismatch`.
* `allele1_rna_n_reads_after_filtering` - the number of reads in the RNAseq HLA allele1 BAM file after the BAM file is filtered to remove reads with more than the input parameter `max_mismatch`.
* `allele2_rna_n_reads_before_filtering` - the number of reads in the RNAseq HLA allele2 BAM file before the BAM file is filtered to remove reads with more than the input parameter `max_mismatch`.
* `allele2_rna_n_reads_after_filtering` - the number of reads in the RNAseq HLA allele2 BAM file after the BAM file is filtered to remove reads with more than the input parameter `max_mismatch`.
* `gene_rpkm` - the gene level RPKM.
* `allele1_rpkm` - the allele1 RPKM
* `allele2_rpkm` - the allele2 RPKM
* `allele1_matched_normal_rpkm` - the RPKM of allele1 in the matched normal sample.
* `allele2_matched_normal_rpkm` - the RPKM of allele2 in the matched normal sample.
* `rna_aib_paired_t_test` - the p-value of the log R t-test.
* `rna_aib_paired_wilcoxon_test` - the p-value of the log R Wilcoxon.
* `rna_aib` - is there RNA AIB in this sample/gene? `TRUE` if `rna_aib_paired_wilcoxon_test < 0.01`.
* `frac_mapping_multi_gene` - the fraction of reads mapping to this gene that also map to other HLA genes.
* `frac_mapping_uniquely` - the fraction of reads that map to both alleles of the given gene.
* `allele1_repression_paired_t_test` - the p-value of the repression t-test for allele 1.
* `allele1_repression_paired_wilcoxon_test` - the p-value of the repression Wilcoxon.
* `allele1_repression_median_tumour_dp` - the median RNA depth across the SNPs in allele 1 in the tumour. 
* `allele1_repression_median_normal_dp` - the median RNA depth across the SNPs in allele 1 in the normal. 
* `allele2_repression_paired_t_test` - the p-value of the repression t-test for allele 2. 
* `allele2_repression_paired_wilcoxon_test` - the p-value of the repression Wilcoxon for allele 2. 
* `allele2_repression_median_tumour_dp` - the median RNA depth across the SNPs in allele 2 in the tumour. 
* `allele2_repression_median_normal_dp` - the median RNA depth across the SNPs in allele 2 in the normal. 
* `allele1_repressed` - is allele 1 repressed in the tumour compared to the normal? `TRUE` if `allele1_repression_paired_wilcoxon_test < 0.01` and `allele1_repression_median_tumour_dp < allele1_repression_median_normal_dp`.
* `allele1_over_expressed` - is allele 1 over-expressed in the tumour compared to the normal? `TRUE` if `allele1_repression_paired_wilcoxon_test < 0.01` and `allele1_repression_median_tumour_dp > allele1_repression_median_normal_dp`.
* `allele2_repressed` - is allele 2 repressed in the tumour compared to the normal? `TRUE` if `allele2_repression_paired_wilcoxon_test < 0.01` and `allele2_repression_median_tumour_dp < allele2_repression_median_normal_dp`.
* `allele2_over_expressed` - is allele 2 over-expressed in the tumour compared to the normal? `TRUE` if `allele2_repression_paired_wilcoxon_test < 0.01` and `allele2_repression_median_tumour_dp > allele2_repression_median_normal_dp`.

QC filters
* `fail_homozygous` - TRUE if the gene is homozygous.
* `fail_n_wxs_snps` - TRUE if the number of filtered DNA SNPs is less than `--min_number_of_snps`
* `fail_expected_depth` - TRUE if the number of filtered SNPs is less than `--min_expected_depth`
* `fail_cn_range` - TRUE if the copy number range is more than `--max_copy_number_range`
* `fail_rna_n_snps` - TRUE if the number of filtered RNA SNPs is less than `--min_number_of_snps`
* `fail_multi_gene_mapping` - TRUE if `frac_mapping_multi_gene` is less than `--max_frac_mapping_multi_gene`
* `fail_unique_allele_mapping` - TRUE if `fail_unique_allele_mapping` is more than `--min_frac_mapping_uniquely`

***

**`novel_splicing_events.csv`**

This contains the novel splicing events, i.e. splicing events that are not defined in the input `mhc.gtf` file. Each row represents a novel splice junction detected in a sample. Columns include:

* `patient` - the patient name, defined in the input inventory file.
* `sample_name` - the sample name, defined in the input inventory file.
* `matched_normal_sample_name` - the matched RNAseq normal sample name, if avaliable, as defined in the input inventory file.
* `gene` - either `HLA-A`, `HLA-B` or `HLA-C`.
* `allele` - the full allele name. 
* `short_allele` - the allele name shortened to two digits.
* `start` - the start position of the splice junction.
* `end` - the end position of the splice junction.
* `n_unique_reads` - The number of uniquely mapping reads supporting the splice junction in the tumour sample. This is the 7th column in splice junction table that is output from STAR (the SJ.out.tab in the STAR documentation).
* `normal_n_unique_reads` - The number of uniquely mapping reads supporting the splice junction in the normal sample. This is the 7th column in splice junction table that is output from STAR.
* `canonical_sj_read_count` - The number of uniquely mapping reads supporting the canonical splice junction in the tumour sample.
* `normal_canonical_sj_read_count` - The number of uniquely mapping reads supporting the canonical splice junction in the normal sample.
* `total_read_count` - The total read count in the tumour sample, i.e. `n_unique_reads` + `canonical_sj_read_count`.
* `normal_total_read_count` - The total read count in the normal sample, i.e. `normal_n_unique_reads + normal_canonical_sj_read_count`.
* `novel_transcript_proportion` - The novel transcript proportion in the tumour sample, i.e. `n_unique_reads/total_read_count`.
* `normal_novel_transcript_proportion` - The novel transcript proportion in the normal sample, i.e. `normal_n_unique_reads/normal_total_read_count`.
* `novel_transcript_proportion_change` - `novel_transcript_proportion - normal_novel_transcript_proportion`
* `novel_transcript_proportion_changev2` - if `novel_transcript_proportion > normal_novel_transcript_proportion`, this is defined as `(tumour_novel_transcript_proportion - normal_novel_transcript_proportion)/tumour_novel_transcript_proportion`. If `novel_transcript_proportion <= normal_novel_transcript_proportion` this is defined as `-1*(normal_novel_transcript_proportion - novel_transcript_proportion)/normal_novel_transcript_proportion`.
* `fisher_pvalue` - the p-value from the Fisher's test checking if the tumour or normal is enriched for the given splice junction.
* `fisher_odds_ratio` - the odds ratio from the Fisher's test checking if the tumour or normal is enriched for the given splice junction.
* `fisher_ci_lower` - the lower value from the confidence inverval from the Fisher's test checking if the tumour or normal is enriched for the given splice junction.
* `fisher_ci_upper` - the upper value from the confidence inverval from the Fisher's test checking if the tumour or normal is enriched for the given splice junction.
* `sj_type` - either 'partial_intron_retention', 'partial_exon_skip', 'complete_exon_skip' or 'other'.
* `exon_intron_name` - the name of the exon/intron that is skipped/retained.
* `premature_stop` - TRUE if the splice junction introduces a premature stop codon into the new sequence.
* `framshift` - TRUE if the splice junction introduces a frameshift into the new sequence.
* `sj_consequence` - either 'Inframe - PTC', 'Frameshift - PTC', 'Frameshift - no PTC' or 'Inframe'.
* `sample_type` - either 'tumour' or 'normal', as defined in the input inventory file.
* `purity` - the purity of the tumour sample, as defined in the input inventory file.
* `purity_scaled_novel_transcript_proportion` - the novel_transcript_proportion scaled for the purity.

Both tumour and normal samples are included as rows. If the `sample_type = normal`, or if `sample_type = tumour` but `matched_normal_sample_name` is empty, then the following columns will be empty: `normal_n_unique_reads`, `normal_canonical_sj_read_count`, `normal_total_read_count`, `normal_novel_transcript_proportion`, `novel_transcript_proportion_change`, `novel_transcript_proportion_changev2`, `fisher_pvalue`, `fisher_odds_ratio`, `fisher_ci_lower`, `fisher_ci_upper`, `purity`, `purity_scaled_novel_transcript_proportion`.

***

**`known_splicing_events.csv`**

This contains the known splicing events, i.e. splicing events that are defined in the input `mhc.gtf` file. Each row represents a known splice junction. If a known splice junction isn't detected in a sample, then it is still included in this table with `n_unique_reads` set to zero.  Columns include:

* `patient` - the patient name, defined in the input inventory file.
* `sample_name` - the sample name, defined in the input inventory file.
* `gene` - either `HLA-A`, `HLA-B` or `HLA-C`.
* `allele` - the full allele name. 
* `short_allele` - the allele name shortened to two digits.
* `start` - the start position of the splice junction.
* `end` - the end position of the splice junction.
* `n_unique_reads` - The number of uniquely mapping reads supporting the splice junction in the tumour sample. This is the 7th column in splice junction table that is output from STAR (the SJ.out.tab in the STAR documentation).
* `known_feature_name` - the name of the intron.
* `sample_type` - either 'tumour' or 'normal', as defined in the input inventory file.
* `purity` - the purity of the tumour sample, as defined in the input inventory file.

***

**`cohort_library_size.csv`**
This contains the library size for every sample in the cohort.

***

**`mosdepth_*.csv`**
This contains the mosdepth output for every allele specific HLA BAM file, including the WES and RNAseq Novoalign BAM files and the RNAseq STAR BAM files. The following columns are included:
* `sample_name` - the sample name, defined in the input inventory file.
* `allele` - the full allele name. 
* `start` - the start position of the feature (intron, exon, UTR)
* `end` - the end position of the feature (intron, exon, UTR)
* `feature_name` - the name of the feature
* `depth` - the mean depth across the given feature.

