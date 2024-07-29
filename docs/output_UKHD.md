# nf-core/rnaseq: Output on the UKHD system

## Introduction

This document outlines the output from the basic RNA-seq pipeline on the UKHD system using default tools and parameters. For the full version with different tools, please refer to the [complete documentation](https://github.com/lucianhu/rnaseq_nextflow/blob/master/docs/output.md).

> **Warning**
> Running Nextflow on Docker is insecure and requires administrative privileges. Using Conda is more secure and allows users to run Nextflow without needing special permissions.

Most of the plots are taken from the MultiQC report generated from the [full-sized test dataset](https://github.com/nf-core/test-datasets/tree/rnaseq#full-test-dataset-origin) for the pipeline using a command similar to the one below:

```bash
$ nextflow run /path/to/rnaseq_nextflow/main.nf -profile conda -params-file nf-params-rna.yaml -work-dir path/to/work_dir
```

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [FastQC](#fastqc) - Raw read QC
  - [TrimGalore](#trimgalore) - Adapter and quality trimming
- [Alignment and quantification](#alignment-and-quantification)
  - [STAR](#star) - Fast spliced aware genome alignment
  - [Salmon](#salmon) - Transcriptome quantification
- [Alignment post-processing](#alignment-post-processing)
  - [SAMtools](#samtools) - Sort and index alignments
  - [picard MarkDuplicates](#picard-markduplicates) - Duplicate read marking
- [Other steps](#other-steps)
  - [StringTie](#stringtie) - Transcript assembly and quantification
  - [featureCounts](#featurecounts) - Read counting relative to gene biotype
  - [BEDTools and bedGraphToBigWig](#bedtools-and-bedgraphtobigwig) - Create bigWig coverage files
- [Quality control](#quality-control)
  - [RSeQC](#rseqc) - Various RNA-seq QC metrics
  - [Qualimap](#qualimap) - Various RNA-seq QC metrics
  - [dupRadar](#dupradar) - Assessment of technical / biological read duplication
  - [DESeq2](#deseq2) - PCA plot and sample pairwise distance heatmap and dendrogram
  - [MultiQC](#multiqc) - Present QC for raw reads, alignment, read counting and sample similiarity
- [Workflow reporting and genomes](#workflow-reporting-and-genomes)
  - [Reference genome files](#reference-genome-files) - Saving reference genome indices/files
  - [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Preprocessing

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots in this directory are generated relative to the raw, input reads. They may contain adapter sequence and regions of low quality. To see how your reads look after adapter and quality trimming please refer to the FastQC reports in the `trimgalore/fastqc/` directory.

```bash
$ fastqc --quiet --threads 6 ${SAMPLE}.read _1.fastq.gz ${SAMPLE}.read_2.fastq.gz
```

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `${SAMPLE}_fastqc.html`: FastQC report containing quality metrics.
  - `${SAMPLE}_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

### TrimGalore

[Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) is a wrapper tool around Cutadapt and FastQC to peform quality and adapter trimming on FastQ files. Cutadapt identifies and removes adapter sequences, primers, poly-A tails, and other unwanted sequences from high-throughput sequencing reads. FastQC, on the other hand, assesses the quality of the sequencing data. It is the default trimming tool used by this pipeline, however you can use fastp instead by specifying the `--trimmer fastp` parameter. You can specify additional options for Trim Galore! via the `--extra_trimgalore_args` parameters.

> **NB:** TrimGalore! will only run using multiple cores if you are able to use more than > 5 and > 6 CPUs for single- and paired-end data, respectively. The total cores available to TrimGalore! will also be capped at 4 (7 and 8 CPUs in total for single- and paired-end data, respectively) because there is no longer a run-time benefit. See [release notes](https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019) and [discussion whilst adding this logic to the nf-core/atacseq pipeline](https://github.com/nf-core/atacseq/pull/65).

```bash
$ trim_galore \
    --fastqc_args '-t 12' \
    --cores 8 \
    --paired \
    --gzip \
    ${SAMPLE}.read_1.fastq.gz \
    ${SAMPLE}.read_2.fastq.gz
```

<details markdown="1">
<summary>Output files</summary>

- `trimgalore/`
  - `${SAMPLE}_val_1.fq.gz` and `${SAMPLE}_val_2.fq.gz` : If `--save_trimmed` is specified, FastQ files **after** adapter trimming will be placed in this directory.
  - `${SAMPLE}_trimming_report.txt`: Log file generated by Trim Galore!.
- `trimgalore/fastqc/`
  - `${SAMPLE}_fastqc.html`: FastQC report containing quality metrics for read 1 (_and read2 if paired-end_) **after** adapter trimming.
  - `${SAMPLE}_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
 
</details>

![MultiQC - cutadapt trimmed sequence length plot](images/mqc_cutadapt_trimmed.png)

## Alignment and quantification

### STAR

[STAR](https://github.com/alexdobin/STAR) is a read aligner designed for splice aware mapping typical of RNA sequencing data. STAR stands for *S*pliced *T*ranscripts *A*lignment to a *R*eference, and has been shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. Using `--aligner star_salmon` is the default alignment and quantification option.

```bash
$ STAR \
    --genomeDir star \
    --readFilesIn input1/${SAMPLE}_val_1.fq.gz input2/${SAMPLE}_val_2.fq.gz \
    --runThreadN 12 \
    --outFileNamePrefix $NAME \
     \
    --sjdbGTFfile genome.gtf \
    --outSAMattrRGline 'ID:NAME' 'SM:NAME'  \
    --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --readFilesCommand zcat --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outSAMstrandField intronMotif
```

<details markdown="1">
<summary>Output files</summary>

- `star_salmon/`
  - `${SAMPLE}.Aligned.out.bam`: If `--save_align_intermeds` is specified the original BAM file containing read alignments to the reference genome will be placed in this directory.
  - `${SAMPLE}.Aligned.toTranscriptome.out.bam`: If `--save_align_intermeds` is specified the original BAM file containing read alignments to the transcriptome will be placed in this directory.
- `star_salmon/log/`
  - `${SAMPLE}.SJ.out.tab`: File containing filtered splice junctions detected after mapping the reads.
  - `${SAMPLE}.Log.final.out`: STAR alignment report containing the mapping results summary.
  - `${SAMPLE}.Log.out` and `${SAMPLE}.Log.progress.out`: STAR log files containing detailed information about the run. Typically only useful for debugging purposes.
- `star_salmon/unmapped/`
  - `${SAMPLE}.fastq.gz`: If `--save_unaligned` is specified, FastQ files containing unmapped reads will be placed in this directory.

The STAR section of the MultiQC report shows a bar plot with alignment rates: good samples should have most reads as _Uniquely mapped_ and few _Unmapped_ reads.

</details>

![MultiQC - STAR alignment scores plot](images/mqc_star.png)

### Salmon

[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) from [Ocean Genomics](https://oceangenomics.com/) allows quantification of reads against an index generated from a reference set of target transcripts. By default, the transcriptome-level BAM files generated by STAR are provided to Salmon for downstream quantification.

```bash
# Prepare the reference files for RSEM (RNA-Seq by Expectation-Maximization)

$ rsem-prepare-reference \
    --gtf genome.gtf \
    --num-threads 12 \  
    rsem/genome.fa \
    rsem/genome

$ cp rsem/genome.transcripts.fa .
            
# Extracts decoy sequences from the genome fasta file
$ grep '^>' genome.fa | cut -d ' ' -f 1 | cut -d $'\t' -f 1 > decoys.txt

# Removes the '>' character from the decoy sequence identifiers
$ sed -i.bak -e 's/>//g' decoys.txt

# Combines the genome and transcriptome fasta files into a gentrome file
$ cat genome.transcripts.fa genome.fa > gentrome.fa

# Indexes the combined gentrome file for Salmon, using the decoy sequences
$ salmon index --threads 6 -t gentrome.fa -d decoys.txt -k 31 -i salmon

# Quantifies transcript abundances using Salmon with alignment-based input and library type ISR
$ salmon quant --geneMap genome.gtf --threads 6 --libType=ISR \
    -t genome.transcripts.fa -a ${SAMPLE}.Aligned.toTranscriptome.out.bam -o ${SAMPLE}

# If the meta_info.json file exists, copy it to the current directory with a modified name
$ if [ -f ${SAMPLE}/aux_info/meta_info.json ]; then
    cp ${SAMPLE}/aux_info/meta_info.json "${SAMPLE}_meta_info.json"
fi

# Generates a gene-to-transcript mapping file
$ tx2gene.py --quant_type salmon --gtf genome.gtf --quants quants \
    --id gene_id --extra gene_name -o tx2gene.tsv
```

The `tximport`  package summarizes Salmon results into matrices for differential analysis, providing count and TPM quantifications at both gene and transcript levels.

```bash    
# Imports Salmon quantification data into R using tximport
$ tximport.r \
    NULL \                      # Placeholder for additional options (no extra options specified here)
    quants \                    # Directory or file path containing Salmon quantification results (e.g., TPM, counts)
    salmon.merged \             # Output directory or prefix where the summarized results will be saved
    salmon \                    # Indicates that Salmon was used for quantification
    tx2gene.tsv                 # File mapping transcript IDs to gene IDs for aggregation of transcript-level data to gene-level
```

The `summarizedexperiment.r` command imports gene expression data into `R` and creates SummarizedExperiment objects (`*.rds`). It uses gene-level count data and TPM values from Salmon, along with a mapping file (`tx2gene.tsv`) that connects transcript IDs to gene IDs. This process helps organize and prepare the data for subsequent analysis, such as differential expression analysis or other downstream tasks.

```bash
# Creates SummarizedExperiment objects from gene counts and TPM data
$ summarizedexperiment.r NULL salmon.merged.gene_counts.tsv salmon.merged.gene_tpm.tsv tx2gene.tsv
$ summarizedexperiment.r NULL salmon.merged.gene_counts_length_scaled.tsv salmon.merged.gene_tpm.tsv tx2gene.tsv
$ summarizedexperiment.r NULL salmon.merged.gene_counts_scaled.tsv salmon.merged.gene_tpm.tsv tx2gene.tsv

# Creates SummarizedExperiment objects from transcript counts and TPM data
$ summarizedexperiment.r NULL salmon.merged.transcript_counts.tsv salmon.merged.transcript_tpm.tsv tx2gene.tsv
```

<details markdown="1">
<summary>Output files</summary>

- `star_salmon/`
  - `star_salmon.merged.gene_counts.tsv`: Matrix of gene-level raw counts across all samples.
  - `star_salmon.gene_tpm.tsv`: Matrix of gene-level TPM values across all samples.
  - `star_salmon.gene_counts.rds`: RDS object that can be loaded in R that contains a [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) container with the TPM (`abundance`), estimated counts (`counts`) and transcript length (`length`) in the assays slot for genes.
  - `star_salmon.merged.gene_counts_scaled.tsv`: Matrix of gene-level library size-scaled counts across all samples.
  - `star_salmon.merged.gene_counts_scaled.rds`: RDS object that can be loaded in R that contains a [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) container with the TPM (`abundance`), estimated library size-scaled counts (`counts`) and transcript length (`length`) in the assays slot for genes.
  - `star_salmon.merged.gene_counts_length_scaled.tsv`: Matrix of gene-level length-scaled counts across all samples.
  - `star_salmon.merged.gene_counts_length_scaled.rds`: RDS object that can be loaded in R that contains a [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) container with the TPM (`abundance`), estimated length-scaled counts (`counts`) and transcript length (`length`) in the assays slot for genes.
  - `star_salmon.merged.transcript_counts.tsv`: Matrix of isoform-level raw counts across all samples.
  - `star_salmon.merged.transcript_tpm.tsv`: Matrix of isoform-level TPM values across all samples.
  - `star_salmon.merged.transcript_counts.rds`: RDS object that can be loaded in R that contains a [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) container with the TPM (`abundance`), estimated isoform-level raw counts (`counts`) and transcript length (`length`) in the assays slot for transcripts.
  - `tx2gene.tsv`: Tab-delimited file containing gene to transcripts IDs mappings.
  </details>

Additional files for Salmon:

<details markdown="1">
<summary>Output files</summary>

- `salmon/${SAMPLE}/`
  - `aux_info/`: Auxiliary info e.g. versions and number of mapped reads.
  - `cmd_info.json`: Information about the Salmon quantification command, version and options.
  - `lib_format_counts.json`: Number of fragments assigned, unassigned and incompatible.
  - `libParams/`: Contains the file `flenDist.txt` for the fragment length distribution.
  - `logs/`: Contains the file `salmon_quant.log` giving a record of Salmon's quantification.
  - `quant.genes.sf`: Salmon _gene_-level quantification of the sample, including feature length, effective length, TPM, and number of reads.
  - `quant.sf`: Salmon _transcript_-level quantification of the sample, including feature length, effective length, TPM, and number of reads.

</details>

![MultiQC - Salmon fragment length distribution plot](images/mqc_salmon.png)

## Alignment post-processing

The pipeline has been written in a way where all the files generated downstream of the alignment are placed in the same directory as specified by `--aligner` e.g. if `--aligner star_salmon` is specified then all the downstream results will be placed in the `star_salmon/` directory.

### SAMtools

The original BAM files generated by the selected alignment algorithm are further processed with [SAMtools](http://samtools.sourceforge.net/) to sort them by coordinate, for indexing, as well as to generate read mapping statistics.

```bash
# Sort BAM File
$ samtools sort -@ 6 -o ${SAMPLE}.sorted.bam -T ${SAMPLE}.sorted ${SAMPLE}.Aligned.out.bam

# Index Sorted BAM File
$ samtools index -@ 1 ${SAMPLE}.sorted.bam

# Generate Statistics
$ samtools stats --threads 1 --reference genome.fa ${SAMPLE}.sorted.bam > ${SAMPLE}.sorted.bam.stats

# Flagstat Summary
$ samtools flagstat --threads 1 ${SAMPLE}.sorted.bam > ${SAMPLE}.sorted.bam.flagstat

# Index Statistics
$ samtools idxstats --threads 0 ${SAMPLE}.sorted.bam > ${SAMPLE}.sorted.bam.idxstats
```

<details markdown="1">
<summary>Output files</summary>

- `star_salmon/`
  - `${SAMPLE}.sorted.bam`: If `--save_align_intermeds` is specified the original coordinate sorted BAM file, which are necessary for many downstream analyses, such as variant calling, visualization in genome browsers, and various statistical analyses.
  - `${SAMPLE}.sorted.bam.bai`: If `--save_align_intermeds` is specified the BAI index file, which enables efficient querying and retrieval of data from the BAM file, which is especially useful for visualization tools like IGV (Integrative Genomics Viewer) and for quickly accessing specific genomic regions.
  - `${SAMPLE}.sorted.bam.csi`: If `--save_align_intermeds --bam_csi_index` is specified the CSI index file, which provides additional indexing capabilities for large datasets, allowing for efficient access to specific regions in very large BAM files.
 
- `<ALIGNER>/samtools_stats/`
  - `${SAMPLE}.sorted.bam.stats`: Detailed statistics on read counts, coverage, and mapping quality, often in a text format.
  - `${SAMPLE}.sorted.bam.flagstat`: Text file with statistics on the alignment process, including counts of mapped vs. unmapped reads and read pairing.
  - `${SAMPLE}.sorted.bam.idxstats`: Text file with per-chromosome or per-reference sequence read counts and lengths.

</details>

![MultiQC - SAMtools alignment scores plot](images/mqc_samtools_mapped.png)

![MultiQC - SAMtools mapped reads per contig plot](images/mqc_samtools_idxstats.png)

### picard MarkDuplicates

If fragmentation occurs before amplification (a common scenario), [picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) is useful for identifying and marking PCR duplicates. However, if fragmentation occurs after amplification, using Unique Molecular Identifiers [UMIs](https://emea.illumina.com/science/sequencing-method-explorer/kits-and-arrays/umi.html) for deduplication is more appropriate. This step **marks** duplicate reads in the alignments to gauge the overall duplication level in your samples but does **NOT** remove them. RNA-Seq data often contain a significant level of true biological duplication, especially from highly expressed genes, where the same fragments are sequenced multiple times.

In RNA-Seq, varying gene expression levels often result in a higher duplication rate compared to DNA-Seq. For highly expressed genes, a high number of duplicate reads is expected and reasonable, whereas less expressed genes should have fewer duplicate reads.

```bash
$ picard \
    -Xmx29491M \
    MarkDuplicates \
    --ASSUME_SORTED true --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp \
    --INPUT ${SAMPLE}.sorted.bam \
    --OUTPUT ${SAMPLE}.markdup.sorted.bam \
    --REFERENCE_SEQUENCE genome.fa \
    --METRICS_FILE ${SAMPLE}.markdup.sorted.MarkDuplicates.metrics.txt
```

<details markdown="1">
<summary>Output files</summary>

- `star_salmon/`
  - `${SAMPLE}.markdup.sorted.bam`: Sorted BAM files are necessary for many downstream analyses, such as variant calling, visualization in genome browsers, and various statistical analyses.
  - `${SAMPLE}.markdup.sorted.bam.bai`: BAI index file for coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM index file and so will be saved by default in the results directory.
  - `${SAMPLE}.markdup.sorted.bam.csi`: CSI index file for coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM index file and so will be saved by default in the results directory. Only generated if `--bam_csi_index` is specified as a parameter.
- `star_salmon/samtools_stats/`
  - SAMtools `${SAMPLE}.markdup.sorted.bam.flagstat`, `${SAMPLE}.markdup.sorted.bam.idxstats` and `${SAMPLE}.markdup.sorted.bam.stats` files generated from the duplicate marked alignment files.
- `star_salmon/picard_metrics/`
  - `${SAMPLE}.markdup.sorted.MarkDuplicates.metrics.txt`: Metrics file from MarkDuplicates.

</details>

![MultiQC - Picard MarkDuplicates metrics plot](images/mqc_picard_markduplicates.png)

## Other steps

### StringTie

[StringTie](https://ccb.jhu.edu/software/stringtie/) is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus. In order to identify differentially expressed genes between experiments, StringTie's output can be processed by specialized software like [Ballgown](https://github.com/alyssafrazee/ballgown), [Cuffdiff](http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/index.html) or other programs ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), etc.).

```bash
$ stringtie \  # Invoke StringTie for transcript assembly and quantification
    ${SAMPLE}.markdup.sorted.bam \  # Input BAM file with aligned RNA-Seq reads (duplicates marked and sorted)
    --rf \  # Assumes a stranded library with fr-firststrand orientation (reverse-forward)
    -G genome.gtf \  # Reference annotation file in GTF format (mouse genome GRCm39)
    -o ${SAMPLE}.transcripts.gtf \  # Output file for assembled transcripts in GTF format
    -A ${SAMPLE}.gene.abundance.txt \  # Output file for gene-level abundance estimates (tab-delimited)
    -C ${SAMPLE}.coverage.gtf \  # Output file for read coverage information in GTF format
    -b ${SAMPLE}.ballgown \  # Output directory for Ballgown input files (*.ctab)
    -p 6 \  # Number of threads to use for computation (6 threads)
    -v \  # Enable verbose mode for detailed logging
    -e  # Estimate expression only for transcripts present in the provided annotation (-G); no novel transcript assembly
```

<details markdown="1">
<summary>Output files</summary>

- `star_salmon/stringtie/`
  - `${SAMPLE}.coverage.gtf`: GTF file containing transcripts that are fully covered by reads.
  - `${SAMPLE}.transcripts.gtf`: GTF file containing all of the assembled transcipts from StringTie.
  - `${SAMPLE}.gene_abundance.txt`: Text file containing gene aboundances and FPKM values.
- `star_salmon/stringtie/<SAMPLE>.ballgown/`: Ballgown output directory.

</details>

### featureCounts

[featureCounts](http://bioinf.wehi.edu.au/featureCounts/) from the [Subread](http://subread.sourceforge.net/) package is a quantification tool used to summarise the mapped read distribution over genomic features such as genes, exons, promotors, gene bodies, genomic bins and chromosomal locations. We can also use featureCounts to count overlaps with different classes of genomic features. This provides an additional QC to check which features are most abundant in the sample, and to highlight potential problems such as rRNA contamination.

```bash
# Quantify Read Counts
$ featureCounts \
    -B -C -g gene_biotype -t exon \  # Count reads mapping to exons, grouped by gene biotype
    -p \  # Process paired-end reads
    -T 6 \  # Use 6 threads for processing
    -a genome.gtf \  # Reference annotation file in GTF format
    -s 2 \  # Strandedness: reverse strand
    -o ${SAMPLE}.featureCounts.txt \  # Output file for count results
    ${SAMPLE}.markdup.sorted.bam  # Input BAM file with aligned RNA-Seq reads

# Extract and Format Counts
$ cut -f 1,7 ${SAMPLE}.featureCounts.txt | tail -n +3 | cat biotypes_header.txt - >> ${SAMPLE}.biotype_counts_mqc.tsv
    # Extract columns 1 and 7 (gene IDs and counts) from featureCounts output,
    # Skip the first two lines (header), 
    # Concatenate with a predefined header file,
    # Append the result to biotype counts file for MultiQC

# Generate RNA Summary
$ mqc_features_stat.py \
    ${SAMPLE}.biotype_counts_mqc.tsv \  # Input file with biotype counts
    -s ${SAMPLE} \  # Sample identifier for naming and labeling
    -f rRNA \  # Focus on rRNA biotype
    -o ${SAMPLE}.biotype_counts_rrna_mqc.tsv  # Output file for rRNA counts summary
```

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/featurecounts/`
  - `*.featureCounts.txt`: featureCounts biotype-level quantification results for each sample.
  - `*.featureCounts.txt.summary`: featureCounts summary file containing overall statistics about the counts.
  - `*_mqc.tsv`: MultiQC custom content files used to plot biotypes in report.

</details>

![MultiQC - featureCounts biotypes plot](images/mqc_featurecounts_biotype.png)

### BEDTools and bedGraphToBigWig

The [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) format is an indexed binary format useful for displaying dense, continuous data in Genome Browsers such as the [UCSC](https://genome.ucsc.edu/cgi-bin/hgTracks) and [IGV](http://software.broadinstitute.org/software/igv/). This mitigates the need to load the much larger BAM files for data visualisation purposes which will be slower and result in memory issues. The bigWig format is also supported by various bioinformatics software for downstream processing such as meta-profile plotting.

<details markdown="1">
<summary>Output files</summary>

- `star_salmon/bigwig/`
  - `*.forward.bigWig`: bigWig coverage file relative to genes on the forward DNA strand.
  - `*.reverse.bigWig`: bigWig coverage file relative to genes on the reverse DNA strand.

</details>

## Quality control

### RSeQC

[RSeQC](<(http://rseqc.sourceforge.net/)>) is a package of scripts designed to evaluate the quality of RNA-seq data. This pipeline runs several, but not all RSeQC scripts. You can tweak the supported scripts you would like to run by adjusting the `--rseqc_modules` parameter which by default will run all of the following: `bam_stat.py`, `inner_distance.py`, `infer_experiment.py`, `junction_annotation.py`, `junction_saturation.py`,`read_distribution.py` and `read_duplication.py`.

The majority of RSeQC scripts generate output files which can be plotted and summarised in the MultiQC report.

#### Infer experiment

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/infer_experiment/`
  - `*.infer_experiment.txt`: File containing fraction of reads mapping to given strandedness configurations.

</details>

This script predicts the "strandedness" of the protocol (i.e. unstranded, sense or antisense) that was used to prepare the sample for sequencing by assessing the orientation in which aligned reads overlay gene features in the reference genome. The strandedness of each sample has to be provided to the pipeline in the input samplesheet (see [usage docs](https://nf-co.re/rnaseq/usage#samplesheet-input)). However, this information is not always available, especially for public datasets. As a result, additional features have been incorporated into this pipeline to auto-detect whether you have provided the correct information in the samplesheet, and if this is not the case then a warning table will be placed at the top of the MultiQC report highlighting the offending samples (see image below). If required, this will allow you to correct the input samplesheet and rerun the pipeline with the accurate strand information. Note, it is important to get this information right because it can affect the final results.

RSeQC documentation: [infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py)

![MultiQC - Strand check table](images/mqc_strand_check.png)

![MultiQC - RSeQC infer experiment plot](images/mqc_rseqc_inferexperiment.png)

#### Read distribution

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/read_distribution/`
  - `*.read_distribution.txt`: File containing fraction of reads mapping to genome feature e.g. CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions etc.

</details>

This tool calculates how mapped reads are distributed over genomic features. A good result for a standard RNA-seq experiments is generally to have as many exonic reads as possible (`CDS_Exons`). A large amount of intronic reads could be indicative of DNA contamination in your sample but may be expected for a total RNA preparation.

RSeQC documentation: [read_distribution.py](http://rseqc.sourceforge.net/#read-distribution-py)

![MultiQC - RSeQC read distribution plot](images/mqc_rseqc_readdistribution.png)

#### Junction annotation

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/junction_annotation/bed/`
  - `*.junction.bed`: BED file containing splice junctions.
  - `*.junction.Interact.bed`: BED file containing interacting splice junctions.
- `<ALIGNER>/rseqc/junction_annotation/log/`
  - `*.junction_annotation.log`: Log file generated by the program.
- `<ALIGNER>/rseqc/junction_annotation/pdf/`
  - `*.splice_events.pdf`: PDF file containing splicing events plot.
  - `*.splice_junction.pdf`: PDF file containing splice junctions plot.
- `<ALIGNER>/rseqc/junction_annotation/rscript/`
  - `*.junction_plot.r`: R script used to generate pdf plots above.
- `<ALIGNER>/rseqc/junction_annotation/xls/`
  - `*.junction.xls`: Excel spreadsheet with junction information.

</details>

Junction annotation compares detected splice junctions to a reference gene model. Splicing annotation is performed in two levels: splice event level and splice junction level.

RSeQC documentation: [junction_annotation.py](http://rseqc.sourceforge.net/#junction-annotation-py)

![MultiQC - RSeQC junction annotation plot](images/mqc_rseqc_junctionannotation.png)

#### Inner distance

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/inner_distance/pdf/`
  - `*.inner_distance_plot.pdf`: PDF file containing inner distance plot.
- `<ALIGNER>/rseqc/inner_distance/rscript/`
  - `*.inner_distance_plot.r`: R script used to generate pdf plot above.
- `<ALIGNER>/rseqc/inner_distance/txt/`
  - `*.inner_distance_freq.txt`: File containing frequency of insert sizes.
  - `*.inner_distance_mean.txt`: File containing mean, median and standard deviation of insert sizes.

</details>

The inner distance script tries to calculate the inner distance between two paired-end reads. It is the distance between the end of read 1 to the start of read 2, and it is sometimes confused with the insert size (see [this blog post](http://thegenomefactory.blogspot.com.au/2013/08/paired-end-read-confusion-library.html) for disambiguation):

This plot will not be generated for single-end data. Very short inner distances are often seen in old or degraded samples (_eg._ FFPE) and values can be negative if the reads overlap consistently.

RSeQC documentation: [inner_distance.py](http://rseqc.sourceforge.net/#inner-distance-py)

![MultiQC - RSeQC inner distance plot](images/mqc_rseqc_innerdistance.png)

#### Junction saturation

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/junction_saturation/pdf/`
  - `*.junctionSaturation_plot.pdf`: PDF file containing junction saturation plot.
- `<ALIGNER>/rseqc/junction_saturation/rscript/`
  - `*.junctionSaturation_plot.r`: R script used to generate pdf plot above.

</details>

This script shows the number of splice sites detected within the data at various levels of subsampling. A sample that reaches a plateau before getting to 100% data indicates that all junctions in the library have been detected, and that further sequencing will not yield any more observations. A good sample should approach such a plateau of _Known junctions_, however, very deep sequencing is typically required to saturate all _Novel Junctions_ in a sample.

RSeQC documentation: [junction_saturation.py](http://rseqc.sourceforge.net/#junction-saturation-py)

![MultiQC - RSeQC junction saturation plot](images/mqc_rseqc_junctionsaturation.png)

#### Read duplication

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/read_duplication/pdf/`
  - `*.DupRate_plot.pdf`: PDF file containing read duplication plot.
- `<ALIGNER>/rseqc/read_duplication/rscript/`
  - `*.DupRate_plot.r`: R script used to generate pdf plot above.
- `<ALIGNER>/rseqc/read_duplication/xls/`
  - `*.pos.DupRate.xls`: Read duplication rate determined from mapping position of read. First column is “occurrence” or duplication times, second column is number of uniquely mapped reads.
  - `*.seq.DupRate.xls`: Read duplication rate determined from sequence of read. First column is “occurrence” or duplication times, second column is number of uniquely mapped reads.

</details>

This plot shows the number of reads (y-axis) with a given number of exact duplicates (x-axis). Most reads in an RNA-seq library should have a low number of exact duplicates. Samples which have many reads with many duplicates (a large area under the curve) may be suffering excessive technical duplication.

RSeQC documentation: [read_duplication.py](http://rseqc.sourceforge.net/#read-duplication-py)

![MultiQC - RSeQC read duplication plot](images/mqc_rseqc_readduplication.png)

#### BAM stat

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/bam_stat/`
  - `*.bam_stat.txt`: Mapping statistics for the BAM file.

</details>

This script gives numerous statistics about the aligned BAM files. A typical output looks as follows:

```txt
#Output (all numbers are read count)
#==================================================
Total records:                                 41465027
QC failed:                                     0
Optical/PCR duplicate:                         0
Non Primary Hits                               8720455
Unmapped reads:                                0

mapq < mapq_cut (non-unique):                  3127757
mapq >= mapq_cut (unique):                     29616815
Read-1:                                        14841738
Read-2:                                        14775077
Reads map to '+':                              14805391
Reads map to '-':                              14811424
Non-splice reads:                              25455360
Splice reads:                                  4161455
Reads mapped in proper pairs:                  21856264
Proper-paired reads map to different chrom:    7648
```

MultiQC plots each of these statistics in a dot plot. Each sample in the project is a dot - hover to see the sample highlighted across all fields.

RSeQC documentation: [bam_stat.py](http://rseqc.sourceforge.net/#bam-stat-py)

### Qualimap

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/qualimap/<SAMPLE>/`
  - `qualimapReport.html`: Qualimap HTML report that can be viewed in a web browser.
  - `rnaseq_qc_results.txt`: Textual results output.
- `<ALIGNER>/qualimap/<SAMPLE>/images_qualimapReport/`: Images required for the HTML report.
- `<ALIGNER>/qualimap/<SAMPLE>/raw_data_qualimapReport/`: Raw data required for the HTML report.
- `<ALIGNER>/qualimap/<SAMPLE>/css/`: CSS files required for the HTML report.

</details>

[Qualimap](http://qualimap.bioinfo.cipf.es/) is a platform-independent application written in Java and R that provides both a Graphical User Interface (GUI) and a command-line interface to facilitate the quality control of alignment sequencing data. Shortly, Qualimap:

- Examines sequencing alignment data according to the features of the mapped reads and their genomic properties.
- Provides an overall view of the data that helps to to the detect biases in the sequencing and/or mapping of the data and eases decision-making for further analysis.

The [Qualimap RNA-seq QC module](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#rna-seq-qc) is used within this pipeline to assess the overall mapping and coverage relative to gene features.

![MultiQC - Qualimap gene coverage plot](images/mqc_qualimap_coverage.png)

![MultiQC - Qualimap genomic origin plot](images/mqc_qualimap_features.png)

### dupRadar

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/dupradar/box_plot/`
  - `*_duprateExpBoxplot.pdf`: PDF file containing box plot for duplicate rate relative to mean expression.
- `<ALIGNER>/dupradar/gene_data/`
  - `*_dupMatrix.txt`: Text file containing duplicate metrics per gene.
- `<ALIGNER>/dupradar/histogram/`
  - `*_expressionHist.pdf`: PDF file containing histogram of reads per kilobase values per gene.
- `<ALIGNER>/dupradar/intercepts_slope/`
  - `*_intercept_slope.txt`: Text file containing intercept slope values.
- `<ALIGNER>/dupradar/scatter_plot/`
  - `*_duprateExpDens.pdf`: PDF file containing typical dupRadar 2D density scatter plot.

See [dupRadar docs](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html) for further information regarding the content of these files.

</details>

[dupRadar](https://www.bioconductor.org/packages/release/bioc/html/dupRadar.html) is a Bioconductor library written in the R programming language. It generates various QC metrics and plots that relate duplication rate with gene expression levels in order to identify experiments with high technical duplication. A good sample with little technical duplication will only show high numbers of duplicates for highly expressed genes. Samples with technical duplication will have high duplication for all genes, irrespective of transcription level.

![dupRadar - Example good and bad experiment plot](images/dupradar_example_plot.png)

> _Credit: [dupRadar documentation](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html)_

### Preseq

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/preseq/`
  - `*.lc_extrap.txt`: Preseq expected future yield file.
- `<ALIGNER>/preseq/log/`
  - `*.command.log`: Standard error output from command.

</details>

The [Preseq](http://smithlabresearch.org/software/preseq/) package is aimed at predicting and estimating the complexity of a genomic sequencing library, equivalent to predicting and estimating the number of redundant reads from a given sequencing depth and how many will be expected from additional sequencing using an initial sequencing experiment. The estimates can then be used to examine the utility of further sequencing, optimize the sequencing depth, or to screen multiple libraries to avoid low complexity samples. A shallow curve indicates that the library has reached complexity saturation and further sequencing would likely not add further unique reads. The dashed line shows a perfectly complex library where total reads = unique reads. Note that these are predictive numbers only, not absolute. The MultiQC plot can sometimes give extreme sequencing depth on the X axis - click and drag from the left side of the plot to zoom in on more realistic numbers.

![MultiQC - Preseq library complexity plot](images/mqc_preseq_plot.png)

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/<ALIGNER>/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools i.e. FastQC, Cutadapt, SortMeRNA, STAR, RSEM, HISAT2, Salmon, SAMtools, Picard, RSeQC, Qualimap, Preseq and featureCounts. Additionally, various custom content has been added to the report to assess the output of dupRadar, DESeq2 and featureCounts biotypes, and to highlight samples failing a mimimum mapping threshold or those that failed to match the strand-specificity provided in the input samplesheet. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Workflow reporting and genomes

### Reference genome files

<details markdown="1">
<summary>Output files</summary>

- `genome/`
  - `*.fa`, `*.gtf`, `*.gff`, `*.bed`, `.tsv`: If the `--save_reference` parameter is provided then all of the genome reference files will be placed in this directory.
- `genome/index/`
  - `star/`: Directory containing STAR indices.
  - `hisat2/`: Directory containing HISAT2 indices.
  - `rsem/`: Directory containing STAR and RSEM indices.
  - `salmon/`: Directory containing Salmon indices.
  - `kallisto/`: Directory containing Kallisto indices.

</details>

A number of genome-specific files are generated by the pipeline because they are required for the downstream processing of the results. If the `--save_reference` parameter is provided then these will be saved in the `genome/` directory. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices so that you can save them somewhere locally. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
