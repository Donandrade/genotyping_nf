# Nextflow Genotyping Pipeline (v2.0)

This pipeline performs a complete genotyping workflow: from raw read quality control to joint variant calling. It is specifically optimized for large-scale execution on the **University of Florida HiPerGator (SLURM)**.

Unlike the [previous Bash-based version](git@github.com:Donandrade/genotyping_pipeline.git), this **Nextflow** implementation automates sample-level and genomic parallelization, ensuring maximum efficiency and reproducibility for bioinformatics research.

---

## 1. Key Improvements (Nextflow vs. Bash)

* **Native Parallelism:** Automatically submits multiple jobs to SLURM; no need for manual job arrays.
* **Genomic Chunking:** Breaks large chromosomes into smaller, manageable chunks to speed up `bcftools merge` and `call`.
* **Fault Tolerance:** Native support for `-resume` (starts from where it left off) and automatic retries for transient HPC errors.
* **Incremental Merging:** Seamlessly integrates new sequence data with historical variant calls using the `--past_calls` flag.

---

## 2. Pipeline Architecture

The workflow consists of the following steps:
1.  **TRIM:** Quality trimming and adapter removal using **Trimmomatic**.
2.  **ALIGN:** **BWA-MEM** alignment, followed by **Samtools** sorting and indexing.
3.  **QC:** BAM statistics and coverage analysis via **Samtools** and **Mosdepth**.
4.  **INDIVIDUAL PILEUP:** Generation of raw per-sample VCFs.
5.  **SPLIT & CHUNK:** Parallel splitting of VCFs into defined genomic regions (chunks).
6.  **MERGE & CALL:** Joint calling of variants across all samples (including historical data if provided).
7.  **CONCATENATE:** Merging chunks back into a single, **genome-wide final VCF**.

---

## 3. Getting Started

## Getting started

Clone this repository into your SLURM-based HPC account before running the pipeline.

```bash
git clone https://github.com/Donandrade/genotyping_nf.git

cd genotyping_nf
```

---

## 4. Repository Organization

The repository is divided into scripts and accessory files (including example FASTQ and sample.tsv files):

```bash
genotyping_nf/
├── examples/               # Example datasets and metadata for testing
│   ├── fastq_set1          # 1,920 simulated paired FASTQ files for workflow testing
│   ├── fastq_set2          # Another 1,920 simulated FASTQ files for workflow testing
│   ├── large_sample_1.tsv  # Sample sheet containing paths for fastq_set1 (1,920 files)
│   ├── large_sample_2.tsv  # Sample sheet containing paths for fastq_set2 (1,920 files)
│   ├── small_sample_1.tsv  # Sample sheet with a reduced set of FASTQ files
│   ├── small_sample_2.tsv  # Sample sheet with a reduced set of FASTQ files
│   ├── small_sample_3.tsv  # Sample sheet with a reduced set of FASTQ files
│   ├── small_sample_4.tsv  # Sample sheet with a reduced set of FASTQ files
│   ├── reference/          # Directory containing a small reference genome and its BWA indexes
│   └── probes.bed          # Example target probes file for the provided datasets
├── main.nf                 # Main Nextflow script containing the variant calling pipeline logic
├── nextflow.config         # Pipeline configuration: defines parameters, I/O paths, and resource allocation (CPU/RAM)
└── run_pipeline.sh         # Shell script to execute the Nextflow pipeline
```

**PLEASE NOTE**

Nextflow provides a highly efficient way to allocate `CPU` and `RAM` resources for each process in the workflow. This means you can customize the specific amount of memory and CPU cores used by each step (e.g., `TRIM`, `ALIGN`, `QC`, `INDIVIDUAL PILEUP`, `SPLIT & CHUNK`, `MERGE` & `CALL` and `CONCATENATE`). These resources are managed and defined within the `nextflow.config` file.

### Prerequisites
* **Nextflow:** (Load via `module load nextflow` on HiPerGator).
* **SLURM Account:** Configured for your specific group (e.g., `munoz`).

### Required Input Files
1.  **`samples.tsv`**: A tab-separated file containing a mandatory header. You can find various examples of these sample sheets in the `examples/` directory (e.g., `large_sample_1.tsv` or `small_sample_1.tsv`):

The following example shows the structure and header of the sample sheet located at `examples/large_sample_1.tsv`

```text
    sample	r1	r2
    sample0001	fastq_set1/sample0001_R1.fq.gz	fastq_set1/sample0001_R2.fq.gz
    sample0002	fastq_set1/sample0002_R1.fq.gz	fastq_set1/sample0002_R2.fq.gz
    sample0003	fastq_set1/sample0003_R1.fq.gz	fastq_set1/sample0003_R2.fq.gz
    sample0004	fastq_set1/sample0004_R1.fq.gz	fastq_set1/sample0004_R2.fq.gz
    sample0005	fastq_set1/sample0005_R1.fq.gz	fastq_set1/sample0005_R2.fq.gz
    sample0006	fastq_set1/sample0006_R1.fq.gz	fastq_set1/sample0006_R2.fq.gz
    sample0007	fastq_set1/sample0007_R1.fq.gz	fastq_set1/sample0007_R2.fq.gz
    sample0008	fastq_set1/sample0008_R1.fq.gz	fastq_set1/sample0008_R2.fq.gz
    sample0009	fastq_set1/sample0009_R1.fq.gz	fastq_set1/sample0009_R2.fq.gz
 ```

2.  **Reference Genome:** A FASTA file indexed with `samtools faidx` and `bwa index` (e.g., in `examples/reference/`).
3.  **Probes (Optional):** A `.bed` file to restrict analysis to specific regions of interest (e.g., `examples/probes.bed`).

---

## 5. Usage

The pipeline is designed to be submitted as a background "manager" job to prevent local timeouts.

### Configuration
Edit `nextflow.config` to adjust:
* SLURM account/QOS.
* Memory and CPU limits per process.
* Tool module versions.

### Execution

This is an example of an execution where probes are targeted during the pileup generation step. Note that existing pileups (see the `MY_CHUNK` parameter below) are merged with the one generated in real-time.

Use the command below to submit the workflow if you have probes and previous pileups to merge with the current run:

```bash
sbatch --mail-user=deandradesilvae@ufl.edu \
    --export=ALL,MY_SAMPLES="examples/small_sample_1.tsv",MY_PROBES="examples/probes.bed",MY_OLD_PILEUPS="false",MY_REF="examples/reference/subgenome_blue.multi.fa",MY_CHUNK=10000 \
    run_pipeline.sh
```

This workflow should work in other scenarios, as in the example below:

1. You want to run the workflow just once, without probes and without previous pileups to merge

```bash
sbatch --mail-user=deandradesilvae@ufl.edu \
    --export=ALL,MY_SAMPLES="examples/small_sample_1.tsv",MY_PROBES="null",MY_OLD_PILEUPS="false",MY_REF="examples/reference/subgenome_blue.multi.fa",MY_CHUNK=10000 \
    run_pipeline.sh
```

If you wanna test the workflow this repository has two simulated FASTQ dataset in `fastq_set1/` and `fastq_set2/` direcoties and the `sample.tsv` exemples in the samples directory

**Variable Descriptions:**

`--mail-user`: Your institutional email address for SLURM job notifications. Providing this allows the scheduler to send you real-time updates regarding job status (e.g., BEGIN, END, or FAIL).

`MY_SAMPLES`: Path to your sample TSV file (see exemples in the `samples/` directory).

`MY_PROBES`: Path to your BED file (or null) - see exemples in the `probes.bed` file.

`MY_OLD_PILEUPS`: Path to the directory containing VCF chunks from a previous run.

`MY_REF`: Path to the index prefix generated by `bwa index`.

**Important:** The chunk size of these historical files must be identical to the value set in the MY_CHUNK parameter to ensure correct merging.

`MY_CHUNK`: Defines the genomic window size (in base pairs) for parallel processing. This value determines how the merged pileups are split into multiple VCF chunks for simultaneous variant calling.
Optimization Note: Avoid using excessively small values (e.g., < 5,000 bp), as creating too many small tasks can overload the cluster's file system and scheduler due to high thread and I/O overhead.

## 6. Directory Structure & Outputs

The pipeline organizes results into the following structure:

```bash
output/
├── 01_trimmed/           # Cleaned FASTQ files and trimming logs
├── 02_bam/               # Sorted and indexed BAM files
├── 04_splited_call/      # Intermediate VCFs split by genomic chunk
├── 04_final_calls/
│   ├── chunks/           # Per-chunk merged and called VCFs (pileup & called)
│   └── global/           # Final 'genome_wide_final.vcf.gz' output
├── pipeline_trace.txt    # Detailed performance metrics (CPU, RAM, Time)
└── execution_report.html # Visual performance and success report
```

