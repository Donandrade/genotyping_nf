# Nextflow Genotyping Pipeline (v2.0)

This pipeline performs a complete genotyping workflow: from raw read quality control to joint variant calling. It is specifically optimized for large-scale execution on the **University of Florida HiPerGator (SLURM)**.

Unlike the previous Bash-based version, this **Nextflow** implementation automates sample-level and genomic parallelization, ensuring maximum efficiency and reproducibility for bioinformatics research.

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

### Prerequisites
* **Nextflow:** (Load via `module load nextflow` on HiPerGator).
* **SLURM Account:** Configured for your specific group (e.g., `munoz`).

### Required Input Files
1.  **`samples.tsv`**: A tab-separated file with a mandatory header:
    ```text
    sample  r1  r2
    ID01    fastq/ID01_R1.fq.gz  fastq/ID01_R2.fq.gz
    ID02    fastq/ID02_R1.fq.gz  fastq/ID02_R2.fq.gz
    ```
2.  **Reference Genome:** A FASTA file indexed with `samtools faidx` and `bwa index`.
3.  **Probes (Optional):** A `.bed` file to restrict analysis to specific regions of interest.

---

## 4. Usage

The pipeline is designed to be submitted as a background "manager" job to prevent local timeouts.

### Configuration
Edit `nextflow.config` to adjust:
* SLURM account/QOS.
* Memory and CPU limits per process.
* Tool module versions.

### Execution

Submit the pipeline using the provided wrapper script:

```bash
sbatch --mail-user=your@email.edu \
       --export=ALL,MY_SAMPLES="samples.tsv",MY_PROBES="probes.bed",MY_OLD_PILEUPS="path/to/old_chunks/",MY_CHUNK=10000 \
       run_pipeline.sh
```

**Variable Descriptions:**

`--mail-user`: Your institutional email address for SLURM job notifications. Providing this allows the scheduler to send you real-time updates regarding job status (e.g., BEGIN, END, or FAIL).

`MY_SAMPLES`: Path to your sample TSV file.

`MY_PROBES`: Path to your BED file (or null).

`MY_OLD_PILEUPS`: Path to the directory containing chunks from a previous run.

`MY_CHUNK`: Defines the genomic window size (in base pairs) for parallel processing. This value determines how the merged pileups are split into multiple VCF chunks for simultaneous variant calling.
Optimization Note: Avoid using excessively small values (e.g., < 5,000 bp), as creating too many small tasks can overload the cluster's file system and scheduler due to high thread and I/O overhead.

## 5. Directory Structure & Outputs

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
