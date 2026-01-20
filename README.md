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
sbatch run_pipeline.sh

Key Command-Line Arguments:

--samples: Path to your sample TSV file.

--probes: Path to BED file (set to null to run genome-wide).

--past_calls: Path to the 04_final_calls/chunks/ directory from a previous run to enable incremental merging.

--chunk_size: Size in base-pairs for genomic splitting (default: 10,000).

```

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
