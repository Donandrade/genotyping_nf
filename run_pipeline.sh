#!/bin/bash
#SBATCH --job-name=nf_blueberry
#SBATCH --output=logs/nf_manager_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=$MY_EMAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=96:00:00
#SBATCH --account=munoz
#SBATCH --qos=munoz-b

# Run this sbatch script as shown below to pass the email via parameter.
## sbatch --mail-user=your-email@ufl.edu run_pipeline.sh

# 1. Preparation
mkdir -p logs
module purge
module load nextflow

# 2. HiPerGator Optimization (Using local /tmp disk)
export NXF_TEMP=$SLURM_TMPDIR
export NXF_WORK="./work"
export NXF_OPTS="-Xms2g -Xmx4g"

echo "========================================================"
echo "Starting Nextflow Pipeline"
echo "Date: $(date)"
echo "Work Dir: $PWD"
echo "========================================================"

# 3. Simplified Execution
# Nextflow will read 'params.ref' and 'params.samples' from nextflow.config

# Use "--probes false" to run the pileup for the entire chromosome
# Use "--past_calls" to set the path for the previous variant calls (chunks)
#    --past_calls 2025_calling/output/04_final_calls/chunks/ \
# nextflow run main.nf \
    # --samples samples/sample_2_last20.tsv \
    # --probes probes.bed \
    # --past_calls 2026_calling/output/04_final_calls/chunks/ \
    # --chunk_size 10000 \
    # -resume \
    # -with-report logs/report.html \
    # -with-timeline logs/timeline.html


nextflow run main.nf \
	--samples $MY_SAMPLES \
	--probes $MY_PROBES \
	--past_calls $MY_OLD_PILEUPS \
	--chunk_size $MY_CHUNK \
	-resume \
        -with-report logs/report.html \
        -with-timeline logs/timeline.html

echo "Finished: $(date)"
