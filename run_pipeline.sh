#!/bin/bash
#SBATCH --job-name=nf_blueberry
#SBATCH --output=logs/nf_manager_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --time=96:00:00
#SBATCH --account=munoz
#SBATCH --qos=munoz

# 1. Preparation
mkdir -p logs
module purge
module load nextflow

# 2. HiPerGator Optimization
export NXF_TEMP=$SLURM_TMPDIR
export NXF_WORK="./work"
export NXF_OPTS="-Xms2g -Xmx4g"

echo "========================================================"
echo "Starting Nextflow Pipeline - Chromosome-based logic"
echo "Date: $(date)"
echo "Work Dir: $PWD"
echo "Samples: $MY_SAMPLES"
echo "Ref: $MY_REF"
echo "========================================================"

# 3. Execution
# Note: --chunk_size foi removido pois a lógica agora é por cromossomo.
# O parâmetro --past_calls agora deve apontar para o diretório /chromosomes/ de rodadas anteriores.

# Exemplo de submissão:
# sbatch --mail-user=seuemail@ufl.edu --export=ALL,MY_SAMPLES="samples.tsv",MY_PROBES="probes.bed",MY_OLD_PILEUPS="false",MY_REF="ref.fa" run_pipeline.sh

nextflow run main.nf \
    --samples "$MY_SAMPLES" \
    --ref "$MY_REF" \
    --probes "$MY_PROBES" \
    --past_calls "$MY_OLD_PILEUPS" \
    -resume \
    -with-report logs/report.html \
    -with-timeline logs/timeline.html

echo "Finished: $(date)"
