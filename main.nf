nextflow.enable.dsl=2

// Definição dos cromossomos (Baseado no seu genotyping.conf)
params.chromosomes = [
  "VaccDscaff1",
  "VaccDscaff2",
  "VaccDscaff4",
  "VaccDscaff6",
  "VaccDscaff7",
  "VaccDscaff11",
  "VaccDscaff12",
  "VaccDscaff13",
  "VaccDscaff17",
  "VaccDscaff20",
  "VaccDscaff21",
  "VaccDscaff22"
]

log.info """
    B L U E B E R R Y   G E N O T Y P I N G   (Nextflow)
    ====================================================
    Samples   : ${params.samples}
    Ref       : ${params.ref}
    Outdir    : ${params.outdir}
    """

process TRIM {
    tag "$sample"
    publishDir "${params.outdir}/trimmomatic", mode: 'copy'

    // O módulo carrega a variável $HPC_TRIMMOMATIC_ADAPTER
    // Definido no nextflow.config: module = 'trimmomatic'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*_paired.fq.gz"), emit: paired
    path "*_unpaired.fq.gz"
    path "*.log"

    script:
    """
    # Usamos \$ para escapar e dizer ao Nextflow: "Isso é uma variável do Bash"
    ADAPTERS="\$HPC_TRIMMOMATIC_ADAPTER/TruSeq3-PE.fa"

    trimmomatic PE -threads ${task.cpus} -phred33 \
        ${reads[0]} ${reads[1]} \
        ${sample}_R1_paired.fq.gz ${sample}_R1_unpaired.fq.gz \
        ${sample}_R2_paired.fq.gz ${sample}_R2_unpaired.fq.gz \
        ILLUMINACLIP:\${ADAPTERS}:2:30:10 \
        SLIDINGWINDOW:4:20 TRAILING:20 MINLEN:50 \
        2> ${sample}.trim.log
    """
}

process ALIGN {
    tag "$sample"
    // Nota: BWA precisa que os índices (.amb, .ann, etc) estejam junto com a ref
    
    input:
    tuple val(sample), path(reads)
    path ref
    path ref_indices // Truque para o Nextflow trazer os índices para a pasta de trabalho

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} -M \
        -R "@RG\\tID:${sample}\\tLB:lib1\\tPL:ILLUMINA\\tPU:10K\\tSM:${sample}" \
        $ref ${reads[0]} ${reads[1]} \
        | samtools view -hb - \
        | samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -

    samtools index ${sample}.sorted.bam
    """
}

process QC_BAM {
    tag "$sample"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path "${sample}.flagstat.txt"
    path "${sample}.stats.txt"

    script:
    """
    samtools flagstat $bam > ${sample}.flagstat.txt
    samtools stats $bam    > ${sample}.stats.txt
    # bam validate removido para economizar tempo, descomente se necessário
    # bam validate --in $bam --so_coord --verbose > ${sample}.validate.txt || true
    """
}

process INDIVIDUAL_CALL {
    tag "$sample"
    publishDir "${params.outdir}/pileup", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref
    path ref_indices

    output:
    path "${sample}.sorted_norm.vcf.gz", emit: vcf
    path "${sample}.sorted_norm.vcf.gz.tbi", emit: tbi

    script:
    """
    # 1. Mpileup
    bcftools mpileup -f $ref \
        --annotate FORMAT/AD,FORMAT/DP --min-MQ 20 \
        $bam -Oz -o raw.vcf.gz

    # 2. Sort & Norm (Exatamente como no seu script bash)
    bcftools sort raw.vcf.gz \
        | bcftools norm -O u --atomize -f $ref \
        | bcftools norm --multiallelics -any -f $ref -O z -o ${sample}.sorted_norm.vcf.gz
    
    tabix -f -p vcf ${sample}.sorted_norm.vcf.gz
    """
}

process MERGE_AND_CALL {
    tag "$region"
    publishDir "${params.outdir}/merge", mode: 'copy'

    input:
    path vcf_files
    path tbi_files
    val region

    output:
    path "merged.*.called.vcf.gz", emit: vcf

    script:
    // Limpeza de caracteres para nome de arquivo seguro
    safe_region = region.replace(':', '_').replace('-', '_')
    """
    # Gera lista de arquivos presentes
    ls -1 *.vcf.gz | sort > vcf_list.txt

    echo "Merging region: ${region}..."

    # Merge usando a região específica (-r) para performance
    bcftools merge -Oz --threads ${task.cpus} -m none \
        --file-list vcf_list.txt \
        -r "${region}" \
        -o merged.${safe_region}.all.vcf.gz

    tabix -f -p vcf merged.${safe_region}.all.vcf.gz

    # Variant Calling
    bcftools call -mv -Oz \
        -o merged.${safe_region}.called.vcf.gz \
        merged.${safe_region}.all.vcf.gz
        
    tabix -f -p vcf merged.${safe_region}.called.vcf.gz
    """
}

workflow {
    // Carrega amostras (Header obrigatório: sample, r1, r2)
    Channel.fromPath(params.samples)
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sample, [file(row.r1), file(row.r2)]) }
        .set { fastq_ch }

    // Define canais de referência e índices para evitar erro do BWA
    ref_file = file(params.ref)
    // Pega todos os arquivos na mesma pasta da ref que começam com o nome da ref
    ref_indices = Channel.fromPath("${params.ref}.*").collect()

    // Processamento
    TRIM(fastq_ch)
    
    ALIGN(TRIM.out.paired, ref_file, ref_indices)
    
    QC_BAM(ALIGN.out.bam)
    
    INDIVIDUAL_CALL(ALIGN.out.bam, ref_file, ref_indices)

    // Coleta todos os VCFs gerados
    vcf_list = INDIVIDUAL_CALL.out.vcf.collect()
    tbi_list = INDIVIDUAL_CALL.out.tbi.collect()
    
    // Canal de Regiões
    regions_ch = Channel.fromList(params.chromosomes)

    // Merge paralelo por região
    MERGE_AND_CALL(vcf_list, tbi_list, regions_ch)
}
