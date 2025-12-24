nextflow.enable.dsl=2

// --- PARÂMETROS ---
//params.chromosomes = ["VaccDscaff1", "VaccDscaff2", "VaccDscaff4", "VaccDscaff6", "VaccDscaff7", "VaccDscaff11", "VaccDscaff12", "VaccDscaff13", "VaccDscaff17", "VaccDscaff20", "VaccDscaff21", "VaccDscaff22"]

log.info """
    B L U E B E R R Y   G E N O T Y P I N G   (Nextflow)
    ====================================================
    Samples   : ${params.samples}
    Ref       : ${params.ref}
    Outdir    : ${params.outdir}
    """

// --- PROCESSOS ---

process TRIM {
    tag "$sample"
    publishDir "${params.outdir}/trimmomatic", mode: 'copy'
    input: tuple val(sample), path(reads)
    output:
        tuple val(sample), path("*_paired.fq.gz"), emit: paired
        path "*.log", emit: log
        path "*_unpaired.fq.gz"
    script:
    """
    ADAPTERS="\$HPC_TRIMMOMATIC_ADAPTER/TruSeq3-PE.fa"
    trimmomatic PE -threads ${task.cpus} -phred33 ${reads[0]} ${reads[1]} \\
        ${sample}_R1_paired.fq.gz ${sample}_R1_unpaired.fq.gz \\
        ${sample}_R2_paired.fq.gz ${sample}_R2_unpaired.fq.gz \\
        ILLUMINACLIP:\${ADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 TRAILING:20 MINLEN:50 \\
        2> ${sample}.trim.log
    """
}

process ALIGN {
    tag "$sample"
    input:
        tuple val(sample), path(reads)
        path ref
        path ref_indices

    output:
        tuple val(sample), path("${sample}.sorted.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} -M $ref ${reads[0]} ${reads[1]} | \\
    samtools view -hb - | \\
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -
    """
}

process ADD_READ_GROUPS {
    tag "$sample"
    input: tuple val(sample), path(bam)
    output: tuple val(sample), path("${sample}.sorted.group.bam"), emit: bam
    script:
    def avail_mem = task.memory ? task.memory.toGiga() - 4 : 4
    """
    java -Xmx${avail_mem}g -jar \$HPC_PICARD_DIR/picard.jar AddOrReplaceReadGroups \\
        INPUT=$bam OUTPUT=${sample}.sorted.group.bam \\
        RGID=$sample RGLB=lib1 RGPL=ILLUMINA RGPU=10K RGSM=$sample
    """
}

process MARK_DUPLICATES {
    tag "$sample"
    publishDir "${params.outdir}/bam_rmdup", mode: 'copy'
    input: tuple val(sample), path(bam)
    output:
        tuple val(sample), path("${sample}.sorted.group.rmdup.bam"), path("${sample}.sorted.group.rmdup.bam.bai"), emit: bam
        path "${sample}.duplicate.metrics", emit: metrics
    script:
    def avail_mem = task.memory ? task.memory.toGiga() - 4 : 8
    """
    java -Xmx${avail_mem}g -jar \$HPC_PICARD_DIR/picard.jar MarkDuplicates \\
        INPUT=$bam OUTPUT=${sample}.sorted.group.rmdup.bam \\
        METRICS_FILE=${sample}.duplicate.metrics ASSUME_SORT_ORDER=coordinate \\
        REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
    samtools index ${sample}.sorted.group.rmdup.bam
    """
}

process QC_BAM {
    tag "$sample"
    publishDir "${params.outdir}/reports", mode: 'copy'
    input: tuple val(sample), path(bam), path(bai)
    output: path "${sample}.*", emit: logs
    script:
    """
    samtools flagstat $bam > ${sample}.flagstat.txt
    samtools stats $bam > ${sample}.stats.txt
    """
}

process COVERAGE_QC {
    tag "$sample"
    input: tuple val(sample), path(bam), path(bai)
    output: path "${sample}.mosdepth.*", emit: logs
    script:
    """
    mosdepth -n -x --threads ${task.cpus} ${sample} ${bam}
    """
}

process INDIVIDUAL_CALL {
    tag "$sample"
    input:
        tuple val(sample), path(bam), path(bai)
        path ref
        path ref_indices

    output:
        path "${sample}.sorted_norm.vcf.gz", emit: vcf
        path "${sample}.sorted_norm.vcf.gz.tbi", emit: tbi

    script:
    """
    bcftools mpileup -f $ref --annotate FORMAT/AD,FORMAT/DP --min-MQ 20 $bam -Oz -o raw.vcf.gz
    bcftools sort raw.vcf.gz | bcftools norm -O u --atomize -f $ref | \\
    bcftools norm --multiallelics -any -f $ref -O z -o ${sample}.sorted_norm.vcf.gz
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
    path "merged.*.called.vcf.gz.tbi", emit: tbi
    path "*.vcf_stats.txt", emit: stats

    script:
    def safe_region = region.replace(':', '_').replace('-', '_')
    """
    # 1. Lista os VCFs recebidos para o merge
    ls -1 *.vcf.gz | sort > vcf_list.txt

    # 2. Merge de todas as amostras nesta região específica
    bcftools merge \\
        --threads ${task.cpus} \\
        -m none \\
        -r "${region}" \\
        --file-list vcf_list.txt \\
        -Oz -o merged.${safe_region}.all.vcf.gz

    # 3. Chamada de variantes (Genotipagem conjunta)
    bcftools call \\
        -mv \\
        --threads ${task.cpus} \\
        -Oz -o merged.${safe_region}.called.vcf.gz \\
        merged.${safe_region}.all.vcf.gz

    # 4. Indexação e Estatísticas
    tabix -p vcf merged.${safe_region}.called.vcf.gz
    bcftools stats merged.${safe_region}.called.vcf.gz > ${safe_region}.vcf_stats.txt

    # Limpeza do arquivo intermediário pesado
    rm merged.${safe_region}.all.vcf.gz
    """
}

process MULTIQC {
    tag "Geral"
    publishDir "${params.outdir}/multiqc_report", mode: 'copy'

    input:
    path logs
    path config

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . -c $config
    """
}

// --- WORKFLOW ---

workflow {
    // 1. Canais de Entrada
    fastq_ch = Channel.fromPath(params.samples)
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sample, [file(row.r1), file(row.r2)]) }

    ref_file = file(params.ref)
    ref_indices = Channel.fromPath("${params.ref}.*").collect()

    // 2. Alinhamento e Processamento
    ALIGN(fastq_ch, ref_file, ref_indices)
    ADD_READ_GROUPS(ALIGN.out.bam)
    MARK_DUPLICATES(ADD_READ_GROUPS.out.bam)

    // 3. Controle de Qualidade
    QC_BAM(MARK_DUPLICATES.out.bam)
    COVERAGE_QC(MARK_DUPLICATES.out.bam)

    // 4. Calling Individual
    INDIVIDUAL_CALL(MARK_DUPLICATES.out.bam, ref_file, ref_indices)

    // 5. Paralelismo por Regiões (Chunks de 5MB)
    regions_ch = Channel.fromPath("regions.txt")
        .splitText()
        .map { it.trim() }

    vcf_list = INDIVIDUAL_CALL.out.vcf.collect()
    tbi_list = INDIVIDUAL_CALL.out.tbi.collect()

    MERGE_AND_CALL(vcf_list, tbi_list, regions_ch)

    // 6. Relatórios Finais
    multiqc_config_file = file("multiqc_config.yaml")

    ch_reports = Channel.empty()
        .mix(MARK_DUPLICATES.out.metrics)
        .mix(QC_BAM.out.logs)
        .mix(COVERAGE_QC.out.logs)
        .mix(MERGE_AND_CALL.out.stats.collect()) // Coleta todas as stats
        .collect()

    MULTIQC(ch_reports, multiqc_config_file)
}
