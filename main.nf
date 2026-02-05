nextflow.enable.dsl=2

// --- FUNÇÕES AUXILIARES ---

def get_chromosomes(fai_path) {
    def chroms = []
    file(fai_path).eachLine { line ->
        def fields = line.split('\t')
        chroms << fields[0]
    }
    return chroms
}

// MELHORIA: Garantir que strings "null" ou "false" sejam tratadas como inválidas
def is_valid(param) {
    return param != null && param.toString() != "null" && param.toString() != "false" && param.toString() != ""
}

// --- PROCESSOS ---

process TRIM {
    tag "$sample"
    input:  tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_R1_paired.fq.gz"), path("${sample}_R2_paired.fq.gz"), emit: paired
        path "${sample}.trim.log", emit: log
        path "*_unpaired.fq.gz"
    script:
    """
    ADAPTERS="\$HPC_TRIMMOMATIC_ADAPTER/TruSeq3-PE.fa"
    trimmomatic PE -threads ${task.cpus} -phred33 \
        ${reads[0]} ${reads[1]} \
        ${sample}_R1_paired.fq.gz ${sample}_R1_unpaired.fq.gz \
        ${sample}_R2_paired.fq.gz ${sample}_R2_unpaired.fq.gz \
        ILLUMINACLIP:\${ADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 TRAILING:20 MINLEN:50 \
        2> ${sample}.trim.log
    """
}

/*
process ALIGN {
    tag "$sample"
    input:  tuple val(sample), path(r1), path(r2); path ref; path ref_indices
    output: tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam
    script:
    """
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sample}\\tLB:lib1\\tPL:ILLUMINA\\tPU:unit1\\tSM:${sample}" $ref $r1 $r2 | \
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -
    samtools index ${sample}.sorted.bam
    """
}
*/

process ALIGN {
    tag "$sample"
    input:  tuple val(sample), path(r1), path(r2); path ref; path ref_indices
    output: 
        tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam
        path "${sample}.stats.txt", emit: stats // NOVO: estatísticas de alinhamento
    script:
    """
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sample}\\tLB:lib1\\tPL:ILLUMINA\\tPU:unit1\\tSM:${sample}" $ref $r1 $r2 | \
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -
    samtools index ${sample}.sorted.bam
    samtools stats -@ ${task.cpus} ${sample}.sorted.bam > ${sample}.stats.txt
    """
}

process INDIVIDUAL_PILEUP {
    tag "$sample"

    input:  tuple val(sample), path(bam), path(bai); path ref; path ref_indices; path probes_file

    output: tuple path("${sample}.raw.sort.norm.vcf.gz"), path("${sample}.raw.sort.norm.vcf.gz.tbi"), emit: vcf

    script:
    // Verificação para evitar a flag -T null
    def has_probes = (probes_file.name != 'EMPTY_FILE' && probes_file.name != 'null')
    def target_opt = has_probes ? "-T ${probes_file}" : ""
    """
    bcftools mpileup ${target_opt} -f ${ref} --annotate FORMAT/AD,FORMAT/DP ${bam} -Oz -o ${sample}.raw.vcf.gz
    tabix -p vcf ${sample}.raw.vcf.gz

    # New sort and normalized block
    bcftools sort ${sample}.raw.vcf.gz \
        | bcftools norm -O u --atomize -f ${ref} \
        | bcftools norm --multiallelics -any -f ${ref} -O z -o ${sample}.raw.sort.norm.vcf.gz

    tabix -f -p vcf ${sample}.raw.sort.norm.vcf.gz

    """
}

process MERGE_AND_CALL_BY_CHROM {
    tag "$chrom"
    publishDir "${params.outdir}/04_final_calls/chromosomes", mode: 'copy'

    input:
        tuple val(chrom), path(vcfs), path(tbis), path(past_vcf), path(past_tbi)
        path ref

    output:
        path "merged.${chrom}.pileup.vcf.gz", emit: vcf_merged
        path "merged.${chrom}.called.vcf.gz", emit: vcf_called
        path "merged.${chrom}.called.vcf.gz.tbi", emit: tbi_called
        path "${chrom}.stats.txt", emit: stats

    script:
    """
    echo "${vcfs.join('\n')}" > vcf_list.txt

    # 1. Merge das novas amostras
    bcftools merge -r "${chrom}" -Oz --threads ${task.cpus} -m none -l vcf_list.txt -o new_samples.vcf.gz
    tabix -f -p vcf new_samples.vcf.gz

    # 2. Lógica de Merge Histórico (CORRIGIDA para evitar conflito de nomes e falta de TBI)
    if [[ -s "${past_vcf}" && "${past_vcf}" != empty* ]]; then
        # Se o TBI histórico não existir no work dir, recria agora
        if [ ! -f "${past_vcf}.tbi" ]; then
            tabix -f -p vcf ${past_vcf}
        fi

        # Merge usando arquivo temporário para não sobrescrever a entrada durante a leitura
        bcftools merge -r "${chrom}" --force-samples -Oz --threads ${task.cpus} -m none \
            new_samples.vcf.gz ${past_vcf} -o total_temp.vcf.gz
        mv total_temp.vcf.gz merged.${chrom}.pileup.vcf.gz
    else
        mv new_samples.vcf.gz merged.${chrom}.pileup.vcf.gz
    fi
    tabix -f -p vcf merged.${chrom}.pileup.vcf.gz

    # 3. Variant Calling
    V_COUNT=\$(bcftools view -H merged.${chrom}.pileup.vcf.gz | head -n 1 | wc -l)
    if [ "\$V_COUNT" -gt 0 ]; then
        bcftools call -m -Oz -o merged.${chrom}.called.vcf.gz merged.${chrom}.pileup.vcf.gz
    else
        bcftools view -h merged.${chrom}.pileup.vcf.gz | bgzip -c > merged.${chrom}.called.vcf.gz
    fi

    tabix -f -p vcf merged.${chrom}.called.vcf.gz
    bcftools stats merged.${chrom}.called.vcf.gz > ${chrom}.stats.txt
    """
}

process CONCATENATE_ALL {
    tag "Final Join"
    publishDir "${params.outdir}/04_final_calls/global", mode: 'copy'
    input:
        path(called_vcfs)
        path(tbis)
    output:
        path "genome_wide_final.vcf.gz", emit: vcf
        path "genome_wide_final.vcf.gz.tbi", emit: tbi

    script:
    """
    bcftools concat -a -Oz -o genome_wide_final.vcf.gz \$(ls *.called.vcf.gz | sort -V)
    tabix -p vcf genome_wide_final.vcf.gz
    """
}


process MULTIQC {
    publishDir "${params.outdir}/05_QC/multiqc", mode: 'copy'
    
    input:
    path all_logs
    path config

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . --config $config
    """
}


// --- WORKFLOW ---

workflow {
    fastq_ch = Channel.fromPath(params.samples).splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sample, [file(row.r1), file(row.r2)]) }

    ref_file    = file(params.ref, checkIfExists: true)
    ref_indices = Channel.fromPath("${params.ref}.*").collect()
    
    // CORREÇÃO: probes_file agora é rigorosamente validado
    probes_file = is_valid(params.probes) ? file(params.probes) : file('EMPTY_FILE')

    TRIM(fastq_ch)
    ALIGN(TRIM.out.paired, ref_file, ref_indices)
    INDIVIDUAL_PILEUP(ALIGN.out.bam, ref_file, ref_indices, probes_file)

    def chromosomes = get_chromosomes("${params.ref}.fai")
    ch_chroms = Channel.fromList(chromosomes)

    ch_to_group = ch_chroms.combine(INDIVIDUAL_PILEUP.out.vcf)
        .map { chrom, vcf, tbi -> tuple(chrom, vcf, tbi) }

    vcf_grouped_by_chrom = ch_to_group.groupTuple(by: 0)

    // 3. Adicionar Lógica do Past VCF (CORRIGIDA contra argumentos nulos)
    ch_merge_input = vcf_grouped_by_chrom.map { chrom, vcfs, tbis ->
        // Inicia com placeholders seguros
        def past_vcf = file("empty_${chrom}.vcf.gz")
        def past_tbi = file("empty_${chrom}.vcf.gz.tbi")

        if (is_valid(params.past_calls)) {
            def vcf_path = file("${params.past_calls}/merged.${chrom}.pileup.vcf.gz")
            if (vcf_path.exists()) {
                past_vcf = vcf_path
                past_tbi = file("${vcf_path}.tbi")
            }
        }

        return tuple(chrom, vcfs, tbis, past_vcf, past_tbi)
    }

    MERGE_AND_CALL_BY_CHROM(ch_merge_input, ref_file)

    vcf_to_concat = MERGE_AND_CALL_BY_CHROM.out.vcf_called.collect()
    tbi_to_concat = MERGE_AND_CALL_BY_CHROM.out.tbi_called.collect()
    CONCATENATE_ALL(vcf_to_concat,tbi_to_concat)

    // Coleta de logs para o MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files
        .mix(TRIM.out.log)             // Survival rate do Trimmomatic
        .mix(ALIGN.out.stats)          // Mapeamento, erros, insert size
        .mix(MERGE_AND_CALL_BY_CHROM.out.stats) // Qualidade de SNPs e Indels
        .collect()

    // Execução do MultiQC
    // Certifique-se de que params.multiqc_config aponta para um arquivo real ou trate-o como opcional
    def mqc_config = params.multiqc_config ? file(params.multiqc_config) : []
    MULTIQC(ch_multiqc_files, mqc_config)

}
