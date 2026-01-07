nextflow.enable.dsl=2

// --- PROCESSOS ---

process TRIM {
    tag "$sample"
    publishDir "${params.outdir}/01_trimmed", mode: 'copy'

    input:
        tuple val(sample), path(reads)

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

process ALIGN {
    tag "$sample"
    publishDir "${params.outdir}/02_bam", mode: 'copy'

    input:
        tuple val(sample), path(r1), path(r2)
        path ref
        path ref_indices

    output:
        tuple val(sample), path("${sample}.sorted.group.bam"), path("${sample}.sorted.group.bam.bai"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} -M \
        -R "@RG\\tID:${sample}\\tLB:lib1\\tPL:ILLUMINA\\tPU:10K\\tSM:${sample}" \
        $ref $r1 $r2 | \
    samtools view -hb - | \
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.group.bam -

    samtools index ${sample}.sorted.group.bam
    """
}

process QC_BAM {
    tag "$sample"
    publishDir "${params.outdir}/05_qc/samtools", mode: 'copy'
    
    input:
        tuple val(sample), path(bam), path(bai)
        
    output:
        path "${sample}.*"

    script:
    """
    samtools flagstat $bam > ${sample}.flagstat.txt
    samtools stats $bam > ${sample}.stats.txt
    """
}

process COVERAGE_QC {
    tag "$sample"
    publishDir "${params.outdir}/05_qc/coverage", mode: 'copy'
    
    input:
        tuple val(sample), path(bam), path(bai)
        
    output:
        path "${sample}.mosdepth.*"

    script:
    """
    mosdepth -n -x --threads ${task.cpus} ${sample} ${bam}
    """
}

process INDIVIDUAL_CALL {
    tag "$sample | $region"
    publishDir "${params.outdir}/03_individual_vcfs", mode: 'copy'

    input:
        tuple val(sample), path(bam), path(bai), val(region)
        path ref
        path ref_indices
        path probes_file

    output:
        tuple val(region), path("${sample}.${safe_region}.vcf.gz"), emit: vcf
        path "${sample}.${safe_region}.vcf.gz.tbi", emit: tbi

    script:
    safe_region = region.replace(':', '_').replace('-', '_')
    def target_opt = params.probes ? "-T ${probes_file}" : "-r ${region}"
    """
    bcftools mpileup ${target_opt} -f ${ref} --annotate FORMAT/AD,FORMAT/DP --min-MQ 0 ${bam} -Oz -o raw.vcf.gz
    bcftools sort raw.vcf.gz | \
    bcftools norm -O u --atomize -f ${ref} | \
    bcftools norm --multiallelics -any -f ${ref} -O z -o ${sample}.${safe_region}.vcf.gz
    tabix -f -p vcf ${sample}.${safe_region}.vcf.gz
    """
}

process MERGE_AND_CALL {
    tag "$region"
    publishDir "${params.outdir}/04_final_calls", mode: 'copy'

    input:
        tuple val(region), path(vcfs)
        path tbis
        path ref

    output:
        path "merged.${safe_region}.all.vcf.gz", emit: vcf_raw
        path "merged.${safe_region}.called.vcf.gz", emit: vcf_called
        path "${safe_region}.vcf_stats.txt", emit: stats

    script:
    safe_region = region.replace(':', '_').replace('-', '_')
    def region_filter = (region == 'target_probes') ? "" : "-r ${region}"

    """
    # Criar lista de arquivos garantida
    for f in ${vcfs}; do
        echo \$f >> vcf_list.txt
    done
    sort -u vcf_list.txt > vcf_list_sorted.txt

    bcftools merge -Oz --threads ${task.cpus} -m none \
        --file-list vcf_list_sorted.txt \
        ${region_filter} \
        -o merged.${safe_region}.all.vcf.gz

    # Mantido sem -v para diagnóstico de conteúdo
    bcftools call -m -Oz -o merged.${safe_region}.called.vcf.gz merged.${safe_region}.all.vcf.gz

    bcftools stats merged.${safe_region}.called.vcf.gz > ${safe_region}.vcf_stats.txt
    tabix -f -p vcf merged.${safe_region}.called.vcf.gz
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
    fastq_ch = Channel.fromPath(params.samples).splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sample, [file(row.r1), file(row.r2)]) }
    
    ref_file = file(params.ref)
    ref_indices = Channel.fromPath("${params.ref}.*").collect()
    mqc_config = file(params.multiqc_config)

    // 1. Processamento Inicial
    TRIM(fastq_ch)
    ALIGN(TRIM.out.paired, ref_file, ref_indices)

    // 2. Controle de Qualidade (QC)
    ch_samtools_logs = QC_BAM(ALIGN.out.bam)
    ch_mosdepth_logs = COVERAGE_QC(ALIGN.out.bam)

    // 3. Seleção de Alvo e Calling
    if (params.probes) {
        ch_regions = Channel.fromList(['target_probes'])
        ch_probes_file = Channel.fromPath(params.probes).collect()
    } else {
        ch_regions = Channel.fromPath(params.regions).splitText().map{ it.trim() }
        ch_probes_file = Channel.value([])
    }

    INDIVIDUAL_CALL(ALIGN.out.bam.combine(ch_regions), ref_file, ref_indices, ch_probes_file)

    // 4. Agrupamento e Merge
    vcf_grouped = INDIVIDUAL_CALL.out.vcf.groupTuple(by: 0)
    tbi_all = INDIVIDUAL_CALL.out.tbi.collect()

    MERGE_AND_CALL(vcf_grouped, tbi_all, ref_file)

    // 5. Relatórios Finais (Coleta todos os logs produzidos)
    ch_multiqc_input = Channel.empty()
        .mix(TRIM.out.log)
        .mix(ch_samtools_logs)
        .mix(ch_mosdepth_logs)
        .mix(MERGE_AND_CALL.out.stats)
        .collect()

    MULTIQC(ch_multiqc_input, mqc_config)
}
