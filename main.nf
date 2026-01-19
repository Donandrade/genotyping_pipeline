nextflow.enable.dsl=2

// --- HELPER FUNCTIONS ---
def create_chunks(fai_path, chunk_size) {
    def chunks = []
    file(fai_path).eachLine { line ->
        def fields = line.split('\t')
        def chrom = fields[0]
        def size = fields[1].toLong()

        if (size <= chunk_size) {
            chunks << "${chrom}"
        } else {
            for (long start = 0; start < size; start += chunk_size) {
                long end = Math.min(start + chunk_size, size)
                chunks << "${chrom}:${start + 1}-${end}"
            }
        }
    }
    return chunks
}

// --- PROCESSES ---

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
        ILLUMINACLIP:\${ADAPTERS}:2:30:10 \
        SLIDINGWINDOW:4:20 \
        TRAILING:20 \
        MINLEN:50 \
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
        tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam
        
    script:
    """
    # Using -R to define a clean Read Group. 
    # This prevents the "2:2:2:" naming snowball by standardizing sample IDs.
    bwa mem -t ${task.cpus} \\
        -R "@RG\\tID:${sample}\\tLB:lib1\\tPL:ILLUMINA\\tPU:unit1\\tSM:${sample}" \\
        $ref $r1 $r2 | \\
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -
    
    samtools index ${sample}.sorted.bam
    """
}


process INDIVIDUAL_PILEUP {
    tag "$sample"
    input:
        tuple val(sample), path(bam), path(bai)
        path ref
        path ref_indices
        path probes_file
    output:
        tuple val(sample), path("${sample}.raw.vcf.gz"), path("${sample}.raw.vcf.gz.tbi"), emit: vcf
    script:
    def target_opt = params.probes ? "-T ${probes_file}" : ""
    """
    bcftools mpileup ${target_opt} -f ${ref} --annotate FORMAT/AD,FORMAT/DP ${bam} -Oz -o ${sample}.raw.vcf.gz
    tabix -p vcf ${sample}.raw.vcf.gz
    """
}

process SPLIT_VCF_TO_CHUNK {
    tag "$sample | $chunk"
    publishDir "${params.outdir}/04_splited_call", mode: 'copy'
    input:
        tuple val(sample), path(vcf), path(tbi), val(chunk)
    output:
        tuple val(chunk), path("${sample}.${safe_chunk}.chunk.vcf.gz"), path("${sample}.${safe_chunk}.chunk.vcf.gz.tbi"), emit: chunk_vcf
    script:
    safe_chunk = chunk.replace(':', '_').replace('-', '_')
    """
    bcftools view -r ${chunk} ${vcf} -Oz -o ${sample}.${safe_chunk}.chunk.vcf.gz
    tabix -p vcf ${sample}.${safe_chunk}.chunk.vcf.gz
    """
}

process MERGE_AND_CALL_BY_CHUNK {
    tag "$chunk"
    publishDir "${params.outdir}/04_final_calls/chunks", mode: 'copy'

    input:
        // O segredo está aqui: o Nextflow renomeia os arquivos no 'stage' do processo
        tuple val(chunk), path(vcfs), path(tbis), path('past.vcf.gz'), path('past.vcf.gz.tbi')
        path ref

    output:
        path "merged.${safe_chunk}.pileup.vcf.gz", emit: vcf_merged
        path "merged.${safe_chunk}.called.vcf.gz", emit: vcf_called
        path "merged.${safe_chunk}.called.vcf.gz.tbi", emit: tbi_called
        path "${safe_chunk}.stats.txt", emit: stats

    script:
    safe_chunk = chunk.replace(':', '_').replace('-', '_')
    """
    # 1. Merge das amostras atuais (2026)
    echo "${vcfs.join('\n')}" > vcf_list.txt
    bcftools merge -Oz --threads ${task.cpus} -l vcf_list.txt -o novas_amostras.vcf.gz
    tabix -f -p vcf novas_amostras.vcf.gz

    # 2. Lógica de Merge Histórico
    # Verificamos se 'past.vcf.gz' é um VCF real (tem conteúdo e não é o dummy vazio)
    if [ -s past.vcf.gz ] && [ "\$(zcat past.vcf.gz | head -c 1 | wc -c)" -ne 0 ]; then
        
        # Garantimos que o índice do passado está atualizado (resolve o erro 255)
        tabix -f -p vcf past.vcf.gz
        
        bcftools merge --force-samples -Oz --threads ${task.cpus} \
            novas_amostras.vcf.gz past.vcf.gz \
            -o tmp_merged.vcf.gz
        mv tmp_merged.vcf.gz merged.${safe_chunk}.pileup.vcf.gz
    else
        # Se for a primeira rodada ou passado inválido, usamos apenas as novas
        mv novas_amostras.vcf.gz merged.${safe_chunk}.pileup.vcf.gz
    fi

    # 3. Indexação do Pileup Final (essencial para o calling)
    tabix -f -p vcf merged.${safe_chunk}.pileup.vcf.gz

    # 4. Variant Calling
    # Mesmo se houver apenas header, o bcftools call lida com isso, mas checamos variantes:
    V_COUNT=\$(bcftools view -H merged.${safe_chunk}.pileup.vcf.gz | head -n 1 | wc -l)

    if [ "\$V_COUNT" -gt 0 ]; then
        bcftools call -m -Oz -o merged.${safe_chunk}.called.vcf.gz merged.${safe_chunk}.pileup.vcf.gz
    else
        echo "Chunk ${chunk} sem variantes. Mantendo apenas header."
        cp merged.${safe_chunk}.pileup.vcf.gz merged.${safe_chunk}.called.vcf.gz
    fi

    tabix -f -p vcf merged.${safe_chunk}.called.vcf.gz
    bcftools stats merged.${safe_chunk}.called.vcf.gz > ${safe_chunk}.stats.txt
    """
}

process CONCATENATE_ALL {
    tag "Final Join"
    publishDir "${params.outdir}/04_final_calls/global", mode: 'copy'
    input: path(called_vcfs)
    output: 
        path "genome_wide_final.vcf.gz"
        path "genome_wide_final.vcf.gz.tbi"

    script:
    """
    # 1. Gera a lista mestre de amostras
    bcftools query -l \$(ls *.called.vcf.gz | head -n 1) > samples_list.txt

    # 2. Padroniza o header e INDEXA cada arquivo novamente
    for f in *.called.vcf.gz; do
        bcftools reheader -s samples_list.txt \$f -o fixed_\$f
        mv fixed_\$f \$f
        # O concat precisa do índice para funcionar com a flag -a
        tabix -p vcf \$f
    done

    # 3. Concatena agora que todos têm headers idênticos e índices novos
    bcftools concat -a -Oz -o genome_wide_final.vcf.gz \$(ls *.called.vcf.gz | sort -V)

    # 4. Índice final do arquivo global
    tabix -p vcf genome_wide_final.vcf.gz
    """
}

// --- WORKFLOW ---

workflow {
    fastq_ch = Channel.fromPath(params.samples).splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sample, [file(row.r1), file(row.r2)]) }

    ref_file = file(params.ref)
    ref_indices = Channel.fromPath("${params.ref}.*").collect()
    probes_file = params.probes ? file(params.probes) : []

    TRIM(fastq_ch)
    ALIGN(TRIM.out.paired, ref_file, ref_indices)
    INDIVIDUAL_PILEUP(ALIGN.out.bam, ref_file, ref_indices, probes_file)

    def chunks_list = create_chunks("${params.ref}.fai", params.chunk_size)
    ch_chunks = Channel.fromList(chunks_list)

    ch_split_input = INDIVIDUAL_PILEUP.out.vcf.combine(ch_chunks)
    SPLIT_VCF_TO_CHUNK(ch_split_input)

    vcf_grouped_by_chunk = SPLIT_VCF_TO_CHUNK.out.chunk_vcf.groupTuple(by: 0)

    ch_combined_input = vcf_grouped_by_chunk.map { chunk, vcfs, tbis ->
        def safe_chunk = chunk.replace(':', '_').replace('-', '_')
        
        // Criamos placeholders locais para quando não houver passado
        def dummy_vcf = file("${workDir}/dummy_${safe_chunk}.vcf")
        if(!dummy_vcf.exists()) dummy_vcf.text = "" // Arquivo vazio de sinalização

        def past_vcf = dummy_vcf
        def past_tbi = dummy_vcf // No NF, podemos passar o mesmo arquivo se não for usado

        if (params.past_calls != "false") {
            def vcf_path = file("${params.past_calls}/merged.${safe_chunk}.pileup.vcf.gz")
            if (vcf_path.exists()) {
                past_vcf = vcf_path
                past_tbi = file("${vcf_path}.tbi")
            }
        }
        return tuple(chunk, vcfs, tbis, past_vcf, past_tbi)
    }

    MERGE_AND_CALL_BY_CHUNK(ch_combined_input, ref_file)

    CONCATENATE_ALL(MERGE_AND_CALL_BY_CHUNK.out.vcf_called.collect())
}
