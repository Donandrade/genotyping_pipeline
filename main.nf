nextflow.enable.dsl=2

// --- HELPER FUNCTIONS ---

/**
 * Generates chunks. 
 * If chromosome < chunk_size: returns "chrom"
 * If chromosome > chunk_size: returns "chrom:start-end"
 */

def create_chunks(fai_path, chunk_size) {
    def chunks = []
    file(fai_path).eachLine { line ->
        def fields = line.split('\t')
        def chrom = fields[0]
        def size = fields[1].toLong()

        if (size <= chunk_size) {
            // Use the whole chromosome if it's smaller than the chunk limit
            chunks << "${chrom}"
        } else {
            // Split into coordinate-based windows
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
    # Define a variável localmente para garantir a expansão correta
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
    output: tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam
    script:
    """
    bwa mem -t ${task.cpus} $ref $r1 $r2 | samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -
    samtools index ${sample}.sorted.bam
    """
}

// 1. Generate Individual Pileup
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

// 2. Split VCF into Chunks (The Scatter Step)
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

// 3. Joint Call per Chunk (The Gather Step)
process MERGE_AND_CALL_BY_CHUNK {
    tag "$chunk"
    publishDir "${params.outdir}/04_final_calls/chunks", mode: 'copy'
    input:
        tuple val(chunk), path(vcfs), path(tbis)
        path ref
    output:
        path "merged.${safe_chunk}.pileup.vcf.gz", emit: vcf_merged
        path "merged.${safe_chunk}.called.vcf.gz", emit: vcf_called
        path "${safe_chunk}.stats.txt", emit: stats
    script:
    safe_chunk = chunk.replace(':', '_').replace('-', '_')
    """
    bcftools merge -Oz --threads ${task.cpus} ${vcfs} -o merged.${safe_chunk}.pileup.vcf.gz
    bcftools call -m -Oz -o merged.${safe_chunk}.called.vcf.gz merged.${safe_chunk}.pileup.vcf.gz
    bcftools stats merged.${safe_chunk}.called.vcf.gz > ${safe_chunk}.stats.txt
    """
}

process CONCATENATE_ALL {
    tag "Final Join"
    publishDir "${params.outdir}/04_final_calls/global", mode: 'copy'
    input: path(called_vcfs)
    output: path "genome_wide_final.vcf.gz"
    script:
    """
    bcftools concat -Oz -o genome_wide_final.vcf.gz \$(ls *.vcf.gz | sort -V)
    tabix -p vcf genome_wide_final.vcf.gz
    """
}

// --- WORKFLOW ---

workflow {
    // Input setup
    fastq_ch = Channel.fromPath(params.samples).splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sample, [file(row.r1), file(row.r2)]) }

    ref_file = file(params.ref)
    ref_indices = Channel.fromPath("${params.ref}.*").collect()
    probes_file = params.probes ? file(params.probes) : []

    // Alignment phase
    TRIM(fastq_ch)
    ALIGN(TRIM.out.paired, ref_file, ref_indices)

    // Step 1: Broad Pileup
    INDIVIDUAL_PILEUP(ALIGN.out.bam, ref_file, ref_indices, probes_file)

    // Step 2: Chunk Generation (with size check)
    def chunks_list = create_chunks("${params.ref}.fai", params.chunk_size)
    ch_chunks = Channel.fromList(chunks_list)

    // Step 3: Scatter VCFs into the defined chunks
    // Every sample VCF is now split into the sub-regions defined in Step 2
    ch_split_input = INDIVIDUAL_PILEUP.out.vcf.combine(ch_chunks)
    SPLIT_VCF_TO_CHUNK(ch_split_input)

    // Step 4: Gather all samples for each specific chunk
    vcf_grouped_by_chunk = SPLIT_VCF_TO_CHUNK.out.chunk_vcf.groupTuple(by: 0)	

    // Step 5: Parallel Joint Calling
    MERGE_AND_CALL_BY_CHUNK(vcf_grouped_by_chunk, ref_file)

    // Step 6: Final Consolidation
    CONCATENATE_ALL(MERGE_AND_CALL_BY_CHUNK.out.vcf_called.collect())
}
