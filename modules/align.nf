process ALIGN {
    publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: "*.bam"
    publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: "*.bam.bai"
    tag { sample_id }

    input:
    tuple val(sample_id), path(reads)
 
    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
    
    script:

    """
    bwa index ${params.ref}
    bwa mem -t $task.cpus ${params.ref} ${reads[0]} ${reads[1]} | samtools view -u -@ $task.cpus - | samtools sort -@ $task.cpus -o ${sample_id}.bam -
    samtools index ${sample_id}.bam ${sample_id}.bam.bai

    """
}

process GETCOVERAGE {
    publishDir "${params.outdir}/coverage_plots/", mode: 'copy', pattern: "*.png"
    tag { sample_id }

    input:
    tuple val(sample_id), path(bam)
 
    output:
    tuple val(sample_id), path("${sample_id}_coverage.tsv"), path("*.png")

    script:

    """
    samtools depth -a -d 0 ${bam} > ${sample_id}_coverage.tsv
    plot_genome_coverage.py -i ${sample_id}_coverage.tsv -d 10 -s ${sample_id} -o ${sample_id}.png 
    plot_genome_coverage.py -i ${sample_id}_coverage.tsv -d 10 -s ${sample_id} -o ${sample_id}_normalized.png -n

    """
}

process GETCONSENSUS {
    publishDir "${params.outdir}/consensus_sequences/", mode: 'copy', pattern: "${sample_id}.fa"
    tag { sample_id }


    input:
    tuple val(sample_id), path(bam)
 
    output:
    tuple val(sample_id), path("${sample_id}.fa")

    script:

    """
    samtools mpileup -aa -A -d 0 -Q 0 ${bam} | ivar consensus -p ${sample_id}

    """
}
