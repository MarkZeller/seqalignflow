#!/usr/bin/env nextflow

log.info"""\

S E Q A L I G N F L O W - P I P E L I N E    
===================================
reference               : ${params.ref}
reads                   : ${params.reads}
output folder           : ${params.output}
adapters for trimming   : ${params.adapt}
"""

include { TRIMREADS } from './modules/preprocessing_reads'
include { ALIGNREADS; GETCOVERAGE; GETCONSENSUS } from './modules/align_reads'


workflow {
    reads = Channel.fromFilePairs( params.reads, checkIfExists: true )

    TRIMREADS( reads ) 
    ALIGNREADS( TRIMREADS.out )
    GETCOVERAGE( ALIGNREADS.out )
    GETCONSENSUS( GETCOVERAGE.out )
}

workflow.onComplete {
    log.info "[AssembleFlow] Pipeline Complete"
}


process TRIMREADS {
    publishDir "${params.outdir}/trimmed_reads/", mode: 'copy', pattern: "*_trimmed.fastq.gz"
    tag { sample_id }

    input:
    tuple val(sample_id), path(reads)
 
    output:
    tuple val(sample_id), path("*_trimmed.fastq.gz")

    script:

    """
    bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=${sample_id}_R1_trimmed.fastq.gz out2=${sample_id}_R2_trimmed.fastq.gz ref=${param.adapters} ktrim=r k=23 mink=4 hdist=1 minlength=50 qtrim=w trimq=20 tpe 

    """
}

process ALIGN {
    publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: "*.bam{.bai}"
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
    publishDir params.output, mode: 'copy', pattern: "*.png"
    tag { sample_id }

    input:
    tuple val(sample_id), path(bam)
 
    output:
    tuple val(sample_id), path("${sample_id}_coverage.tsv"), path("*.png")

    script:

    """
    samtools depth -a -d 0 ${sample_id}.bam > ${sample_id}_coverage.tsv
    python3 plot_genome_coverage.py -i ${sample_id}_coverage.tsv -d 10 -s ${sample_id} -o ${sample_id}.png 
    python3 plot_genome_coverage.py -i ${sample_id}_coverage.tsv -d 10 -s ${sample_id} -o ${sample_id}_normalized.png -n

    """
}

process GETCONSENSUS {
    publishDir "${params.outdir}/consensus_sequences/", mode: 'copy', pattern: "${sample_id}*"
    tag { sample_id }


    input:
    tuple val(sample_id), path(bam)
 
    output:
    tuple val(sample_id), path("${sample_id}*")

    script:

    """
    samtools mpileup -aa -A -d 0 -Q 0 ${bam} | ivar consensus -p ${sample_id}

    """
}

workflow.onComplete {
    log.info "[SeqAlignFlow] Pipeline Complete"
}

