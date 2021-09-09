#!/usr/bin/env nextflow
 
params.reads = './*_R{1,2}.fastq.gz'
params.ref = ''
params.output = ''
params.cpus = 8
params.mem = '16GB'

log.info"""\

S E Q A L I G N F L O W - P I P E L I N E    
===================================
reference               : ${params.ref}
reads                   : ${params.reads}
output folder           : ${params.output}
adapters for trimming   : ${params.adapt}
"""

Channel 
    .fromFilePairs( params.reads )
    .ifEmpty { error "Oops! Cannot find any file matching: ${params.reads}"  }
    .set { read_pairs }

adapters = file(params.adapt)
reference = file(params.ref)

process deduplicateReads {

    input:
    set pair_id, file(reads) from read_pairs
 
    output:
    set pair_id, '*deduplicated.fastq.gz' into deduplicated_reads

    script:

    """
    clumpify.sh in1=${reads[0]} in2=${reads[1]} out1=R1_deduplicated.fastq.gz out2=R2_deduplicated.fastq.gz dedupe subs=0 reorder
    """
}

process mergeReads {

    input:
    set pair_id, file(reads) from deduplicated_reads
 
    output:
    file 'merged.fastq.gz' into merged_reads
    set pair_id, '*unmerged.fastq.gz' into unmerged_reads


    script:

    """
    bbmerge.sh in1=${reads[0]} in2=${reads[1]} out=merged.fastq.gz outu1=R1_unmerged.fastq.gz outu2=R2_unmerged.fastq.gz vstrict

    """
}

process trimReads {

    input:
    file(merged) from merged_reads
    set pair_id, file(reads) from unmerged_reads
 
    output:
    set pair_id, '*_trimmed.fastq.gz' into trimmed_reads_paired_for_alignment
    file 'trimmed.fastq.gz' into trimmed_reads_single_for_alignment

    script:

    """
    bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=R1_trimmed.fastq.gz out2=R2_trimmed.fastq.gz ref=${adapters} ktrim=r k=23 mink=4 hdist=1 minlength=50 qtrim=w trimq=20 tpe tbo 
    bbduk.sh in=${merged} out=trimmed.fastq.gz ref=${adapters} ktrim=r k=23 mink=4 hdist=1 minlength=50 qtrim=w trimq=20 tpe tbo 

    """
}

process alignReads {
    publishDir params.output, mode: 'copy'

    input:
    set pair_id, file(reads) from trimmed_reads_paired_for_alignment
    file(single_reads) from trimmed_reads_single_for_alignment
 
    output:
    file 'alignment.bam'
    file 'alignment.bai'

    script:

    """
    bwa index ${reference}
    bwa mem ${reference} ${reads[0]} ${reads[1]} -t ${task.cpus} | samtools view -u -@ $task.cpus - | samtools sort -@ $task.cpus -o paired.bam -
    bwa mem ${reference} ${single_reads} -t $task.cpus | samtools view -u -@ $task.cpus - | samtools sort -@ $task.cpus -o single.bam -
    samtools merge alignment.bam paired.bam single.bam 
    samtools index alignment.bam alignment.bai

    """
}

workflow.onComplete {
    log.info "[SeqAlignFlow] Pipeline Complete"
}
