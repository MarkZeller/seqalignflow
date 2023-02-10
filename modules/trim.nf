process TRIM {
    publishDir "${params.outdir}/trimmed_reads/", mode: 'copy', pattern: "*_trimmed.fastq.gz"
    tag { sample_id }

    input:
    tuple val(sample_id), path(reads)
 
    output:
    tuple val(sample_id), path("*_trimmed.fastq.gz")

    script:

    """
    bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=${sample_id}_R1_trimmed.fastq.gz out2=${sample_id}_R2_trimmed.fastq.gz ref=${params.adapt} ktrim=r k=23 mink=4 hdist=1 minlength=50 qtrim=w trimq=20 tpe 

    """
}
