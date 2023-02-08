#!/usr/bin/env nextflow

log.info"""\

S E Q A L I G N F L O W - P I P E L I N E    
===================================
reference               : ${params.ref}
reads                   : ${params.reads}
output folder           : ${params.outdir}
adapters for trimming   : ${params.adapt}
"""

include { TRIMREADS } from './modules/trim_reads'
include { ALIGNREADS; GETCOVERAGE; GETCONSENSUS } from './modules/align_reads'

workflow {
    reads = Channel.fromFilePairs( params.reads, checkIfExists: true )

    TRIMREADS( reads ) 
    ALIGNREADS( TRIMREADS.out )
    GETCOVERAGE( ALIGNREADS.out )
    GETCONSENSUS( ALIGNREADS.out )
}

workflow.onComplete {
    log.info "[SeqAlignFlow] Pipeline Complete"
}

