#!/usr/bin/env nextflow

log.info"""\

S E Q A L I G N F L O W - P I P E L I N E    
===================================
reference               : ${params.ref}
reads                   : ${params.reads}
output folder           : ${params.outdir}
adapters for trimming   : ${params.adapt}
"""

include { TRIM } from './modules/trim'
include { ALIGN; GETCOVERAGE; GETCONSENSUS } from './modules/align'

workflow {
    reads = Channel.fromFilePairs( params.reads, checkIfExists: true )

    TRIM( reads ) 
    ALIGN( TRIM.out )
    GETCOVERAGE( ALIGN.out )
    GETCONSENSUS( ALIGN.out )
}

workflow.onComplete {
    log.info "[SeqAlignFlow] Pipeline Complete"
}

