# SeqAlignFlow

SeqAlignFlow is a pipeline that aligns reads against a reference and generates a bam file, consensus sequence and coverage plots. SeqAlignFlow is specifically developed for metagenomic sequencing datasets.

### Requirements
* Linux or MacOS
* bbtools (https://jgi.doe.gov/data-and-tools/bbtools/)
* bwa
* samtools
* Seaborn

### Contents

The entire workflow has been implemented in NextFlow. An overview of each step in assembleflow:

* Trim reads
* Align reads to reference
* Generate (normailized) coverage plot
* Generate consensus sequence

### Usage

SeqAlignFlow uses the `nextflow.config` file to provide parameters to the pipeline. To run the pipeline simply upodate the `nextflow.config` file:
Areas to configure: 
* reads: directory containing paired-end reads
* adapt: adapter sequences in fasta format
* ref: reference sequence in fasta format (only 1 reference is allowed at this moment)
* outdir: output directory

### Running
Once your workflow.sh has been configured you can start the workflow by simply running:
`nextflow run main.nf`