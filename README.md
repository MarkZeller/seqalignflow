#MetaSeqAlignFlow

MetaSeqAlignFlow is a straight-forward pipeline that aligns reads against a reference and generates a bam file for metagenomic sequencing datasets. MetaSeqAlignFlow requires bbtools (https://jgi.doe.gov/data-and-tools/bbtools/), bwa (https://sourceforge.net/projects/bio-bwa/), samtools (http://samtools.sourceforge.net/), and iVar (https://github.com/andersen-lab/ivar) to be available in your path. Specify a reference fasta file and the output directory, go to the directory containing your fastq files and run:

```
nextflow run main.nf
```

MetaSeqAlignFlow contains the following processes:
- deduplicate reads
- merge reads
- trim reads
- align reads
- combine bam files
- generate consensus
