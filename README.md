#SeqAlignFlow

SeqAlignFlow is a straight-forward pipeline that aligns reads against a reference and generates a bam file for each sample. SeqAlignFlow requires bbtools (https://jgi.doe.gov/data-and-tools/bbtools/).
Specify the location of the reference fasta and the output folder, go to the directory containing your reads and run

```
nextflow run main.nf
```

SeqAlignFlow contains the following steps:
- deduplication
- merge reads
- trim reads
- align reads
- combine bam files
