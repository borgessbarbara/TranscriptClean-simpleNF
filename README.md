# TranscriptClean-simpleNF
Simplified TranscriptClean implementation into Nextflow.

## Mandatory inputs

- SAM file (`.sam`, input with `--sam`)
- FASTA file (`.fasta` or `.fa`, input with `--fasta`)

## Basic Usage 

```
nextflow run transcriptclean.nf \\ 
--sam path/to/file.sam \\
--fasta path/to/file.fasta \\
```

## Optional inputs and params

- Tabular `.tab`  file with splice junctions outputed from STAR alignment for splice junction correction (set `--sj_correction true` and input with `--splice_junctions`)
- GTF for splice junction extraction before splice junction correction (set `--sj_correction true` and input with `--gtf`)
- VCF for variant-aware correction (set `--variant_aware true` and input with `--vcf`)
