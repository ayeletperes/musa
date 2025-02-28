The pre-processing nextflow pipeline for the R24 macaque sequences.

### The library protocol is as such:

The 5' of read 1 is comprised of (i) a given number of random bases for a sample (4-8) (ii) followed by an Illumina index  (iii) followed by the RACE primer sequence.

The UMI in the header is comprised of the 4-8 random bases. The TRIM comprises of (i),(ii), and (iii) which is retained in the header after trimming the sequence.



### Inputs:

1. Two fastq file of paired sequencing
2. sample_name
3. run_FastQC (yes/no)



### Outputs:

1. {sampleName}.fasta
2. log tab file for each steps
3. report for each steps
4. pipeline_statistic.csv - table of pass and fail reads for each step.


### Pipeline container:

The pipeline uses two docker images:

1. immcantation/suite:4.4.0
2. biocontainers/fastqc:v0.11.9_cv8


### How to run the pre-processing pipeline:

```bash
nextflow run main.nf --reads '{sampleName}_{R1,R2}.fastq.gz' -w {work_dir} --outdir {out_dir} --sample_name {sampleName} --run_FastQC yes --nproc {nproc} -with-docker -resume
```

### Sequence processing steps:


**1. FastQC**

> fastqc

**2. Paired-end assembly**

> AssemblePairs align

**3. Quality control**

> FilterSeq




