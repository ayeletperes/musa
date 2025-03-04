# VDJbase repertoire annotation and downstream analysis for the light chain R24 macaque sequences


The nextflow pipeline performs anotation and downstrean analysis of AIRR-seq.

The pipeline can be devided into seven main componenet:

**1. Initial repertoire alignment and annotation based reference set**

> In this section, the repertoire sequences are annotated using IgBlast and MakeDb (chango) against the supplied reference set and collapse tham.

**2. Undocumented allele inference**

> In this section, inference for undocumented V allele (not found in the initial reference set) is performed using TIgGER.

**3. Second repertoire alignment and annotation with the discovered undocumented alleles.**

> In this section, in case undocumented alleles were inferred the repertoire is re-aligned and annotated with the additional alleles.

**4. Pre genotype inference**

> In this section, Sequences were required to have no mutations within the V region, a single V allele assignment, and alignment starting at position one of the V germline.

**5. Genotype inference**

> In this section, genotypes are inferred for each of the IG calls: V, D, and J using TIgGER Bayesian inferernce tool, and a personal genotype reference set is created.

**6. Third repertoire alignment and annotation with the personal reference set.**

> In this section, the repertoire sequences are annotated using IgBlast and MakeDb (chango) against the personal genotype reference set from step 5.

**7. OGRDB statistics.**

> the reperotire statistics are deduced using ogrdbstats.


### Input files:

1. An AIRR-seq in fasta format.
2. heavy_chain no
3. chain - (IGLV/IGKV)
4. ndm_chain - (IGL/IGK)
5. Personal reference set files for IG\*V, IGHD, and IG*J alleles in fasta format.
6. Optimized thresholds files for IG\*V, IG*J alleles.


### Output:

1. The aligned repertoire with the personal reference set.
2. The aligned repertoire with the init reference set.
3. A genotype reports
4. A log files
5. An ogrdb statistics report for the init aligned repertoire
6. An ogrdb statistics report for the final aligned repertoire
7. pipeline_statistic.csv - table of pass and fail reads for some of the steps.

### Docker images: 

The pipeline uses two docker images:

1. peresay/suite
2. williamlees/ogrdbstats



### How to run the annotation pipeline:

```bash
nextflow run main.nf --airr_seq {sampleName}.fasta --v_germline_file V_gapped.fasta --d_germline D.fasta --j_germline J.fasta
--v_optimized_thresholds IGHV_optimized_thresholds.tsv --j_optimized_thresholds IGHJ_optimized_thresholds.tsv 
--chain IGLV\IGKV --ndm_chain IGL\IGK --heavy_chain yes -w {work_dir} --outdir {out_dir} --sample_name {sampleName} --nproc {nproc} -with-docker -resume
```
