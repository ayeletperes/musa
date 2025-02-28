$HOSTNAME = ""
params.outdir = 'results'  

// global param
params.mate = "pair"
params.mate2 = "single"


// Process Parameters:
// Process Parameters for Assemble_pairs_assemble_pairs:
params.Assemble_pairs_assemble_pairs.method = "align" 
params.Assemble_pairs_assemble_pairs.coord = "illumina" 
params.Assemble_pairs_assemble_pairs.rc = "tail"
params.Assemble_pairs_assemble_pairs.head_fields_R1 = "UMI TRIM"
params.Assemble_pairs_assemble_pairs.head_fields_R2 = ""
params.Assemble_pairs_assemble_pairs.failed = "false"
params.Assemble_pairs_assemble_pairs.fasta = "false"
params.Assemble_pairs_assemble_pairs.nproc = "${params.nproc}"
params.Assemble_pairs_assemble_pairs.alpha = 0.00001
params.Assemble_pairs_assemble_pairs.maxerror = 0.3
params.Assemble_pairs_assemble_pairs.minlen = 8
params.Assemble_pairs_assemble_pairs.maxlen = 1000
params.Assemble_pairs_assemble_pairs.scanrev = "false"
params.Assemble_pairs_assemble_pairs.minident = 0.5
params.Assemble_pairs_assemble_pairs.evalue = 0.00001
params.Assemble_pairs_assemble_pairs.maxhits = 100
params.Assemble_pairs_assemble_pairs.fill = "false"
params.Assemble_pairs_assemble_pairs.aligner = "blastn"
params.Assemble_pairs_assemble_pairs.// align_exec =  ""
params.Assemble_pairs_assemble_pairs.// dbexec =  ""   
params.Assemble_pairs_assemble_pairs.gap = 0
params.Assemble_pairs_assemble_pairs.usearch_version = "11.0.667"
params.Assemble_pairs_assemble_pairs.assemble_reference =  ''

// Process Parameters for Assemble_pairs_parse_log_AP:
params.Assemble_pairs_parse_log_AP.field_to_parse = "ID REFID LENGTH OVERLAP GAP ERROR IDENTITY PVALUE EVALUE1 EVALUE2"

// Process Parameters for Filter_Sequence_Quality_filter_seq_quality:
params.Filter_Sequence_Quality_filter_seq_quality.method = "quality"
params.Filter_Sequence_Quality_filter_seq_quality.nproc = "${params.nproc}"
params.Filter_Sequence_Quality_filter_seq_quality.q = "20"
params.Filter_Sequence_Quality_filter_seq_quality.n_length = "35"
params.Filter_Sequence_Quality_filter_seq_quality.n_missing = "10"
params.Filter_Sequence_Quality_filter_seq_quality.fasta = "false"


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.mate2){params.mate2 = ""} 

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_6_reads_g_11}
 } else {  
	g_6_reads_g_11 = Channel.empty()
 }

Channel.value(params.mate).into{g_7_mate_g_11;g_7_mate_g_22;g_7_mate_g1_15;g_7_mate_g1_19;g_7_mate_g1_12}
Channel.value(params.mate2).into{g_10_mate_g2_7;g_10_mate_g2_5;g_10_mate_g2_0}


process unizp {

input:
 set val(name),file(reads) from g_6_reads_g_11
 val mate from g_7_mate_g_11

output:
 set val(name),file("*.fastq")  into g_11_reads0_g_22, g_11_reads0_g_27, g_11_reads0_g1_12

script:

readArray = reads.toString().split(' ')	
R1 = readArray.grep(~/.*R1.*/)[0]
R2 = readArray.grep(~/.*R2.*/)[0]

"""
case "$R1" in
*.gz | *.tgz ) 
        gunzip -c $R1 > R1.fastq
        ;;
*)
        cp $R1 ./R1.fastq
        echo "$R1 not gzipped"
        ;;
esac

case "$R2" in
*.gz | *.tgz ) 
        gunzip -c $R2 > R2.fastq
        ;;
*)
        cp $R2 ./R2.fastq
        echo "$R2 not gzipped"
        ;;
esac
"""
}


process preprocessing_meta_data {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*json$/) "meta_data/$filename"}
input:
 set val(name), file(files) from g_11_reads0_g_27

output:
 file "*json"  into g_27_outputFile00


"""
#!/usr/bin/env Rscript

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)


json_data <- list(
  sample = list(
    data_processing = list(
      preprocessing = list(
        paired_reads_assembly = "align",
        quality_thresholds = "20",
        software_versions = list(
            AssemblePairs = "0.7.1",
            FilterSeq = "0.7.1"
    	  )
	   )
	 )
  )
)

# Convert to JSON string without enclosing scalar values in arrays
json_string <- toJSON(json_data, pretty = TRUE, auto_unbox = TRUE)

# Write the JSON string to a file
writeLines(json_string, "pre_processed_metadata.json")

"""
}

//* params.run_FastQC =  "no"  //* @dropdown @options:"yes","no"
if (params.run_FastQC == "no") { println "INFO: FastQC will be skipped"}


process FastQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.(html|zip)$/) "fastqc/$filename"}
input:
 set val(name), file(reads) from g_11_reads0_g_22
 val mate from g_7_mate_g_22

output:
 file '*.{html,zip}'  into g_22_FastQCout00

errorStrategy 'retry'
maxRetries 5

script:
nameAll = reads.toString()
if (nameAll.contains('.gz')) {
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    file =  nameAll 
    runGzip = ''
}
"""
if [ "${params.run_FastQC}" == "yes" ]; then
    ${runGzip}
    fastqc ${file} 
else
    touch process.skiped.html
fi
"""
}


process Assemble_pairs_assemble_pairs {

input:
 set val(name),file(reads) from g_11_reads0_g1_12
 val mate from g_7_mate_g1_12

output:
 set val(name),file("*_assemble-pass.f*")  into g1_12_reads0_g2_0
 set val(name),file("AP_*")  into g1_12_logFile1_g1_15
 set val(name),file("*_assemble-fail.f*") optional true  into g1_12_reads_failed22
 set val(name),file("out*")  into g1_12_logFile3_g23_0

script:
method = params.Assemble_pairs_assemble_pairs.method
coord = params.Assemble_pairs_assemble_pairs.coord
rc = params.Assemble_pairs_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_assemble_pairs.failed
fasta = params.Assemble_pairs_assemble_pairs.fasta
nproc = params.Assemble_pairs_assemble_pairs.nproc
alpha = params.Assemble_pairs_assemble_pairs.alpha
maxerror = params.Assemble_pairs_assemble_pairs.maxerror
minlen = params.Assemble_pairs_assemble_pairs.minlen
maxlen = params.Assemble_pairs_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_assemble_pairs.scanrev
minident = params.Assemble_pairs_assemble_pairs.minident
evalue = params.Assemble_pairs_assemble_pairs.evalue
maxhits = params.Assemble_pairs_assemble_pairs.maxhits
fill = params.Assemble_pairs_assemble_pairs.fill
aligner = params.Assemble_pairs_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_assemble_pairs.// dbexec
gap = params.Assemble_pairs_assemble_pairs.gap
usearch_version = params.Assemble_pairs_assemble_pairs.usearch_version
assemble_reference = params.Assemble_pairs_assemble_pairs.assemble_reference
head_seqeunce_file = params.Assemble_pairs_assemble_pairs.head_seqeunce_file
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	if(R1.contains("."+head_seqeunce_file)){
		R1 = readArray[0]
		R2 = readArray[1]
	}else{
		R2 = readArray[0]
		R1 = readArray[1]
	}
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi

	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc}  2>&1 | tee out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process Assemble_pairs_parse_log_AP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "reports/$filename"}
input:
 set val(name),file(log_file) from g1_12_logFile1_g1_15
 val mate from g_7_mate_g1_15

output:
 set val(name),file("*.tab")  into g1_15_logFile0_g1_19, g1_15_logFile0_g1_25

script:
field_to_parse = params.Assemble_pairs_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Assemble_pairs_report_assemble_pairs {

input:
 set val(name),file(log_files) from g1_15_logFile0_g1_19
 val matee from g_7_mate_g1_19

output:
 file "*.rmd"  into g1_19_rMarkdown0_g1_25



shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')
	assemble = readArray[0]
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("assemble_length", "Histogram showing the distribution assembled sequence lengths in 
	                            nucleotides for the Align step (top) and Reference step (bottom).")
	figures("assemble_overlap", "Histogram showing the distribution of overlapping nucleotides between 
	                             mate-pairs for the Align step (top) and Reference step (bottom).
	                             Negative values for overlap indicate non-overlapping mate-pairs
	                             with the negative value being the number of gap characters between
	                             the ends of the two mate-pairs.")
	figures("assemble_error", "Histograms showing the distribution of paired-end assembly error 
	                           rates for the Align step (top) and identity to the reference germline 
	                           for the Reference step (bottom).")
	figures("assemble_pvalue", "Histograms showing the distribution of significance scores for 
	                            paired-end assemblies. P-values for the Align mode are shown in the top
	                            panel. E-values from the Reference step's alignment against the 
	                            germline sequences are shown in the bottom panel for both input files
	                            separately.")
	```
	
	```{r, echo=FALSE, warning=FALSE}
	assemble_log <- loadLogTable(file.path(".", "!{assemble}"))
	
	# Subset to align and reference logs
	align_fields <- c("ERROR", "PVALUE")
	ref_fields <- c("REFID", "GAP", "EVALUE1", "EVALUE2", "IDENTITY")
	align_log <- assemble_log[!is.na(assemble_log$ERROR), !(names(assemble_log) %in% ref_fields)]
	ref_log <- assemble_log[!is.na(assemble_log$REFID), !(names(assemble_log) %in% align_fields)]
	
	# Build log set
	assemble_list <- list()
	if (nrow(align_log) > 0) { assemble_list[["Align"]] <- align_log }
	if (nrow(ref_log) > 0) { assemble_list[["Reference"]] <- ref_log }
	plot_titles <- names(assemble_list)
	```
	
	# Paired-End Assembly
	
	Assembly of paired-end reads is performed using the AssemblePairs tool which 
	determines the read overlap in two steps. First, de novo assembly is attempted 
	using an exhaustive approach to identify all possible overlaps between the 
	two reads with alignment error rates and p-values below user-defined thresholds. 
	This method is denoted as the `Align` method in the following figures. 
	Second, those reads failing the first stage of de novo assembly are then 
	mapped to the V-region reference sequences to create a full length sequence, 
	padding with Ns, for any amplicons that have insufficient overlap for 
	de novo assembly. This second stage is referred to as the `Reference` step in the
	figures below.
	
	## Assembled sequence lengths
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="length", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_length")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="overlap", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_overlap")`
	
	## Alignment error rates and significance
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="error", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="pvalue", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```

	`r figures("assemble_pvalue")`

	EOF
	
	open OUT, ">AP_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}
}


process Assemble_pairs_presto_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "reports/$filename"}
input:
 file rmk from g1_19_rMarkdown0_g1_25
 file log_file from g1_15_logFile0_g1_25

output:
 file "*.html" optional true  into g1_25_outputFileHTML00
 file "*csv" optional true  into g1_25_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Filter_Sequence_Quality_filter_seq_quality {

input:
 set val(name),file(reads) from g1_12_reads0_g2_0
 val mate from g_10_mate_g2_0

output:
 set val(name), file("*_${method}-pass.fast*")  into g2_0_reads0_g_29
 set val(name), file("FS_*")  into g2_0_logFile1_g2_5
 set val(name), file("*_${method}-fail.fast*") optional true  into g2_0_reads22
 set val(name),file("out*") optional true  into g2_0_logFile3_g23_0

script:
method = params.Filter_Sequence_Quality_filter_seq_quality.method
nproc = params.Filter_Sequence_Quality_filter_seq_quality.nproc
q = params.Filter_Sequence_Quality_filter_seq_quality.q
n_length = params.Filter_Sequence_Quality_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Quality_filter_seq_quality.n_missing
fasta = params.Filter_Sequence_Quality_filter_seq_quality.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = "-q ${q}"
		n_length = ""
		n_missing = ""
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} >> out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} >> out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} >> out_${R1}_FS.log
	"""
}


}


process maccac_fastq_fasta {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.fasta$/) "reads/$filename"}
input:
 set val(name),  file(reads) from g2_0_reads0_g_29

output:
 set val(name),  file("*.fasta")  into g_29_fastaFile00

script:
	
readArray_reads = reads.toString().split(' ')[0]	

"""

sed -n '1~4s/^@/>/p;2~4p' ${readArray_reads} > ${name}.fasta

 
"""
}


process Filter_Sequence_Quality_parse_log_FS {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*table.tab$/) "reports/$filename"}
input:
 set val(name), file(log_file) from g2_0_logFile1_g2_5
 val mate from g_10_mate_g2_5

output:
 set val(name), file("*table.tab")  into g2_5_logFile0_g2_7, g2_5_logFile0_g2_16

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}


process Filter_Sequence_Quality_report_filter_Seq_Quality {

input:
 val mate from g_10_mate_g2_7
 set val(name), file(log_files) from g2_5_logFile0_g2_7

output:
 file "*.rmd"  into g2_7_rMarkdown0_g2_16


shell:

if(mate=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	quality_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, quality_log_2, titles=plot_titles, sizing="figure")
	```
	
	`r figures("quality")`
		
	EOF
	
	open OUT, ">FSQ_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read")#params$quality_titles
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1],
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, titles=plot_titles[1], sizing="figure")
	```
	
	`r figures("quality")`
	
	EOF
	
	open OUT, ">FSQ_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Filter_Sequence_Quality_presto_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "reports/$filename"}
input:
 file rmk from g2_7_rMarkdown0_g2_16
 file log_file from g2_5_logFile0_g2_16

output:
 file "*.html" optional true  into g2_16_outputFileHTML00
 file "*csv" optional true  into g2_16_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process macca_report_pre_proceesing_pipline_cat_all_file {

input:
 set val(name), file(log_file) from g1_12_logFile3_g23_0
 set val(name), file(log_file) from g2_0_logFile3_g23_0

output:
 set val(name), file("all_out_file.log")  into g23_0_logFile0_g23_9

script:
readArray = log_file.toString()

"""

echo $readArray
cat out* >> all_out_file.log
"""

}


process macca_report_pre_proceesing_pipline_report_pipeline_macca {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.csv$/) "pipeline_statistics/$filename"}
input:
 set val(name), file(log_files) from g23_0_logFile0_g23_9

output:
 file "*.csv"  into g23_9_csvFile00


script:

readArray = log_files.toString().split(' ')
R1 = readArray[0]

"""
#!/usr/bin/env Rscript 


library(prestor)
library(knitr)
library(captioner)


console_log <- loadConsoleLog("$R1")

count_df <- plotConsoleLog(console_log, sizing="figure")

df<-count_df[,c("task","total", "pass", "fail")]



write.csv(df,"pipeline_statistics.csv") 


"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
