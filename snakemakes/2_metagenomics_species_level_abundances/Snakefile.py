## JSB 
## This snakemake runs all steps for species-level metagenomics data, (kraken/bracken), aligns reads to reference genome for PHLAME, and removes human reads for upload to SRA 
# kraken executable is in env snakemake_env.yml

#Reminders:
# Put your own email address in cluster.slurm.json so that you get emails about failures. No one else wants your slurm emails.

''' VARIABLES '''
#USER defined variables (in theory do not need to be touched)
#spls = "allsamples_noprimateprojectormock_SepidermidisATCC12228.csv"
spls = "MetagenomicsSampleNames.csv"
import sys
SCRIPTS_DIRECTORY = "./scripts"
sys.path.insert(0, SCRIPTS_DIRECTORY)
import basic_snakemake_functions as bsf # generic functions to read samples.csv etc.

SCRIPTS_DIRECTORY = "/nfs/tamilab001/c3ddb-scratch-mit_lieberman/scripts"
sys.path.insert(0, SCRIPTS_DIRECTORY)


from read_move_link_samplesCSV import *
from itertools import compress



''' PRE-SNAKEMAKE '''
# define couple of lists from samples.csv
# Format: Path,Sample,ReferenceGenome,ProviderName,Subject
[PATH_ls, SAMPLE_ls, REF_Genome_ls, PROVIDER_ls, CLADEID_ls] = read_samplesCSV(spls)
# Write sample_info.csv for each sample
split_samplesCSV(PATH_ls,SAMPLE_ls,REF_Genome_ls,PROVIDER_ls,CLADEID_ls)


''' SNAKEMAKE '''
rule all:
	input:
		expand("kraken_filtered/{sampleID}/R1_nohuman.fq",sampleID=SAMPLE_ls),
		expand("kraken_filtered/{sampleID}/R2_nohuman.fq",sampleID=SAMPLE_ls),
	
		
rule dedupe:
	# spades uses fwd and rev reads
	input:
		fq1="/nfs/tamilab001/c3ddb-scratch-mit_lieberman/projects/jsb_cuti/all_samples/AllSteps/1-data_processed/{sampleID}/filt1.fq.gz",
		fq2="/nfs/tamilab001/c3ddb-scratch-mit_lieberman/projects/jsb_cuti/all_samples/AllSteps/1-data_processed/{sampleID}/filt2.fq.gz",
	params:
		interm="reads_deduped/{sampleID}/deduped_concatenated.fq",
	output:
		fq1="reads_deduped/{sampleID}/R1.fq",
		fq2="reads_deduped/{sampleID}/R2.fq",
	shell:
		"""
		dedupe.sh threads=10 --overwrite=true --tossbrokenreads in1={input.fq1} in2={input.fq2} out={params.interm}
		reformat.sh in={params.interm} out1={output.fq1} out2={output.fq2}
		"""

rule kraken2:
	#quality assessment based only on fwd 
	input:
		fq1 = rules.dedupe.output.fq1,
		fq2 = rules.dedupe.output.fq2,
	output:
		report="new_krakenreps/{sampleID}_krakenRep.txt",
		seq_results="new_kraken_seqs/{sampleID}.kraken",
	shell:
		"kraken2 --minimum-hit-groups 3 --confidence .1 --paired --threads 10 "
		"--db /orcd/nese/tami/001/databases/krakendb_plus_pf {input.fq1} {input.fq2} "
		"--report {output.report} > {output.seq_results} "

rule extract_reads:
	input:
		rep = rules.kraken2.output.seq_results,
		fq1 = rules.dedupe.output.fq1,
		fq2 = rules.dedupe.output.fq2,
	output:
		fq1o="kraken_filtered/{sampleID}/R1_nohuman.fq",
		fq2o="kraken_filtered/{sampleID}/R2_nohuman.fq",
	shell:
		"""
		python3 extract_kraken_reads.py --fastq-output --exclude -t 9606 -k {input.rep} -s1 {input.fq1} -s2 {input.fq2} -o {output.fq1o} -o2 {output.fq2o}
		"""

rule brackenS:
	input:
		report1 = rules.kraken2.output.report,
		report2 = rules.kraken2_nofilters.output.report,
	output:
		bracken_rep1="Kmers_new/bracken/{sampleID}.S.bracken",
		bracken_rep2="Kmers_new/bracken_nofilter/{sampleID}.S.bracken",
		bracken_rep1_reads="Kmers_new/bracken/{sampleID}.S.reads.bracken",
		bracken_rep2_reads="Kmers_new/bracken_nofilter/{sampleID}.S.reads.bracken",
	shell:
		"""
		bracken -d /orcd/nese/tami/001/databases/krakendb_plus_pf -i {input.report1} -o {output.bracken_rep1} -l S -w {output.bracken_rep1_reads}
		bracken -d /orcd/nese/tami/001/databases/krakendb_plus_pf -i {input.report2} -o {output.bracken_rep2} -l S -w {output.bracken_rep2_reads}
		"""

rule brackenG:
	input:
		report1 = rules.kraken2.output.report,
		report2 = rules.kraken2_nofilters.output.report,
	output:
		bracken_rep1="Kmers_new/bracken/{sampleID}.G.bracken",
		bracken_rep2="Kmers_new/bracken_nofilter/{sampleID}.G.bracken",
		bracken_rep1_reads="Kmers_new/bracken/{sampleID}.G.reads.bracken",
		bracken_rep2_reads="Kmers_new/bracken_nofilter/{sampleID}.G.reads.bracken",
	shell:
		"""
		bracken -d /orcd/nese/tami/001/databases/krakendb_plus_pf -i {input.report1} -o {output.bracken_rep1} -l G -w {output.bracken_rep1_reads}
		bracken -d /orcd/nese/tami/001/databases/krakendb_plus_pf -i {input.report2} -o {output.bracken_rep2} -l G -w {output.bracken_rep2_reads}
		"""

