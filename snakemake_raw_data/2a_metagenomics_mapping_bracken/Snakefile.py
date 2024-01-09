## JSB all steps master Snakemake 2020
## Assembly, Mapping, & Kraken

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
import mergeFQ_runSPAdes as msp # functions for spades file prep rule
from itertools import compress



''' PRE-SNAKEMAKE '''
# define couple of lists from samples.csv
# Format: Path,Sample,ReferenceGenome,ProviderName,Subject
[PATH_ls, SAMPLE_ls, REF_Genome_ls, PROVIDER_ls, CLADEID_ls] = read_samplesCSV(spls)
# Write sample_info.csv for each sample
split_samplesCSV(PATH_ls,SAMPLE_ls,REF_Genome_ls,PROVIDER_ls,CLADEID_ls)


# set(list_patient) provides only unique clade IDs


''' SNAKEMAKE '''
rule all:
	input:
		# expand("data/{sampleID}/R1.fq.gz",sampleID=set(SAMPLE_ls)),
		# expand("data/{sampleID}/R2.fq.gz",sampleID=set(SAMPLE_ls)),
		# expand("Mapping/3-bowtie2/{sampleID}_ref_{references}_aligned.sorted.bam",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
		# expand("Mapping/3-bowtie2/{sampleID}_ref_{references}_aligned.sam",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
		# expand("Mapping/3-bowtie2/{sampleIDs}_ref_{references}_cleanup_done.txt",zip, sampleIDs=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
		# expand("Mapping/4-vcf/{sampleID}_ref_{references}_aligned.sorted.strain.variant.vcf.gz",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
		# expand("Mapping/6-diversity/{sampleID}_ref_{references}_aligned.sorted.strain.variant.diversity.mat",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
		# expand("Mapping/5-quals/{sampleID}_ref_{references}_aligned.sorted.strain.variant.quals.mat",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
		# expand("data/{sampleID}/deduped_1.fq.gz",sampleID=SAMPLE_ls),
		# expand("data/{sampleID}/deduped_2.fq.gz",sampleID=SAMPLE_ls),
		expand("Kmers/bracken/{sampleID}.S.bracken",sampleID=SAMPLE_ls),
		expand("Kmers/bracken/{sampleID}.G.bracken",sampleID=SAMPLE_ls),

#PART 0: prepare filtered, clean FASTQ samples

rule make_data_links:
	# NOTE: All raw data needs to be names fastq.gz. No fq! The links will be names fq though.
	input:
		sample_info_csv="data/{sampleID}/sample_info.csv",
	output:
		# Recommend using symbolic links to your likely many different input files
		fq1o="data/{sampleID}/R1.fq.gz",
		fq2o="data/{sampleID}/R2.fq.gz",
	run:
		# get stuff out of mini csv file
		with open(input.sample_info_csv,'r') as f:
			this_sample_info = f.readline() # only one line to read
		this_sample_info = this_sample_info.strip('#').split(',')
		path = this_sample_info[0] # remember python indexing starts at 0
		paths = path.split(' ')
		sample = this_sample_info[1]
		provider = this_sample_info[3]
		# make links
		bsf.cp_append_files_beta(paths, sample)

		
rule dedupe:
	# spades uses fwd and rev reads
	input:
		fq1=rules.make_data_links.output.fq1o,
		fq2=rules.make_data_links.output.fq2o,
	output:
		fqio1="data/{sampleID}/deduped1.fq.gz",
		fqio2="data/{sampleID}/deduped2.fq.gz",
	shell:
		"dedupe.sh -Xmx100g in1={input.fq1} in2={input.fq2} out1={output.fqio1} out2={output.fqio2}"


rule kraken2:
	#quality assessment based only on fwd 
	input:
		fq1 = rules.dedupe.output.fq1o,
		fq2 = rules.dedupe.output.fq2o,
	output:
		report="Kmers/kraken2/{sampleID}_krakenRep.txt",
		seq_results="0-tmp/{sampleID}_krakSeq.txt.gz",
	shell:
		"kraken2 --minimum-hit-groups 3 --confidence .1 --paired --gzip-compressed --threads 20 "
		"--db /nfs/tamilab001/c3ddb-scratch-mit_lieberman/tools/databases/kraken2 {input.fq1} {input.fq2} "
		"--report {output.report} |gzip > {output.seq_results} "

rule brackenS:
	input:
		report = rules.kraken2.output.report,
	output:
		bracken_rep="Kmers/bracken/{sampleID}.S.bracken",
	shell:
		"bracken -d /nfs/tamilab001/c3ddb-scratch-mit_lieberman/tools/databases/kraken2 -i {input.report} -o {output.bracken_rep} -l S"

rule brackenG:
	input:
		report = rules.kraken2.output.report,
	output:
		bracken_rep="Kmers/bracken/{sampleID}.G.bracken",
	shell:
		"bracken -d /nfs/tamilab001/c3ddb-scratch-mit_lieberman/tools/databases/kraken2 -i {input.report} -o {output.bracken_rep} -l G"

