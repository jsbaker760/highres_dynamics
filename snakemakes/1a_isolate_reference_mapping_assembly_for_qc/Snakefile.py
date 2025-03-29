## JSB
## This snakemake takes each samples and aligns them 
## To the reference genomes Pacnes_C1 and SepidermidisATCC12228

''' VARIABLES '''
# each row is one sample/reference pair. each isolate will have two rows
# one for Sepidermidis_ATCC12228 and Pacnes_C1

spls = "isolate_samples.csv"

import sys
SCRIPTS_DIRECTORY = "./scripts"
sys.path.insert(0, SCRIPTS_DIRECTORY)
SCRIPTS_DIRECTORY = "/scratch/mit_lieberman/scripts"
sys.path.insert(0, SCRIPTS_DIRECTORY)


from read_move_link_samplesCSV import *
import basic_snakemake_functions as bsf # generic functions to read samples.csv etc.
import mergeFQ_runSPAdes as msp # functions for spades file prep rule
from itertools import compress

''' PRE-SNAKEMAKE '''
# define couple of lists from samples.csv
# Format: Path,Sample,ReferenceGenome,ProviderName,Subject
[PATH_ls, SAMPLE_ls, REF_Genome_ls, PROVIDER_ls, CLADEID_ls] = read_samplesCSV(spls)
# Write sample_info.csv for each sample
# split_samplesCSV(PATH_ls,SAMPLE_ls,REF_Genome_ls,PROVIDER_ls,CLADEID_ls)


''' SNAKEMAKE '''
rule all:
	input:
		# Just through spades
		expand("Assembly/isolates/spades/{sampleID}/contigs.fasta",sampleID=SAMPLE_ls),
		# Just through prokka
		expand("Assembly/isolates/annotation/{sampleID}/prokka_out.faa",sampleID=SAMPLE_ls),
		# information about each assembly 
		expand("Assembly/isolates/assemblystats/{sampleID}/sumStats_assembly_annotation.tsv",sampleID=SAMPLE_ls),
		#  all mapping steps Through all steps (part 2) # #		
		#through mapping
		expand("Mapping/4-vcf/{sampleID}_ref_{references}_aligned.sorted.strain.variant.vcf.gz",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
		expand("Mapping/6-diversity/{sampleID}_ref_{references}_aligned.sorted.strain.variant.diversity.mat",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
		expand("Mapping/5-quals/{sampleID}_ref_{references}_aligned.sorted.strain.variant.quals.mat",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),


#PART 0: prepare filtered, clean FASTQ samples

rule make_data_links:
	# NOTE: All raw data needs to be names fastq.gz. No fq! The links will be names fq though.
	input:
		sample_info_csv="data/{sampleID}/sample_info.csv",
	output:
		# Recommend using symbolic links to your likely many different input files
		fq1="data/{sampleID}/R1.fq.gz",
		fq2="data/{sampleID}/R2.fq.gz",
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
		if len(paths)>1:
			cp_append_files(paths, sample, provider)
		else:
			makelink(path, sample, provider)

rule cutadapt:
	input:
		# Recommend using symbolic links to your likely many different input files
		fq1 = rules.make_data_links.output.fq1,
		fq2 = rules.make_data_links.output.fq2,
	output:
		fq1o="0-tmp/{sampleID}_R1_trim.fq.gz",
		fq2o="0-tmp/{sampleID}_R2_trim.fq.gz",
	log:
		log="logs/cutadapt_{sampleID}.txt",
	conda:
		"envs/cutadapt.yaml"
	shell:
		"cutadapt -a CTGTCTCTTAT -o {output.fq1o} {input.fq1} 1> {log};"
		"cutadapt -a CTGTCTCTTAT -o {output.fq2o} {input.fq2} 1>> {log};"

rule sickle2050:
	input:
		fq1o = rules.cutadapt.output.fq1o,
		fq2o = rules.cutadapt.output.fq2o,
	output:
		fq1o="1-data_processed/{sampleID}/filt1.fq.gz",
		fq2o="1-data_processed/{sampleID}/filt2.fq.gz",
		fqSo="1-data_processed/{sampleID}/filt_sgls.fq.gz",
	log:
		log="logs/sickle2050_{sampleID}.txt",
	shell:
		"module add c3ddb/sickle/1.33 ;"
		"sickle pe -g -f {input.fq1o} -r {input.fq2o} -t sanger -o {output.fq1o} -p {output.fq2o} -s {output.fqSo} -q 20 -l 50 -x -n"

# PART 1: map to reference genome

rule refGenome_index: 
    input:
        fasta="/refgenomedir/{reference}/genome.fasta"
    params:
        "/refgenomedir/{reference}/genome_bowtie2",
    output:
        bowtie2idx="/refgenomedir/{reference}/genome_bowtie2.1.bt2"
    conda:
        "envs/bowtie2.yaml"
    shell:
        "bowtie2-build -q {input.fasta} {params} "

rule bowtie2:
    input:
        fq1=ancient("1-data_processed/{sampleID}/filt1.fq.gz"),
        fq2=ancient("1-data_processed/{sampleID}/filt2.fq.gz"),
        bowtie2idx=ancient(rules.refGenome_index.output.bowtie2idx) # put here, so rule botie2 only executed after rule refGenome_index done
    params:
        refGenome="/refgenomedir/{reference}/genome_bowtie2",
    output:
        samA="Mapping/3-bowtie2/{sampleID}_ref_{reference}_aligned.sam",
    log:
        log="logs/bowtie2_{sampleID}_ref_{reference}.txt",
    conda:
        "envs/bowtie2.yaml"
    shell:
        # 8 threads coded into json
        "bowtie2 --very-sensitive --threads 8 -X 2000 --no-mixed --dovetail -x {params.refGenome} -1 {input.fq1} -2 {input.fq2} -S {output.samA} 2> {log} "
        
rule sam2bam:
    input:
        samA=ancient(rules.bowtie2.output.samA),
    output:
        bamA="Mapping/3-bowtie2/{sampleID}_ref_{reference}_aligned.sorted.bam",
    conda:
        "envs/samtools15_bcftools12.yaml"
    shell:
        # 8 threads coded into json
        " samtools view -bS {input.samA} | samtools sort - -o {output.bamA} ;"
        " samtools index {output.bamA} ;"
        " rm {input.samA} ;"

rule mpileup2vcf:
    input:
        bamA=ancient(rules.sam2bam.output.bamA),
        ref="/refgenomedir/{reference}/genome.fasta"
    output:
        pileup="Mapping/4-vcf/{sampleID}_ref_{reference}_aligned.sorted.pileup",
        variants="Mapping/4-vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.vcf.gz",
        vcf_strain="Mapping/4-vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.vcf.gz",
    params:
        vcf_raw="Mapping/4-vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.gz",
    conda:
        "envs/samtools15_bcftools12.yaml"
    shell:
        " samtools faidx {input.ref} ; "
        " samtools mpileup -q30 -x -s -O -d3000 -f {input.ref} {input.bamA} > {output.pileup} ;" 
        " samtools mpileup -q30 -t SP -d3000 -vf {input.ref} {input.bamA} > {params.vcf_raw} ;"
        " bcftools call -c -Oz -o {output.vcf_strain} {params.vcf_raw} ;"
        " bcftools view -Oz -v snps -q .75 {output.vcf_strain} > {output.variants} ;"
        " tabix -p vcf {output.variants} ;"
        " rm {params.vcf_raw}"

# strain.vcf ==> vcf_to_quals.m ==> quals.mat
rule vcf2quals:
    input:
        vcf_strain = ancient(rules.mpileup2vcf.output.vcf_strain),
    params:
        refGenomeDir="/refgenomedir/{reference}/" 
    output:
        file_quals = "Mapping/5-quals/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.quals.mat", 
    shell:
        """
        module add mit/matlab/2015b; 
        matlab -r "path('{SCRIPTS_DIRECTORY}',path); vcf_to_quals_snakemake( '{input.vcf_strain}', '{output.file_quals}', '{params.refGenomeDir}' )"
        """

# strain.pileup ==> pileup_to_diversity.m ==> diversity.mat
rule pileup2diversity_matrix:
    input:
        pileup = ancient(rules.mpileup2vcf.output.pileup),
    params:
        refGenomeDir="/refgenomedir/{reference}/", 
    output:
        file_diversity = "Mapping/6-diversity/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.diversity.mat",
        file_coverage = "Mapping/6-diversity/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.coverage.mat",
    shell:
        """
        module add mit/matlab/2015b; 
        matlab -r "path('{SCRIPTS_DIRECTORY}',path); pileup_to_diversity_matrix_snakemake( '{input.pileup}', '{output.file_diversity}', '{output.file_coverage}', '{params.refGenomeDir}' )" ;
        rm {input.pileup};
        """
# Part 2: Assembling each isolate into contigs for filtering with CDHIT

rule spades:
	input:
		fastq1=rules.sickle2050.output.fq1o,
		fastq2=rules.sickle2050.output.fq2o,
	params:
		threads=16,
		outdir="Assembly/isolates/spades/{sampleID}"
	conda:
		"envs/spades.yaml"
	output:
		fasta="Assembly/isolates/spades/{sampleID}/contigs.fasta", # produced by spades''
	shell:
		"spades.py --isolate -k	21,33,55,77,99,127	--phred-offset	33	-t {params.threads} -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir}"

rule prokka:
	# prokka annotations; not currently using --proteins option to specify trusted gbk file
	input:
		rules.spades.output.fasta,
	params:
		outdir="Assembly/isolates/annotation/{sampleID}",
	threads: 16
	output:
		txt="Assembly/isolates/annotation/{sampleID}/prokka_out.txt",
		faa="Assembly/isolates/annotation/{sampleID}/prokka_out.faa",
	conda:
		"envs/prokka.yaml"
	shell:
		"prokka --compliant --force --cpus {threads} --outdir {params.outdir} --prefix prokka_out {input} ; conda deactivate"

# samples we were able to assemble
rule sumstats:
	input:
		infile="samples.csv",
		fastq=rules.concatenate_spades_input.output.fastq1,
		contig=rules.spades.output.fasta,
		assembly=rules.prokka.output.txt,
	output:
		"Assembly/isolates/assemblystats/{cladeID}/sumStats_assembly_annotation.tsv",
	conda: 
		"envs/py_for_snakemake.yaml",
	shell:
		"python3 scripts/getSumStats_SPAdes_prokka_v2.py -s {input.infile} -p {wildcards.cladeID} -f {input.fastq} -c {input.contig} -a {input.assembly} -o {output} ; conda deactivate " 

rule merge_sumstats:
	input:
		expand("Assembly/isolates/assemblystats/{cladeID}/sumStats_assembly_annotation.tsv",cladeID=CLADES_ls),
	output:
		"Assembly/5-assemblystats/all_clades_sumStats_assembly_annotation.tsv",
	shell:
		"cat {input} |sed -n '3~2!p' > {output}"

## cd-hit is generated by concatenating all prokkas together with references (supplement) 
## into the file combined_proteins.fa and running the command:
## cd-hit -c .90 -i combined_proteins.fa -o CDhit_outp90

