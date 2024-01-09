## JSB all steps master Snakemake 2020
## Assembly, Mapping, & Kraken

#Reminders:
# Put your own email address in cluster.slurm.json so that you get emails about failures. No one else wants your slurm emails.

''' VARIABLES '''
#USER defined variables (in theory do not need to be touched)
#spls = "allsamples_noprimateprojectormock_SepidermidisATCC12228.csv"
spls = "cacnes_cluster_assembly_mapping.csv"
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


CLADES_ls = set(CLADEID_ls)


# set(list_patient) provides only unique clade IDs


''' SNAKEMAKE '''
rule all:
    input:
        # expand("data/{sampleID}/sample_info.csv",sampleID=SAMPLE_ls)
        # Only data links #
        #expand("data/{sampleID}/R1.fq.gz",sampleID=SAMPLE_ls),
        #expand("data/{sampleID}/R2.fq.gz",sampleID=SAMPLE_ls),
        # Through all steps (part 1) #
        #"0-tmp/SubjectsPassKraken.txt",
        #expand("1-kraken2/{sampleID}_krakenRep.txt",sampleID=SAMPLE_ls),
        #expand("0-tmp/{sampleID}_krakSeq.txt.gz",sampleID=SAMPLE_ls),
        #expand("3-bracken/{sampleID}.bracken",sampleID=SAMPLE_ls),
        #expand("0-tmp/{sampleID}_R1_trim.fq.gz",sampleID=SAMPLE_ls),
        #expand("0-tmp/{sampleID}_R2_trim.fq.gz",sampleID=SAMPLE_ls),
        # expand("1-data_processed/{sampleID}/filt1.fq.gz",sampleID=SAMPLE_ls),
        # expand("1-data_processed/{sampleID}/filt2.fq.gz",sampleID=SAMPLE_ls),
        # expand("1-data_processed/{sampleID}/filt_sgls.fq.gz",sampleID=SAMPLE_ls),
        # Just through spades
        # expand("Assembly/spadesk33/{sampleID}/contigs.fasta",sampleID=SAMPLE_ls),
        # expand("Assembly/spadesk55/{sampleID}/contigs.fasta",sampleID=SAMPLE_ls),
        # Just through prokka
        # expand("Assembly/annotation_filtered/{sampleID}/prokka_out.faa",sampleID=SAMPLE_ls),
        # Through all steps (part 2) #
        # "Assembly/assemblystats_filtered/sumStats_assembly_annotation.tsv",
        # "Assembly/orthologinfo_filtered/annotation_orthologs.tsv",
        #"Assembly/orthologinfo/annotation_orthologs.tsv"
        #expand("Assembly/sepi_orthologinfo/annotation_orthologs.tsv"),
        #expand("Assembly/cacnes_orthologinfo/annotation_orthologs.tsv"),
        #through mapping
        # expand("Mapping/3-bowtie2/{sampleID}_ref_{references}_aligned.sam",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
        # expand("Mapping/4-vcf/{sampleID}_ref_{references}_aligned.sorted.strain.variant.vcf.gz",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
        # expand("Mapping/6-diversity/{sampleID}_ref_{references}_aligned.sorted.strain.variant.diversity.mat",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
        # expand("Mapping/5-quals/{sampleID}_ref_{references}_aligned.sorted.strain.variant.quals.mat",zip, sampleID=SAMPLE_ls, references=REF_Genome_ls),#,zip, sampleID=subset[0], references=subset[1]),
        #expand("1-data_processed/{sampleID}/filt1.fq.gz",sampleID=SAMPLE_ls),
        #expand("1-data_processed/{sampleID}/filt2.fq.gz",sampleID=SAMPLE_ls),
        #expand("1-data_processed/{sampleID}/filt_sgls.fq.gz",sampleID=SAMPLE_ls),
        #expand("1-data_processed/{sampleID}/k14.cat.fq.gz.msh",sampleID=SAMPLE_ls),
        #expand("1-data_processed/{sampleID}/k19.cat.fq.gz.msh",sampleID=SAMPLE_ls),
        #expand("3-bracken/{sampleID}.S.bracken",sampleID=SAMPLE_ls),
        #expand("3-bracken/{sampleID}.S8.bracken",sampleID=SAMPLE_ls),
        #expand("3-bracken/{sampleID}.S9.bracken",sampleID=SAMPLE_ls),
        #expand("3-bracken/{sampleID}.S.uncracked.bracken",sampleID=SAMPLE_ls),
        #expand("3-bracken/{sampleID}.G.uncracked.bracken",sampleID=SAMPLE_ls),
        #expand("3-bracken/{sampleID}.S9.uncracked.bracken",sampleID=SAMPLE_ls),
        #expand("SKA/{sampleID}15comparisons.tsv",sampleID=SAMPLE_ls),
        # expand("Assembly/clade_assemblies/cacnes/{reference}/prokka_out.faa",reference=set(REF_Genome_ls)),
        "Assembly/sepi_orthologinfo_filtered/annotation_orthologs.tsv",
        "Assembly/cacnes_orthologinfo_filtered/annotation_orthologs.tsv"


# returns sample names of samples for the clade specified by arguement
def get_cladefq(wildcards):
      SID=get_clade_samples(wildcards.reference)
      fq=expand("1-data_processed/{sampleID}/filt{n}.fq.gz",sampleID=SID,n=[1,2])
      return fq
      
#returns sample names of samples for the clade specified by arguement
def get_clade_samples(reference):
      is_clade = [int(i == reference) for i in REF_Genome_ls]
      samples_clade = list(compress(SAMPLE_ls,is_clade))
      return samples_clade
      
rule concatenate_spades_input:
    #calls functions from mergeFQ_runSPAdes.py which concatenate num_reads from fastq files
    #where num_reads is the number of reads in the smallest fastq files. Thhis maximizes the number of reads
    #used in the assembly without biasing ancestor calls towards high-coverage samples.
    input:
        clade_fq=get_cladefq,
    output:
        fastq1="0-tmp/in1_spades_{reference}.fq.gz", # produced by py script, and used by spades for assemebly
        fastq2="0-tmp/in2_spades_{reference}.fq.gz", # 
    run:
        clade_ls = get_clade_samples(wildcards.reference)
        file_names = msp.build_sample_file_list(clade_ls)
        outfileLs = msp.merge_fq(file_names,wildcards.reference)

rule spades_clade_cacnes:
    input:
      fastq1=rules.concatenate_spades_input.output.fastq1,
      fastq2=rules.concatenate_spades_input.output.fastq2,
    params:
      outdir="Assembly/clade_assemblies/cacnes/{reference}"
    conda:
      "envs/spades.yaml"
    threads: 16
    output:
      fasta="Assembly/clade_assemblies/cacnes/{reference}/contigs.fasta", # produced by spades''
    shell:
      "spades.py -m 500 -k 21,33,55,77 --phred-offset 33 --careful -t {threads} -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir}"

rule refGenome_index: 
    input:
        fasta=ancient(rules.spades_clade_cacnes.output.fasta),
    params:
        "Assembly/clade_assemblies/cacnes/{reference}/genome_bowtie2",
    output:
        bowtie2idx="Assembly/clade_assemblies/cacnes/{reference}/genome_bowtie2.1.bt2"
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
        refGenome="Assembly/clade_assemblies/cacnes/{reference}/genome_bowtie2",
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
        ref="Assembly/clade_assemblies/cacnes/{reference}/contigs.fasta",
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
        refGenomeDir="/scratch/mit_lieberman/projects/jsb_cuti/all_samples/AllSteps/Assembly/clade_assemblies/cacnes/{reference}/", 
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
        refGenomeDir="/scratch/mit_lieberman/projects/jsb_cuti/all_samples/AllSteps/Assembly/clade_assemblies/cacnes/{reference}/", 
    output:
        file_diversity = "Mapping/6-diversity/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.diversity.mat",
        file_coverage = "Mapping/6-diversity/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.coverage.mat",
    shell:
        """
        module add mit/matlab/2015b; 
        matlab -r "path('{SCRIPTS_DIRECTORY}',path); pileup_to_diversity_matrix_snakemake( '{input.pileup}', '{output.file_diversity}', '{output.file_coverage}', '{params.refGenomeDir}' )" ;
        rm {input.pileup};
        """

