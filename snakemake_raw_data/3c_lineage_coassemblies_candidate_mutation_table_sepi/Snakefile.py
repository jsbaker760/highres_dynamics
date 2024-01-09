#########################################
# LIEBERMAN LAB SNAKEFILE FOR CASE STEP #
#########################################

# REMINDERS:
# Put your own email address in cluster.slurm.json and in myjob.slurm so that you get emails about failures. 

# Version History:
# # 2022.03.16: JSB: added mods to allow multiple different candidate_mutation_tables.mat to be generated
# # 2020.06.24: cleaned up json and log file names
# # 2019.02.21, Arolyn: changed combine_positions so that it takes into account if a sample is an outgroup
# # 2019.02.16, FMK added python commands to create matlab input files, and tweaked snakemakeslurm.sh to resubmit failed jobs, and some other minor tuning.
# # 2019.02.13, Arolyn: updated json to request more memory for variants2positions rule and candidate_mutation_table rule.
# # 2019.02.12, Arolyn: rule candidate_mutation_table now calling make_candidate_mutation_table_snakemake_indels.m instead of make_candidate_mutation_table_snakemake.m
# # 2018.12, Felix/Arolyn: Original Snakefile from the lab hackathon

# Version Notes:
# This is the matlab version of the case step. Includes indel_counter in candidate_mutation_table.mat.
# Use python version of the case step to generate coverage matrix.

''' PRE-SNAKEMAKE '''

## User defined variables (in theory do not need to be touched)
maxFQ = -30 # purity threshold (from mapping quality) for including position

import sys, subprocess

SCRIPTS_DIRECTORY = "./scripts"
sys.path.insert(0, SCRIPTS_DIRECTORY)
SCRIPTS_DIRECTORY = "/scratch/mit_lieberman/scripts"
sys.path.insert(0, SCRIPTS_DIRECTORY)

spls = "samples_case.csv"
#from read_samplesCSV import read_samples_caseCSV

from read_move_link_samplesCSV import *
import basic_snakemake_functions as bsf # generic functions to read samples.csv etc.
import mergeFQ_runSPAdes as msp # functions for spades file prep rule
from itertools import compress
import sys

# Define couple of lists from samples_case.csv
# # Path = points to snakemake directory used for mapping step (NOT raw data!)
# # Sample = sample name; must match mapping step
# # ReferenceGenome = name of reference genome directory
# # Outgroup = 1 or 0 depending on if sample is or is not an outgroup
# # Clade = the number of the cluster to which the sample belongs
[PATH_ls, SAMPLE_ls, REF_Genome_ls, OUTGROUP_ls, CLADEID_ls] = read_samplesCSV(spls)
# # The unique clades
CLADES_ls = set(CLADEID_ls)


''' SNAKEMAKE '''
rule all:
    input:
        # # Data links only # #
        # expand("data/vcf/{sampleID}_ref_{reference}_outgroup{outgroup}.vcf.gz",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
        # expand("data/qual/{sampleID}_ref_{reference}_outgroup{outgroup}.quals.mat",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
        # expand("data/diversity/{sampleID}_ref_{reference}_outgroup{outgroup}.diversity.mat",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
        # # temp_pos only
        # expand("1-temp_pos/{sampleID}_ref_{reference}_outgroup{outgroup}_positions.mat",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
        # # Everything # #
        expand("2-candidate_mutation_table/clade_{cladeID}_candidate_mutation_table.mat",cladeID=CLADES_ls),
        # expand("1-temp_pos/clade_{cladeID}_string_sampleID_names.txt",cladeID=CLADES_ls),
        # expand("1-temp_pos/clade_{cladeID}_string_outgroup_bool.txt",cladeID=CLADES_ls),
        # expand("1-temp_pos/clade_{cladeID}_string_file_other_p_to_consider.txt",cladeID=CLADES_ls),
        # expand("1-temp_pos/clade_{cladeID}_string_diversity_mat.txt",cladeID=CLADES_ls),
        # expand("1-temp_pos/clade_{cladeID}_string_qual_mat.txt",cladeID=CLADES_ls),



''' HELPER FUNCTIONS for switching from sample-wise rules to reference-wise rules'''

def get_clade_wildcards(cladeID):
    is_clade = [int(i == cladeID) for i in CLADEID_ls]
    sampleID_clade = list(compress(SAMPLE_ls,is_clade))
    reference_clade = list(compress(REF_Genome_ls,is_clade))
    outgroup_clade = list(compress(OUTGROUP_ls,is_clade))
    return sampleID_clade,reference_clade,outgroup_clade
    
def get_sampleID_names(wildcards):  
    sampleID_clade,_,_ = get_clade_wildcards(wildcards.cladeID)
    return sampleID_clade

def get_outgroup_bool(wildcards):  
    _,_,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
    return outgroup_clade

def get_mat_positions_prep(wildcards):
    sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
    mat_positions_prep=expand("1-temp_pos/{sampleID}_ref_{reference}_outgroup{outgroup}_positions.mat",zip,sampleID=sampleID_clade, reference=reference_clade, outgroup=outgroup_clade)
    return mat_positions_prep

def get_diversity_mat(wildcards):
    sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
    diversity_mat = expand("data/diversity/{sampleID}_ref_{reference}_outgroup{outgroup}.diversity.mat",zip,sampleID=sampleID_clade, reference=reference_clade, outgroup=outgroup_clade)
    return diversity_mat   

def get_quals_mat(wildcards):
    sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
    quals_mat = expand("data/qual/{sampleID}_ref_{reference}_outgroup{outgroup}.quals.mat",zip,sampleID=sampleID_clade, reference=reference_clade, outgroup=outgroup_clade)
    return quals_mat 

rule build_data_links:
    output:
        vcf_links = expand("data/vcf/{sampleID}_ref_{reference}_outgroup{outgroup}.vcf.gz",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
        qual_mat_links = expand("data/qual/{sampleID}_ref_{reference}_outgroup{outgroup}.quals.mat",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
        div_mat_links = expand("data/diversity/{sampleID}_ref_{reference}_outgroup{outgroup}.diversity.mat",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
    run:
        import subprocess
        subprocess.run( "rm -fr data/ " ,shell=True) # clean it up prior run
        subprocess.run( "mkdir -p data/vcf/ data/qual/ data/diversity/ " ,shell=True)
        for i in range(len(SAMPLE_ls)):
            subprocess.run( "ln -fs -T " + PATH_ls[i] + "/6-diversity/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "_aligned.sorted.strain.variant.diversity.mat data/diversity/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "_outgroup" + OUTGROUP_ls[i] + ".diversity.mat" ,shell=True)
            subprocess.run( "ln -fs -T " + PATH_ls[i] + "/5-quals/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "_aligned.sorted.strain.variant.quals.mat data/qual/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "_outgroup" + OUTGROUP_ls[i] + ".quals.mat" ,shell=True)
            subprocess.run( "ln -fs -T " + PATH_ls[i] + "/4-vcf/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "_aligned.sorted.strain.variant.vcf.gz data/vcf/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "_outgroup" + OUTGROUP_ls[i] + ".vcf.gz " ,shell=True)

rule variants2positions:
    input:
        variants = "data/vcf/{sampleID}_ref_{reference}_outgroup{outgroup}.vcf.gz",
    params:
        REF_GENOME_DIRECTORY = "/scratch/mit_lieberman/reference_genomes/sepidermidis_assemblies_72/{reference}/",
        outgroup_tag = 0, # boolean (0==ingroup or 1==outgroup)
    output:
        mat_positions = "1-temp_pos/{sampleID}_ref_{reference}_outgroup{outgroup}_positions.mat",
    shell:
        """
        module add mit/matlab/2015b; 
        matlab -r "path('{SCRIPTS_DIRECTORY}',path); generate_positions_single_sample_snakemake( '{input.variants}', '{output.mat_positions}', {maxFQ}, '{params.REF_GENOME_DIRECTORY}',{params.outgroup_tag}  ) "  
        """

rule string_sampleID_names:
    params:
        sampleID_names = get_sampleID_names,
    output:
        string_sampleID_names = "1-temp_pos/clade_{cladeID}_string_sampleID_names.txt",
    run:
        with open( output.string_sampleID_names, "w") as f: 
            print(*params.sampleID_names, sep=" ", file=f)

# build input for candidate_mutation_table
rule string_outgroup_bool:
    params:
      outgroup_bool = get_outgroup_bool,
    output:
      string_outgroup_bool = "1-temp_pos/clade_{cladeID}_string_outgroup_bool.txt",
    run:
      with open( output.string_outgroup_bool ,"w") as f: 
        print(*params.outgroup_bool, sep=" ", file=f)

rule combine_positions_prep:
    input:
        mat_positions = get_mat_positions_prep,
    output:
        string_input_p_positions = "1-temp_pos/clade_{cladeID}_string_file_other_p_to_consider.txt",
    run:
        with open( output.string_input_p_positions ,"w") as f: 
            print(*input.mat_positions, sep=" ", file=f)

rule combine_positions:
    input:
        string_input_p_positions = "1-temp_pos/clade_{cladeID}_string_file_other_p_to_consider.txt", #"1-temp_pos/clade_{cladeID}_string_file_other_p_to_consider.txt",
        string_outgroup_bool = "1-temp_pos/clade_{cladeID}_string_outgroup_bool.txt",       #"1-temp_pos/clade_{cladeID}_string_outgroup_bool.txt",
    params:
        REF_GENOME_DIRECTORY = "/scratch/mit_lieberman/reference_genomes/sepidermidis_assemblies_72/clade_{cladeID}/",
        file_other_p_to_consider = "add_positions/{cladeID}_ther_positions.mat",
    output:
        mat_positions = "1-temp_pos/clade_{cladeID}_allpositions.mat",
    shell:
        """
        module add mit/matlab/2015b; 
        matlab -r "path('{SCRIPTS_DIRECTORY}',path); combine_positions_snakemake_out( '{input.string_input_p_positions}', '{params.file_other_p_to_consider}', '{output.mat_positions}', '{input.string_outgroup_bool}', '{params.REF_GENOME_DIRECTORY}', {maxFQ} )" 
        """

# build input for candidate_mutation_table
rule string_diversity_mat:
    input:
        diversity_mat = get_diversity_mat,
    output:
        string_diversity_mat = "1-temp_pos/clade_{cladeID}_string_diversity_mat.txt",
    run:
        with open( output.string_diversity_mat ,"w") as f: 
            print(*input.diversity_mat, sep=" ", file=f)

# build input for candidate_mutation_table
rule string_quals_mat:
    input:
      quals_mat = get_quals_mat,
    output:
      string_qual_mat = "1-temp_pos/clade_{cladeID}_string_qual_mat.txt",
    run:
      with open( output.string_qual_mat ,"w") as f: 
        print(*input.quals_mat, sep=" ", file=f)

# build input for candidate_mutation_table


rule candidate_mutation_table:
    input:
        mat_positions = rules.combine_positions.output.mat_positions, # "1-temp_pos/clade_{cladeID}_allpositions.mat",
        string_diversity_mat = rules.string_diversity_mat.output.string_diversity_mat, # "1-temp_pos/clade_{cladeID}_string_diversity_mat.txt",
        string_qual_mat = rules.string_quals_mat.output.string_qual_mat, # "1-temp_pos/clade_{cladeID}_string_qual_mat.txt",
        string_sampleID_names = rules.string_sampleID_names.output.string_sampleID_names, # "1-temp_pos/clade_{cladeID}_string_sampleID_names.txt",
        string_outgroup_bool = rules.string_outgroup_bool.output.string_outgroup_bool, # "1-temp_pos/clade_{cladeID}_string_outgroup_bool.txt",
    output:
        candidate_mutation_table = "2-candidate_mutation_table/clade_{cladeID}_candidate_mutation_table.mat",
    shell:
        """
        module add mit/matlab/2015b; 
        matlab -r "path('{SCRIPTS_DIRECTORY}',path); make_candidate_mutation_table_snakemake_indels( '{input.mat_positions}', '{input.string_sampleID_names}', '{input.string_outgroup_bool}', '{input.string_qual_mat}', '{input.string_diversity_mat}', '{output.candidate_mutation_table}' )" 
        """


