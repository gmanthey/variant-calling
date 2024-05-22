import re

configfile: "config.yml"
localrules: filtered, index_bam, index_reference, merge_trimmed, index_raw_vcf, retain_list

def get_chromosomes():
    chromosomes = {}
    with open(config["chromosome_file"], 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if len(line) == 1:
                chromosomes[line[0]] = line
            else:
                chromosomes[line[0]] = line[1:]

    return chromosomes

def get_individuals():  
    individuals = dict()
    with open(config["individual_file"], 'r') as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue
            line = line.split()
            if line[0] not in individuals:
                individuals[line[0]] = []
            individuals[line[0]].append(line[1])

    return individuals

rule filtered:
    input:
        expand('{vcf_dir}/genome.IF-GF-MM2-RM.vcf.gz', vcf_dir = config["vcf_dir"]),
        expand('{vcf_dir}/genome.IF-GF-MM2-RM.vcf.gz.csi', vcf_dir = config["vcf_dir"])

include: "subworkflows/call.snakefile"

include: "subworkflows/filter.snakefile"