import os
import pandas as pd

def raw_vcf_individual(wildcards):
    ro_vcf_dirs = [config["ro_ind_vcf_dir"]] if isinstance(config.get("ro_ind_vcf_dir", []), list)  else config.get("ro_ind_vcf_dir", [])
    for ro_vcf_dir in ro_vcf_dirs:
        vcf_ro = f"{ro_vcf_dir}/{wildcards.individual}.raw.vcf.gz"
        if os.path.exists(vcf_ro):
            vcf_base = vcf_ro
            break
    else:
        vcf_base = f"{config['vcf_dir']}/individuals/{wildcards.individual}.raw.vcf.gz"
    return [vcf_base, vcf_base + ".csi"]

rule filter_ind_quality:
    input:
        raw_vcf_individual
    output:
        temp(expand("{vcf_dir}/individuals/{{individual}}.QUAL.vcf.gz", vcf_dir = config["vcf_dir"]))
    log: expand("{logs}/{{individual}}/filter_qual.log", logs=config["log_dir"])
    run:
        shell(f"bcftools view --threads {{threads}} -e 'QUAL < 20' -Oz -o {{output}} {input[0]} > {{log}} 2>&1")

def individual_vcfs(wildcards):
    individuals = get_individuals(include_outgroup=False)

    vcf_list = expand("{vcf_dir}/individuals/{individual}.QUAL.vcf.gz", vcf_dir = config["vcf_dir"], individual = individuals.keys())
    return vcf_list

def individual_vcf_indices(wildcards):
    individuals = get_individuals(include_outgroup=False)

    vcf_list = expand("{vcf_dir}/individuals/{individual}.QUAL.vcf.gz.csi", vcf_dir = config["vcf_dir"], individual = individuals.keys())
    return vcf_list

def chromosome_filtered_vcfs(wildcards):
    chromosomes = get_chromosomes()
    return expand("{vcf_dir}/chromosomes/{chromosome}.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"], chromosome = chromosomes.keys())

rule merge_vcf_into_chromosomes:
    input:
        vcf=individual_vcfs,
        index=individual_vcf_indices
    output:
        expand("{vcf_dir}/chromosomes/{{chromosome}}.allsites.vcf.gz", vcf_dir = config["vcf_dir"])
    threads: 2
    log: expand("{logs}/{{chromosome}}/merge.log", logs=config["log_dir"])
    run:
        chromosomes = ','.join(get_chromosomes()[wildcards.chromosome])
        shell(f"bcftools merge --threads {{threads}} -r {chromosomes} -Oz -o {{output}} {{input.vcf}} > {{log}} 2>&1")

rule filter_qual_depth_missing_rpbz:
    input:
        expand("{vcf_dir}/chromosomes/{{chromosome}}.allsites.vcf.gz", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/chromosomes/{{chromosome}}.IF.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/{{chromosome}}/filter_qual_depth_missing.log", logs=config["log_dir"])
    run:
        sampn = int(shell("bcftools query -l {input} | wc -l", read=True))
        avgdp = shell("bcftools query -f '%DP\n' {input} | datamash median 1 | datamash round 1", read=True)
        if avgdp == '':
            shell(f"cp {input} {output}")
            logout = open(log[0], 'w')
            logout.write('Empty vcf file, forwarding.')
            logout.close()
            return
        
        avgdp = int(avgdp)
        dphi = 2 * avgdp

        shell(f"""bcftools view --types snps --threads {{threads}} -e "INFO/DP > {dphi} || INFO/DP < {sampn} || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o {output} {input} > {log} 2>&1""")

rule filter_gf:
    input:
        expand("{vcf_dir}/chromosomes/{{chromosome}}.IF.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/chromosomes/{{chromosome}}.IF.vcf.gz.csi", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/chromosomes/{{chromosome}}.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/{{chromosome}}/filter_lcs.log", logs=config["log_dir"])
    shell:
        """bcftools +setGT -Oz -o {output} {input[0]} -- -t q -i "FMT/DP < {config[min_depth]}" -n "./." > {log} 2>&1"""

rule merge_chromosome_vcf:
    input:
        chromosome_filtered_vcfs
    output:
        expand("{vcf}/genome.IF-GF.vcf.gz", vcf = config["vcf_dir"])
    log: expand("{logs}/merge_all.log", logs=config["log_dir"])
    shell:
        "bcftools concat -Oz -o {output} {input} > {log} 2>&1"

rule sample_stats:
    input:
        expand("{vcf_dir}/genome.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/genome.IF-GF.vcf.gz.csi", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/sample.stats", vcf_dir = config["vcf_dir"])
    shell:
        """echo -e "ID\tnREF\tnALT\tnHET\tnTs\tnTv\tavgDP\tSingletons\tMissing_Sites\tproportion_Missing" > {output}
        bcftools stats --threads {threads} -S- {input[0]} | grep 'PSC' | grep -v '#' | tr ' ' '_' | awk '{{OFS="\t"}}{{print $3,$4,$5,$6,$7,$8,$10,$11,$14,$14/($4+$5+$6+$14)}}' >> {output}
        """

rule retain_list:
    input:
        expand("{vcf_dir}/sample.stats", vcf_dir = config["vcf_dir"])
    output:
        "results/retain.list"
    run:
        sample_stats = pd.read_csv(input[0], sep='\t')
        individuals = sample_stats[sample_stats['proportion_Missing'] < config['max_missingness_individual']]
        individuals.to_csv(output[0], index=False, header=False, columns=['ID'])
    

rule filter_genotype_missing_ind:
    input:
        expand("{vcf_dir}/genome.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/genome.IF-GF.vcf.gz.csi", vcf_dir = config["vcf_dir"]),
        expand("results/retain.list", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/genome.IF-GF-MM2.vcf.gz", vcf_dir = config["vcf_dir"])
    log: 
        expand("{logs}/filter_genotype_missing_samples.log", logs=config["log_dir"]),
        expand("{logs}/filter_genotype_missing_min.log", logs=config["log_dir"])
    shell:
        """bcftools view --threads {threads} --samples-file {input[2]} --force-samples -Ou {input[0]} 2> {log[0]} | bcftools view --min-ac 1 --threads {threads} -i 'F_MISSING<{config[max_missingness_site]}' -Oz -o {output} > {log[1]} 2>&1""" 

rule join_outgroup:
    input:
        full_vcf=expand("{vcf_dir}/genome.IF-GF-MM2.vcf.gz", vcf_dir = config["vcf_dir"]),
        full_vcf_index=expand("{vcf_dir}/genome.IF-GF-MM2.vcf.gz.csi", vcf_dir = config["vcf_dir"]),
        outgroup_vcfs=expand("{vcf_dir}/outgroup/{individual}.raw.vcf.gz", vcf_dir = config["vcf_dir"], individual=config.get('outgroup_individuals', [])),
        outgroup_vcf_indices=expand("{vcf_dir}/outgroup/{individual}.raw.vcf.gz.csi", vcf_dir = config["vcf_dir"], individual=config.get('outgroup_individuals', [])),
    output:
        expand("{vcf_dir}/genome.IF-GF-MM2-OG.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/join_outgroup.log", logs=config["log_dir"])
    shell:
        """bcftools merge --force-single --threads {threads} -Oz -o {output} {input.full_vcf} {input.outgroup_vcfs} > {log} 2>&1"""

rule sam_index_reference:
    input:
        config["genome"]
    output:
        "results/genome/genome.fai"
    log: expand("{logs}/sam_index_reference.log", logs=config["log_dir"])
    shell:
        "samtools faidx {input} --fai-idx {output} > {log} 2>&1"

rule filter_repeats:
    input:
        expand("{vcf_dir}/genome.IF-GF-MM2-OG.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/genome.IF-GF-MM2-OG.vcf.gz.csi", vcf_dir = config["vcf_dir"]),
        "results/genome/genome.fai"
    output:
        expand("{vcf_dir}/genome.IF-GF-MM2-OG-RM.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/filter_repeats.log", logs=config["log_dir"])
    shell:
        """bcftools view --threads {threads} -T <(bedtools complement -i {config[repeat_bed]} -g {input[2]}) -Oz -o {output} {input[0]} > {log} 2>&1"""

