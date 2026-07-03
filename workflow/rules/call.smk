import os

def individual_bam(wildcards):
    return get_bam_file(wildcards.individual)

rule individual_file:
    input:
        config["individual_file"]
    output:
        temp("results/{group}/{individual}.txt")
    run:
        with open(output[0], "w") as f:
            f.write(wildcards.individual)

rule call:
    input:
        individual_bam,
        config["genome"],
    output:
        temp(f"{config['vcf_dir']}/individuals/{{individual}}.raw.unnamed.vcf.gz")
    threads: 2
    log: 
        expand("{logs}/{{individual}}/mpileup.log", logs=config["log_dir"]),
        expand("{logs}/{{individual}}/call.log", logs=config["log_dir"])
    shell:
        "bcftools mpileup --threads {threads} -q 20 -Q 20 -Ou -s {wildcards.individual} --ignore-RG -f {input[2]} {input[0]} -a \"AD,ADF,ADR,DP,SP\" 2> {log[0]} | bcftools call --threads {threads} -a \"GQ\" --ploidy {config[ploidy]} -m -Oz -o {output} > {log[1]} 2>&1"

rule rename_individual:
    input:
        expand("{vcf_dir}/{{group}}/{{individual}}.raw.unnamed.vcf.gz", vcf_dir=config["vcf_dir"]),
        "results/{group}/{individual}.txt"
    output:
        expand("{vcf_dir}/{{group}}/{{individual}}.raw.vcf.gz", vcf_dir=config["vcf_dir"])
    log:
        expand("{logs}/{{individual}}/{{group}}_rename.log", logs=config["log_dir"])
    threads: 2
    shell:
        "bcftools reheader --threads {threads} -s {input[1]} {input[0]} -o {output} > {log} 2>&1"

rule outgroup_variant_file:
    input:
        expand("{vcf_dir}/genome.IF-GF-MM2.vcf.gz", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/outgroup/regions.txt", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/outgroup_variant_file.log", logs=config["log_dir"])
    shell:
        """bcftools query -f '%CHROM\t%POS\n' {input} > {output} 2> {log}"""

rule call_outgroup:
    input:
        config["genome"],
        expand("{vcf_dir}/outgroup/regions.txt", vcf_dir = config["vcf_dir"]),
        individual_bam,
    output:
        temp(f"{config['vcf_dir']}/outgroup/{{individual}}.raw.unnamed.vcf.gz")
    threads: 2
    log: 
        expand("{logs}/{{individual}}/mpileup.log", logs=config["log_dir"]),
        expand("{logs}/{{individual}}/call.log", logs=config["log_dir"])
    shell:
        "bcftools mpileup --threads {threads} -Q 20 -R {input[1]} -Ou -s {wildcards.individual} --ignore-RG -f {input[0]} {input[2]} -a \"AD,ADF,ADR,DP,SP\" 2> {log[0]} | bcftools call --threads {threads} -a \"GQ\" --ploidy {config[ploidy]} -m -Oz -o {output} > {log[1]} 2>&1"

