rule call:
    input:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config['bam_dir'], extension = config['final_bam_extension']),
        config["genome"],
        expand("{bam_dir}/{{individual}}{extension}.bam.bai", bam_dir = config['bam_dir'], extension = config['final_bam_extension']) if exists(expand("{bam_dir}/{{individual}}{extension}.bam.bai", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])) else expand("{bam_dir}/{{individual}}{extension}.bai", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    output:
        expand("{vcf_dir}/{{individual}}.raw.vcf.gz", vcf_dir = config["vcf_dir"])
    threads: 2
    log: 
        expand("{logs}/{{individual}}/mpielup.log", logs=config["log_dir"]),
        expand("{logs}/{{individual}}/call.log", logs=config["log_dir"])
    shell:
        "bcftools mpileup --threads {threads} -q 20 -Q 20 -C 50 -Ou -s {wildcards.individual} -f {input[1]} {input[0]} -a \"AD,ADF,ADR,DP,SP\" 2> {log[0]} | bcftools call --threads {threads} --ploidy 2 -m -Oz -o {output} > {log[1]} 2>&1"


include: "index.snakefile"