


def individual_fasta(wildcards):
    individuals = get_individuals()
    fasta_list = expand("{fasta_dir}/{individual}.fasta", bam_dir = config["consensus_fasta_dir"], individual = individuals.keys())
    return fasta_list 

rule consensus:
    input:
        individual_fasta

rule consensus_individual:
    input:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config["bam_dir"], extension = config["final_bam_extension"]),
        expand("{bam_dir}/{{individual}}{extension}.bam.bai", bam_dir = config["bam_dir"], extension = config["final_bam_extension"])
    output:
        expand("{fasta_dir}/{{individual}}.fasta", fasta_dir = config["consensus_fasta_dir"])
    log: expand("{logs}/{{individual}}/consensus.log", logs=config["log_dir"])
    threads: 2
    shell:
        "samtools consensus -@ {threads} -o {output} {input[0]} > {log} 2>&1"