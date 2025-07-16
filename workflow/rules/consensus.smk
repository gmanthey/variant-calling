from os import path


def individual_fasta(wildcards):
    individuals = get_individuals()
    fasta_list = expand("{fasta_dir}/individuals/{individual}.fasta", fasta_dir = config["consensus_fasta_dir"], individual = individuals.keys())
    return fasta_list 

def genome_sequences():
    sequences = []
    with open(config["genome"], 'r') as f:
        for line in f:
            if line.startswith('>'):
                sequences.append(line.strip().split()[0][1:])

    return sequences

rule consensus:
    input:
        expand("{fasta_dir}/combined/{sequence}.fasta", fasta_dir = config["consensus_fasta_dir"], sequence = genome_sequences())

rule consensus_individual:
    input:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config["bam_dir"], extension = config["final_bam_extension"]),
        expand("{bam_dir}/{{individual}}{extension}.bam.bai", bam_dir = config["bam_dir"], extension = config["final_bam_extension"])
    output:
        expand("{fasta_dir}/individuals/{{individual}}.fasta", fasta_dir = config["consensus_fasta_dir"])
    log: expand("{logs}/{{individual}}/consensus.log", logs=config["log_dir"])
    threads: 2
    shell:
        "samtools consensus -@ {threads} -o {output} {input[0]} > {log} 2>&1"



rule consensus_combined:
    input:
        individual_fasta
    output:
        expand("{fasta_dir}/combined/{sequence}.fasta", fasta_dir = config["consensus_fasta_dir"], sequence = genome_sequences())
    run:
        output_files = {path.splitext(path.basename(o))[0]: open(o, 'w') for o in output}
        curr_output = None
        for file in input:
            ind_id = path.splitext(path.basename(file))[0]
            with open(file, 'r') as in_file:
                for line in in_file:
                    if line.startswith('>'):
                        seqname = line.strip().split()[0][1:]
                        curr_output = output_files[seqname]
                        curr_output.write(">" + ind_id + "\n")
                    else:
                        curr_output.write(line)

        for output_file in output_files.values():
            output_file.close()