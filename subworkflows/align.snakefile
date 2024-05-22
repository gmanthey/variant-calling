
def individual_bams(wildcards):
    individuals = get_individuals()
    bam_list = expand("{bam_dir}/{individual}{extension}.bam", bam_dir = config["bam_dir"], individual = individuals.keys())
    bam_index_list = expand("{bam_dir}/{individual}{extension}.bai", bam_dir = config["bam_dir"], individual = individuals.keys())
    return bam_list + bam_index_list

def get_reads(wildcards):
    individual_reads = get_individuals()
    all_reads = sum(individual_reads.values(), start=[])
    reads = sorted([read for read in all_reads if read.startswith(wildcards.id)])

    return expand("{fastq_dir}/{read}", fastq_dir = config["raw_fastq_dir"], read = reads)

rule bams:
    input:
        individual_bams


rule index_reference:
    input: 
        config["genome"]
    output:
        "output/genome.pac"
    params:
        "output/genome"
    log: expand("{logs}/index_reference.log", logs=config["log_dir"])
    shell:
        "bwa-mem2 index -p {params[0]} {input[0]} > {log} 2>&1"



rule trim_paired_reads:
    input:
        get_reads
    output:
        temp(expand("{fastq_trimmed_dir}/{{id}}_R{read}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], read=[1, 2]))
    log: expand("{logs}/{{id}}/trim.log", logs=config["log_dir"])
    threads: 4
    resources:
        mem_mb = 40000
    params:
        memory = "40G"
    shell:
        "bbduk.sh t={threads} -Xmx{params.memory} overwrite=true in={input[0]} in2={input[1]} out={output[0]} out2={output[1]} ref={config[adapters]} ktrim=r k=23 mink=25 hdist=1 tpe tbo > {log[0]} 2>&1"

def individual_trimmed(wildcards):
    individuals = get_individuals()

    individual_reads = individuals[wildcards.individual]

    individual_ids = {}

    for read in individual_reads:
        read_id_pos = re.search(r"R[12]", read).start()
        individual_id = read[:read_id_pos]
        if individual_id not in individual_ids:
            individual_ids[individual_id] = []
        individual_ids[individual_id].append(read)
    
    for id, reads in individual_ids.items():
        if len(reads) > 2:
            raise ValueError(f"Individual {id} has more than 2 reads ({reads})")

    return sorted(expand("{fastq_trimmed_dir}/{id}_R{{read}}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], id = individual_ids.keys()))

rule merge_trimmed:
    input:
        individual_trimmed
    output:
        expand("{fastq_trimmed_dir}/{{individual}}_R{{read}}.trimmed.all.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"])
    shell:
        "cat {input} > {output}"

rule align:
    input:
        "output/genome.pac",
        expand("{fastq_trimmed_dir}/{{individual}}_R{read}.trimmed.all.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], read=[1, 2])
    params:
        genome_idx = "output/genome",
        memory = "8G"
    output:
        temp(expand("{bam_dir}/{{individual}}.sorted.bam", bam_dir = config["bam_dir"]))
    threads: 8
    resources:
        mem_mb = 100000
    log: 
        expand("{logs}/{{individual}}/bwa.log", logs=config["log_dir"]),
        expand("{logs}/{{individual}}/fixmate.log", logs=config["log_dir"]),
        expand("{logs}/{{individual}}/sort.log", logs=config["log_dir"])
    shell:
        "bwa-mem2 mem -t {threads} {params.genome_idx} {input[1]} {input[2]} 2> {log[0]} | samtools fixmate -@ {threads} -m - - 2> {log[1]} | samtools sort -@ {threads} -m {params.memory} -o {output} > {log[2]} 2>&1"

rule markdup:
    input:
        expand("{bam_dir}/{{individual}}.sorted.bam", bam_dir = config['bam_dir'])
    output:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    log: expand("{logs}/{{individual}}/markdup.log", logs=config["log_dir"])
    threads: 4
    shell:
        "samtools markdup -@ {threads} -d {config[optical_dup_dist]} -S -r {input} {output} > {log} 2>&1"

rule index_bam:
    input:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    output:
        expand("{bam_dir}/{{individual}}{extension}.bai", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    log: expand("{logs}/{{individual}}/index.log", logs=config["log_dir"])
    shell:
        "samtools index -@ {threads} {input} {output} > {log} 2>&1"
