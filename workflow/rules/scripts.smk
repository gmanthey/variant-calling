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
            if len(line) == 1:
                individuals[line[0]].append(line[0])
            else:
                individuals[line[0]].append(line[1])

    return individuals

def get_bam_file(individual):
    if os.path.exists(f"{config['ro_bam_dir']}/{individual}{config['final_bam_extension']}.bam"):
        bam_prefix = f"{config['ro_bam_dir']}/{individual}{config['final_bam_extension']}"
    else:
        bam_prefix = f"{config['bam_dir']}/{individual}{config['final_bam_extension']}"

    bam_file = bam_prefix + ".bam"

    if os.path.exists(bam_prefix + ".bam") and os.path.exists(bam_prefix + ".bai"):
            return [bam_file, bam_prefix + ".bai"]

    return [bam_file, bam_prefix + ".bam.bai"]

def get_summary_files(wildcards):
    individuals = get_individuals()

    temp_file_ids = []
    for individual in individuals:
        fastq_files_per_individual = individuals[individual]

        run_ids = {}

        for fastq_file in fastq_files_per_individual:
            if fastq_file.endswith('gz'):
                file = gzip.open(os.path.join(config['raw_fastq_dir'], fastq_file), 'rt')
            else:
                file = open(os.path.join(config['raw_fastq_dir'], fastq_file), 'r')
            first_read = file.readline().strip().split()[0]
            file.close()
            first_read = first_read[1:].split('/')[0] # BGI names (and possibly others) have a /1 or /2 at the end
            if first_read not in run_ids:
                run_ids[first_read] = []
            run_ids[first_read].append(os.path.basename(fastq_file))

        for read_id, run_fastq_files in run_ids.items():
            if len(run_fastq_files) != 2:
                raise ValueError(f"Read id {read_id} does not have 2 files associate with it. (Found {', '.join(run_fastq_files)})")
    
        temp_file_ids.extend(['_'.join(sorted(run_fastq_files)) for run_fastq_files in run_ids.values()])

    return expand("results/fastq_stats/{run_id}.summary", run_id = temp_file_ids)