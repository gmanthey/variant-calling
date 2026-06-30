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

def get_individuals(include_outgroup = False):
    individuals = dict()
    with open(config["individual_file"], 'r') as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue
            line = line.split()
            if not include_outgroup and line[0] in config.get('outgroup_individuals', []):
                continue
            if line[0] not in individuals:
                individuals[line[0]] = []
            if len(line) > 1:
                individuals[line[0]].append(line[1])

    return individuals

def get_bam_file(individual):
    ro_bam_dirs = config.get("ro_bam_dir", [])
    ro_bam_dirs = ro_bam_dirs if isinstance(ro_bam_dirs, list)  else [ro_bam_dirs]
    for ro_bam_dir in ro_bam_dirs:
        if os.path.exists(f"{ro_bam_dir}/{individual}.bam"):
            bam_prefix = f"{ro_bam_dir}/{individual}"
            break
    else:
        bam_prefix = f"{config['bam_dir']}/{individual}"

    bam_file = bam_prefix + ".bam"

    if os.path.exists(bam_prefix + ".bam") and os.path.exists(bam_prefix + ".bai"):
        return [bam_file, bam_prefix + ".bai"]

    return [bam_file, bam_prefix + ".bam.bai"]
