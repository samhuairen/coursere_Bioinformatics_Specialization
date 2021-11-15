def compute_n_length(contig_length_list, ratio):
    """return the N50 length of assembly contigs"""
    # check the contig length list,
    for contig in contig_length_list:
        if not isinstance(contig, int):
            break
            return 'The contig list should be a int list'
    total_length = sum(contig_length_list)
    sorted_list = sorted(contig_length_list, reverse=True)
    comulative_length = 0
    for contig in sorted_list:
        comulative_length += contig
        if comulative_length >= total_length * ratio:
            break
    return contig


def compute_ng50(contig_list, genome_length):
    sorted_contig_list = sorted(contig_list, reverse=True)
    cumulative_length = 0
    for item in sorted_contig_list:
        cumulative_length += item
        if cumulative_length > genome_length * 0.5:
            break
    return item
    

def compute_genome_length(genome_fasta):
    fasta_dict = {}
    with open(genome_fasta, 'r') as f:
        for line in f:
            if line.startswith(">"):
                genome_id = line[1:]
                fasta_dict[genome_id] = ""
            else:
                fasta_dict[genome_id] += line.strip()
    total_length = 0
    for key, value in fasta_dict.items():
        total_length += len(value)
    return total_length

    
def summarize_assembly(assembly_file):
    """out put N50 lenght, the number of long contigs which exceed 1000bp, and the total lenght of contigs"""
    contig_names = []
    contig_length_list = []
    coverage = []
    summary = open(assembly_file, 'r')
    for line in summary:
        if line.startswith("#"):
            continue
        line_text = line.strip().split("\t")
        contig_names.append(line_text[0])
        contig_length_list.append(int(line_text[1]))
        coverage.append(float(line_text[2]))
    return contig_length_list


contig_list1 = summarize_assembly('contig_stats_25.tabular')
contig_list2 = summarize_assembly('contig_stats_55.tabular')
contig_list3 = summarize_assembly('contig_stats_85.tabular')

n50_length_25 = compute_n50_length(contig_list1)
n50_length_55 = compute_n50_length(contig_list2)
n50_length_85 = compute_n50_length(contig_list3)


number_of_long_contigs1 = len([contig for contig in contig_list1 if contig >= 1000])
number_of_long_contigs2 = len([contig for contig in contig_list2 if contig >= 1000])
number_of_long_contigs3 = len([contig for contig in contig_list3 if contig >= 1000])

total_length_long1 = sum([contig for contig in contig_list1 if contig >= 1000])
total_length_long2 = sum([contig for contig in contig_list2 if contig >= 1000])
total_length_long3 = sum([contig for contig in contig_list3 if contig >= 1000])

genome_length = compute_genome_length('staph_genome.fasta')

cons1 = [20, 20, 30, 30, 60, 60, 80, 100, 200]