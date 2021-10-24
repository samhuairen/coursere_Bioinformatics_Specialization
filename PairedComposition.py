def paired_omposition(text, k, d):
    paired_composition = []
    for i in range(0, len(text) - 2*k):
        pattern1 = text[i:i+k]
        pattern2 = text[i+k+d: i+k+d+k]
        paired = (pattern1+"|"+pattern2)
        paired_composition.append(paired)
        paired_composition.sort()
    return paired_composition


def paired_deGruign(gapped_patterns):
    from collections import defaultdict
    adjecency_dict = defaultdict(list)
    for pattern in gapped_patterns:
        part1_pattern = pattern.split("|")[0]
        part2_pattern = pattern.split("|")[1]
        prefix = part1_pattern[:-1]+"|"+part2_pattern[:-1]
        suffex = part1_pattern[1:]+"|"+part2_pattern[1:]
        adjecency_dict[prefix].append(suffex)
    return adjecency_dict




gapped_patterns = ['AAT|CCA','ATG|CAT','ATG|GAT','CAT|GGA','CCA|GGG',
                   'GCC|TGG','GGA|GTT','GGG|TGT','TAA|GCC','TGC|ATG','TGG|ATG']


text2 = 'TAATGCCATGGGATGTT'
paired = paired_omposition(text2, 3,2)
results = ""
for item in paired:
    item = item.split("|")
    results += "("+item[0]+'|'+item[1]+")"+" "
results


