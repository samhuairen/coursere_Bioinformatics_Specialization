import  re
def reverse_complement(dna_seq):
    """reverse complement of dna string"""
    base_dict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    complement_seq = ''
    for base in dna_seq:
        complement_seq += base_dict[base]
    return complement_seq[::-1]


def transcript(dna_seq):
    trans_dict = {
        'A': 'A',
        'T': 'U',
        'C': 'C',
        'G': 'G'
    }
    transcript_seq = ''
    for base in dna_seq:
        transcript_seq += trans_dict[base]
    return transcript_seq

    
def translation(transcript):
    """translate RNA to protein"""
    codon_table = {
        'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGU': 'S',
        'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I',
        'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAU': 'H',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAU': 'D',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'UAA': '*', 'UAC': 'Y', 'UAG': '*', 'UAU': 'Y',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UGA': '*', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C',
        'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F'
    }
    
    protein_seq = ''
    for i in range(0, len(transcript) - len(transcript) % 3, 3):
        codon = transcript[i:i + 3]
        if codon in codon_table.keys():
            amino_acid = codon_table[codon]
            protein_seq += amino_acid
    return protein_seq


def finding_substrings_peptide_encoding_genome(dna_seq, sub_peptide):
    out_strings = []
    frames = [0, 1, 2]
    len_substring = len(sub_peptide)
    for frame in frames:
        dna_frame = dna_seq[frame::]
        rna_frame = transcript(dna_frame)
        protein_frame = translation(rna_frame)
        finding_peptide = re.finditer(sub_peptide, protein_frame)
        for item in finding_peptide:
            start = item.start()
            out_strings.append(dna_frame[start*3: start*3+len_substring*3])
    for frame in frames:
        res_dna_frame = reverse_complement(dna_seq)[frame::]
        res_rna_frame = transcript(res_dna_frame)
        res_protein_frame = translation(res_rna_frame)
        finding_peptide2 = re.finditer(sub_peptide, res_protein_frame)
        for item in finding_peptide2:
            start = item.start()
            out_strings.append(reverse_complement(res_dna_frame[start*3: start*3+len_substring*3]))
    return out_strings
    
    
        
with open('Bacillus_brevis2.txt', 'r') as f:
    dna_seq = f.readlines()
    dna_seq2 = "".join([seq.strip() for seq in dna_seq])
    peptide_seq = 'VKLFPWFNQY'
    resu = finding_substrings_peptide_encoding_genome(dna_seq2, peptide_seq)


mRNA1 = 'CCGAGGACCGAAAUCAAC'
mRNA2 = 'CCCAGUACCGAAAUUAAC'
mRNA3 = 'CCCAGGACUGAGAUCAAU'
mRNA4 = 'CCCCGUACGGAGAUGAAA'

frames = [0, 1, 2]
for frame in frames:
    mRNA = mRNA4[frame:]
    protein = translation(mRNA)
    print(protein)
    


