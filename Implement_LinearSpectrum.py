Alphabet = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

aminoAcidMass = {
    'G': 57,
    'A': 71,
    'S': 87,
    'P': 97,
    'V': 99,
    'T': 101,
    'C': 103,
    'I': 113,
    'L': 113,
    'N': 114,
    'D': 115,
    'K': 128,
    'Q': 128,
    'E': 129,
    'M': 131,
    'H': 137,
    'F': 147,
    'R': 156,
    'Y': 163,
    'W': 186
}


def linear_spectrum(pepetides, alphabet, amino_acid_mass):
    """ return the linear spectrum of Peptide
    """
    len_peptides = len(pepetides)
    prefix_mass = [0]
    amino_mass = [aminoAcidMass[amino] for amino in pepetides]
    for i in range(0, len_peptides):
        prefix_mass.append(sum(amino_mass[0:i+1]))
    linear_spectrum = [0]
    for i in range(0, len(prefix_mass)):
        for j in range(i+1, len(prefix_mass)):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(linear_spectrum)






