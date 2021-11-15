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


def linear_spectrum(pepetides):
    """ return the linear spectrum of Peptide
    """
    len_peptides = len(pepetides)
    prefix_mass = [0]
    amino_mass = [aminoAcidMass[amino] for amino in pepetides]
    for i in range(0, len_peptides):
        prefix_mass.append(sum(amino_mass[0:i + 1]))
    linear_spectrum = [0]
    for i in range(0, len(prefix_mass)):
        for j in range(i + 1, len(prefix_mass)):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(linear_spectrum)


def cyclo_spectrum(peptides):
    """ return the cylic theoretical peptides spectrum """
    prefix_pepetides_without_last_pepetide = linear_spectrum(peptides[:-1])
    total_mass = linear_spectrum(peptides)[-1]
    difference = [total_mass - item for item in prefix_pepetides_without_last_pepetide]
    cylic_spectrum = prefix_pepetides_without_last_pepetide + difference
    return sorted(cylic_spectrum)


p1 = 'TLAM'
p2 = 'MAIT'
p3 = 'TMLA'
p4 = 'TMIA'
p5 = 'MTAI'
p6 = 'TAIM'

ts = [0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]
ts1 = cyclo_spectrum(p1)
ts2 = cyclo_spectrum(p2)
ts3 = cyclo_spectrum(p3)
ts4 = cyclo_spectrum(p4)
ts5 = cyclo_spectrum(p5)
ts6 = cyclo_spectrum(p6)
