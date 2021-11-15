
aminoAcidMass = ['57', '71', '87', '97', '99', '101', '103', '113', '114', '115', '128', '129',
                 '131', '137', '147', '156', '163', '186']

def linear_spectrum(aalist):
    """ return the linear spectrum of Peptide
    """
    len_peptides = len(aalist)
    prefix_mass = [0]
    for i in range(0, len_peptides):
        prefix_mass.append(sum(aalist[0:i + 1]))
    linear_spectrum = [0]
    for i in range(0, len(prefix_mass)):
        for j in range(i + 1, len(prefix_mass)):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(linear_spectrum)


def cyclo_spectrum(aalist):
    """ return the cylic theoretical peptides spectrum """
    prefix_pepetides_without_last_pepetide = linear_spectrum(aalist[:-1])
    total_mass = linear_spectrum(aalist)[-1]
    difference = [total_mass - item for item in prefix_pepetides_without_last_pepetide]
    cylic_spectrum = prefix_pepetides_without_last_pepetide + difference
    return sorted(cylic_spectrum)


def mass(peptides):
    """ return the sum of mass of pepetides"""
    total_mass = 0
    for aa in peptides.split("-"):
        total_mass += int(aa)
    return total_mass


def parent_mass(spectrum):
    """ return the largest mass value of the specturm"""
    return max(spectrum)


def expand(candidate_pepetides):
    """ all return peptides which length is one longer than the candidate pepetides"""
    expand_candidate = []
    for peptide in candidate_pepetides:
        if peptide == '':
            for amino in aminoAcidMass:
                expand_candidate.append(amino)
        else:
            for amino in aminoAcidMass:
                expand_candidate.append(peptide+"-"+amino)
    return expand_candidate


def consistent(peptides, spectrum):
    aa_list = [int(item) for item in peptides.split("-")]
    linear_spectrum_aa = linear_spectrum(aa_list)
    linear_spectrum_dict = {item: linear_spectrum_aa.count(item) for item in set(linear_spectrum_aa)}
    spectrum_dict = {item: spectrum.count(item) for item in set(spectrum)}
    for key, value in linear_spectrum_dict.items():
        if value > spectrum_dict.get(key, 0):
            return False
    return True


def consistentCycloSpectrum(peptides, spectrum):
    aalist = [int(aa) for aa in peptides.split("-")]
    if cyclo_spectrum(aalist) == spectrum:
        return True
    else:
        return False


def cyclo_peptide_sequencing(specturm):
    """ return the pepetides which theoretical specturm equals input spectrum"""
    candidate_pepetides = {''}
    final_pepetides = []
    while candidate_pepetides:
        candidate_pepetides = expand(candidate_pepetides)
        for peptide in list(candidate_pepetides):
            if mass(peptide) == parent_mass(specturm):
                if consistentCycloSpectrum(peptide, specturm) and peptide not in final_pepetides:
                    final_pepetides.append(peptide)
                candidate_pepetides.remove(peptide)
            elif not consistent(peptide, specturm):
                candidate_pepetides.remove(peptide)
    return final_pepetides


aa_string = '0 71 97 97 97 99 99 99 103 113 129 131 170 184 194 196 198 198 200 216 226 228 260 269 283 287 293 295 297 313 325 329 357 357 368 382 384 386 392 394 424 442 454 454 460 465 481 483 485 491 513 523 551 553 557 562 573 578 582 584 612 622 644 650 652 654 670 675 681 681 693 711 741 743 749 751 753 767 778 778 806 810 822 838 840 842 848 852 866 875 907 909 919 935 937 937 939 941 951 965 1004 1006 1022 1032 1036 1036 1036 1038 1038 1038 1064 1135'

spectrum_example = [int(item) for item in aa_string.split(" ")]
resu = cyclo_peptide_sequencing(spectrum_example)









