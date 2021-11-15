
aminoAcidMass = ['57', '71', '87', '97', '99', '101', '103', '113', '114', '115', '128', '129',
                 '131', '137', '147', '156', '163', '186']


def mass(peptides):
    """ return the sum of mass of pepetides"""
    total_mass = 0
    for aa in peptides.split("-"):
        total_mass += int(aa)
    return total_mass


def parent_mass(spectrum):
    """ return the largest mass value of the specturm"""
    return max(spectrum)


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


def scoring(peptides, experimental_speactrum):
    aalist = [int(item) for item in peptides.split("-")]
    cyclic_pepetides = cyclo_spectrum(aalist)
    score = []
    all_masses = set(cyclic_pepetides + experimental_speactrum)
    a_times = {}
    b_times = {}
    for item in all_masses:
        a_times[item] = cyclic_pepetides.count(item)
        b_times[item] = experimental_speactrum.count(item)
    for key in a_times.keys():
        if a_times[key] == b_times[key]:
            score.append(a_times[key])
        if a_times[key] > b_times[key]:
            if b_times[key] == 0:
                continue
            else:
                score.append(b_times[key])
        if a_times[key] < b_times[key]:
            if a_times[key] == 0:
                continue
            else:
                score.append(a_times[key])
    return sum(score)


def expand(candidate_pepetides):
    """ all return peptides which length is one longer than the candidate pepetides"""
    expand_candidate = set()
    for peptide in candidate_pepetides:
        if peptide == '':
            for amino in aminoAcidMass:
                expand_candidate.add(amino)
        else:
            for amino in aminoAcidMass:
                expand_candidate.add(peptide+"-"+amino)
    return expand_candidate


def linear_score(pepetides, spectrum):
    """return the score between a linear specturm and experimental spectrum"""
    aalist = [int(item) for item in pepetides.split("-")]
    linear_spect = linear_spectrum(aalist)
    score = []
    all_masses = set(linear_spect + spectrum)
    a_times = {}
    b_times = {}
    for item in all_masses:
        a_times[item] = linear_spect.count(item)
        b_times[item] = spectrum.count(item)
    for key in a_times.keys():
        if a_times[key] == b_times[key]:
            score.append(a_times[key])
        if a_times[key] > b_times[key]:
            if b_times[key] == 0:
                continue
            else:
                score.append(b_times[key])
        if a_times[key] < b_times[key]:
            if a_times[key] == 0:
                continue
            else:
                score.append(a_times[key])
    
    return sum(score)


def trim(lead_board, spectrum, N):
    """ trim the lead_board pepetides according to the linear_score between each pepetides
    in lead_broad and the spectrum, left N number of pepetides in lead_board"""
    l = len(lead_board)
    pepetide_score_dict = {}
    for item in lead_board:
        linear_scores = linear_score(item, spectrum)
        pepetide_score_dict[item] = linear_scores
    peptide_score_dict_sorted = sorted(pepetide_score_dict.items(), key=lambda x: x[1], reverse=True)
    lead_board_sorted = [item[0] for item in peptide_score_dict_sorted]
    linear_scores_sorted = [item[1] for item in peptide_score_dict_sorted]
    for j in range(N, l):
        if linear_scores_sorted[j] < linear_scores_sorted[N - 1]:
            return lead_board_sorted[:j]
    return lead_board


def leaderboardCyclopeptideSequencing(spectrum, N):
    leader_board = {''}
    leader_peptide = '0'
    while leader_board:
        leader_board = expand(leader_board)
        for peptide in list(leader_board):
            if mass(peptide) == parent_mass(spectrum):
                if scoring(peptide, spectrum) > scoring(leader_peptide, spectrum):
                    leader_peptide = peptide
            elif mass(peptide) > parent_mass(spectrum):
                leader_board.remove(peptide)
        leader_board = trim(leader_board, spectrum, N)
    return leader_peptide


aa_strings = '0 71 71 87 101 113 113 114 114 114 115 128 128 128 129 129 131 131 147 147 147 156 163 172 184 184 185 186 186 215 218 241 242 243 245 258 261 261 262 269 275 275 278 285 287 291 294 298 301 312 315 332 340 342 346 358 362 365 371 372 376 389 392 392 402 403 404 408 422 426 441 444 445 448 455 463 471 471 473 475 479 490 493 493 500 507 517 518 526 535 536 549 550 555 558 562 573 576 584 586 589 600 604 606 607 620 621 627 627 640 646 647 649 655 663 663 670 678 683 686 693 704 713 714 720 720 733 733 734 735 736 742 754 756 768 775 784 791 791 794 796 798 813 825 827 833 834 842 847 856 862 864 865 867 867 867 882 882 885 889 898 911 919 928 938 942 953 954 955 957 961 978 981 981 982 984 989 993 996 998 999 1003 1011 1026 1028 1048 1057 1066 1066 1067 1071 1075 1082 1085 1085 1096 1097 1102 1112 1117 1125 1127 1128 1134 1140 1141 1156 1168 1173 1179 1180 1185 1186 1189 1194 1196 1198 1199 1204 1204 1213 1229 1230 1243 1254 1259 1269 1269 1269 1274 1282 1297 1300 1301 1308 1309 1313 1317 1318 1320 1325 1326 1327 1333 1335 1340 1346 1357 1360 1374 1382 1383 1388 1402 1425 1428 1429 1437 1438 1440 1447 1447 1448 1449 1453 1454 1460 1461 1471 1474 1474 1482 1483 1488 1489 1511 1519 1530 1538 1560 1561 1566 1567 1575 1575 1578 1588 1589 1595 1596 1600 1601 1602 1602 1609 1611 1612 1620 1621 1624 1647 1661 1666 1667 1675 1689 1692 1703 1709 1714 1716 1722 1723 1724 1729 1731 1732 1736 1740 1741 1748 1749 1752 1767 1775 1780 1780 1780 1790 1795 1806 1819 1820 1836 1845 1845 1850 1851 1853 1855 1860 1863 1864 1869 1870 1876 1881 1893 1908 1909 1915 1921 1922 1924 1932 1937 1947 1952 1953 1964 1964 1967 1974 1978 1982 1983 1983 1992 2001 2021 2023 2038 2046 2050 2051 2053 2056 2060 2065 2067 2068 2071 2088 2092 2094 2095 2096 2107 2111 2121 2130 2138 2151 2160 2164 2167 2167 2182 2182 2182 2184 2185 2187 2193 2202 2207 2215 2216 2222 2224 2236 2251 2253 2255 2258 2258 2265 2274 2281 2293 2295 2307 2313 2314 2315 2316 2316 2329 2329 2335 2336 2345 2356 2363 2366 2371 2379 2386 2386 2394 2400 2402 2403 2409 2422 2422 2428 2429 2442 2443 2445 2449 2460 2463 2465 2473 2476 2487 2491 2494 2499 2500 2513 2514 2523 2531 2532 2542 2549 2556 2556 2559 2570 2574 2576 2578 2578 2586 2594 2601 2604 2605 2608 2623 2627 2641 2645 2646 2647 2657 2657 2660 2673 2677 2678 2684 2687 2691 2703 2707 2709 2717 2734 2737 2748 2751 2755 2758 2762 2764 2771 2774 2774 2780 2787 2788 2788 2791 2804 2806 2807 2808 2831 2834 2863 2863 2864 2865 2865 2877 2886 2893 2902 2902 2902 2918 2918 2920 2920 2921 2921 2921 2934 2935 2935 2935 2936 2936 2948 2962 2978 2978 3049'
aa_list = [int(item) for item in aa_strings.split(" ")]
N = 202
resu = leaderboardCyclopeptideSequencing(aa_list, N)