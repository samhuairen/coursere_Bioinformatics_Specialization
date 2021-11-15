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


def scoring(peptides, experimental_speactrum):
    cyclic_pepetides = cyclo_spectrum(peptides)
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


def linear_score(pepetides, spectrum):
    """return the score between a linear specturm and experimental spectrum"""
    linear_spect = linear_spectrum(pepetides)
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
    for j in range(0, l):
        peptide = lead_board[j]
        linear_scores = linear_score(peptide, spectrum)
        pepetide_score_dict[peptide] = linear_scores
    peptide_score_dict_sorted = sorted(pepetide_score_dict.items(), key=lambda x: x[1], reverse=True)
    lead_board_sorted = [item[0] for item in peptide_score_dict_sorted]
    linear_scores_sorted = [item[1] for item in peptide_score_dict_sorted]
    for j in range(N, l):
        if linear_scores_sorted[j] < linear_scores_sorted[N-1]:
            return lead_board_sorted[:j]
    return lead_board
    
lead_board_string = 'YVWHNHFVAEELAIGTNRPPRRWKMFIHVDWWSAHLGIWTHMCFMIQ GLPRHAYHQKSFVSQWYYSTKLYLVIKHSRSKMIGFIGHKMIIYVCW FRGRQGVAHVIEVYFIDRGCTQFFVAPTKLWDKNDYAAQYAEWDAHD DMWQRMKVQKHEMPIIFDGRVNKMGDWKPHWMHVTYPQKTVKWYVYP NQNGWSGHPGMWVGMNALQCAQDHDSAPPVTLICSTWWYTINDIYQA RMYEVQFVVWPMKYGSDCRNSEPGWTQVKSVKGAGVETDHHGINRTQ CIGWWATFGGKGGRRDAGPPKTLFWRVESEYIAHPEENGYKSPLDWV PHKNDFPAPCQADHISEKNIKAEFATIETYCYGGLMFEQPTMWTCQI EYWFWILEMMSHTEKCERDNSHDDSIPVHLMCSCGEPLRDHQKNTGR CSWNCQDEPETIALNPFTNAKHSECNIKLRLDHLNMADTMAKWCKTF LLVIKATEPPRYWFCVTGGFKMHKDLDRHEHGVRPKTYGLSPMYYVP AVSWQVWTIMQGMCRMIENVLANMAANNCHMVCFIAMFGKGFNTKWN VEWHEGGEMVYNSSNIKGHHLHCSNAYNQRSHEAMFWTVCQSWCDGQ MKAHPWVGMVVQELCEHGAKHFFSNLRNYGSMESKPWETSRHSTRNA HSMDTVRQVAHRRGPREREHKIHDMSNEMYCQLKNQRKQHLAEHHHC'
experimental_spectrum = '0 71 87 87 87 97 97 99 101 101 103 103 103 103 103 103 113 113 113 113 113 113 113 113 114 114 114 114 115 115 128 128 128 128 129 129 129 131 131 131 137 147 147 156 156 156 163 163 200 200 200 200 200 202 202 208 210 214 215 216 216 216 216 216 216 225 226 227 228 229 230 232 234 234 241 242 242 242 243 243 244 245 250 256 259 259 266 268 269 270 276 278 301 310 312 312 313 313 313 313 316 316 328 328 328 329 329 329 331 331 333 339 339 345 347 347 347 354 356 358 363 363 363 363 369 370 371 372 372 373 383 390 391 396 398 399 413 414 415 415 416 425 425 426 426 427 429 430 432 443 445 446 450 450 453 457 457 459 460 460 467 470 476 478 482 483 484 486 487 491 494 502 502 503 510 511 519 526 526 528 528 528 528 529 530 533 538 539 543 545 546 554 554 558 563 563 565 566 571 571 577 578 579 581 584 585 592 596 598 606 614 615 615 616 617 623 629 639 639 639 641 641 641 642 643 646 646 649 650 657 657 659 667 667 668 675 676 676 679 680 682 691 692 697 699 699 702 710 712 712 719 721 726 730 730 730 738 742 744 745 752 752 752 754 759 760 762 762 766 770 771 777 778 779 788 789 793 795 796 796 797 798 806 820 823 825 826 830 833 835 839 839 841 841 848 851 855 855 857 858 859 865 865 867 873 873 877 877 880 890 891 891 892 893 897 898 899 900 909 910 917 926 927 934 939 949 952 952 953 954 954 954 954 955 961 962 962 967 968 968 978 980 986 986 991 993 993 996 998 1004 1004 1006 1012 1013 1025 1030 1038 1040 1040 1040 1041 1045 1048 1056 1057 1062 1063 1064 1065 1065 1067 1067 1068 1068 1071 1075 1080 1089 1090 1093 1093 1094 1099 1099 1099 1106 1107 1110 1126 1127 1135 1135 1139 1143 1149 1155 1159 1162 1167 1168 1168 1169 1169 1169 1169 1170 1178 1180 1181 1187 1191 1193 1195 1202 1202 1206 1207 1208 1208 1212 1212 1213 1218 1220 1230 1238 1238 1240 1255 1256 1256 1262 1263 1266 1266 1268 1268 1281 1281 1282 1283 1283 1290 1290 1299 1308 1309 1315 1315 1315 1316 1319 1322 1322 1323 1326 1330 1331 1343 1343 1351 1353 1358 1369 1369 1369 1369 1369 1371 1376 1376 1380 1381 1384 1391 1396 1401 1403 1403 1411 1412 1412 1418 1418 1419 1422 1428 1432 1436 1437 1439 1444 1446 1456 1468 1471 1472 1474 1478 1482 1482 1483 1485 1490 1493 1497 1498 1500 1504 1506 1507 1515 1515 1518 1519 1525 1531 1531 1532 1532 1539 1545 1547 1550 1559 1567 1569 1569 1577 1578 1585 1591 1593 1594 1595 1596 1596 1601 1602 1603 1611 1618 1619 1621 1621 1628 1628 1629 1632 1632 1632 1644 1644 1646 1648 1653 1662 1682 1690 1692 1693 1694 1695 1697 1697 1698 1705 1706 1706 1714 1714 1715 1716 1716 1716 1723 1731 1732 1734 1735 1745 1745 1745 1746 1747 1749 1757 1759 1775 1794 1796 1801 1809 1809 1809 1810 1818 1819 1819 1819 1820 1825 1826 1827 1829 1829 1831 1844 1845 1845 1846 1848 1848 1851 1858 1858 1859 1860 1860 1877 1878 1897 1914 1922 1922 1923 1929 1932 1932 1934 1938 1938 1940 1945 1945 1948 1948 1957 1957 1957 1958 1958 1959 1960 1961 1961 1965 1965 1972 1974 1990 2007 2025 2026 2026 2041 2042 2045 2045 2047 2048 2052 2053 2053 2058 2058 2060 2060 2061 2061 2061 2075 2078 2085 2085 2086 2088 2088 2088 2103 2108 2121 2129 2139 2148 2156 2156 2160 2161 2161 2166 2170 2171 2172 2173 2173 2174 2174 2181 2189 2198 2199 2202 2206 2207 2208 2214 2216 2216 2217 2234 2242 2261 2269 2269 2270 2271 2274 2274 2276 2277 2285 2285 2288 2292 2301 2303 2308 2311 2312 2321 2328 2329 2330 2330 2331 2337 2344 2355 2370 2372 2373 2373 2374 2386 2388 2389 2398 2402 2405 2414 2415 2416 2416 2417 2421 2424 2425 2442 2444 2457 2458 2459 2468 2468 2471 2476 2484 2485 2486 2487 2499 2501 2511 2519 2528 2530 2530 2531 2535 2536 2545 2545 2555 2570 2571 2571 2571 2572 2581 2584 2586 2587 2589 2590 2597 2598 2602 2624 2639 2642 2642 2644 2657 2659 2659 2666 2668 2673 2684 2684 2684 2686 2698 2699 2699 2700 2700 2700 2701 2702 2727 2737 2745 2752 2755 2756 2758 2770 2771 2787 2787 2787 2787 2799 2800 2802 2812 2813 2813 2813 2813 2815 2829 2831 2857 2858 2859 2874 2883 2884 2884 2900 2900 2900 2902 2908 2913 2914 2915 2915 2916 2918 2918 2926 2926 2944 2944 2958 2986 2987 2997 3001 3005 3012 3013 3013 3013 3015 3016 3022 3028 3029 3029 3031 3033 3037 3047 3047 3057 3071 3100 3100 3114 3115 3116 3125 3126 3126 3128 3133 3134 3141 3146 3147 3148 3150 3152 3160 3168 3176 3185 3201 3213 3227 3228 3229 3231 3238 3239 3244 3247 3247 3249 3260 3263 3265 3277 3278 3283 3296 3314 3316 3326 3341 3341 3341 3341 3342 3348 3363 3364 3368 3376 3376 3378 3387 3391 3396 3410 3411 3428 3429 3438 3444 3454 3461 3463 3472 3477 3479 3489 3492 3499 3504 3524 3524 3525 3539 3541 3543 3551 3557 3559 3564 3575 3576 3578 3590 3593 3607 3627 3630 3638 3652 3654 3654 3654 3655 3677 3679 3680 3692 3704 3704 3704 3706 3739 3741 3757 3767 3767 3767 3777 3783 3791 3792 3793 3801 3811 3817 3835 3854 3867 3870 3880 3897 3898 3904 3906 3906 3906 3914 3920 3922 3939 3964 3967 3982 3983 3993 4020 4026 4026 4037 4043 4045 4053 4067 4067 4070 4079 4095 4108 4130 4139 4140 4151 4173 4174 4182 4182 4192 4198 4222 4226 4230 4245 4261 4269 4287 4295 4295 4295 4302 4345 4354 4359 4376 4382 4382 4389 4398 4416 4416 4451 4458 4469 4490 4503 4504 4510 4529 4545 4554 4561 4597 4607 4616 4618 4618 4632 4658 4694 4710 4717 4732 4744 4745 4745 4747 4797 4832 4841 4860 4861 4873 4873 4944 4947 4960 4960 4970 4974 5057 5061 5073 5075 5107 5160 5172 5189 5236 5275 5286 5323 5389 5438 5552'
N = 5
lead_board = [str(item) for item in lead_board_string.split(" ")]
spectrum = [int(item) for item in experimental_spectrum.split(" ")]
resu = trim(lead_board, spectrum, N)

Spectrum = [0, 71, 178, 202, 202, 202, 333, 333, 333, 404, 507, 507]
Peptide = 'MAMA'
Peptide2 = 'PEEP'
Spectrum2 = [int(item) for item in '0 97 97 129 129 194 203 226 226 258 323 323 323 355 403 452'.split(" ")]
score = scoring(Peptide, Spectrum)
score2 = linear_score(Peptide2, Spectrum2)

