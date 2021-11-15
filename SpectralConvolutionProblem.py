def convolution(specturm):
    """return all positive difference in the spectrum"""
    specturm_sorted = sorted(specturm)
    specturm2 = specturm_sorted[:-1]
    convolution_list = []
    for item1 in specturm_sorted:
        for item2 in specturm2:
            convolution_list.append(item1-item2)
    convolution_result = filter(lambda x: x>0, convolution_list)
    return list(convolution_result)

aa_strings = '0 71 71 71 87 99 115 128 129 147 163 170 186 200 215 216 218 218 234 243 262 285 287 289 314 330 333 333 344 347 381 401 404 413 415 418 434 448 452 459 480 500 505 519 530 530 551 551 562 576 581 601 622 629 633 647 663 666 668 677 680 700 734 737 748 748 751 767 792 794 796 819 838 847 863 863 865 866 881 895 911 918 934 952 953 966 982 994 1010 1010 1010 1081'
spectrums = [int(item) for item in aa_strings.split(" ")]
resu = convolution(spectrums)

0 86 160 234 308 320 382
Strings = '0 86 160 234 308 320 382'
Spectrums = [int(item) for item in Strings.split(" ")]
resu = convolution(Spectrums)

resu_dict = {}
for item in resu:
    resu_dict[item] = resu.count(item)
