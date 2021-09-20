result_test = open('result_test.txt', 'r')
content = result_test.readlines()
position = []
sequences = open('test_motif_sequence.fasta', 'w')
for line in content:
    position.append(line.strip().split("/")[1].split(" ")[0])
    sequences.write(line.strip().split("/")[1].split(" ")[3]+"\n")
sequences.close()

meme_resu = open("meme.results.txt", 'r')
position_meme = []
for line in meme_resu.readlines():
    line_content = line.strip().split("\t")
    start_pos = line_content[3]
    position_meme.append(start_pos)
    
meme_resu.close()