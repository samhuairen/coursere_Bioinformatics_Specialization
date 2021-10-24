def composition_kmer(text, k):
    kmers = []
    for i in range(0, len(text)-k+1):
        kmer = text[i:i+k]
        kmers.append(kmer)
    return kmers


example_data = open("dataset_197_3-2.txt", 'r')
k = int(example_data.readline().strip())
text = example_data.readline().strip()
example_data.close()

out = open("results_composition.txt", 'w')
kmers = composition_kmer(text, k)
for kmer in kmers:
    out.write(kmer+"\n")
out.close()
