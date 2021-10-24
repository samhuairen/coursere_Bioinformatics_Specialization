def spell_path_genome(strings):
    path = ''
    for string in strings[:-1]:
        path += string[0]
    return path+strings[-1]


data = open("dataset_198_3-2.txt", "r")
strings = data.read().strip().split('\n')
data.close()

path = spell_path_genome(strings)
print(path)


    