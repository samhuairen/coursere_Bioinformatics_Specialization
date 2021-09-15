profile = {
    'A': [0.2, 0.2, 0, 0, 0, 0, 0.9, 0.1, 0.1, 0.1, 0.3, 0],
    'T': [0.1, 0.6, 0, 0, 0, 0, 0, 0.4, 0.1, 0.2, 0.4, 0.6],
    'C': [0, 0, 1, 1, 0.9, 0.9, 0.1, 0, 0, 0, 0, 0],
    'G': [0.7, 0.2, 0, 0, 0.1, 0.1, 0, 0.5, 0.8, 0.7, 0.3, 0.4]
}

import math

entropy_init = 0
for key, value in profile.items():
    for p in value:
        if p > 0:
            entropy_init += p * math.log(p, 2)

entropy = entropy_init * (-1)