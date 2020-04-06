import numpy as np


positions = []
with open('methylation_frequency.tsv') as f:
    f.readline()
    for line in f:
        positions.append(line.split()[1])
dim = len(positions)

summary = {}
with open('methylation_calls.tsv') as f:
    f.readline()
    current_read = ''
    for line in f:
        row = line.split()
        if row[4] != current_read:
            if current_read != '':
                summary[current_read] = vector
            n = 0
            current_read = row[4]
            vector = np.zeros(dim)
        else:
            vector[n] = row[5]
            n += 1


def nonzeros(vector):
    n = 0
    for i in vector:
        if i != 0:
            n += 1
    return n

print(dim)
for i in summary.values():
    temp = nonzeros(i)
    if temp > 800:
        print(temp)
