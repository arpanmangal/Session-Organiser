import random

print(1)
print(10)
print(6)
print(10)
print(2)

arr = []
for i in range(400):
    arr.append([])
    for j in range(400):
        arr[i].append(0)

for i in range(400):
    for j in range(i+1):
        if (i == j):
            arr[i][j] = 0
        else:
            arr[i][j] = random.random()
            arr[j][i] = arr[i][j]

for i in range(400):
    for j in range(399):
        print (".1f" % arr[i][j], end='')
    print ("%.1f" % arr[i][399], end='')
    print()