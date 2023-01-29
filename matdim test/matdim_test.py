import numpy as np

starting = np.random.randint(5, size=(4, 4, 4))

output = np.zeros((2, 2, 2))


for i in range(2):
    for j in range(2):
        for k in range(2):
            output[i, j, k] = np.mean(starting[2*i:2*(i+1), 2*j:2*(j+1), 2*k:2*(k+1)])

print(starting)
print(output)


