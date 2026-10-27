import numpy as np

cont=0.0
cont2=0.0
matriz=np.loadtxt("matriz.txt",delimiter=',')
matriz2=np.loadtxt("matriz2.txt", delimiter=',')
N,M= matriz.shape
for j in range(N):
    for k in range(M):
        if matriz[j][k]==1:
            cont+=1
print(cont)
for j in range(N):
    for k in range(M):
        if matriz2[j][k]==1:
            cont2+=1
#print(cont2)
