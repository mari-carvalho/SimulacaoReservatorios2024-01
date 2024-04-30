# Fluxo Linear com Condição de Contorno Pressão-Pressão (Dirchlet)

# Importando as Bibliotecas:

import numpy as np 
import math as mt
from fluxo_linerar_dirchlet import solucaoPressPress
import matplotlib.pyplot as plt

# Dados de Entrada:

p0 = 19000000
pw = 9000000
qw = 0.01
mi = 0.001
k = 9.869233e-14
h = 10 
phi = 0.2
c = 2.04e-9
L = 20
A = 30 

x = np.linspace(0,L,1000) # de 0 a 10, com 1000 elementos
t = np.linspace(0,100,10) # de 0 a 1000, com 10 elementos

# não importa o tamanho de cada vetor, são independentes e formam uma matriz de i linhas e j colunas 

p_dirchlet_matriz = np.zeros((len(t), len(x))) # matriz de 10 linhas e 1000 colunas 

solucao_dirchlet = solucaoPressPress(p0, pw, phi, mi, k, L, c) # chamando a classe 

for i in range(len(t)): # vai percorrer todas as linhas, todos os tempos
    p_dirchlet_matriz[i, :] = solucao_dirchlet.PressPress(x, t[i]) # chamando a função # para cada valor de tempo (linhas), vai calcular um vetor de pressões com os valores de posição x, a linha da vez (tempo) por todas as colunas (posições x)
print(p_dirchlet_matriz)

for i in range(len(t)):
    plt.plot(x, p_dirchlet_matriz[i, :], linestyle='-') # vai plotar todo o vetor de x por cada um dos vetores de pressão gerados, na linha da vez (tempo) vai plotar todas as colunas (vetor de posições x) por todas as linhas (vetor de pressões gerado)
plt.xlabel('x [m]')
plt.ylabel('Pressão [kgf/cm²]')
plt.grid(True)
plt.show()
