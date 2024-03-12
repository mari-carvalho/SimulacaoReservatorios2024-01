# Fluxo Linear com Condição de Contorno Pressão-Pressão (Dirchlet)

# Importando Bibliotecas:

import numpy as np 
import math as mt
from test_aula import solucaoPressPress
import matplotlib.pyplot as plt

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

x = np.linspace(0,L,1000) # de 0 a 10, com 10 elementos
t = np.linspace(0,1000,10) # de 1 a 10, com 10 elementos

p_dirchlet_matriz = np.zeros((len(t), len(x)))

solucao_dirchlet = solucaoPressPress(p0, pw, phi, mi, k, L, c)

for i in range(len(t)):
    p_dirchlet_matriz[i, :] = solucao_dirchlet.PressPress(x, t[i])
print(p_dirchlet_matriz)

for i in range(len(t)):
    plt.plot(x, p_dirchlet_matriz[i, :], linestyle='-')
plt.xlabel('x [m]')
plt.ylabel('Pressão [kgf/cm²]')
plt.grid(True)
plt.show()
