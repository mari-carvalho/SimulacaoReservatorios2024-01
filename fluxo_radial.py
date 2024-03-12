# Fluxo Radial - Exponencial Integral

# Importando Bibliotecas:

import numpy as np 
import math as mt 
import scipy.special
import scipy 
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

r = np.linspace(0,L,1000) # de 0 a 10, 10 elementos 
print(r)
t = np.linspace(1,100,10) # de 1 a 10, 10 elementos 
print(t)

p_radial_matriz = np.zeros((len(t), len(r)))
print(p_radial_matriz)

# Criando a Função para a Solução Analítica:

def calculate_radial(po:float, qw:float, mi:float, k:float, h:float, phi:float, c:float, r:float, t:float) -> float:


    x = (phi*mi*c*(r[j]**2))/(4*k*t[i])
    p_radial = po - ((qw*mi)/(4*mt.pi*k*h)) * -scipy.special.expi(-x)

    return p_radial 

# Chamando as Soluções:

for i in range(len(t)): # t
    for j in range(len(r)): # r
        if i == 0:
            p_radial_matriz[i,j] = p0
        else:
            p_radial_matriz[i,j] = calculate_radial(p0, qw, mi, k, h, phi, c, r, t)
print(p_radial_matriz)

# Plotagem dos Resultados da Solução de Fluxo Radial:

for i in range(len(t)): # matriz de pressões
    plt.plot(r, p_radial_matriz[i, :], linestyle='-')
plt.xlabel('x [m]')
plt.ylabel('Pressão [kgf/cm²]')
plt.grid(True)
plt.show()