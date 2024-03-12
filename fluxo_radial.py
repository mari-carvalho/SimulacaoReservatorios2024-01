# Fluxo Radial - Exponencial Integral

# Importando Bibliotecas:

import numpy as np 
import math as mt 
import scipy.special

# Dados de Entrada:

k = 100 # md
mi = 3 # cp
phi = 0.20  
c = 130*10**-6 # (kgf/cm²)-1
p0 = 70 # kgf/cm²
pw = 110 # kgf/cm²
qw = 35 # m³std/d
L = 200 # m 
A = 2000 # m²
h = 10 # m

r = np.linspace(0,10,10) # de 0 a 10, 10 elementos 
print(r)
t = np.linspace(1,10,10) # de 1 a 10, 10 elementos 
print(t)

p_radial_matriz = np.zeros((len(r), len(t)))
print(p_radial_matriz)

# Criando a Função para a Solução Analítica:

def calculate_radial(po:float, qw:float, mi:float, k:float, h:float, phi:float, c:float, r:float, t:float) -> float:


    x = (phi*mi*c*(r[j]**2))/(4*k*t[i])
    p_radial = po - ((qw*mi)/(4*mt.pi*k*h)) * scipy.special.expn(1, x)

    return p_radial 

# Chamando as Soluções:

for i in range(len(t)): # t
    for j in range(len(r)): # r
        if i == 0 and j == 0:
            p_radial_matriz[i,j] = p0
        elif i == 10 and j == 10:
            p_radial_matriz[i,j] = pw
        else:
            p_radial_matriz[i,j] = calculate_radial(p0, qw, mi, k, h, phi, c, r, t)
print(p_radial_matriz)