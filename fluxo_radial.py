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
po = 70 # kgf/cm²
pw = 110 # kgf/cm²
qw = 35 # m³std/d
L = 200 # m 
A = 2000 # m²
h = 10 # m

r = np.linspace(0,10,10)
print(r)
t = np.linspace(1,10,10)
print(t)

p_radial_list = []
x_list = []

# Criando a Função para a Solução Analítica:

def calculate_x(phi:float, mi:float, c:float, r:float, k:float, t:float) -> float:

    x = (phi*mi*c*(r[i]**2))/(4*k*t[j])

    return x 

def calculate_radial(po:float, qw:float, mi:float, k:float, h:float, phi:float, c:float, r:float, t:float, x_list:float) -> float:

    p_radial = po - ((qw*mi)/(4*mt.pi*k*h)) * scipy.special.expn(1, x)

    return p_radial 

# Chamando as Soluções:

for i in range(1,10): # r
    for j in range(1,10): # t 

        x = calculate_x(phi, mi, c, r, k, t)
        print('valor de x', x)
        p_radial = calculate_radial(po, qw, mi, k, h, phi, c, r, t, x)
        p_radial_list.append(p_radial)

print(p_radial_list)