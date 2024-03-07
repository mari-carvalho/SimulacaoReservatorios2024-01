# Fluxo Linear com Condição de Contorno Fluxo-Pressão (Neumann)

# Importando as Bibliotecas:

import numpy as np 
import math as mt 
from scipy.special import erfc
import matplotlib.pyplot as plt 

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

x = np.linspace(0,10,10)
print(x)
t = np.linspace(1,10,10)
print(t)

p_neumann_list = [] 

# Criando a Função para a Solução Analítica:

def calculate_difusividade(k:float, phi:float, mi:float, c:float) -> float:

    difusividade = k/(phi*mi*c)

    return difusividade

def calculate_neumann(po:float, qw:float, mi:float, L:float, k:float, A:float, difusividade:float, x:float, t:float) -> float:

    p_neumann = po - ((qw*mi*L)/(k*A)) * ((mt.sqrt((4*difusividade*t[i])/(mt.pi*(L**2)))) * (mt.exp((-x[j]**2)/(4*difusividade*t[i]))) - (x[j]/L)*(erfc((x[j])/(mt.sqrt(4*difusividade*t[i])))))

    return p_neumann

# Chamando as Soluções:

difusividade = calculate_difusividade(k, phi, mi, c)
print(difusividade)

for i in range(0,10): # t
    for j in range(0,10): # x

        p_neumann = calculate_neumann(po, qw, mi, L, k, A, difusividade, x, t)
        p_neumann_list.append(p_neumann)
    print('Pressão de Neumann', p_neumann_list)




        

