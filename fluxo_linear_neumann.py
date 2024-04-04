# Fluxo Linear com Condição de Contorno Fluxo-Pressão (Neumann)

# Importando as Bibliotecas:

import numpy as np 
import math as mt 
from scipy.special import erfc
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

x = np.linspace(0,L,1000) # de 0 a 10, com 10 elementos
t = np.linspace(0,100,10) # de 1 a 10, com 10 elementos

p_neumann_matriz = np.zeros((len(t), len(x)))

# Criando a Função para a Solução Analítica:

def calculate_difusividade(k:float, phi:float, mi:float, c:float) -> float:

    difusividade = k/(phi*mi*c)

    return difusividade

def calculate_neumann(p0:float, qw:float, mi:float, L:float, k:float, A:float, difusividade:float, x:float, t:float) -> float:

    p_neumann = p0 - qw*mi*L/(k*A)*(np.sqrt(4*difusividade*t[i]/(np.pi*L**2))*mt.exp(-x[j]**2/(4*difusividade*t[i]))-x[j]/L*erfc(x[j]/np.sqrt(4*difusividade*t[i])))
                
    return p_neumann

# Chamando as Soluções:

difusividade = calculate_difusividade(k, phi, mi, c)
print('difusividade', difusividade)

for i in range(len(t)): # t para cada linha vai calcular um vetor de pressões de acordo com a posição x da coluna 
    for j in range(len(x)): # x para cada posição de x (coluna) vai calcular um valor de pressão e ao final vai resultar em todo um vetor de pressões 
        if i == 0:
            p_neumann_matriz[i,j] = p0 
        else:
            p_neumann_matriz[i,j] = calculate_neumann(p0, qw, mi, L, k, A, difusividade, x, t)

# Plotagem dos Resultados da Solução de Fluxo Linear - Neumann:

for i in range(len(t)): # matriz de pressões
    plt.plot(x, p_neumann_matriz[i, :], linestyle='-')
plt.xlabel('x [m]')
plt.ylabel('Pressão [kgf/cm²]')
plt.grid(True)
plt.show()


        

