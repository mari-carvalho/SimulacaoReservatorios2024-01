# Formulação Explícita BTCS:

import numpy as np 
import matplotlib.pyplot as plt 
import sympy as sp 
from solver_gauss_seidel import gauss_seidel


# Definindo as Variáveis de Entrada:

p0 = 19000000
pw = 9000000
qw = 0.01
q0 = 100
entrada = pw 
saida = p0
mi = 0.001
k = 9.869233e-14
h = 10 
phi = 0.2
c = 2.04e-9
L = 10
A = 30  
x0 = 0
xf = L
t0 = 0
tf = 100
h_t = 0.25
h_x = 0.5

n_x = (xf-x0)/(h_x)
n_t = (tf-t0)/(h_t)

x = np.zeros(int(n_x)+1) # de 0 ao tamanho do reservatório com 10 elementos na malha 
t = np.zeros(int(n_t)+1) # de 0 a 10 segundos com 10 elementos
p = np.zeros((int(n_x)+1,int(n_x)+1))

# Alimentando os vetores:
for i in range(len(x)):
    if i == 0:
        x[i] = x0
    elif i == 1:
        x[i] = i*(h_x/2)
    elif i == len(x):
        x[i] = L 
    elif i == len(x)-1:
        x[i] = x[i-1] + (h_x/2)
    else:
        x[i] = x[i-1] + h_x
        
for i in range(len(t)):
    if i == 0:
        t[i] = t0
    elif i == len(t):
        t[i] = tf
    else:
        t[i] = i*h_t      

print('x', x)
print('t', t)
print('p', p)

def calculate_eta(k:float, phi:float, mi:float, c:float) -> float:
    eta = k/(phi*mi*c)

    return eta 

eta = calculate_eta(k,phi,mi,c)

def calculate_rx(h_t:float, h_x:float) -> float:
    rx = (h_t)/(h_x**2)

    return rx 

rx = calculate_rx(h_t, h_x)

# Criando o método MDF de BTCS:

# Critérios de Scarvorought, 1966:
n = 12 
Eppara = 0.5*10**12-n

# Número Máximo de Interações:
maxit = 1000

ai = -eta*rx 
bi = 1 + 2*rx*eta
an = -(4/3)*rx*eta
b1 = 1 + 4*rx*eta

p_coeficientes = np.zeros((int(n_x), int(n_x)))
p_old = np.ones(int(n_x))*p0 # vai atualizar cada linha da matriz 
p_solucoes = np.zeros((int(n_t)+1, int(n_x)))
d = np.zeros(int(n_x)) # vai guardar os valores de p no tempo anterior mais 8/3*eta*rx*(Pw ou P0)
p_solucoes[h, :]  = p0

for i in range(len(x)): # variando a linha
    if i == 0:
        p_coeficientes[i,0] = b1
        p_coeficientes[i,1] = an
        d[i] = p_old[i] + 8/3*rx*eta*pw
    elif i == len(x)-1: # o último, N+1
        p_coeficientes[i,len(x)-2] = an 
        p_coeficientes[i,len(x)-1] = b1
        d[i] = p_old[i] + 8/3*rx*eta*p0
    else:
        print(i)
        p_coeficientes[i,i-1] = ai # linha 1, coluna 0 (i-1)
        p_coeficientes[i,i] = bi
        p_coeficientes[i,i+1] = ai
        d[i] = p_old[i] # condição central é 0 

    x0 = p_old # os primeiros valores de chute inicial vão ser os valores de p calculadas no tempo anterior 
    p_new = gauss_seidel(p_coeficientes,d,x0,Eppara,maxit)
    p_old = p_new # atualiza a matriz, coloca o vetor de pressões calculado no tempo anterior (p_new) em p_old 
    p_solucoes[h, :] = p_new # vai guardar na matriz de solucoes todos os vetores de pressão calculados nos tempos 
            
print(p_solucoes)


