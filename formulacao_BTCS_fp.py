# Formulação Explícita BTCS:

import numpy as np 
import matplotlib.pyplot as plt 
import sympy as sp 

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

ai = -eta*rx 
bi = 1 + 2*rx*eta
an = -rx*eta
b1 = 1 + rx*eta
bn = eta*rx
ci = 1 - rx*eta


for i in range(len(x)): # variando a linha
    if i == 0:
        p[i,0] = b1
        p[i,1] = an
    elif i == len(x)-1: # o último, N+1
        p[i,len(x)-2] = bn 
        p[i,len(x)-1] = ci
    else:
        print(i)
        p[i,i-1] = ai # linha 1, coluna 0 (i-1)
        p[i,i] = bi
        p[i,i+1] = ai
            
print(p)

