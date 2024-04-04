# Formulação Explícita FTCS:

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
L = 20
A = 30  
x0 = 0
xf = L
x_elementos = 10
t0 = 0
tf = 10
t_elementos = 10
delta_t = (tf-t0)/t_elementos # é o passo, delta t 
delta_x = (xf-x0)/x_elementos
print('delta_x', delta_x)
print('delta_t', delta_t)

x = np.zeros(x_elementos) # de 0 ao tamanho do reservatório com 10 elementos na malha 
t = np.zeros(t_elementos) # de 0 a 10 segundos com 10 elementos
p = np.zeros((x_elementos,t_elementos))
print('x', x)
print('t', t)
print('p', p)

def calculate_eta(k:float, phi:float, mi:float, c:float) -> float:
    eta = k/(phi*mi*c)

    return eta 

eta = calculate_eta(k,phi,mi,c)

def calculate_rx(delta_t:float, delta_x:float) -> float:
    rx = (delta_t)/(delta_x**2)

    return rx 

rx = calculate_rx(delta_t, delta_x)

# Criando o método MDF de FTCS:

for i in range(t_elementos):
    for j in range(x_elementos):
        if entrada == pw and saida == p0:
            if j == 0:
                p[j, i+delta_t] = pw
            elif j == 1: # bloco 1
                p[j, i+delta_t] = (8/3)*eta*rx*pw + (1-(4*eta*rx))*p[1,i] + (4/3)*eta*rx*p[2,i]
            elif j == len(x): # bloco N 
                p[j, i+delta_t] = (4/3)*eta*rx*p[j-1,i] + (1-4*eta*rx)*p[j,i] + (8/3)*eta*rx*p0
            elif j == len(x) + 1: # N+1
                p[j, i+delta_t] = p0
            else: # blocos interiores 
                p[j, i+delta_t] = eta*rx*p[j-1,i] + (1-2*eta*rx)*p[j,i] + eta*rx*p[j+1,i]
        if entrada == qw and saida == p0:
            if j == 0:
                p[j, i+delta_t] = qw
            elif j == 1: # bloco 1
                p[j, i+delta_t] = eta*rx*p[2,i] + (1-eta*rx)*p[1,i] + eta*rx*qw*((mi*delta_x)/(k*A))
            elif j == len(x): #bloco N
                p[j, i+delta_t] = (4/3)*eta*rx*p[j-1,i] + (1-4*eta*rx)*p[j,i] + (8/3)*eta*rx*p0
            elif j == len(x) + 1: # N+1
                p[j, i+delta_t] = p0
            else: # blocos interiores
                p[j, i+delta_t] = eta*rx*p[j-1,i] + (1-2*eta*rx)*p[j,i] + eta*rx*p[j+1,i]
        if entrada == qw and saida == q0:
            if j == 0:
                p[j, i+delta_t] = qw
            elif j == 1: # bloco 1
                p[j, i+delta_t] = eta*rx*p[2,j] + (1-eta*rx)*p[1,j] + eta*rx*qw*((mi*delta_x)/(k*A))
            elif j == len(x): #bloco N
                p[j, i+delta_t] = -eta*rx*p[j-1,i] + (1+eta*rx)*p[j,i] - eta*rx*q0*((mi*delta_x)/(k*A))
            elif j == len(x) + 1: # N+1
                p[j, i+delta_t] = p0
            else: # blocos interiores
                p[j, i+delta_t] = eta*rx*p[j-1,i] + (1-2*eta*rx)*p[j,i] + eta*rx*p[j+1,i]
        x[j+1] = x[j] + delta_x
    t[i+1] = t[i] + delta_t

            

