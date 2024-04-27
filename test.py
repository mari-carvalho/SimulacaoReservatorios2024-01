import numpy as np
import matplotlib.pyplot as plt 
import sympy as sp 
from solver_gauss_seidel import gauss_seidel

p0 = 19000000
pw = 9000000
qw = 0.01
q0 = 100
cc = 'pp'
mi = 0.001
k = 9.869233e-14
h = 10 
phi = 0.2
c = 2.04e-9
L = 20
A = 30  
x0 = 0
xf = L
t0 = 0
tf = 100
h_t = 0.25
h_x = 0.5

def calculate():
        
        n_x = (xf-x0)/(h_x)
        n_t = (tf-t0)/(h_t)
        print('n_x',n_x)
        print('n_t',n_t)

        x = np.zeros(int(n_x)+1) # de 0 ao tamanho do reservatório com 10 elementos na malha 
        t = np.zeros(int(n_t)+1) # de 0 a 10 segundos com 10 elementos

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
        Eppara = 0.5*10**-n

        # Número Máximo de Interações:
        maxit = 1000

        ai = +1/2*eta*rx 
        bi = 1/2 - rx*eta
        an = +(2/3)*rx*eta
        b1 = 1/2 - 2*rx*eta

        d_coeficientes = np.zeros((int(n_x)-1, int(n_t)-1))

        for j in range(len(t)): # 0 a 4 (tamanho de t), tempo 4 elemento 5; 1 a 4 (tamanho de t, mesmo for), tempo 4 elemento 4; precisa de mais um elemento 
            d_coeficientes = np.zeros((int(n_x)-1, int(n_x)-1))
            for i in range(len(d_coeficientes)): # variando a linha
                if i == 0:
                    d_coeficientes[i,0] = b1
                    d_coeficientes[i,1] = an
                elif i == len(d_coeficientes)-1: # o último, N
                    d_coeficientes[i,len(d_coeficientes)-2] = an
                    d_coeficientes[i,len(d_coeficientes)-1] = b1
                else:
                    d_coeficientes[i,i-1] = ai # linha 1, coluna 0 (i-1)
                    d_coeficientes[i,i] = bi
                    d_coeficientes[i,i+1] = ai
        print(d_coeficientes)

calculate = calculate()