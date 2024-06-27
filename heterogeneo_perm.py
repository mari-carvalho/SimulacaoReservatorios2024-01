# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 16:27:51 2024

@author: 03950025081
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import math as mt

def Gauss_Seidel(A, b, x0, Eppara, maxit):
    ne = len(b)
    x = np.zeros(ne) if x0 is None else np.array(x0)
    iter = 0
    Epest = np.linspace(100,100,ne)

    while np.max(Epest) >= Eppara and iter <= maxit:
        x_old = np.copy(x)

        for i in range(ne):
            sum1 = np.dot(A[i, :i], x[:i])
            sum2 = np.dot(A[i, i + 1:], x_old[i + 1:])
            x[i] = (b[i] - sum1 - sum2) / A[i, i]

        # Critério de parada
        Epest = np.abs((x - x_old) / x) * 100

        iter += 1

    return x

# Critério de Scarvorought, 1966
n = 12 # Números de algarismos significativos
Eppara = 0.5*10**(2-n) # Termo relativos

# Número Máximo de Iterações
maxit = 1000

# Propriedades do Material - Cobre
p0 = 19000000
pwf = 9000000
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
rw = 10 
kw = 50


Lx = 20  # cm
Ly = 20  # cm

# Dados Iniciais
tempo_maximo = 1000  # segundos
Po = 19000000  # °C
Pn = 0  # °C
Ps = 0  # °C
Pw = 0  # °C
Pe = 0  # °C

# Parâmetros de simulação
nx = 4
ny = 4
N = nx * ny
nt = 1000

# Cálculos Iniciais
dx = Lx / (nx - 1)
dy = Ly / (ny - 1)
dt = tempo_maximo / nt

def calculate_beta(phi: float, mi: float, c: float) -> float:
  beta = 1 / (phi * mi * c)

  return beta

beta = calculate_beta(phi, mi, c)

def calculate_rx(dt: float, dx: float) -> float:
  rx = (dt) / (dx ** 2)

  return rx

rx = calculate_rx(dt, dx)

def calculate_ry(dt: float, dy: float) -> float:
  ry = (dt) / (dy ** 2)

  return ry

ry = calculate_ry(dt, dy)

# Solução da Equação do Calor - MDFI
inicio = time.time()

# Condição Inicial
P_old = np.full((nx, ny), Po)
P_new = np.full((nx, ny), Po)
Pold = np.full(N,Po)

# Inicializando as Matrizes
A = np.zeros((N, N))
B = np.zeros(N)
tempo = 0
h = 0

# Matriz de Permeabilidades Equivalentes:

# Definindo os índices da matriz 2D para o vetor 1D
ind = np.arange(N).reshape(nx,ny)

while tempo < tempo_maximo:
    h += 1

    # Parte central
    caminho_arquivo = 'mapa_perm.xlsx'
    df = pd.read_excel(caminho_arquivo)
    perm = df.values
    print(perm)

    nx = len(perm)-1
    for i in range(len(perm)):
        for j in range(len(perm)):
            if j != 0 and j != nx and i != 0 and i != nx:
                print(perm[i,j])
                print(perm[i+1,j])
                ki_mais1 = 2/((1/perm[i,j]) + (1/perm[i+1,j]))
                print(perm[i,j+1])
                kj_mais1 = 2/((1/perm[i,j]) + (1/perm[i,j+1]))
                print(perm[i-1,j])
                ki_menos1 = 2/((1/perm[i,j]) + (1/perm[i-1,j]))
                print(perm[i,j-1])
                kj_menos1 = 2/((1/perm[i,j]) + (1/perm[i,j-1]))

    for i in range(1,ny-1):
      for m in range(i*nx+1,(i+1)*nx-1):

        Ap = 1 + beta*rx*(ki_menos1+ki_mais1) + beta*ry*(kj_menos1+kj_mais1)
        Aw = -rx*beta*ki_menos1
        Ae = -rx*beta*ki_mais1
        As = -ry*beta*kj_mais1
        An = -ry*beta*kj_menos1
        poco = nx/2
        if i == poco and m == poco:
            req = 0.5612*dx # Van Pollen:
            S = Pold[m] + ((2 * mt.pi * kw) / (phi * mi * c)) * ((pwf - Pold[m]) / (dx * dy * (mt.log(req / rw))))
        else:
            S = Pold[m]

        A[m,m] = Ap
        A[m,m-1] = Aw
        A[m,m+1] = Ae
        A[m,m+nx] = As
        A[m,m-nx] = An
        B[m] = S
    
    # Canto Superior Esquerdo
    #perm = np.array([[2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2]])
    nx = len(perm)-1
    for j in range(len(perm)):
        for i in range(len(perm)):
            if j == 0 and i == 0:
                print(perm[i,j])
                print(perm[i+1,j])
                ki_mais1 = 2/((1/perm[i,j]) + (1/perm[i+1,j]))
                print(perm[i,j+1])
                kj_mais1 = 2/((1/perm[i,j]) + (1/perm[i,j+1]))
    m = 0
    Ap = 1 + beta*rx*ki_mais1 + beta*ry*kj_mais1
    Ae = -beta*rx*ki_mais1
    As = -beta*ry*kj_mais1
    S = Pold[m]

    A[m,m] = Ap
    A[m,m+1] = Ae
    A[m,m+nx] = As
    B[m] = S

    # Fronteira Norte
    #perm = np.array([[2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2]])
    nx = len(perm)-1
    for i in range(len(perm)):
        for j in range(len(perm)):
            if i == 0 and j != 0 and j != nx:
                print(perm[i,j])
                print(perm[i+1,j])
                kj_mais1 = 2/((1/perm[i,j]) + (1/perm[i+1,j]))

    for m in range(1,nx-1):
      Ap = 1 + beta*ry*kj_mais1
      Aw = -rx*beta
      Ae = -rx*beta
      As = -beta*ry*kj_mais1
      S = Pold[m]

      A[m,m] = Ap
      A[m,m-1] = Aw
      A[m,m+1] = Ae
      A[m,m+nx] = As
      B[m] = S

    # Canto Superior Direito 
    #perm = np.array([[2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2]])
    nx = len(perm)-1
    for i in range(len(perm)):
        for j in range(len(perm)):
            if i == 0 and j == nx:
                print(perm[i,j])
                print(perm[i,j-1])
                ki_menos1 = 2/((1/perm[i,j]) + (1/perm[i,j-1]))
                print(perm[i+1,j])
                kj_mais1 = 2/((1/perm[i,j]) + (1/perm[i+1,j]))
    m = nx-1
    Ap = 1 + beta*rx*ki_menos1 + beta*ry*kj_mais1
    Aw = -beta*rx*ki_menos1
    As = -beta*ry*kj_mais1
    S = Pold[m]

    A[m,m] = Ap
    A[m,m-1] = Aw
    A[m,m+nx] = As
    B[m] = S

    # Fronteira Oeste
    #perm = np.array([[2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2]])
    nx = len(perm)-1
    for i in range(len(perm)):
        for j in range(len(perm)):
            if j == 0 and i != 0 and i != nx:
                print(perm[i,j])
                print(perm[i,j+1])
                ki_mais1 = 2/((1/perm[i,j]) + (1/perm[i,j+1]))
                
    for m in range(nx,(ny-2)*nx+1,nx):
      Ap = 1 + beta*rx*ki_mais1
      Ae = -beta*rx*ki_mais1
      An = -ry*beta
      As = -ry*beta
      S = Pold[m]

      A[m,m] = Ap
      A[m,m+1] = Ae
      A[m,m+nx] = As
      A[m,m-nx] = An
      B[m] = S

    # Fronteira Leste
    #perm = np.array([[2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2]])
    nx = len(perm)-1
    for i in range(len(perm)):
        for j in range(len(perm)):
            if j == nx and i != 0 and i != nx:
                print(perm[i,j])
                print(perm[i,j-1])
                ki_menos1 = 2/((1/perm[i,j]) + (1/perm[i,j-1]))
                
    for m in range(2*nx-1,(ny-1)*nx,nx):
      Ap = 1 + beta*rx*ki_menos1
      Aw = -beta*rx*ki_menos1
      An = -ry*beta
      As = -ry*beta
      S = Pold[m]

      A[m,m] = Ap
      A[m,m-1] = Aw
      A[m,m-nx] = An
      A[m,m+nx] = As
      B[m] = S

    # Canto Inferior Esquerdo
    #perm = np.array([[2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2]])
    nx = len(perm)-1
    for i in range(len(perm)):
        for j in range(len(perm)):
            if i == nx and j == 0:
                print(perm[i,j])
                print(perm[i,j+1])
                ki_mais1 = 2/((1/perm[i,j]) + (1/perm[i,j+1]))
                print(perm[i-1,j])
                kj_menos1 = 2/((1/perm[i,j]) + (1/perm[i-1,j]))
    m = (ny-1)*nx
    Ap = 1 + beta*rx*ki_mais1 + beta*ry*kj_menos1
    Ae = -beta*rx*ki_mais1
    An = -beta*ry*kj_menos1
    S = Pold[m]

    A[m,m] = Ap
    A[m,m+1] = Ae
    A[m,m-nx] = An
    B[m] = S

    # Fronteira Sul
    #perm = np.array([[2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2]])
    nx = len(perm)-1
    for i in range(len(perm)):
        for j in range(len(perm)):
            if i == nx and j != 0 and j != nx:
                print(perm[i,j])
                print(perm[i-1,j])
                kj_menos1 = 2/((1/perm[i,j]) + (1/perm[i-1,j]))
                
    for m in range((ny-1)*nx+1,ny*nx-1):
      Ap = 1 + beta*ry*kj_menos1
      Aw = -rx*beta
      Ae = -rx*beta
      An = -beta*ry*kj_menos1
      S = Pold[m]

      A[m,m] = Ap
      A[m,m-1] = Aw
      A[m,m+1] = Ae
      A[m,m-nx] = An
      B[m] = S

    # Canto Inferior Direito
    #perm = np.array([[2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 2, 2]])
    nx = len(perm)-1
    for i in range(len(perm)):
        for j in range(len(perm)):
            if i == nx and j == nx:
                print(perm[i,j])
                print(perm[i,j-1])
                ki_menos1 = 2/((1/perm[i,j]) + (1/perm[i,j-1]))
                print(perm[i-1,j])
                kj_menos1 = 2/((1/perm[i,j]) + (1/perm[i-1,j]))
    m = ny*nx-1
    Ap = 1 + beta*rx*ki_menos1 + beta*ry*kj_menos1
    Aw = -beta*rx*ki_menos1
    An = -beta*ry*kj_menos1
    S = Pold[m]

    A[m,m] = Ap
    A[m,m-1] = Aw
    A[m,m-nx] = An
    B[m] = S

    # Solução do sistema Linear
    x0 = np.ones(N)
    P = Gauss_Seidel(A, B, x0, Eppara, maxit)

    P_new = np.zeros((nx, ny), dtype=np.float64)
    # Extração do valor da variável T
    for i in range(nx):
        for j in range(ny):
            P_new[i, j] = P[ind[i, j]]

    Pold = P.copy()
    P_old = P_new.copy()
    tempo += dt


# Plot
plt.figure()
plt.imshow(perm, extent=[0, Lx, 0, Ly], origin='lower', aspect='auto', cmap='jet')
plt.colorbar()
plt.xlabel('X [m]')
plt.ylabel('Y [m]')
plt.title('Distribuição de Pressão')
plt.show()

# Plot
plt.figure()
plt.imshow(P_new, extent=[0, Lx, 0, Ly], origin='lower', aspect='auto', cmap='jet')
plt.colorbar()
plt.xlabel('X [m]')
plt.ylabel('Y [m]')
plt.title('Distribuição de Pressão')
plt.show()




print(f'Tempo de simulação: {time.time() - inicio:.2f} segundos')
