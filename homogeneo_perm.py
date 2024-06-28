# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 11:02:51 2024

@author: 03950025081
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.animation import FuncAnimation
import time
import pandas as pd
import math as mt
from pandas import DataFrame as Df


def plot_animation_map_2d(grid: dict):
    times_key = list(grid.keys())
    # data_keys = list(grid.values())
    # frames = [i for i in range(len(times_key))]

    # Converta os dados do DataFrame para tipo numérico e trate valores nulos
    # for key in times_key:
    #     grid[key] = grid[key].apply(pd.to_numeric, errors='coerce').fillna(0)

    # ------------------------------------------------------------------------------------------------------------------
    def update(frame):
        if frame == 0:
            pass
        else:
            dataframe = grid[times_key[frame]]

            plt.cla()  # Limpa o eixo atual para atualizar o gráfico

            plt.imshow(dataframe, cmap='rainbow', interpolation='bicubic', origin='upper',
                       extent=(0, dataframe.index.max(), 0, dataframe.columns.max()),
                       norm=colors.Normalize(vmin=a, vmax=b)
                       )

            plt.xlabel('Comprimento x (m)')
            plt.ylabel('Comprimento y (m)')
            plt.title("Mapa de Pressão")
            plt.tight_layout()

    # Configuração do gráfico
    fig, ax = plt.subplots()
    df_map = grid[times_key[0]]
    a = grid[times_key[-1]].min().min()
    b = grid[times_key[0]].max().max()
    map_to_plot = plt.imshow(df_map, cmap='rainbow', interpolation='bicubic', origin='lower',
                             extent=(0, df_map.index.max(), 0, df_map.columns.max()),
                             norm=colors.Normalize(vmin=a, vmax=b))
    fig.colorbar(mappable=map_to_plot, ax=ax)
    plt.xlabel('Comprimento x (m)')
    plt.ylabel('Comprimento y (m)')
    plt.title("Mapa de Pressão")
    plt.tight_layout()
    ani = FuncAnimation(fig, update, frames=len(times_key) - 1, interval=500)  # Intervalo de 1000ms entre frames

    # Salvar a animação como GIF
    ani.save(f'animacao_map_homo.gif', writer='pillow', fps=60)  # 1 frame por segundo
    plt.close()


def calculate_beta(phi: float, mi: float, c: float) -> float:
    beta = 1 / (phi * mi * c)
    return beta


def calculate_rx(dt: float, dx: float) -> float:
    rx = (dt) / (dx ** 2)
    return rx


def calculate_ry(dt: float, dy: float) -> float:
    ry = (dt) / (dy ** 2)
    return ry


def Gauss_Seidel(A, b, x0, Eppara, maxit):
    ne = len(b)
    x = np.zeros(ne) if x0 is None else np.array(x0)
    iter = 0
    Epest = np.linspace(100, 100, ne)

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
n = 12  # Números de algarismos significativos
Eppara = 0.5 * 10 ** (2 - n)  # Termo relativos

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
rw = 0.1


Lx = 20  # cm
Ly = 20  # cm

# Dados Iniciais
tempo_maximo = 10000  # segundos
Po = 19000000  # °C
Pn = 0  # °C
Ps = 0  # °C
Pw = 0  # °C
Pe = 0  # °C

# Parâmetros de simulação
nx = 20
ny = 20
N = nx * ny
nt = 100

# Cálculos Iniciais
dx = Lx / (nx - 1)
dy = Ly / (ny - 1)
dt = tempo_maximo / nt

beta = calculate_beta(phi, mi, c)
rx = calculate_rx(dt, dx)
ry = calculate_ry(dt, dy)

# Solução da Equação do Calor - MDFI
inicio = time.time()

# Condição Inicial
P_old = np.full((nx, ny), Po)
P_new = np.full((nx, ny), Po)
Pold = np.full(N, Po)

# Inicializando as Matrizes
A = np.zeros((N, N))
B = np.zeros(N)
tempo = 0
h = 0

# Matriz de Permeabilidades Equivalentes:

# Definindo os índices da matriz 2D para o vetor 1D
ind = np.arange(N).reshape(nx, ny)

caminho_arquivo = 'mapa_perm_homogeneo.xlsx'
df = pd.read_excel(caminho_arquivo)
perm = df.values
print(perm)

poco = Lx / 2

results = {}

while tempo < tempo_maximo:
    h += 1
    pp = 0
    # Parte central
    px = len(perm) - 1
    for i in range(len(perm)):
        for j in range(len(perm)):
            pp += 1

            if j != 0 and j != px and i != 0 and i != px:
                print(perm[i, j])
                print(perm[i + 1, j])
                ki_mais1 = 2 / ((1 / perm[i, j]) + (1 / perm[i + 1, j])) * 9.869233e-16
                print(perm[i, j + 1])
                kj_mais1 = 2 / ((1 / perm[i, j]) + (1 / perm[i, j + 1])) * 9.869233e-16
                print(perm[i - 1, j])
                ki_menos1 = 2 / ((1 / perm[i, j]) + (1 / perm[i - 1, j])) * 9.869233e-16
                print(perm[i, j - 1])
                kj_menos1 = 2 / ((1 / perm[i, j]) + (1 / perm[i, j - 1])) * 9.869233e-16

                if i == poco and j == poco: # poço no centro
                    req = 0.5612 * dx
                    Ap = 1 + beta * rx * (ki_menos1 + ki_mais1) + beta * ry * (kj_menos1 + kj_mais1) + ((2 * mt.pi * perm[i, j]) / (phi * mi * c * dx * dy * (mt.log(req / rw))))
                    S = Pold[j] + ((2 * mt.pi * perm[i, j]) / (phi * mi * c * dx * dy * (mt.log(req / rw)))) * pwf
                else:
                    Ap = 1 + beta * rx * (ki_menos1 + ki_mais1) + beta * ry * (kj_menos1 + kj_mais1)
                    S = Pold[j]
                Aw = -rx * beta * ki_menos1
                Ae = -rx * beta * ki_mais1
                As = -ry * beta * kj_mais1
                An = -ry * beta * kj_menos1

                A[pp-1, pp-1] = Ap
                A[pp-1, pp-1 - 1] = Aw
                A[pp-1, pp-1 + 1] = Ae
                A[pp-1, pp-1 + nx] = As
                A[pp-1, pp-1 - nx] = An
                B[pp-1] = S

            # Canto Superior Esquerdo ----------------------------------------------------------------------------------
            elif j == 0 and i == 0:
                print(perm[i, j])
                print(perm[i + 1, j])
                ki_mais1 = 2 / ((1 / perm[i, j]) + (1 / perm[i + 1, j])) * 9.869233e-16
                print(perm[i, j + 1])
                kj_mais1 = 2 / ((1 / perm[i, j]) + (1 / perm[i, j + 1])) * 9.869233e-16
                m = 0
                Ap = 1 + beta * rx * ki_mais1 + beta * ry * kj_mais1
                Ae = -beta * rx * ki_mais1
                As = -beta * ry * kj_mais1
                S = Pold[m]

                A[m, m] = Ap
                A[m, m + 1] = Ae
                A[m, m + nx] = As
                B[m] = S

            # Fronteira Norte ------------------------------------------------------------------------------------------
            elif i == 0 and j != 0 and j != px:
                print(perm[i, j])
                print(perm[i + 1, j])
                kj_mais1 = 2 / ((1 / perm[i, j]) + (1 / perm[i + 1, j])) * 9.869233e-16

                Ap = 1 + beta * ry * kj_mais1
                # Aw = -rx * beta
                # Ae = -rx * beta
                As = -beta * ry * kj_mais1
                S = Pold[j]

                A[pp-1, pp-1] = Ap
                # A[pp-1, pp-1 - 1] = Aw
                # A[pp-1, pp-1 + 1] = Ae
                A[pp-1, pp-1 + nx] = As
                B[pp-1] = S

            # Canto Superior Direito -----------------------------------------------------------------------------------
            elif i == 0 and j == px:
                print(perm[i, j])
                print(perm[i, j - 1])
                ki_menos1 = 2 / ((1 / perm[i, j]) + (1 / perm[i, j - 1])) * 9.869233e-16
                print(perm[i + 1, j])
                kj_mais1 = 2 / ((1 / perm[i, j]) + (1 / perm[i + 1, j])) * 9.869233e-16
                m = nx - 1
                Ap = 1 + beta * rx * ki_menos1 + beta * ry * kj_mais1
                Aw = -beta * rx * ki_menos1
                As = -beta * ry * kj_mais1
                S = Pold[m]

                A[m, m] = Ap
                A[m, m - 1] = Aw
                A[m, m + nx] = As
                B[m] = S

            # Fronteira Oeste ------------------------------------------------------------------------------------------
            elif j == 0 and i != 0 and i != px:
                print(perm[i, j])
                print(perm[i, j + 1])
                ki_mais1 = 2 / ((1 / perm[i, j]) + (1 / perm[i, j + 1])) * 9.869233e-16

                Ap = 1 + beta * rx * ki_mais1
                Ae = -beta * rx * ki_mais1
                S = Pold[pp-1]

                A[pp-1, pp-1] = Ap
                A[pp-1, pp-1 + 1] = Ae
                B[pp-1] = S

            # Fronteira Leste ------------------------------------------------------------------------------------------
            elif j == px and i != 0 and i != px:
                print(perm[i, j])
                print(perm[i, j - 1])
                ki_menos1 = 2 / ((1 / perm[i, j]) + (1 / perm[i, j - 1])) * 9.869233e-16

                Ap = 1 + beta * rx * ki_menos1
                Aw = -beta * rx * ki_menos1
                S = Pold[pp - 1]

                A[pp - 1, pp - 1] = Ap
                A[pp - 1, pp - 1 - 1] = Aw
                B[pp - 1] = S

            # Canto Inferior Esquerdo ----------------------------------------------------------------------------------
            elif i == px and j == 0:
                print(perm[i, j])
                print(perm[i, j + 1])
                ki_mais1 = 2 / ((1 / perm[i, j]) + (1 / perm[i, j + 1])) * 9.869233e-16
                print(perm[i - 1, j])
                kj_menos1 = 2 / ((1 / perm[i, j]) + (1 / perm[i - 1, j])) * 9.869233e-16

                m = (ny - 1) * nx

                Ap = 1 + beta * rx * ki_mais1 + beta * ry * kj_menos1
                Ae = -beta * rx * ki_mais1
                An = -beta * ry * kj_menos1
                S = Pold[m]

                A[m, m] = Ap
                A[m, m + 1] = Ae
                A[m, m - nx] = An
                B[m] = S

            # Fronteira Sul --------------------------------------------------------------------------------------------
            elif i == px and j != 0 and j != px:
                print(perm[i, j])
                print(perm[i - 1, j])
                kj_menos1 = 2 / ((1 / perm[i, j]) + (1 / perm[i - 1, j])) * 9.869233e-16

                Ap = 1 + beta * ry * kj_menos1
                An = -beta * ry * kj_menos1
                S = Pold[pp - 1]

                A[pp - 1, pp - 1] = Ap
                A[pp - 1, pp - 1 - nx] = An
                B[pp - 1] = S

            # Canto Inferior Direito
            elif i == px and j == px:
                print(perm[i, j])
                print(perm[i, j - 1])
                ki_menos1 = 2 / ((1 / perm[i, j]) + (1 / perm[i, j - 1])) * 9.869233e-16
                print(perm[i - 1, j])
                kj_menos1 = 2 / ((1 / perm[i, j]) + (1 / perm[i - 1, j])) * 9.869233e-16

                m = ny * nx - 1
                Ap = 1 + beta * rx * ki_menos1 + beta * ry * kj_menos1
                Aw = -beta * rx * ki_menos1
                An = -beta * ry * kj_menos1
                S = Pold[m]

                A[m, m] = Ap
                A[m, m - 1] = Aw
                A[m, m - nx] = An
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

    results[h] = Df(P_new)


plot_animation_map_2d(grid=results)


# Plot
plt.figure()
plt.imshow(perm, extent=[0, Lx, 0, Ly], origin='upper', aspect='auto', cmap='jet')
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
