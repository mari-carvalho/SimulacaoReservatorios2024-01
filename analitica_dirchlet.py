import numpy as np 
import math as mt
import matplotlib.pyplot as plt
from fluxo_linerar_dirchlet import solucaoPressPress
from scipy.special import erfc

class analiticas():
    def calculate_analitica_dirchlet(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,i,j,n_t,n_x, variancia):

        x = np.zeros(int(n_x)+1) # de 0 ao tamanho do reservatório com 10 elementos na malha
        t = np.zeros(int(n_t)+1) # de 0 a 10 segundos com 10 elementos
        p = np.zeros((int(n_t)+1,int(n_x)+1))
        tam = len(x)
        print('tam_', tam)

        if variancia == 'tempo':
            v = 'Steps de Tempo'
        elif variancia == 'malha':
            v = 'Malha'

        h_t = i
        h_x = j

        # Alimentando os vetores:
        for i in range(len(x)):
            if i == 0:
                x[i] = x0
            elif i == 1:
                x[i] = i*(h_x/2)
            elif i == len(x):
                x[i] = L
            elif i == len(x)-1:
                print(i)
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

        # não importa o tamanho de cada vetor, são independentes e formam uma matriz de i linhas e j colunas

        p_dirchlet_matriz = np.zeros((len(t), len(x))) # matriz de 10 linhas e 1000 colunas
        tam2 = len(p_dirchlet_matriz)
        print(tam2)

        solucao_dirchlet = solucaoPressPress(p0, pw, phi, mi, k, L, c) # chamando a classe

        for i in range(len(t)): # vai percorrer todas as linhas, todos os tempos
            for j in range(len(x)):
                if i == 0:
                    p_dirchlet_matriz[i,j] = p0
                else:
                    p_dirchlet_matriz[i, j] = solucao_dirchlet.PressPress(x[j], t[i]) # chamando a função # para cada valor de tempo (linhas), vai calcular um vetor de pressões com os valores de posição x, a linha da vez (tempo) por todas as colunas (posições x)
        print(p_dirchlet_matriz)

        # Plotagem:
        time = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        for i in range(len(t)):
            if t[i] in time:
                plt.plot(x, p_dirchlet_matriz[i, :], linestyle='-', label=f't = {t[i]}')

        legend_label = f'{v} {n_x: .2f}' if variancia == "malha" else f'{v} {n_t: .2f}'
        plt.legend(labels=[legend_label])
        plt.title('Solução Analítica - Dirchlet')
        plt.xlabel('Comprimento [m]')
        plt.ylabel('Pressão [Pa]')
        plt.grid()
        plt.show()

        #Plot 3D BTCS_pp
        X, T = np.meshgrid(x, t)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, T, p_dirchlet_matriz[:, :], cmap='viridis')
        ax.set_xlabel('Comprimento [m]')
        ax.set_ylabel('Tempo [s]')
        ax.set_zlabel('Pressão(x,y) [Pa]')
        ax.set_title('Solução Analítica - Dirchlet')
        fig.text(0.02, 0.02, legend_label, color='black', ha='left')
        plt.show()


        return x, t, p_dirchlet_matriz

    def calculate_analitica_neumann(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,i,j,n_t,n_x,variancia):

        x = np.zeros(int(n_x) + 1)  # de 0 ao tamanho do reservatório com 10 elementos na malha
        t = np.zeros(int(n_t) + 1)  # de 0 a 10 segundos com 10 elementos
        p = np.zeros((int(n_t) + 1, int(n_x) + 1))
        tam = len(x)
        print('tam_', tam)

        if variancia == 'tempo':
            v = 'Steps de Tempo'
        elif variancia == 'malha':
            v = 'Malha'

        h_t = i
        h_x = j

        # Alimentando os vetores:
        for i in range(len(x)):
            if i == 0:
                x[i] = x0
            elif i == 1:
                x[i] = i * (h_x / 2)
            elif i == len(x):
                x[i] = L
            elif i == len(x) - 1:
                print(i)
                x[i] = x[i - 1] + (h_x / 2)
            else:
                x[i] = x[i - 1] + h_x

        for i in range(len(t)):
            if i == 0:
                t[i] = t0
            elif i == len(t):
                t[i] = tf
            else:
                t[i] = i * h_t
        print('x', x)
        print('t', t)
        print('p', p)
        tam1 = len(x)
        print(tam1)
        # não importa o tamanho de cada vetor, são independentes e formam uma matriz de i linhas e j colunas

        p_neumann_matriz = np.zeros((len(t), len(x)))  # matriz de 10 linhas e 1000 colunas

        h_t = i
        h_x = j

        def calculate_difusividade(k: float, phi: float, mi: float, c: float) -> float:
            difusividade = k / (phi * mi * c)

            return difusividade
        def calculate_neumann(p0: float, qw: float, mi: float, L: float, k: float, A: float, difusividade: float, x: float,
                              t: float) -> float:
            p_neumann = p0 - qw * mi * L / (k * A) * (
                        np.sqrt(4 * difusividade * t[i] / (np.pi * L ** 2)) * mt.exp(-x[j] ** 2 / (4 * difusividade * t[i])) -
                        x[j] / L * erfc(x[j] / np.sqrt(4 * difusividade * t[i])))

            return p_neumann

        # Chamando as Soluções:

        difusividade = calculate_difusividade(k, phi, mi, c)
        print('difusividade', difusividade)

        for i in range(len(t)):  # t para cada linha vai calcular um vetor de pressões de acordo com a posição x da coluna
            for j in range(
                    len(x)):  # x para cada posição de x (coluna) vai calcular um valor de pressão e ao final vai resultar em todo um vetor de pressões
                if i == 0:
                    p_neumann_matriz[i, j] = p0
                else:
                    p_neumann_matriz[i, j] = calculate_neumann(p0, qw, mi, L, k, A, difusividade, x, t)

        # Plotagem:
        time = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        for i in range(len(t)):
            if t[i] in time:
                plt.plot(x, p_neumann_matriz[i, :], linestyle='-', label=f't = {t[i]}')

        legend_label = f'{v} {n_x: .2f}' if variancia == "malha" else f'{v} {n_t: .2f}'
        plt.legend(labels=[legend_label])
        plt.title('Solução Analítica - Neumann')
        plt.xlabel('Comprimento [m]')
        plt.ylabel('Pressão [Pa]')
        plt.grid()
        plt.show()

        #Plot 3D BTCS_pp
        X, T = np.meshgrid(x, t)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, T, p_neumann_matriz[:, :], cmap='viridis')
        ax.set_xlabel('Comprimento [m]')
        ax.set_ylabel('Tempo [s]')
        ax.set_zlabel('Pressão(x,y) [Pa]')
        ax.set_title('Solução Analítica - Neumann')
        fig.text(0.02, 0.02, legend_label, color='black', ha='left')
        plt.show()

        return x, t, p_neumann_matriz
