import numpy as np 
import math as mt
import matplotlib.pyplot as plt
from fluxo_linerar_dirchlet import solucaoPressPress

def calculate_analitica_dirchlet(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,i,j,n_t,n_x):

    x = np.zeros(int(n_x)+1) # de 0 ao tamanho do reservatório com 10 elementos na malha 
    t = np.zeros(int(n_t)+1) # de 0 a 10 segundos com 10 elementos
    p = np.zeros((int(n_t)+1,int(n_x)+1))
    tam = len(x)
    print('tam_', tam)

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
    print('x', x)
    print('t', t)
    print('p', p)
    tam1 = len(x)
    print(tam1)
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

    #for i in range(len(t)):
    #    plt.plot(x, p_dirchlet_matriz[i, :], linestyle='-') # vai plotar todo o vetor de x por cada um dos vetores de pressão gerados, na linha da vez (tempo) vai plotar todas as colunas (vetor de posições x) por todas as linhas (vetor de pressões gerado)
    #plt.xlabel('x [m]')
    #plt.ylabel('Pressão [kgf/cm²]')
    #plt.grid(True)
    #plt.show()

    # Plotagem:
    time = [0,10,20,30,40,50,60,70,80,90,100]
    for i in range(len(t)):
        if t[i] in time:
            plt.plot(x, p_dirchlet_matriz[i, :], linestyle='-', label=f't = {t[i]}')

    plt.legend()
    plt.title('Solução Analítica - Dirchlet')
    plt.xlabel('Comprimento (m)')
    plt.ylabel('Pressão (psia)')
    plt.grid()
    plt.show()

    return x, t, p_dirchlet_matriz
