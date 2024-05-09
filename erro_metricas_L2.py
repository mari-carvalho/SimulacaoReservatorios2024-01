# Função para analisar inclinação das retas do decaimento do Tempo (Linear ou Quadrática)

import numpy as np
import matplotlib.pyplot as plt
from solvers import solvers
from analitica_dirchlet import calculate_analitica_dirchlet
from FTCS import FTCS
from BTCS import BTCS
from CN import CN
from prettytable import PrettyTable
class erros():
    def calculate_erros_tempo():

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

        h_t = [0.2, 0.1, 0.05]
        h_x = 0.8
        j = h_x

        def calculate_n_t(tf,t0,i):

            n_t = (tf - t0) / (i)
            return n_t

        n_x = (xf - x0) / (h_x)

        def calculate_h_t_ex():
            x_ex_array = []
            t_ex_array = []
            p_ex_array = []
            n_t_array = []
            for i in h_t:

                n_t = calculate_n_t(tf, t0, i)
                x_ex, t_ex, p_ex = calculate_analitica_dirchlet(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,i,j,n_t,n_x)

                x_ex_array.append(x_ex)
                p_ex_array.append(p_ex)
                t_ex_array.append(t_ex)
                n_t_array.append(n_t)

            return x_ex_array, t_ex_array, p_ex_array, n_t_array


        x_ex_array, t_ex_array, p_ex_array, n_t_array = calculate_h_t_ex()


        def calculate_h_t_calc():
            x_calc_array = []
            t_calc_array = []
            p_calc_array = []
            n_t_array = []
            for i in h_t:

                n_t = calculate_n_t(tf, t0, i)
                x_calc, t_calc, p_calc = BTCS.calculate_BTCS_pp_gs(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x)

                x_calc_array.append(x_calc)
                p_calc_array.append(p_calc)
                t_calc_array.append(t_calc)
                n_t_array.append(n_t)

            return x_calc_array, t_calc_array, p_calc_array, n_t_array


        x_calc_array, t_calc_array, p_calc_array, n_t_array = calculate_h_t_calc()

        # Cálculo do Erro:

        # Norma L2
        L2_list = []
        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            sum = 0
            y_calc = p_calc[50] # selecinar a coluna
            y_ex = p_ex[50] # selecinar a coluna
            for k in range(len(y_ex)): # acesso a cada linha
                sum = sum + ((abs((y_ex[k]-y_calc[k])/(y_ex[k])))**2)
            L2 = np.sqrt((1/(n_x**2))*sum)
            L2_list.append(L2)
        print('L2', L2_list)

        L2_log_list = np.log(L2_list)

        # Norma E_inf
        E_inf_depois_list = []

        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[50]
            y_ex = p_ex[50]
            E_inf_antes_list = []
            for k in range(len(y_ex)):
                E_inf_antes = abs((y_ex[k] - y_calc[k]))
                E_inf_antes_list.append(E_inf_antes)
            E_inf_depois = max(E_inf_antes_list)
            E_inf_depois_list.append(E_inf_depois)

        E_inf_depois_log_list = np.log(E_inf_depois_list)

        # Norma E_rel
        err_rel_total_list = []

        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[50]
            y_ex = p_ex[50]
            err_rel_list = []
            sum = 0
            for k in range(len(y_ex)):
                err_rel = abs((y_ex[k] - y_calc[k])/(y_ex[k]))
                err_rel_list.append(err_rel)
            for j in range(len(err_rel_list)):
                sum = sum + err_rel_list[j]
            err_rel_total = 1/n_x * sum
            err_rel_total_list.append(err_rel_total)

        err_rel_total_log_list = np.log(err_rel_total_list)

        h_t_log_list = np.log(h_t)

        # Plotagem:
        plt.plot(h_t_log_list, L2_list, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
        plt.title('Norma Euclidiana - Norma L2')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup t$ [s]')
        plt.ylabel('L2')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
        coeffs = np.polyfit(h_t_log_list, L2_log_list, grau) # ajusta a linha de tendência aos dados e retorna os coeficientes do polinômio
        tendencia = np.poly1d(coeffs) # cria um polinômio a partir dos coeficientes retornados por polyfit
        x_tendencia = np.linspace(min(h_t_log_list), max(h_t_log_list), 100) # cria um conjunto de pontos para suavizar a linha de tendência
        plt.plot(x_tendencia, tendencia(x_tendencia), color='green', linestyle='dashed', label='Linha de Tendência Linear')
        plt.plot(h_t_log_list, L2_log_list, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
        plt.title('Norma Euclidiana - Norma L2')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup t$ [s]')
        plt.ylabel('ln (L2)')
        plt.show()

        # Plotagem:
        plt.plot(h_t_log_list, E_inf_depois_list, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
        plt.title('Erro Absoluto Máximo - Norma E$ \infty$')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup t$ [s]')
        plt.ylabel('E$ \infty$')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
        coeffs = np.polyfit(h_t_log_list, E_inf_depois_log_list, grau) # ajusta a linha de tendência aos dados e retorna os coeficientes do polinômio
        tendencia = np.poly1d(coeffs) # cria um polinômio a partir dos coeficientes retornados por polyfit
        x_tendencia = np.linspace(min(h_t_log_list), max(h_t_log_list), 100) # cria um conjunto de pontos para suavizar a linha de tendência
        plt.plot(x_tendencia, tendencia(x_tendencia), color='green', linestyle='dashed', label='Linha de Tendência Linear')
        plt.plot(h_t_log_list, E_inf_depois_log_list, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
        plt.title('Erro Absoluto Máximo - Norma E$ \infty$')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup t$ [s]')
        plt.ylabel('ln (E$ \infty$)')
        plt.show()

        # Plotagem:
        plt.plot(h_t_log_list, err_rel_total_list, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
        plt.title('Erro Relativo - Norma L1')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup t$ [s]')
        plt.ylabel('L1')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
        coeffs = np.polyfit(h_t_log_list, err_rel_total_log_list, grau) # ajusta a linha de tendência aos dados e retorna os coeficientes do polinômio
        tendencia = np.poly1d(coeffs) # cria um polinômio a partir dos coeficientes retornados por polyfit
        x_tendencia = np.linspace(min(h_t_log_list), max(h_t_log_list), 100) # cria um conjunto de pontos para suavizar a linha de tendência
        plt.plot(x_tendencia, tendencia(x_tendencia), color='green', linestyle='dashed', label='Linha de Tendência Linear')
        plt.plot(h_t_log_list, err_rel_total_log_list, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
        plt.title('Erro Relativo - Norma L1')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup t$ [s]')
        plt.ylabel('ln (L1)')
        plt.show()

        # Tabelas:
        tabela = PrettyTable(['N° de Blocos', 'Steps de Tempo', 'Norma Euclidiana - L2', 'Erro Absoluto Máximo - E_inf', 'Erro Relativo - Norma L1'])

        for n_t_val, L2_val, E_inf_val, err_rel_val in zip(n_t_array, L2_list, E_inf_depois_list, err_rel_total_list):
            rounded_L2_val = round(L2_val, 3)
            rounded_E_inf_val = round(E_inf_val, 3)
            rounded_err_rel_val = round(err_rel_val, 3)
            rounded_n_t_val = round(n_t_val, 3)
            tabela.add_row([n_x, rounded_n_t_val, rounded_L2_val, rounded_E_inf_val, rounded_err_rel_val])

        print(tabela)

    def calculate_erros_malha():

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

        h_x = [5, 4, 3, 2, 1]
        h_t = 0.8
        i = h_t

        def calculate_n_x(xf, x0, j):

            n_x = (xf - x0) / (j)
            return n_x

        n_t = (tf - t0) / (h_t)

        def calculate_h_t_ex():
            x_ex_array = []
            t_ex_array = []
            p_ex_array = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_ex, t_ex, p_ex = calculate_analitica_dirchlet(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0,
                                                                tf, i, j, n_t, n_x)

                x_ex_array.append(x_ex)
                p_ex_array.append(p_ex)
                t_ex_array.append(t_ex)
                n_x_array.append(n_x)

            return x_ex_array, t_ex_array, p_ex_array, n_x_array

        x_ex_array, t_ex_array, p_ex_array, n_x_array = calculate_h_t_ex()

        def calculate_h_t_calc():
            x_calc_array = []
            t_calc_array = []
            p_calc_array = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_calc, t_calc, p_calc = BTCS.calculate_BTCS_pp_gs(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0,
                                                                tf, i, j, n_t, n_x)

                x_calc_array.append(x_calc)
                p_calc_array.append(p_calc)
                t_calc_array.append(t_calc)
                n_x_array.append(n_x)

            return x_calc_array, t_calc_array, p_calc_array, n_x_array

        x_calc_array, t_calc_array, p_calc_array, n_x_array = calculate_h_t_calc()

        # Cálculo do Erro:

        # Norma L2
        L2_list = []
        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            sum = 0
            y_calc = p_calc[:, 3]  # selecinar a coluna
            y_ex = p_ex[:, 3]  # selecinar a coluna
            n_x = n_x_array[i]
            for k in range(len(y_ex)):  # acesso a cada linha
                sum = sum + ((abs((y_ex[k] - y_calc[k]) / (y_ex[k]))) ** 2)
            L2 = np.sqrt((1 / (n_x ** 2)) * sum)
            L2_list.append(L2)
        print('L2', L2_list)

        L2_log_list = np.log(L2_list)

        # Norma E_inf
        E_inf_depois_list = []

        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[:, 3]
            y_ex = p_ex[:, 3]
            E_inf_antes_list = []
            n_x = n_x_array[i]
            for k in range(len(y_ex)):
                E_inf_antes = abs((y_ex[k] - y_calc[k]))
                E_inf_antes_list.append(E_inf_antes)
            E_inf_depois = max(E_inf_antes_list)
            E_inf_depois_list.append(E_inf_depois)

        E_inf_depois_log_list = np.log(E_inf_depois_list)

        # Norma E_rel
        err_rel_total_list = []

        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[:, 3]
            y_ex = p_ex[:, 3]
            err_rel_list = []
            sum = 0
            n_x = n_x_array[i]
            for k in range(len(y_ex)):
                err_rel = abs((y_ex[k] - y_calc[k]) / (y_ex[k]))
                err_rel_list.append(err_rel)
            for j in range(len(err_rel_list)):
                sum = sum + err_rel_list[j]
            err_rel_total = 1 / n_x * sum
            err_rel_total_list.append(err_rel_total)

        err_rel_total_log_list = np.log(err_rel_total_list)

        h_x_log_list = []
        # Log de h_t
        for i in range(len(h_x)):
            h_x_novo = (h_x[i]) ** 2
            h_x_novo2 = np.log(h_x_novo)
            h_x_log_list.append(h_x_novo2)
        print('h_t_log', h_x_log_list)

        # Plotagem:
        plt.plot(h_x_log_list, L2_list, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Norma Euclidiana - Norma L2')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [s]')
        plt.ylabel('L2')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 2
        coeffs = np.polyfit(h_x_log_list, L2_log_list,
                            grau)  # ajusta a linha de tendência aos dados e retorna os coeficientes do polinômio
        tendencia = np.poly1d(coeffs)  # cria um polinômio a partir dos coeficientes retornados por polyfit
        x_tendencia = np.linspace(min(h_x_log_list), max(h_x_log_list),
                                  100)  # cria um conjunto de pontos para suavizar a linha de tendência
        plt.plot(x_tendencia, tendencia(x_tendencia), color='green', linestyle='dashed',
                 label='Linha de Tendência Linear')
        plt.plot(h_x_log_list, L2_log_list, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Norma Euclidiana - Norma L2')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [s]')
        plt.ylabel('ln (L2)')
        plt.show()

        # Plotagem:
        plt.plot(h_x_log_list, E_inf_depois_list, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Erro Absoluto Máximo - Norma E$ \infty$')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [s]')
        plt.ylabel('E$ \infty$')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 2
        coeffs = np.polyfit(h_x_log_list, E_inf_depois_log_list,
                            grau)  # ajusta a linha de tendência aos dados e retorna os coeficientes do polinômio
        tendencia = np.poly1d(coeffs)  # cria um polinômio a partir dos coeficientes retornados por polyfit
        x_tendencia = np.linspace(min(h_x_log_list), max(h_x_log_list),
                                  100)  # cria um conjunto de pontos para suavizar a linha de tendência
        plt.plot(x_tendencia, tendencia(x_tendencia), color='green', linestyle='dashed',
                 label='Linha de Tendência Linear')
        plt.plot(h_x_log_list, E_inf_depois_log_list, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Erro Absoluto Máximo - Norma E$ \infty$')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [s]')
        plt.ylabel('ln (E$ \infty$)')
        plt.show()

        # Plotagem:
        plt.plot(h_x_log_list, err_rel_total_list, linestyle='none', marker='o', color="#FF007F",label='Erro Analítica/Explícita')
        plt.title('Erro Relativo - Norma L1')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [s]')
        plt.ylabel('L1')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 2
        coeffs = np.polyfit(h_x_log_list, err_rel_total_log_list,
                            grau)  # ajusta a linha de tendência aos dados e retorna os coeficientes do polinômio
        tendencia = np.poly1d(coeffs)  # cria um polinômio a partir dos coeficientes retornados por polyfit
        x_tendencia = np.linspace(min(h_x_log_list), max(h_x_log_list),
                                  100)  # cria um conjunto de pontos para suavizar a linha de tendência
        plt.plot(x_tendencia, tendencia(x_tendencia), color='green', linestyle='dashed',
                 label='Linha de Tendência Linear')
        plt.plot(h_x_log_list, err_rel_total_log_list, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Erro Relativo - Norma L1')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [s]')
        plt.ylabel('ln (L1)')
        plt.show()

        # Tabelas:
        tabela = PrettyTable(['Steps de Tempo', 'N° de Blocos', 'Norma Euclidiana - L2', 'Erro Absoluto Máximo - E_inf', 'Erro Relativo - Norma L1'])

        for n_x_val, L2_val, E_inf_val, err_rel_val in zip(n_x_array, L2_list, E_inf_depois_list, err_rel_total_list):
            rounded_L2_val = round(L2_val, 3)
            rounded_E_inf_val = round(E_inf_val, 3)
            rounded_err_rel_val = round(err_rel_val, 3)
            rounded_n_x_val = round(n_x_val, 3)
            tabela.add_row([n_t, rounded_n_x_val, rounded_L2_val, rounded_E_inf_val, rounded_err_rel_val])

        print(tabela)








