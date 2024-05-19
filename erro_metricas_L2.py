# Função para analisar inclinação das retas do decaimento do Tempo (Linear ou Quadrática)

import numpy as np
import matplotlib.pyplot as plt
from solvers import solvers
from analitica_dirchlet import analiticas
from FTCS import FTCS
from BTCS import BTCS
from CN import CN
from prettytable import PrettyTable

class erros_pp_gs():
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
        variancia = 'tempo'

        h_t = [0.8, 0.7042253521, 0.6024096386, 0.5, 0.4]
        h_x = 0.2
        j = h_x

        def calculate_n_t(tf,t0,i):

            n_t = (tf - t0) / (i)
            return n_t

        n_x = (xf - x0) / (h_x)

        def calculate_h_t_ex(variancia):
            x_ex_array = []
            t_ex_array = []
            p_ex_array = []
            n_t_array = []
            for i in h_t:

                n_t = calculate_n_t(tf, t0, i)
                x_ex, t_ex, p_ex = analiticas.calculate_analitica_dirchlet(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,i,j,n_t,n_x,variancia)

                x_ex_array.append(x_ex)
                p_ex_array.append(p_ex)
                t_ex_array.append(t_ex)
                n_t_array.append(n_t)

            return x_ex_array, t_ex_array, p_ex_array, n_t_array


        x_ex_array, t_ex_array, p_ex_array, n_t_array = calculate_h_t_ex(variancia)


        def calculate_h_t_calc(variancia):
            x_calc_array = []
            t_calc_array = []
            p_calc_array = []
            n_t_array = []
            for i in h_t:

                n_t = calculate_n_t(tf, t0, i)
                x_calc, t_calc, p_calc = BTCS.calculate_BTCS_pp_gs(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array.append(x_calc)
                p_calc_array.append(p_calc)
                t_calc_array.append(t_calc)
                n_t_array.append(n_t)

            return x_calc_array, t_calc_array, p_calc_array, n_t_array


        x_calc_array, t_calc_array, p_calc_array, n_t_array = calculate_h_t_calc(variancia)

        # Cálculo do Erro:

        # Norma L2
        L2_list_tempo_pp_gs = []
        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            sum = 0
            y_calc = p_calc[:, 83] # selecinar a coluna
            y_ex = p_ex[:, 83] # selecinar a coluna
            for k in range(len(y_ex)): # acesso a cada linha
                sum = sum + ((abs((y_ex[k]-y_calc[k])/(y_ex[k])))**2)
            L2 = np.sqrt((1/(n_x**2))*sum)
            L2_list_tempo_pp_gs.append(L2)
        print('L2', L2_list_tempo_pp_gs)

        L2_log_list = np.log(L2_list_tempo_pp_gs)

        # Norma E_inf
        E_inf_depois_list_tempo_pp_gs = []

        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[:, 83]
            y_ex = p_ex[:, 83]
            E_inf_antes_list = []
            for k in range(len(y_ex)):
                E_inf_antes = abs((y_ex[k] - y_calc[k]))
                E_inf_antes_list.append(E_inf_antes)
            E_inf_depois = max(E_inf_antes_list)
            E_inf_depois_list_tempo_pp_gs.append(E_inf_depois)

        E_inf_depois_log_list = np.log(E_inf_depois_list_tempo_pp_gs)

        # Norma E_rel
        err_rel_total_list_tempo_pp_gs = []

        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[:, 83]
            y_ex = p_ex[:, 83]
            err_rel_list = []
            sum = 0
            for k in range(len(y_ex)):
                err_rel = abs((y_ex[k] - y_calc[k])/(y_ex[k]))
                err_rel_list.append(err_rel)
            for j in range(len(err_rel_list)):
                sum = sum + err_rel_list[j]
            err_rel_total = 1/n_x * sum
            err_rel_total_list_tempo_pp_gs.append(err_rel_total)

        err_rel_total_log_list = np.log(err_rel_total_list_tempo_pp_gs)

        h_t_log_list = np.log(h_t)

        # Plotagem:
        plt.plot(h_t, L2_list_tempo_pp_gs, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s])')
        plt.ylabel('ln (L2)')
        plt.show()

        # Plotagem:
        plt.plot(h_t, E_inf_depois_list_tempo_pp_gs, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s])')
        plt.ylabel('ln (E$ \infty$)')
        plt.show()

        # Plotagem:
        plt.plot(h_t, err_rel_total_list_tempo_pp_gs, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s])')
        plt.ylabel('ln (L1)')
        plt.show()

        # Tabelas:
        tabela = PrettyTable(['N° de Blocos', 'Steps de Tempo', 'Norma Euclidiana - L2', 'Erro Absoluto Máximo - E_inf', 'Erro Relativo - Norma L1'])

        for n_t_val, L2_val, E_inf_val, err_rel_val in zip(n_t_array, L2_list_tempo_pp_gs, E_inf_depois_list_tempo_pp_gs, err_rel_total_list_tempo_pp_gs):
            rounded_L2_val = round(L2_val, 6)
            rounded_E_inf_val = round(E_inf_val, 6)
            rounded_err_rel_val = round(err_rel_val, 6)
            rounded_n_t_val = round(n_t_val, 6)
            tabela.add_row([n_x, rounded_n_t_val, rounded_L2_val, rounded_E_inf_val, rounded_err_rel_val])

        print(tabela)

        return L2_list_tempo_pp_gs, E_inf_depois_list_tempo_pp_gs, err_rel_total_list_tempo_pp_gs, n_t_array, n_x

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
        variancia = 'malha'

        h_x = [1, 0.5, 0.25, 0.125, 0.0625, 0.5]
        h_t = 0.1
        i = h_t


        def calculate_n_x(xf, x0, j):

            n_x = (xf - x0) / (j)
            return n_x

        n_t = (tf - t0) / (h_t)

        def calculate_h_t_ex(variancia):
            x_ex_array = []
            t_ex_array = []
            p_ex_array = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_ex, t_ex, p_ex = analiticas.calculate_analitica_dirchlet(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0,
                                                                tf, i, j, n_t, n_x, variancia)

                x_ex_array.append(x_ex)
                p_ex_array.append(p_ex)
                t_ex_array.append(t_ex)
                n_x_array.append(n_x)

            return x_ex_array, t_ex_array, p_ex_array, n_x_array

        x_ex_array, t_ex_array, p_ex_array, n_x_array = calculate_h_t_ex(variancia)

        def calculate_h_t_calc(variancia):
            x_calc_array = []
            t_calc_array = []
            p_calc_array = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_calc, t_calc, p_calc = BTCS.calculate_BTCS_pp_gs(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0,
                                                                tf, i, j, n_t, n_x, variancia)

                x_calc_array.append(x_calc)
                p_calc_array.append(p_calc)
                t_calc_array.append(t_calc)
                n_x_array.append(n_x)

            return x_calc_array, t_calc_array, p_calc_array, n_x_array

        x_calc_array, t_calc_array, p_calc_array, n_x_array = calculate_h_t_calc(variancia)

        # Cálculo do Erro:

        # Norma L2
        L2_list_malha_pp_gs = []
        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            sum = 0
            y_calc = p_calc[904]  # selecinar a coluna
            y_ex = p_ex[904]  # selecinar a coluna
            n_x = n_x_array[i]
            for k in range(len(y_ex)):  # acesso a cada linha
                sum = sum + ((abs((y_ex[k] - y_calc[k]) / (y_ex[k]))) ** 2)
            L2 = np.sqrt((1 / (n_x ** 2)) * sum)
            L2_list_malha_pp_gs.append(L2)
        print('L2', L2_list_malha_pp_gs)

        L2_log_list = np.log(L2_list_malha_pp_gs)

        # Norma E_inf
        E_inf_depois_list_malha_pp_gs = []

        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[904]
            y_ex = p_ex[904]
            E_inf_antes_list = []
            n_x = n_x_array[i]
            for k in range(len(y_ex)):
                E_inf_antes = abs((y_ex[k] - y_calc[k]))
                E_inf_antes_list.append(E_inf_antes)
            E_inf_depois = max(E_inf_antes_list)
            E_inf_depois_list_malha_pp_gs.append(E_inf_depois)

        E_inf_depois_log_list = np.log(E_inf_depois_list_malha_pp_gs)

        # Norma E_rel
        err_rel_total_list_malha_pp_gs = []

        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[904] # mariana
            y_ex = p_ex[904]
            err_rel_list = []
            sum = 0
            n_x = n_x_array[i]
            for k in range(len(y_ex)):
                err_rel = abs((y_ex[k] - y_calc[k]) / (y_ex[k]))
                err_rel_list.append(err_rel)
            for j in range(len(err_rel_list)):
                sum = sum + err_rel_list[j]
            err_rel_total = 1 / n_x * sum
            err_rel_total_list_malha_pp_gs.append(err_rel_total)

        err_rel_total_log_list = np.log(err_rel_total_list_malha_pp_gs)

        h_x_log_list = []
        # Log de h_t
        for i in range(len(h_x)):
            h_x_novo = (h_x[i])
            h_x_novo2 = np.log(h_x_novo)
            h_x_log_list.append(h_x_novo2)
        print('h_t_log', h_x_log_list)

        # Plotagem:
        plt.plot(h_x, L2_list_malha_pp_gs, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Norma Euclidiana - Norma L2')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [m]')
        plt.ylabel('L2')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [s])')
        plt.ylabel('ln (L2)')
        plt.show()

        # Plotagem:
        plt.plot(h_x, E_inf_depois_list_malha_pp_gs, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Erro Absoluto Máximo - Norma E$ \infty$')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [m]')
        plt.ylabel('E$ \infty$')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (E$ \infty$)')
        plt.show()

        # Plotagem:
        plt.plot(h_x, err_rel_total_list_malha_pp_gs, linestyle='none', marker='o', color="#FF007F",label='Erro Analítica/Explícita')
        plt.title('Erro Relativo - Norma L1')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [m]')
        plt.ylabel('L1')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (L1)')
        plt.show()

        # Tabelas:
        tabela = PrettyTable(['Steps de Tempo', 'N° de Blocos', 'Norma Euclidiana - L2', 'Erro Absoluto Máximo - E_inf', 'Erro Relativo - Norma L1'])

        for n_x_val, L2_val, E_inf_val, err_rel_val in zip(n_x_array, L2_list_malha_pp_gs, E_inf_depois_list_malha_pp_gs, err_rel_total_list_malha_pp_gs):
            rounded_L2_val = round(L2_val, 6)
            rounded_E_inf_val = round(E_inf_val, 6)
            rounded_err_rel_val = round(err_rel_val, 6)
            rounded_n_x_val = round(n_x_val, 6)
            tabela.add_row([n_t, rounded_n_x_val, rounded_L2_val, rounded_E_inf_val, rounded_err_rel_val])

        print(tabela)

        return L2_list_malha_pp_gs, E_inf_depois_list_malha_pp_gs, err_rel_total_list_malha_pp_gs, n_x_array, n_t

class erros_fp_gs():
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
        A = 200
        x0 = 0
        xf = L
        t0 = 0
        tf = 100
        variancia = 'tempo'

        h_t = [0.8, 0.7042253521, 0.6024096386, 0.5, 0.4]
        h_x = 0.2
        j = h_x

        def calculate_n_t(tf,t0,i):

            n_t = (tf - t0) / (i)
            return n_t

        n_x = (xf - x0) / (h_x)

        def calculate_h_t_ex(variancia):
            x_ex_array = []
            t_ex_array = []
            p_ex_array = []
            n_t_array = []
            for i in h_t:

                n_t = calculate_n_t(tf, t0, i)
                x_ex, t_ex, p_ex = analiticas.calculate_analitica_neumann(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,i,j,n_t,n_x, variancia)

                x_ex_array.append(x_ex)
                p_ex_array.append(p_ex)
                t_ex_array.append(t_ex)
                n_t_array.append(n_t)

            return x_ex_array, t_ex_array, p_ex_array, n_t_array


        x_ex_array, t_ex_array, p_ex_array, n_t_array = calculate_h_t_ex(variancia)


        def calculate_h_t_calc(variancia):
            x_calc_array = []
            t_calc_array = []
            p_calc_array = []
            n_t_array = []
            for i in h_t:

                n_t = calculate_n_t(tf, t0, i)
                x_calc, t_calc, p_calc = BTCS.calculate_BTCS_fp_gs(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array.append(x_calc)
                p_calc_array.append(p_calc)
                t_calc_array.append(t_calc)
                n_t_array.append(n_t)

            return x_calc_array, t_calc_array, p_calc_array, n_t_array


        x_calc_array, t_calc_array, p_calc_array, n_t_array = calculate_h_t_calc(variancia)

        # Cálculo do Erro:

        # Norma L2
        L2_list_tempo_fp_gs = []
        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            sum = 0
            y_calc = p_calc[:, 12] # selecinar a coluna
            y_ex = p_ex[:, 12] # selecinar a coluna
            for k in range(len(y_ex)): # acesso a cada linha
                sum = sum + ((abs((y_ex[k]-y_calc[k])/(y_ex[k])))**2)
            L2 = np.sqrt((1/(n_x**2))*sum)
            L2_list_tempo_fp_gs.append(L2)
        print('L2', L2_list_tempo_fp_gs)

        L2_log_list = np.log(L2_list_tempo_fp_gs)

        # Norma E_inf
        E_inf_depois_list_tempo_fp_gs = []

        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[:, 12]
            y_ex = p_ex[:, 12]
            E_inf_antes_list = []
            for k in range(len(y_ex)):
                E_inf_antes = abs((y_ex[k] - y_calc[k]))
                E_inf_antes_list.append(E_inf_antes)
            E_inf_depois = max(E_inf_antes_list)
            E_inf_depois_list_tempo_fp_gs.append(E_inf_depois)

        E_inf_depois_log_list = np.log(E_inf_depois_list_tempo_fp_gs)

        # Norma E_rel
        err_rel_total_list_tempo_fp_gs = []

        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[:, 12]
            y_ex = p_ex[:, 12]
            err_rel_list = []
            sum = 0
            for k in range(len(y_ex)):
                err_rel = abs((y_ex[k] - y_calc[k])/(y_ex[k]))
                err_rel_list.append(err_rel)
            for j in range(len(err_rel_list)):
                sum = sum + err_rel_list[j]
            err_rel_total = 1/n_x * sum
            err_rel_total_list_tempo_fp_gs.append(err_rel_total)

        err_rel_total_log_list = np.log(err_rel_total_list_tempo_fp_gs)

        h_t_log_list = np.log(h_t)

        # Plotagem:
        plt.plot(h_t, L2_list_tempo_fp_gs, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.plot(h_t_log_list, E_inf_depois_list_tempo_fp_gs, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
        plt.title('Erro Absoluto Máximo - Norma E$ \infty$')
        plt.legend()
        plt.xlabel(r' ln ($\bigtriangleup t$ [s])')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s])')
        plt.ylabel('ln (E$ \infty$)')
        plt.show()

        # Plotagem:
        plt.plot(h_t_log_list, err_rel_total_list_tempo_fp_gs, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s])')
        plt.ylabel('ln (L1)')
        plt.show()

        # Tabelas:
        tabela = PrettyTable(['N° de Blocos', 'Steps de Tempo', 'Norma Euclidiana - L2', 'Erro Absoluto Máximo - E_inf', 'Erro Relativo - Norma L1'])

        for n_t_val, L2_val, E_inf_val, err_rel_val in zip(n_t_array, L2_list_tempo_fp_gs, E_inf_depois_list_tempo_fp_gs, err_rel_total_list_tempo_fp_gs):
            rounded_L2_val = round(L2_val, 6)
            rounded_E_inf_val = round(E_inf_val, 6)
            rounded_err_rel_val = round(err_rel_val, 6)
            rounded_n_t_val = round(n_t_val, 6)
            tabela.add_row([n_x, rounded_n_t_val, rounded_L2_val, rounded_E_inf_val, rounded_err_rel_val])

        print(tabela)

        return L2_list_tempo_fp_gs, E_inf_depois_list_tempo_fp_gs, err_rel_total_list_tempo_fp_gs, n_t_array, n_x

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
        variancia = 'malha'

        h_x = [1, 0.5, 0.25, 0.125, 0.0625, 0.5]
        h_t = 0.1
        i = h_t

        def calculate_n_x(xf, x0, j):

            n_x = (xf - x0) / (j)
            return n_x

        n_t = (tf - t0) / (h_t)

        def calculate_h_t_ex(variancia):
            x_ex_array = []
            t_ex_array = []
            p_ex_array = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_ex, t_ex, p_ex = analiticas.calculate_analitica_neumann(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0,
                                                                tf, i, j, n_t, n_x, variancia)

                x_ex_array.append(x_ex)
                p_ex_array.append(p_ex)
                t_ex_array.append(t_ex)
                n_x_array.append(n_x)

            return x_ex_array, t_ex_array, p_ex_array, n_x_array

        x_ex_array, t_ex_array, p_ex_array, n_x_array = calculate_h_t_ex(variancia)

        def calculate_h_t_calc(variancia):
            x_calc_array = []
            t_calc_array = []
            p_calc_array = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_calc, t_calc, p_calc = BTCS.calculate_BTCS_fp_gs(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0,
                                                                tf, i, j, n_t, n_x, variancia)

                x_calc_array.append(x_calc)
                p_calc_array.append(p_calc)
                t_calc_array.append(t_calc)
                n_x_array.append(n_x)

            return x_calc_array, t_calc_array, p_calc_array, n_x_array

        x_calc_array, t_calc_array, p_calc_array, n_x_array = calculate_h_t_calc(variancia)

        # Cálculo do Erro:

        # Norma L2
        L2_list_malha_fp_gs = []
        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            sum = 0
            y_calc = p_calc[250]  # selecinar a coluna
            y_ex = p_ex[250]  # selecinar a coluna
            n_x = n_x_array[i]
            for k in range(len(y_ex)):  # acesso a cada linha
                sum = sum + ((abs((y_ex[k] - y_calc[k]) / (y_ex[k]))) ** 2)
            L2 = np.sqrt((1 / (n_x ** 2)) * sum)
            L2_list_malha_fp_gs.append(L2)
        print('L2', L2_list_malha_fp_gs)

        L2_log_list = np.log(L2_list_malha_fp_gs)

        # Norma E_inf
        E_inf_depois_list_malha_fp_gs = []

        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[250]
            y_ex = p_ex[250]
            E_inf_antes_list = []
            n_x = n_x_array[i]
            for k in range(len(y_ex)):
                E_inf_antes = abs((y_ex[k] - y_calc[k]))
                E_inf_antes_list.append(E_inf_antes)
            E_inf_depois = max(E_inf_antes_list)
            E_inf_depois_list_malha_fp_gs.append(E_inf_depois)

        E_inf_depois_log_list = np.log(E_inf_depois_list_malha_fp_gs)

        # Norma E_rel
        err_rel_total_list_malha_fp_gs = []

        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[250]
            y_ex = p_ex[250]
            err_rel_list = []
            sum = 0
            n_x = n_x_array[i]
            for k in range(len(y_ex)):
                err_rel = abs((y_ex[k] - y_calc[k]) / (y_ex[k]))
                err_rel_list.append(err_rel)
            for j in range(len(err_rel_list)):
                sum = sum + err_rel_list[j]
            err_rel_total = 1 / n_x * sum
            err_rel_total_list_malha_fp_gs.append(err_rel_total)

        err_rel_total_log_list = np.log(err_rel_total_list_malha_fp_gs)

        h_x_log_list = []
        # Log de h_t
        for i in range(len(h_x)):
            h_x_novo = (h_x[i]) ** 2
            h_x_novo2 = np.log(h_x_novo)
            h_x_log_list.append(h_x_novo2)
        print('h_t_log', h_x_log_list)

        # Plotagem:
        plt.plot(h_x, L2_list_malha_fp_gs, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Norma Euclidiana - Norma L2')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [m]')
        plt.ylabel('L2')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (L2)')
        plt.show()

        # Plotagem:
        plt.plot(h_x, E_inf_depois_list_malha_fp_gs, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Erro Absoluto Máximo - Norma E$ \infty$')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [m]')
        plt.ylabel('E$ \infty$')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (E$ \infty$)')
        plt.show()

        # Plotagem:
        plt.plot(h_x, err_rel_total_list_malha_fp_gs, linestyle='none', marker='o', color="#FF007F",label='Erro Analítica/Explícita')
        plt.title('Erro Relativo - Norma L1')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [m]')
        plt.ylabel('L1')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (L1)')
        plt.show()

        # Tabelas:
        tabela = PrettyTable(['Steps de Tempo', 'N° de Blocos', 'Norma Euclidiana - L2', 'Erro Absoluto Máximo - E_inf', 'Erro Relativo - Norma L1'])

        for n_x_val, L2_val, E_inf_val, err_rel_val in zip(n_x_array, L2_list_malha_fp_gs, E_inf_depois_list_malha_fp_gs, err_rel_total_list_malha_fp_gs):
            rounded_L2_val = round(L2_val, 6)
            rounded_E_inf_val = round(E_inf_val, 6)
            rounded_err_rel_val = round(err_rel_val, 6)
            rounded_n_x_val = round(n_x_val, 6)
            tabela.add_row([n_t, rounded_n_x_val, rounded_L2_val, rounded_E_inf_val, rounded_err_rel_val])

        print(tabela)

        return L2_list_malha_fp_gs, E_inf_depois_list_malha_fp_gs, err_rel_total_list_malha_fp_gs, n_x_array, n_t

#________________Análise do Solver Scipy________________

class erros_pp_solv():
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
        A = 200
        x0 = 0
        xf = L
        t0 = 0
        tf = 100
        variancia = 'tempo'

        h_t = [0.8, 0.7042253521, 0.6024096386, 0.5, 0.4]
        h_x = 0.2
        j = h_x
        def calculate_n_t(tf,t0,i):

            n_t = (tf - t0) / (i)
            return n_t

        n_x = (xf - x0) / (h_x)

        def calculate_h_t_ex(variancia):
            x_ex_array = []
            t_ex_array = []
            p_ex_array = []
            n_t_array = []
            for i in h_t:

                n_t = calculate_n_t(tf, t0, i)
                x_ex, t_ex, p_ex = analiticas.calculate_analitica_dirchlet(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,i,j,n_t,n_x, variancia)

                x_ex_array.append(x_ex)
                p_ex_array.append(p_ex)
                t_ex_array.append(t_ex)
                n_t_array.append(n_t)

            return x_ex_array, t_ex_array, p_ex_array, n_t_array


        x_ex_array, t_ex_array, p_ex_array, n_t_array = calculate_h_t_ex(variancia)


        def calculate_h_t_calc(variancia):
            x_calc_array = []
            t_calc_array = []
            p_calc_array = []
            n_t_array = []
            for i in h_t:

                n_t = calculate_n_t(tf, t0, i)
                x_calc, t_calc, p_calc = BTCS.calculate_BTCS_pp_solv(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array.append(x_calc)
                p_calc_array.append(p_calc)
                t_calc_array.append(t_calc)
                n_t_array.append(n_t)

            return x_calc_array, t_calc_array, p_calc_array, n_t_array


        x_calc_array, t_calc_array, p_calc_array, n_t_array = calculate_h_t_calc(variancia)

        # Cálculo do Erro:

        # Norma L2
        L2_list_tempo_pp_solv = []
        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            sum = 0
            y_calc = p_calc[:, 12] # selecinar a coluna
            y_ex = p_ex[:, 12] # selecinar a coluna
            for k in range(len(y_ex)): # acesso a cada linha
                sum = sum + ((abs((y_ex[k]-y_calc[k])/(y_ex[k])))**2)
            L2 = np.sqrt((1/(n_x**2))*sum)
            L2_list_tempo_pp_solv.append(L2)
        print('L2', L2_list_tempo_pp_solv)

        L2_log_list = np.log(L2_list_tempo_pp_solv)

        # Norma E_inf
        E_inf_depois_list_tempo_pp_solv = []

        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[:, 12]
            y_ex = p_ex[:, 12]
            E_inf_antes_list = []
            for k in range(len(y_ex)):
                E_inf_antes = abs((y_ex[k] - y_calc[k]))
                E_inf_antes_list.append(E_inf_antes)
            E_inf_depois = max(E_inf_antes_list)
            E_inf_depois_list_tempo_pp_solv.append(E_inf_depois)

        E_inf_depois_log_list = np.log(E_inf_depois_list_tempo_pp_solv)

        # Norma E_rel
        err_rel_total_list_tempo_pp_solv = []

        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[:, 12]
            y_ex = p_ex[:, 12]
            err_rel_list = []
            sum = 0
            for k in range(len(y_ex)):
                err_rel = abs((y_ex[k] - y_calc[k])/(y_ex[k]))
                err_rel_list.append(err_rel)
            for j in range(len(err_rel_list)):
                sum = sum + err_rel_list[j]
            err_rel_total = 1/n_x * sum
            err_rel_total_list_tempo_pp_solv.append(err_rel_total)

        err_rel_total_log_list = np.log(err_rel_total_list_tempo_pp_solv)

        h_t_log_list = np.log(h_t)

        # Plotagem:
        plt.plot(h_t, L2_list_tempo_pp_solv, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s]')
        plt.ylabel('ln (L2)')
        plt.show()

        # Plotagem:
        plt.plot(h_t, E_inf_depois_list_tempo_pp_solv, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s]')
        plt.ylabel('ln (E$ \infty$)')
        plt.show()

        # Plotagem:
        plt.plot(h_t, err_rel_total_list_tempo_pp_solv, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s]')
        plt.ylabel('ln (L1)')
        plt.show()

        # Tabelas:
        tabela = PrettyTable(['N° de Blocos', 'Steps de Tempo', 'Norma Euclidiana - L2', 'Erro Absoluto Máximo - E_inf', 'Erro Relativo - Norma L1'])

        for n_t_val, L2_val, E_inf_val, err_rel_val in zip(n_t_array, L2_list_tempo_pp_solv, E_inf_depois_list_tempo_pp_solv, err_rel_total_list_tempo_pp_solv):
            rounded_L2_val = round(L2_val, 6)
            rounded_E_inf_val = round(E_inf_val, 6)
            rounded_err_rel_val = round(err_rel_val, 6)
            rounded_n_t_val = round(n_t_val, 6)
            tabela.add_row([n_x, rounded_n_t_val, rounded_L2_val, rounded_E_inf_val, rounded_err_rel_val])

        print(tabela)

        return L2_list_tempo_pp_solv, E_inf_depois_list_tempo_pp_solv, err_rel_total_list_tempo_pp_solv, n_t_array, n_x

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
        variancia = 'malha'

        h_x = [2, 1.333333333, 1, 0.8, 0.6666666667, 0.571428571, 0.5]
        h_t = 0.1
        i = h_t

        def calculate_n_x(xf, x0, j):

            n_x = (xf - x0) / (j)
            return n_x

        n_t = (tf - t0) / (h_t)

        def calculate_h_t_ex(variancia):
            x_ex_array = []
            t_ex_array = []
            p_ex_array = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_ex, t_ex, p_ex = analiticas.calculate_analitica_dirchlet(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0,
                                                                tf, i, j, n_t, n_x, variancia)

                x_ex_array.append(x_ex)
                p_ex_array.append(p_ex)
                t_ex_array.append(t_ex)
                n_x_array.append(n_x)

            return x_ex_array, t_ex_array, p_ex_array, n_x_array

        x_ex_array, t_ex_array, p_ex_array, n_x_array = calculate_h_t_ex(variancia)

        def calculate_h_t_calc(variancia):
            x_calc_array = []
            t_calc_array = []
            p_calc_array = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_calc, t_calc, p_calc = BTCS.calculate_BTCS_pp_solv(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0,
                                                                tf, i, j, n_t, n_x, variancia)

                x_calc_array.append(x_calc)
                p_calc_array.append(p_calc)
                t_calc_array.append(t_calc)
                n_x_array.append(n_x)

            return x_calc_array, t_calc_array, p_calc_array, n_x_array

        x_calc_array, t_calc_array, p_calc_array, n_x_array = calculate_h_t_calc(variancia)

        # Cálculo do Erro:

        # Norma L2
        L2_list_malha_pp_solv = []
        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            sum = 0
            y_calc = p_calc[250]  # selecinar a coluna
            y_ex = p_ex[250]  # selecinar a coluna
            n_x = n_x_array[i]
            for k in range(len(y_ex)):  # acesso a cada linha
                sum = sum + ((abs((y_ex[k] - y_calc[k]) / (y_ex[k]))) ** 2)
            L2 = np.sqrt((1 / (n_x ** 2)) * sum)
            L2_list_malha_pp_solv.append(L2)
        print('L2', L2_list_malha_pp_solv)

        L2_log_list = np.log(L2_list_malha_pp_solv)

        # Norma E_inf
        E_inf_depois_list_malha_pp_solv = []

        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[250]
            y_ex = p_ex[250]
            E_inf_antes_list = []
            n_x = n_x_array[i]
            for k in range(len(y_ex)):
                E_inf_antes = abs((y_ex[k] - y_calc[k]))
                E_inf_antes_list.append(E_inf_antes)
            E_inf_depois = max(E_inf_antes_list)
            E_inf_depois_list_malha_pp_solv.append(E_inf_depois)

        E_inf_depois_log_list = np.log(E_inf_depois_list_malha_pp_solv)

        # Norma E_rel
        err_rel_total_list_malha_pp_solv = []

        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[250]
            y_ex = p_ex[250]
            err_rel_list = []
            sum = 0
            n_x = n_x_array[i]
            for k in range(len(y_ex)):
                err_rel = abs((y_ex[k] - y_calc[k]) / (y_ex[k]))
                err_rel_list.append(err_rel)
            for j in range(len(err_rel_list)):
                sum = sum + err_rel_list[j]
            err_rel_total = 1 / n_x * sum
            err_rel_total_list_malha_pp_solv.append(err_rel_total)

        err_rel_total_log_list = np.log(err_rel_total_list_malha_pp_solv)

        h_x_log_list = []
        # Log de h_t
        for i in range(len(h_x)):
            h_x_novo = (h_x[i]) ** 2
            h_x_novo2 = np.log(h_x_novo)
            h_x_log_list.append(h_x_novo2)
        print('h_t_log', h_x_log_list)

        # Plotagem:
        plt.plot(h_x, L2_list_malha_pp_solv, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Norma Euclidiana - Norma L2')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [m]')
        plt.ylabel('L2')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (L2)')
        plt.show()

        # Plotagem:
        plt.plot(h_x, E_inf_depois_list_malha_pp_solv, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Erro Absoluto Máximo - Norma E$ \infty$')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [s]')
        plt.ylabel('E$ \infty$')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (E$ \infty$)')
        plt.show()

        # Plotagem:
        plt.plot(h_x, err_rel_total_list_malha_pp_solv, linestyle='none', marker='o', color="#FF007F",label='Erro Analítica/Explícita')
        plt.title('Erro Relativo - Norma L1')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [s]')
        plt.ylabel('L1')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (L1)')
        plt.show()

        # Tabelas:
        tabela = PrettyTable(['Steps de Tempo', 'N° de Blocos', 'Norma Euclidiana - L2', 'Erro Absoluto Máximo - E_inf', 'Erro Relativo - Norma L1'])

        for n_x_val, L2_val, E_inf_val, err_rel_val in zip(n_x_array, L2_list_malha_pp_solv, E_inf_depois_list_malha_pp_solv, err_rel_total_list_malha_pp_solv):
            rounded_L2_val = round(L2_val, 6)
            rounded_E_inf_val = round(E_inf_val, 6)
            rounded_err_rel_val = round(err_rel_val, 6)
            rounded_n_x_val = round(n_x_val, 6)
            tabela.add_row([n_t, rounded_n_x_val, rounded_L2_val, rounded_E_inf_val, rounded_err_rel_val])

        print(tabela)

        return L2_list_malha_pp_solv, E_inf_depois_list_malha_pp_solv, err_rel_total_list_malha_pp_solv, n_x_array, n_t

class erros_fp_solv():
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
        A = 200
        x0 = 0
        xf = L
        t0 = 0
        tf = 100
        variancia = 'tempo'

        h_t = [0.8, 0.7042253521, 0.6024096386, 0.5, 0.4]
        h_x = 0.2
        j = h_x

        def calculate_n_t(tf,t0,i):

            n_t = (tf - t0) / (i)
            return n_t

        n_x = (xf - x0) / (h_x)

        def calculate_h_t_ex(variancia):
            x_ex_array = []
            t_ex_array = []
            p_ex_array = []
            n_t_array = []
            for i in h_t:

                n_t = calculate_n_t(tf, t0, i)
                x_ex, t_ex, p_ex = analiticas.calculate_analitica_neumann(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,i,j,n_t,n_x, variancia)

                x_ex_array.append(x_ex)
                p_ex_array.append(p_ex)
                t_ex_array.append(t_ex)
                n_t_array.append(n_t)

            return x_ex_array, t_ex_array, p_ex_array, n_t_array


        x_ex_array, t_ex_array, p_ex_array, n_t_array = calculate_h_t_ex(variancia)


        def calculate_h_t_calc(variancia):
            x_calc_array = []
            t_calc_array = []
            p_calc_array = []
            n_t_array = []
            for i in h_t:

                n_t = calculate_n_t(tf, t0, i)
                x_calc, t_calc, p_calc = BTCS.calculate_BTCS_fp_solv(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array.append(x_calc)
                p_calc_array.append(p_calc)
                t_calc_array.append(t_calc)
                n_t_array.append(n_t)

            return x_calc_array, t_calc_array, p_calc_array, n_t_array


        x_calc_array, t_calc_array, p_calc_array, n_t_array = calculate_h_t_calc(variancia)

        # Cálculo do Erro:

        # Norma L2
        L2_list_tempo_fp_solv = []
        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            sum = 0
            y_calc = p_calc[:, 12] # selecinar a coluna
            y_ex = p_ex[:, 12] # selecinar a coluna
            for k in range(len(y_ex)): # acesso a cada linha
                sum = sum + ((abs((y_ex[k]-y_calc[k])/(y_ex[k])))**2)
            L2 = np.sqrt((1/(n_x**2))*sum)
            L2_list_tempo_fp_solv.append(L2)
        print('L2', L2_list_tempo_fp_solv)

        L2_log_list = np.log(L2_list_tempo_fp_solv)

        # Norma E_inf
        E_inf_depois_list_tempo_fp_solv = []

        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[:, 12]
            y_ex = p_ex[:, 12]
            E_inf_antes_list = []
            for k in range(len(y_ex)):
                E_inf_antes = abs((y_ex[k] - y_calc[k]))
                E_inf_antes_list.append(E_inf_antes)
            E_inf_depois = max(E_inf_antes_list)
            E_inf_depois_list_tempo_fp_solv.append(E_inf_depois)

        E_inf_depois_log_list = np.log(E_inf_depois_list_tempo_fp_solv)

        # Norma E_rel
        err_rel_total_list_tempo_fp_solv = []

        for i in range(len(p_calc_array)): # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[:, 12]
            y_ex = p_ex[:, 12]
            err_rel_list = []
            sum = 0
            for k in range(len(y_ex)):
                err_rel = abs((y_ex[k] - y_calc[k])/(y_ex[k]))
                err_rel_list.append(err_rel)
            for j in range(len(err_rel_list)):
                sum = sum + err_rel_list[j]
            err_rel_total = 1/n_x * sum
            err_rel_total_list_tempo_fp_solv.append(err_rel_total)

        err_rel_total_log_list = np.log(err_rel_total_list_tempo_fp_solv)

        h_t_log_list = np.log(h_t)

        # Plotagem:
        plt.plot(h_t, L2_list_tempo_fp_solv, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s])')
        plt.ylabel('ln (L2)')
        plt.show()

        # Plotagem:
        plt.plot(h_t, E_inf_depois_list_tempo_fp_solv, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s])')
        plt.ylabel('ln (E$ \infty$)')
        plt.show()

        # Plotagem:
        plt.plot(h_t, err_rel_total_list_tempo_fp_solv, linestyle='none', marker='o', color="#FF007F", label='Erro Analítica/Explícita')
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
        plt.xlabel(r' ln ($\bigtriangleup t$ [s])')
        plt.ylabel('ln (L1)')
        plt.show()

        # Tabelas:
        tabela = PrettyTable(['N° de Blocos', 'Steps de Tempo', 'Norma Euclidiana - L2', 'Erro Absoluto Máximo - E_inf', 'Erro Relativo - Norma L1'])

        for n_t_val, L2_val, E_inf_val, err_rel_val in zip(n_t_array, L2_list_tempo_fp_solv, E_inf_depois_list_tempo_fp_solv, err_rel_total_list_tempo_fp_solv):
            rounded_L2_val = round(L2_val, 6)
            rounded_E_inf_val = round(E_inf_val, 6)
            rounded_err_rel_val = round(err_rel_val, 6)
            rounded_n_t_val = round(n_t_val, 6)
            tabela.add_row([n_x, rounded_n_t_val, rounded_L2_val, rounded_E_inf_val, rounded_err_rel_val])

        print(tabela)

        return L2_list_tempo_fp_solv, E_inf_depois_list_tempo_fp_solv, err_rel_total_list_tempo_fp_solv, n_t_array, n_x

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
        variancia = 'malha'

        h_x = [2, 1.333333333, 1, 0.8, 0.6666666667, 0.571428571, 0.5]
        h_t = 0.1
        i = h_t

        def calculate_n_x(xf, x0, j):

            n_x = (xf - x0) / (j)
            return n_x

        n_t = (tf - t0) / (h_t)

        def calculate_h_t_ex(variancia):
            x_ex_array = []
            t_ex_array = []
            p_ex_array = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_ex, t_ex, p_ex = analiticas.calculate_analitica_neumann(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0,
                                                                tf, i, j, n_t, n_x, variancia)

                x_ex_array.append(x_ex)
                p_ex_array.append(p_ex)
                t_ex_array.append(t_ex)
                n_x_array.append(n_x)

            return x_ex_array, t_ex_array, p_ex_array, n_x_array

        x_ex_array, t_ex_array, p_ex_array, n_x_array = calculate_h_t_ex(variancia)

        def calculate_h_t_calc(variancia):
            x_calc_array = []
            t_calc_array = []
            p_calc_array = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_calc, t_calc, p_calc = BTCS.calculate_BTCS_fp_solv(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0,
                                                                tf, i, j, n_t, n_x, variancia)

                x_calc_array.append(x_calc)
                p_calc_array.append(p_calc)
                t_calc_array.append(t_calc)
                n_x_array.append(n_x)

            return x_calc_array, t_calc_array, p_calc_array, n_x_array

        x_calc_array, t_calc_array, p_calc_array, n_x_array = calculate_h_t_calc(variancia)

        # Cálculo do Erro:

        # Norma L2
        L2_list_malha_fp_solv = []
        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            sum = 0
            y_calc = p_calc[250]  # selecinar a coluna
            y_ex = p_ex[250]  # selecinar a coluna
            n_x = n_x_array[i]
            for k in range(len(y_ex)):  # acesso a cada linha
                sum = sum + ((abs((y_ex[k] - y_calc[k]) / (y_ex[k]))) ** 2)
            L2 = np.sqrt((1 / (n_x ** 2)) * sum)
            L2_list_malha_fp_solv.append(L2)
        print('L2', L2_list_malha_fp_solv)

        L2_log_list = np.log(L2_list_malha_fp_solv)

        # Norma E_inf
        E_inf_depois_list_malha_fp_solv = []

        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[250]
            y_ex = p_ex[250]
            E_inf_antes_list = []
            n_x = n_x_array[i]
            for k in range(len(y_ex)):
                E_inf_antes = abs((y_ex[k] - y_calc[k]))
                E_inf_antes_list.append(E_inf_antes)
            E_inf_depois = max(E_inf_antes_list)
            E_inf_depois_list_malha_fp_solv.append(E_inf_depois)

        E_inf_depois_log_list = np.log(E_inf_depois_list_malha_fp_solv)

        # Norma E_rel
        err_rel_total_list_malha_fp_solv = []

        for i in range(len(p_calc_array)):  # acesso a matriz menor
            p_calc = p_calc_array[i]
            p_ex = p_ex_array[i]
            y_calc = p_calc[250]
            y_ex = p_ex[250]
            err_rel_list = []
            sum = 0
            n_x = n_x_array[i]
            for k in range(len(y_ex)):
                err_rel = abs((y_ex[k] - y_calc[k]) / (y_ex[k]))
                err_rel_list.append(err_rel)
            for j in range(len(err_rel_list)):
                sum = sum + err_rel_list[j]
            err_rel_total = 1 / n_x * sum
            err_rel_total_list_malha_fp_solv.append(err_rel_total)

        err_rel_total_log_list = np.log(err_rel_total_list_malha_fp_solv)

        h_x_log_list = []
        # Log de h_t
        for i in range(len(h_x)):
            h_x_novo = (h_x[i]) ** 2
            h_x_novo2 = np.log(h_x_novo)
            h_x_log_list.append(h_x_novo2)
        print('h_t_log', h_x_log_list)

        # Plotagem:
        plt.plot(h_x, L2_list_malha_fp_solv, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Norma Euclidiana - Norma L2')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [m]')
        plt.ylabel('L2')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (L2)')
        plt.show()

        # Plotagem:
        plt.plot(h_x, E_inf_depois_list_malha_fp_solv, linestyle='none', marker='o', color="#FF007F",
                 label='Erro Analítica/Explícita')
        plt.title('Erro Absoluto Máximo - Norma E$ \infty$')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [s]')
        plt.ylabel('E$ \infty$')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (E$ \infty$)')
        plt.show()

        # Plotagem:
        plt.plot(h_x, err_rel_total_list_malha_fp_solv, linestyle='none', marker='o', color="#FF007F",label='Erro Analítica/Explícita')
        plt.title('Erro Relativo - Norma L1')
        plt.legend()
        plt.xlabel(r'$\bigtriangleup x$ [s]')
        plt.ylabel('L1')
        plt.show()

        # Ajustando uma linha de tendência nos Gráficos LogxLog:
        grau = 1
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
        plt.xlabel(r' ln ($\bigtriangleup x$ [m])')
        plt.ylabel('ln (L1)')
        plt.show()

        # Tabelas:
        tabela = PrettyTable(['Steps de Tempo', 'N° de Blocos', 'Norma Euclidiana - L2', 'Erro Absoluto Máximo - E_inf', 'Erro Relativo - Norma L1'])

        for n_x_val, L2_val, E_inf_val, err_rel_val in zip(n_x_array, L2_list_malha_fp_solv, E_inf_depois_list_malha_fp_solv, err_rel_total_list_malha_fp_solv):
            rounded_L2_val = round(L2_val, 6)
            rounded_E_inf_val = round(E_inf_val, 6)
            rounded_err_rel_val = round(err_rel_val, 6)
            rounded_n_x_val = round(n_x_val, 6)
            tabela.add_row([n_t, rounded_n_x_val, rounded_L2_val, rounded_E_inf_val, rounded_err_rel_val])

        print(tabela)

        return L2_list_malha_fp_solv, E_inf_depois_list_malha_fp_solv, err_rel_total_list_malha_fp_solv, n_x_array, n_t


#L2_list_tempo_pp_gs, E_inf_depois_list_tempo_pp_gs, err_rel_total_list_tempo_pp_gs, n_t_array, n_x = erros_pp_gs.calculate_erros_tempo()
L2_list_malha_pp_gs, E_inf_depois_list_malha_pp_gs, err_rel_total_list_malha_pp_gs, n_x_array, n_t = erros_pp_gs.calculate_erros_malha()
#L2_list_tempo_fp_gs, E_inf_depois_list_tempo_fp_gs, err_rel_total_list_tempo_fp_gs, n_t_array, n_x = erros_fp_gs.calculate_erros_tempo()
#L2_list_malha_fp_gs, E_inf_depois_list_malha_fp_gs, err_rel_total_list_malha_fp_gs, n_x_array, n_t = erros_fp_gs.calculate_erros_malha()
#L2_list_tempo_pp_solv, E_inf_depois_list_tempo_pp_solv, err_rel_total_list_tempo_pp_solv, n_t_array, n_x = erros_pp_solv.calculate_erros_tempo()
#L2_list_tempo_fp_solv, E_inf_depois_list_tempo_fp_solv, err_rel_total_list_tempo_fp_solv, n_t_array, n_x = erros_fp_solv.calculate_erros_tempo()
#L2_list_malha_pp_solv, E_inf_depois_list_malha_pp_solv, err_rel_total_list_malha_pp_solv, n_x_array, n_t = erros_pp_solv.calculate_erros_malha()
#L2_list_malha_fp_solv, E_inf_depois_list_malha_fp_solv, err_rel_total_list_malha_fp_solv, n_x_array, n_t = erros_fp_solv.calculate_erros_malha()

# Tabelas Tempo:
tabela = PrettyTable(['N° de Blocos', 'Steps de Tempo', 'L2 - Dirchlet', 'L2 - Neumann', 'EAM - Dirchlet', 'EAM - Neumann',
                      'L1 - Dirchlet', 'L1 - Neumann'])

for n_t_val, L2_pp_gs_val, L2_fp_gs_val, E_inf_pp_gs_val,  E_inf_fp_gs_val, err_rel_pp_gs_val, err_rel_fp_gs_val in zip(n_t_array, L2_list_tempo_pp_gs,
                                                   L2_list_tempo_fp_gs, E_inf_depois_list_tempo_pp_gs, E_inf_depois_list_tempo_fp_gs,
                                                   err_rel_total_list_tempo_pp_gs, err_rel_total_list_tempo_fp_gs):
    rounded_L2_pp_gs_val = round(L2_pp_gs_val, 6)
    rounded_E_inf_pp_gs_val = round(E_inf_pp_gs_val, 6)
    rounded_err_rel_pp_gs_val = round(err_rel_pp_gs_val, 6)
    rounded_L2_fp_gs_val = round(L2_fp_gs_val, 6)
    rounded_E_inf_fp_gs_val = round(E_inf_fp_gs_val, 6)
    rounded_err_rel_fp_gs_val= round(err_rel_fp_gs_val, 6)
    rounded_n_t_val = round(n_t_val, 6)
    tabela.add_row([n_x, rounded_n_t_val, rounded_L2_pp_gs_val, rounded_L2_fp_gs_val, rounded_E_inf_pp_gs_val, rounded_E_inf_fp_gs_val, rounded_err_rel_pp_gs_val, rounded_err_rel_fp_gs_val])

print(tabela)

# Tabelas Malha:
tabela = PrettyTable(['Steps de Tempo', 'N° de Blocos', 'L2 - Dirchlet', 'L2 - Neumann', 'EAM - Dirchlet', 'EAM - Neumann',
                      'L1 - Dirchlet', 'L1 - Neumann'])

for n_x_val, L2_pp_gs_val, L2_fp_gs_val, E_inf_pp_gs_val,  E_inf_fp_gs_val, err_rel_pp_gs_val, err_rel_fp_gs_val in zip(n_x_array, L2_list_malha_pp_gs,
                                                   L2_list_malha_fp_gs, E_inf_depois_list_malha_pp_gs, E_inf_depois_list_malha_fp_gs,
                                                   err_rel_total_list_malha_pp_gs, err_rel_total_list_malha_fp_gs):
    rounded_L2_pp_gs_val = round(L2_pp_gs_val, 6)
    rounded_E_inf_pp_gs_val = round(E_inf_pp_gs_val, 6)
    rounded_err_rel_pp_gs_val = round(err_rel_pp_gs_val, 6)
    rounded_L2_fp_gs_val = round(L2_fp_gs_val, 6)
    rounded_E_inf_fp_gs_val = round(E_inf_fp_gs_val, 6)
    rounded_err_rel_fp_gs_val= round(err_rel_fp_gs_val, 6)
    rounded_n_x_val = round(n_x_val, 6)
    tabela.add_row([n_t, rounded_n_x_val, rounded_L2_pp_gs_val, rounded_L2_fp_gs_val, rounded_E_inf_pp_gs_val, rounded_E_inf_fp_gs_val, rounded_err_rel_pp_gs_val, rounded_err_rel_fp_gs_val])

print(tabela)

# Tabelas Tempo:
tabela = PrettyTable(['N° de Blocos', 'Steps de Tempo', 'L2 - Dirchlet - GS', 'L2 - Dirchlet - Scipy', 'L2 - Neumann - GS', 'L2 - Neumann - Scipy'])

for n_t_val, L2_pp_gs_val, L2_pp_solv_val, L2_fp_gs_val, L2_fp_solv_val in zip(n_t_array, L2_list_tempo_pp_gs, L2_list_tempo_pp_solv,
                                                   L2_list_tempo_fp_gs, L2_list_tempo_fp_solv):
    rounded_L2_pp_gs_val = round(L2_pp_gs_val, 6)
    rounded_L2_pp_solv_val = round(L2_pp_solv_val, 6)
    rounded_L2_fp_gs_val = round(L2_fp_gs_val, 6)
    rounded_L2_fp_solv_val = round(L2_fp_solv_val, 6)
    rounded_n_t_val = round(n_t_val, 6)
    tabela.add_row([n_x, rounded_n_t_val, rounded_L2_pp_gs_val, rounded_L2_pp_solv_val, rounded_L2_fp_gs_val, rounded_L2_fp_solv_val])

print(tabela)

# Tabelas Malha:
tabela = PrettyTable(['Steps de Tempo', 'N° de Blocos', 'L2 - Dirchlet - GS', 'L2 - Dirchlet - Scipy', 'L2 - Neumann - GS', 'L2 - Neumann - Scipy'])

for n_x_val, L2_pp_gs_val, L2_pp_solv_val, L2_fp_gs_val, L2_fp_solv_val in zip(n_x_array, L2_list_malha_pp_gs, L2_list_malha_pp_solv,
                                                   L2_list_malha_fp_gs, L2_list_malha_fp_solv):
    rounded_L2_pp_gs_val = round(L2_pp_gs_val, 6)
    rounded_L2_pp_solv_val = round(L2_pp_solv_val, 6)
    rounded_L2_fp_gs_val = round(L2_fp_gs_val, 6)
    rounded_L2_fp_solv_val = round(L2_fp_solv_val, 6)
    rounded_n_x_val = round(n_x_val, 6)
    tabela.add_row([n_t, rounded_n_x_val, rounded_L2_pp_gs_val, rounded_L2_pp_solv_val, rounded_L2_fp_gs_val, rounded_L2_fp_solv_val])

print(tabela)

