from analitica_dirchlet import analiticas
from BTCS import BTCS
from FTCS import FTCS
from CN import CN
import matplotlib.pyplot as plt

class comparacao_metodos():
    def calculate_comparacao_metodos():

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
        variancia = 'malha'

        h_x = [0.4, 0.3, 0.2, 0.1]
        h_t = 0.01
        i = h_t

        def calculate_n_x(xf,x0,j):

            n_x = (xf - x0) / (j)
            return n_x

        n_t = (tf - t0) / (h_t)

        def calculate_h_x_ex(variancia):
            x_ex_array_analitica = []
            t_ex_array_analitica = []
            P_ex_array_analitica = []
            n_x_array_analitica = []
            for j in h_x:

                n_x = calculate_n_x(xf, x0, j)
                x_ex_analitica, t_ex_analitica, P_ex_analitica = analiticas.calculate_analitica_dirchlet(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_ex_array_analitica.append(x_ex_analitica)
                P_ex_array_analitica.append(P_ex_analitica)
                t_ex_array_analitica.append(t_ex_analitica)
                n_x_array_analitica.append(n_x)

            return x_ex_array_analitica, t_ex_array_analitica, P_ex_array_analitica, n_x_array_analitica


        def calculate_h_x_calc_ftcs(variancia):
            x_calc_array_ftcs = []
            t_calc_array_ftcs = []
            P_calc_array_ftcs = []
            n_x_array_ftcs = []
            for j in h_x:

                n_x = calculate_n_x(xf, x0, j)
                x_calc_ftcs, t_calc_ftcs, T_calc_ftcs = FTCS.calculate_FTCS_pp(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array_ftcs.append(x_calc_ftcs)
                P_calc_array_ftcs.append(P_calc_ftcs)
                t_calc_array_ftcs.append(t_calc_ftcs)
                n_x_array_ftcs.append(n_x)

            return x_calc_array_ftcs, t_calc_array_ftcs, P_calc_array_ftcs, n_x_array_ftcs


        def calculate_h_x_calc_btcs(variancia):
            x_calc_array_btcs = []
            t_calc_array_btcs= []
            P_calc_array_btcs = []
            n_x_array_btcs = []
            for j in h_x:

                n_x = calculate_n_x(xf, x0, j)
                x_calc_btcs, t_calc_btcs, T_calc_btcs = BTCS.calculate_BTCS_pp(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array_btcs.append(x_calc_btcs)
                P_calc_array_btcs.append(T_calc_btcs)
                t_calc_array_btcs.append(t_calc_btcs)
                n_x_array_btcs.append(n_x)

            return x_calc_array_btcs, t_calc_array_btcs, P_calc_array_btcs, n_x_array_btcs


        def calculate_h_x_calc_cn(variancia):
            x_calc_array_cn = []
            t_calc_array_cn  = []
            P_calc_array_cn  = []
            n_x_array_cn  = []
            for j in h_x:

                n_x = calculate_n_x(xf, x0, j)
                x_calc_cn , t_calc_cn , T_calc_cn  = CN.calculate_CN_pp(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array_cn .append(x_calc_cn )
                P_calc_array_cn .append(T_calc_cn )
                t_calc_array_cn .append(t_calc_cn )
                n_x_array_cn .append(n_x)

            return x_calc_array_cn , t_calc_array_cn , P_calc_array_cn , n_x_array_cn


        x_ex_array_analitica, t_ex_array_analitica, P_ex_array_analitica, n_t_array_analitica = calculate_h_x_ex(
            variancia)
        x_calc_array_ftcs, P_calc_array_ftcs, P_calc_array_ftcs, n_t_array_ftcs = calculate_h_x_calc_ftcs(variancia)
        x_calc_array_btcs, P_calc_array_btcs, P_calc_array_btcs, n_t_array_btcs = calculate_h_x_calc_btcs(variancia)
        x_calc_array_cn, P_calc_array_cn, P_calc_array_cn, n_t_array_cn = calculate_h_x_calc_cn(variancia)

        for i in range(len(P_calc_array_cn)):
            P_ex_analitica = P_ex_array_analitica[i]
            P_calc_ftcs = P_calc_array_ftcs[i]
            P_calc_btcs = P_calc_array_btcs[i]
            P_calc_cn = P_calc_array_cn[i]
            y_ex_analitica = P_ex_analitica[1000]
            y_calc_ftcs = P_calc_ftcs[1000]
            y_calc_btcs = P_calc_btcs[1000]
            y_calc_cn = P_calc_cn[1000]
            x_calc_cn = x_calc_array_cn[i]
            rosa_vibrante = (255 / 255, 20 / 255, 147 / 255)
            roxo_vibrante = (138 / 255, 43 / 255, 226 / 255)
            plt.plot(x_calc_cn, y_ex_analitica, linestyle='dashed', color=rosa_vibrante, label='Analítica')
            plt.plot(x_calc_cn, y_calc_ftcs, linestyle='-', color=roxo_vibrante, label='Explícita')
            plt.plot(x_calc_cn, y_calc_btcs, linestyle='-', color='blue', label='Implícita')
            plt.plot(x_calc_cn, y_calc_cn, linestyle='-', color='green', label='Crank Nicolson')

            plt.legend()
            plt.xlabel('Comprimento [cm]')
            plt.ylabel('Temperatura [°C]')
            plt.title('Comparação dos Métodos de Solução')
            plt.grid()
            plt.show()



calc = comparacao_metodos.calculate_comparacao_metodos()