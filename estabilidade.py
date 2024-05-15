# Análise da Estabilidade para o Método Explícito

import numpy as np
import matplotlib.pyplot as plt
from solvers import solvers
from analitica_dirchlet import analiticas
from FTCS import FTCS
from BTCS import BTCS
from CN import CN
from prettytable import PrettyTable

class estabilidade():
    def calculate_estabilidade_h_t():

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

        h_t = [0.2, 0.1, 0.05, 0.005]
        h_x = 0.5
        j = h_x

        tempos_totais = []

        def calculate_n_t(tf, t0, i):
            n_t = (tf - t0) / (i)
            return n_t

        n_x = (xf - x0) / (h_x)

        rxn = []

        def calculate_h_t_calc_pp(variancia):
            x_calc_array_pp = []
            t_calc_array_pp = []
            p_calc_array_pp = []
            rxn_calc_array_pp = []
            n_t_array = []
            for i in h_t:
                n_t = calculate_n_t(tf, t0, i)
                x_calc_pp, t_calc_pp, p_calc_pp, rxn_calc_pp = FTCS.calculate_FTCS_pp(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array_pp.append(x_calc_pp)
                p_calc_array_pp.append(p_calc_pp)
                t_calc_array_pp.append(t_calc_pp)
                n_t_array.append(n_t)
                rxn_calc_array_pp.append(rxn_calc_pp)


            return x_calc_array_pp, t_calc_array_pp, p_calc_array_pp, rxn_calc_array_pp, n_t_array

        def calculate_h_t_calc_fp(variancia):
            x_calc_array_fp = []
            t_calc_array_fp = []
            p_calc_array_fp = []
            rxn_calc_array_fp = []
            n_t_array = []
            for i in h_t:
                n_t = calculate_n_t(tf, t0, i)
                x_calc_fp, t_calc_fp, p_calc_fp, rxn_calc_fp = FTCS.calculate_FTCS_fp(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array_fp.append(x_calc_fp)
                p_calc_array_fp.append(p_calc_fp)
                t_calc_array_fp.append(t_calc_fp)
                n_t_array.append(n_t)
                rxn_calc_array_fp.append(rxn_calc_fp)

            return x_calc_array_fp, t_calc_array_fp, p_calc_array_fp, rxn_calc_array_fp, n_t_array

        x_calc_array_pp, t_calc_array_pp, p_calc_array_pp, rxn_calc_array_pp, n_t_array = calculate_h_t_calc_pp(variancia)
        rxn.append(rxn_calc_array_pp)
        x_calc_array_fp, t_calc_array_fp, p_calc_array_fp, rxn_calc_array_fp, n_t_array = calculate_h_t_calc_fp(variancia)
        rxn.append(rxn_calc_array_fp)

        conditions = ['Dirchlet', 'Neumann']
        print('rxn', rxn)

        # Tabelas:
        tabela = PrettyTable(['Nº de Blocos', 'Steps de Tempo', 'Condição de Contorno', 'Estabilidade'])
        for cond_idx, condition in enumerate(conditions):
            rx_n = rxn[cond_idx]
            for n_t_idx, delta_t in enumerate(n_t_array):
                rxn_value = rx_n[n_t_idx]
                rounded_value = round(rxn_value, 3)
                rxn_value_str = str(rounded_value)
                rounded_value_n_x = round(n_x, 3)
                rounded_value_n_t = round(delta_t, 3)
                tabela.add_row([rounded_value_n_x, rounded_value_n_t, condition, rxn_value_str])

        print(tabela)

        return rxn_calc_array_pp, rxn_calc_array_fp, n_t_array

    def calculate_estabilidade_h_x():

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

        tempos_totais = []

        def calculate_n_x(xf, x0, j):
            n_x = (xf - x0) / (j)
            return n_x

        n_t = (tf - t0) / (h_t)

        rxn = []

        def calculate_h_x_calc_pp(variancia):
            x_calc_array_pp = []
            t_calc_array_pp = []
            p_calc_array_pp = []
            rxn_calc_array_pp = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_calc_pp, t_calc_pp, p_calc_pp, rxn_calc_pp = FTCS.calculate_FTCS_pp(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array_pp.append(x_calc_pp)
                p_calc_array_pp.append(p_calc_pp)
                t_calc_array_pp.append(t_calc_pp)
                n_x_array.append(n_x)
                rxn_calc_array_pp.append(rxn_calc_pp)


            return x_calc_array_pp, t_calc_array_pp, p_calc_array_pp, rxn_calc_array_pp, n_x_array

        def calculate_h_x_calc_fp(variancia):
            x_calc_array_fp = []
            t_calc_array_fp = []
            p_calc_array_fp = []
            rxn_calc_array_fp = []
            n_x_array = []
            for j in h_x:
                n_x = calculate_n_x(xf, x0, j)
                x_calc_fp, t_calc_fp, p_calc_fp, rxn_calc_fp = FTCS.calculate_FTCS_fp(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x, variancia)

                x_calc_array_fp.append(x_calc_fp)
                p_calc_array_fp.append(p_calc_fp)
                t_calc_array_fp.append(t_calc_fp)
                n_x_array.append(n_x)
                rxn_calc_array_fp.append(rxn_calc_fp)

            return x_calc_array_fp, t_calc_array_fp, p_calc_array_fp, rxn_calc_array_fp, n_x_array

        x_calc_array_pp, t_calc_array_pp, p_calc_array_pp, rxn_calc_array_pp, n_x_array = calculate_h_x_calc_pp(variancia)
        rxn.append(rxn_calc_array_pp)
        x_calc_array_fp, t_calc_array_fp, p_calc_array_fp, rxn_calc_array_fp, n_x_array = calculate_h_x_calc_fp(variancia)
        rxn.append(rxn_calc_array_fp)

        conditions = ['Dirchlet', 'Neumann']
        print('rxn', rxn)

        # Tabelas:
        tabela = PrettyTable(['Steps de Tempo', 'N° de Blocos', 'Condição de Contorno', 'Estabilidade'])
        for cond_idx, condition in enumerate(conditions):
            rx_n = rxn[cond_idx]
            for n_x_idx, delta_x in enumerate(n_x_array):
                rxn_value = rx_n[n_x_idx]
                rounded_value = round(rxn_value, 3)
                rxn_value_str = str(rounded_value)
                rounded_value_n_t = round(n_t, 3)
                rounded_value_n_x = round(delta_x, 3)
                tabela.add_row([rounded_value_n_t, rounded_value_n_x, condition, rxn_value_str])

        print(tabela)

calc_estabilidade = estabilidade.calculate_estabilidade_h_t()
calc_estabilidade = estabilidade.calculate_estabilidade_h_x()

