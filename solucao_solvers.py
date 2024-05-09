# Análise do Tempo Computacional para cada Solver

import numpy as np
import matplotlib.pyplot as plt
from solvers import solvers
from analitica_dirchlet import calculate_analitica_dirchlet
from FTCS import FTCS
from BTCS import BTCS
from CN import CN
from prettytable import PrettyTable
import time

class tempo_computacional():
    def calculate_tempo_computacional():

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

        h_t = [0.2, 0.1]
        h_x = 0.8
        j = h_x

        tempos_totais = []

        def calculate_n_t(tf, t0, i):
            n_t = (tf - t0) / (i)
            return n_t


        n_x = (xf - x0) / (h_x)

        def calculate_h_t_calc_gs():
            x_calc_array_pp_gs = []
            t_calc_array_pp_gs = []
            p_calc_array_pp_gs = []
            tempo_total_array_pp_gs = []
            for i in h_t:
                inicio_gs = time.time()
                n_t = calculate_n_t(tf, t0, i)
                x_calc_pp_gs, t_calc_pp_gs, p_calc_pp_gs = BTCS.calculate_BTCS_pp_gs(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x)
                fim_gs = time.time()
                tempo_total_gs = fim_gs - inicio_gs

                x_calc_array_pp_gs.append(x_calc_pp_gs)
                p_calc_array_pp_gs.append(p_calc_pp_gs)
                t_calc_array_pp_gs.append(t_calc_pp_gs)

                tempo_total_array_pp_gs.append(tempo_total_gs)
            return x_calc_array_pp_gs, t_calc_array_pp_gs, p_calc_array_pp_gs, tempo_total_array_pp_gs

        def calculate_h_t_calc_tdma():
            x_calc_array_pp_tdma = []
            t_calc_array_pp_tdma = []
            p_calc_array_pp_tdma = []
            tempo_total_array_pp_tdma = []
            for i in h_t:
                inicio_tdma = time.time()
                n_t = calculate_n_t(tf, t0, i)
                x_calc_pp_tdma, t_calc_pp_tdma, p_calc_pp_tdma = BTCS.calculate_BTCS_pp_tdma(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x)
                fim_tdma = time.time()
                tempo_total_tdma = fim_tdma - inicio_tdma

                x_calc_array_pp_tdma.append(x_calc_pp_tdma)
                p_calc_array_pp_tdma.append(p_calc_pp_tdma)

                tempo_total_array_pp_tdma.append(tempo_total_tdma)
            return x_calc_array_pp_tdma, t_calc_array_pp_tdma, p_calc_array_pp_tdma, tempo_total_array_pp_tdma

        def calculate_h_t_calc_jac():
            x_calc_array_pp_jac = []
            t_calc_array_pp_jac = []
            p_calc_array_pp_jac = []
            tempo_total_array_pp_jac= []
            for i in h_t:
                inicio_jac = time.time()
                n_t = calculate_n_t(tf, t0, i)
                x_calc_pp_jac, t_calc_pp_jac, p_calc_pp_jac = BTCS.calculate_BTCS_pp_jac(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x)
                fim_jac = time.time()
                tempo_total_jac = fim_jac - inicio_jac

                x_calc_array_pp_jac.append(x_calc_pp_jac)
                p_calc_array_pp_jac.append(p_calc_pp_jac)
                t_calc_array_pp_jac.append(t_calc_pp_jac)

                tempo_total_array_pp_jac.append(tempo_total_jac)
            return x_calc_array_pp_jac, t_calc_array_pp_jac, p_calc_array_pp_jac, tempo_total_array_pp_jac

        def calculate_h_t_calc_gsr():
            x_calc_array_pp_gsr = []
            t_calc_array_pp_gsr= []
            p_calc_array_pp_gsr= []
            tempo_total_array_pp_gsr = []
            for i in h_t:
                inicio_gsr = time.time()
                n_t = calculate_n_t(tf, t0, i)
                x_calc_pp_gsr, t_calc_pp_gsr, p_calc_pp_gsr = BTCS.calculate_BTCS_pp_gsr(p0, pw, qw, q0, cc, mi, k, h, phi, c, L, A, x0, xf, t0, tf, i, j, n_t, n_x)
                fim_jac = time.time()
                tempo_total_gsr = fim_jac - inicio_gsr

                x_calc_array_pp_gsr.append(x_calc_pp_gsr)
                p_calc_array_pp_gsr.append(p_calc_pp_gsr)
                t_calc_array_pp_gsr.append(t_calc_pp_gsr)

                tempo_total_array_pp_gsr.append(tempo_total_gsr)
            return x_calc_array_pp_gsr, t_calc_array_pp_gsr, p_calc_array_pp_gsr, tempo_total_array_pp_gsr

        x_calc_array_pp_gs, t_calc_array_pp_gs, p_calc_array_pp_gs, tempo_total_array_pp_gs = calculate_h_t_calc_gs()
        tempos_totais.append(tempo_total_array_pp_gs)
        x_calc_array_pp_tdma, t_calc_array_pp_tdma, p_calc_array_pp_tdma, tempo_total_array_pp_tdma = calculate_h_t_calc_tdma()
        tempos_totais.append(tempo_total_array_pp_tdma)
        x_calc_array_pp_jac, t_calc_array_pp_jac, p_calc_array_pp_jac, tempo_total_array_pp_jac = calculate_h_t_calc_jac()
        tempos_totais.append(tempo_total_array_pp_jac)
        x_calc_array_pp_gsr, t_calc_array_pp_gsr, p_calc_array_pp_gsr, tempo_total_array_pp_gsr = calculate_h_t_calc_gsr()
        tempos_totais.append(tempo_total_array_pp_gsr)

        solvers = ['Guass Seidel', 'TDMA', 'Jacobi', 'Guass Seidel Relaxamento']
        print('tempos', tempos_totais)

        # Tabelas:
        tabela = PrettyTable(['Solver', 'delta t', 'Tempo Computacional [s]'])
        for solver_idx, solver in enumerate(solvers): # enumerate adiciona um contador à iteração
            tempos_solver = tempos_totais[solver_idx] # obtida a lista de tempos correspondente ao solver da iteração
            for h_t_idx, delta_t in enumerate(h_t): # enumerate adiciona um contador à iteração
                tempo_total = tempos_solver[h_t_idx] # pega o tempo correspondente ao delta t
                rounded_value = round(tempo_total, 3)
                tempo_total_str = str(rounded_value)
                tabela.add_row([solver, delta_t, tempo_total_str])

        print(tabela)