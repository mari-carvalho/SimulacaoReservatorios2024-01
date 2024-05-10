import numpy as np
import matplotlib.pyplot as plt
from solvers import solvers
from FTCS import FTCS
from BTCS import BTCS
from CN import CN
from erro_metricas_L2 import erros_pp_gs
from erro_metricas_L2 import erros_fp_gs
from solucao_solvers import tempo_computacional_pp
from solucao_solvers import tempo_computacional_fp
from estabilidade import estabilidade

# Definindo as Variáveis de Entrada:

#calc_FTCS = FTCS.calculate_FTCS_pp()
#calc_erros_tempo = erros_pp_gs.calculate_erros_tempo()
#calc_erros_tempo = erros_fp_gs.calculate_erros_tempo()
#calc_erros_malha = erros.calculate_erros_malha()
calc_tempo_comp = tempo_computacional_fp.calculate_tempo_computacional_h_x()
#calc_tempo_comp = tempo_computacional_fp.calculate_tempo_computacional_h_t()
#calc_estabilidade = estabilidade.calculate_estabilidade()