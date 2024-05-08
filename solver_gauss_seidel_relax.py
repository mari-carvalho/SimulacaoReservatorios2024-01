import numpy as np

def gauss_seidel_relax(A, b, x0, Eppara, maxit, Lambda):
    ne = len(b)
    x = np.zeros(ne) if x0 is None else np.array(x0)

    iter = 0
    Epest = np.linspace(100,100,ne)

    while np.max(Epest) >= Eppara and iter <= maxit:
        x_old = np.copy(x)

        for i in range(ne):
            sum1 = np.dot(A[i, :i], x[:i])
            sum2 = np.dot(A[i, i+1:], x_old[i +1:])
            x[i] = (b[i] - sum1 - sum2) / A[i, i]

        Epest = np.abs((x - x_old) / x) * 100

        iter += 1

        # Relaxamento
        x = Lambda*x + (1-Lambda)*x_old

    return x
