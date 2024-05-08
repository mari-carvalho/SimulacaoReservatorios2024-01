import numpy as np

def jacobi(A, b, x0, Eppara, maxit):
    ne = len(b)
    x = np.zeros(ne) if x0 is None else np.array(x0)

    iter = 0
    Epest = np.linspace(100,100,ne)
    while np.max(Epest) >= Eppara and iter <= maxit:
        x_old = np.copy(x)
        for i in range(ne):
            sum = np.dot(A[i, :i], x_old[:i]) + np.dot(A[i, i+1:], x_old[i + 1:])
            x[i] = (b[i] - sum)/A[i, i]

        Epest = np.abs((x - x_old)/x)*100

        iter += 1

    return x