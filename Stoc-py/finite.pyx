cimport numpy as np
import numpy as np

def finite_difference(float r, float sigma, int K, int S_0,
                     int F_max, float dS, float T):
    """ Arbitrary Parameters
    r     : Tasa de interes libre de riesgo
    sigma : Volatilidad
    K     : Precio de ejercicio

    Initial Conditions
    S_0   : Valor de accion en tiempo actual
    S_max : Cota superior de la accion
    dS    : Tamaño de intervalos del precio
    M     : Intervalos del precio
    T     : Tiempo de muestreo
    N     : Intervalos del tiempo
    dT    : Tamaños de intervalos del tiempo
    """

    cdef int S_max = S_0 * F_max
    cdef Py_ssize_t M = int(S_max/dS)
    cdef float dT_tmp = (dS/(sigma*S_max))**2
    cdef Py_ssize_t N = int(np.ceil(T/dT_tmp))
    cdef float dT = T/N

    # Ponderación
    cdef np.ndarray[double] a, b, c, J
    a = np.zeros(M+1, dtype=np.float)
    b = np.zeros(M+1, dtype=np.float)
    c = np.zeros(M+1, dtype=np.float)
    J = np.arange(M+1, dtype=np.float)

    a = (dT/(1+r*dT)) * ((sigma**2 * J**2)/2 - (r * J/2))
    b = (dT/(1+r*dT)) * ((1/dT) - sigma**2 * J**2)
    c = (dT/(1+r*dT)) * ((sigma**2 * J**2)/2 + (r * J/2))

    cdef np.ndarray[double, ndim=2] F = np.zeros((N+1, M+1))
    # Condiciones de frontera

    F[:, 0] = max(-K, 0)
    F[N, :] = (J*dS - K).clip(min=0)
    F[:, M] = max(S_max - K, 0)

    # Llenar Malla
    cdef int i, j
    for i in range(N-1, -1, -1):
        for j in range(1, M):
            F[i, j] = a[j] * F[i+1, j-1] + b[j] * F[i+1, j] \
                      + c[j] * F[i+1, j+1]

    P = int(S_0/dS)
    FF = F[0, P]
    return F, FF
