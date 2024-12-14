import numpy as np
from scipy.integrate import quad


# Define basis functions
def phi(x, n=1, a=1):
    return pow(x/a, n)


def d2_phi_dx2(x, i):
    """Second derivative of basis function: x^i"""
    if i < 2:
        return 0  # Second derivative of x^i for i < 2 is zero
    return i * (i - 1) * x**(i - 2)


def get_K(L, EI, GJ, m, N):  # N = number of basis functions
    # Construct stiffness and mass matrices
    K_bending = np.zeros((N, N))
    K_torsion = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            # Bending stiffness matrix
            K_bending[i, j] = quad(lambda x: EI * d2_phi_dx2(x, i) * d2_phi_dx2(x, j), 0, L)[0]
            # Torsional stiffness matrix
            K_torsion[i, j] = quad(lambda x: GJ * phi(x, i) * phi(x, j), 0, L)[0]

    return K_bending, K_torsion
