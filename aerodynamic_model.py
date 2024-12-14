import numpy as np
from scipy.integrate import quad

# Wing properties


# Returns the lift, drag, and circulation distribution
def get_lift_and_drag(b, V, alpha, rho, N, y, Lambda):
    # Initialize circulation
    Gamma = np.zeros(N)

    # Iterative solution
    for k in range(100):  # Convergence loop
        w = np.zeros(N)
        for i in range(N):
            def integrand(y_prime):
                if y[i] == y_prime:  # Avoid singularity
                    return 0
                return Gamma[np.abs(y - y_prime).argmin()] / (y[i] - y_prime)

            w[i] = quad(integrand, -b, b)[0] / (4 * np.pi)

        alpha_eff = alpha - w / V
        Gamma_new = np.pi * b * V * alpha_eff / (1 + (Lambda / (np.pi * b)) ** 2)

        if np.linalg.norm(Gamma_new - Gamma) < 1e-6:
            break
        Gamma = Gamma_new

    # Compute lift and drag
    L = np.sum(rho * V * Gamma * (2 * b / N))  # Total lift
    D_ind = np.sum(rho * Gamma * w * (2 * b / N))  # Induced drag

    return L, D_ind, Gamma

