import numpy as np
from scipy.integrate import dblquad



class RayleighRitzMethod:
    '''
    Uses Rayleigh-Ritz method to find beam properties in beam theory. So far the implementation finds frequencies.
    Steps in Rayleigh-Ritz:
    1) Assume w(x,t) = psi(x)q(t).
    2) Choose psi(x), such as (x/a)^2. We base it off of a familiar function we can use to estimate the shape.
       - To get a more accurate estimation, satisfy as many boundary conditions for the 1st function.
       -
    3) "Shove" it into Euler-Lagrange, where we use the oscillation mode to find KE and U
    4) This gives us a neat second-order ODE, which we can use to solve for frequency
    '''

    # The basis functions are in the form (a-x)^(2+n) * x^(2+n). This is so that the boundary conditions for the panel
    # being clamped on all ends are met. The higher the power, the more boundary conditions are satisfied.
    @staticmethod
    def get_basis_functions(a, N):
        phi = []
        for i in range(N):
            phi_i = lambda x: pow(a - x, i + 2) * pow(x, i + 2)
            phi.append(phi_i)
        return phi

    @staticmethod
    def get_second_derivative(a, N):
        d2phi = []
        for k in range(N):
            d2phi_i = lambda x: (k + 2) * (a - x) ** k * x ** k * ((4 * k + 6) * x ** 2 + (-4 * a * k - 6 * a) * x + a ** 2 * k + a ** 2)
            d2phi.append(d2phi_i)
        return d2phi

    @staticmethod
    def solve_eigenvalue_problem(A1, A2):
        A2_inv = np.linalg.inv(A2)
        problem = A2_inv @ A1
        eigenvalues, eigenvectors = np.linalg.eig(problem)
        return eigenvalues, eigenvectors

    def __init__(self, N, L, m, nu):
        self.L = L
        self.phi = RayleighRitzMethod.get_basis_functions(a=self.L, N=N)
        self.d2phi = RayleighRitzMethod.get_second_derivative(a=self.L, N=N)
        self.N = N
        self.m = m
        self.nu = nu

    def __call__(self, E, h):
        # The solution to the Euler-Lagrangian is the determinant of (K - w^2 * M).
        # E and h are used because the panel has more dimensions than a beam, so a beam requires Young's modulus.
        K = self.get_generalized_spring_constant_matrix(E, h)
        M = self.get_generalized_mass_matrix()
        eigenvalues, eigenvectors = RayleighRitzMethod.solve_eigenvalue_problem(A1=K, A2=M)
        frequencies = np.sqrt(eigenvalues)
        return frequencies, eigenvectors

    def get_generalized_mass_matrix(self):
        M = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                integrand = lambda y, x: (self.phi[i](x) * self.phi[i](y)) * (self.phi[j](x) * self.phi[j](y))
                M[i, j], _ = self.m * dblquad(integrand, 0, self.L, 0, self.L)
        return M

    def get_generalized_spring_constant_matrix(self, E, h):
        B = (E * h**3) / (12 * (1 - self.nu**2))
        K = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                K[i, j] =
        return K






