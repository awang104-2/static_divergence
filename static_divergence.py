from rayleigh_ritz import get_K
from aerodynamic_model import get_lift_and_drag
from scipy.linalg import eig
import numpy as np

# Example parameters for Rayleigh-Ritz method
L = 10.0  # Wing length (m)
EI = 1e4  # Bending stiffness (N*m^2)
GJ = 5e3  # Torsional stiffness (N*m^2)
m = 10.0  # Mass per unit length (kg/m)
N = 5

# Example parameters for aerodynamic model
b = 10  # Semi-span
V = 50  # Freestream velocity
alpha = 5 * np.pi / 180  # Geometric angle of attack (radians)
rho = 1.225  # Air density (kg/m^3)
N_panels = 50  # Number of panels
y = np.linspace(-b, b, N)  # Spanwise locations
Lambda = 0  # Sweep angle (radians)

# Find the stiffness and lift, drag, and circulation distribution
K_bending, K_torsion = get_K(N)
lift, drag, gamma = get_lift_and_drag(b=b, V=V, alpha=alpha, rho=rho, N=N_panels, y=y, Lambda=Lambda)

# Combine the stiffness matrices
K = K_bending + K_torsion  # Total structural stiffness matrix
A = np.outer(gamma, gamma) / (2 * b / N_panels)  # Simplified aerodynamic stiffness matrix

# Solve for divergence speed
q_critical = eig(K, A, left=False, right=False).min()  # Smallest eigenvalue gives critical dynamic pressure
V_divergence = np.sqrt(2 * q_critical / rho)  # Calculate divergence speed
