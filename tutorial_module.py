from petsc4py import PETSc
import numpy as np
from math import log, sqrt, exp, erf


def Ncdf(x):
    """
    Standard normal cumulative distribution function
    """
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


def exact_black_scholes_call(S, K, T, r, sigma):
    """
    Computes the exact Black-Scholes price of a European call option
    """
    d1 = (log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return S * Ncdf(d1) - K * exp(-r * T) * Ncdf(d2)


def run_simulation(S_max=400.0, K=100.0, T=1.0, r=0.05, sigma=0.2, N=400, M=400):
    """
    Solve the Black-Scholes PDE for a European call option using
    an implicit finite-difference scheme and PETSc

    Parameters
    ----------
    S_max : float
        Maximum stock price in the spatial grid
    K : float
        Strike price
    T : float
        Time to maturity
    r : float
        Risk-free interest rate
    sigma : float
        Volatility
    N : int
        Number of spatial steps
    M : int
        Number of time steps

    Returns
    -------
    S : numpy.ndarray
        Stock price grid
    V : numpy.ndarray
        Numerical option values at t = 0
    price : float
        Numerical option price at S = K
    exact : float
        Exact Black-Scholes call price at S = K
    error : float
        Absolute error between numerical and exact prices
    """
    dS = S_max / N
    dt = T / M

    S = np.linspace(0.0, S_max, N + 1)

    # Initial condition in tau = 0
    V = np.maximum(S - K, 0.0)

    # Build A = I - dt * L
    A = PETSc.Mat().createAIJ([N + 1, N + 1], nnz=3)
    A.setUp()

    for i in range(1, N):
        Si = S[i]
        alpha = 0.5 * sigma**2 * Si**2
        beta = r * Si

        a = alpha / dS**2 - beta / (2.0 * dS)
        b = -2.0 * alpha / dS**2 - r
        c = alpha / dS**2 + beta / (2.0 * dS)

        A.setValue(i, i - 1, -dt * a)
        A.setValue(i, i,     1.0 - dt * b)
        A.setValue(i, i + 1, -dt * c)

    # Dirichlet boundary rows
    A.setValue(0, 0, 1.0)
    A.setValue(N, N, 1.0)

    A.assemblyBegin()
    A.assemblyEnd()

    # PETSc vectors
    bvec = PETSc.Vec().createSeq(N + 1)
    xvec = PETSc.Vec().createSeq(N + 1)

    # Solver
    ksp = PETSc.KSP().create()
    ksp.setOperators(A)
    ksp.setType("preonly")
    ksp.getPC().setType("lu")
    ksp.setFromOptions()

    # Time stepping
    for n in range(M):
        tau_new = (n + 1) * dt

        bvec.setValues(range(N + 1), V)
        bvec[0] = 0.0
        bvec[N] = S_max - K * np.exp(-r * tau_new)

        bvec.assemblyBegin()
        bvec.assemblyEnd()

        ksp.solve(bvec, xvec)

        reason = ksp.getConvergedReason()
        if reason <= 0:
            raise RuntimeError(f"KSP failed at step {n}, reason = {reason}")

        V = xvec.getArray().copy()

    price = np.interp(K, S, V)
    exact = exact_black_scholes_call(K, K, T, r, sigma)
    error = abs(price - exact)

    return S, V, price, exact, error


if __name__ == "__main__":
    S, V, price, exact, error = run_simulation()
    print("Option price at S=K:", price)
    print("Exact Black-Scholes price:", exact)
    print("Absolute error:", error)