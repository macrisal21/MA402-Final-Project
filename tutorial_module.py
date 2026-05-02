from petsc4py import PETSc
import numpy as np
from math import log, sqrt, exp, erf


def Ncdf(x):
    """
    Computes the cumulative distribution function (CDF) of the standard normal distribution.

    Parameters
    ----------
    x : float
        The input value where the standard normal CDF is evaluated.

    Returns
    -------
    float
        The probability P(Z <= x), where Z is a standard normal random variable.
    """
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


def exact_black_scholes_call(S, K, T, r, sigma):
    """
    Computes the exact Black-Scholes price of a European call option.

    Parameters
    ----------
    S : float
        Current stock price.
    K : float
        Strike price of the option.
    T : float
        Time to maturity.
    r : float
        Risk-free interest rate.
    sigma : float
        Volatility of the underlying asset.

    Returns
    -------
    float
        The exact Black-Scholes price of the European call option.
    """
    if S <= 0:
        return 0.0
    d1 = (log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return S * Ncdf(d1) - K * exp(-r * T) * Ncdf(d2)

def run_simulation(S_max=400.0, K=100.0, T=1.0, r=0.05, sigma=0.2, N=400, M=400):
    """
    Solves the Black-Scholes PDE for a European call option using an implicit
    finite-difference scheme and PETSc's KSP linear solver.

    Parameters
    ----------
    S_max : float, optional
        Maximum stock price used to truncate the computational domain.
    K : float, optional
        Strike price of the option.
    T : float, optional
        Time to maturity.
    r : float, optional
        Risk-free interest rate.
    sigma : float, optional
        Volatility of the underlying asset.
    N : int, optional
        Number of spatial grid intervals.
    M : int, optional
        Number of time steps.

    Returns
    -------
    S : numpy.ndarray
        Grid of stock prices.
    V : numpy.ndarray
        Numerical option values on the stock-price grid at the initial time.
    price : float
        Numerical option price interpolated at S = K.
    exact : float
        Exact Black-Scholes price at S = K.
    error : float
        Absolute error between the numerical and exact prices.
    K : float
        Strike price used in the simulation.
    T : float
        Time to maturity used in the simulation.
    r : float
        Risk-free interest rate used in the simulation.
    sigma : float
        Volatility used in the simulation.
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
    """
    Sets the preconditioner type for the KSP solver to LU factorization.

    Parameters
    ----------
    pc : PETSc.PC
        The preconditioner object associated with the KSP solver.

    Returns
    -------
    None
    """

    ksp.setFromOptions()
    """
    Configures the KSP solver using runtime PETSc options.

    Parameters
    ----------
    ksp : PETSc.KSP
        The Krylov subspace solver object.

    Returns
    -------
    None
    """

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
        """
        Returns the reason why the KSP solver terminated.

        Parameters
        ----------
        ksp : PETSc.KSP
            The Krylov subspace solver object.

        Returns
        -------
        int
            An integer indicating the convergence reason. Positive values indicate
            convergence, while negative values indicate divergence or failure.
        """

        if reason <= 0:
            raise RuntimeError(f"KSP failed at step {n}, reason = {reason}")

        V = xvec.getArray().copy()

    price = np.interp(K, S, V)
    exact = exact_black_scholes_call(K, K, T, r, sigma)
    error = abs(price - exact)

    return S, V, price, exact, error, K, T, r, sigma


if __name__ == "__main__":
    S, V, price, exact, error, K, T, r, sigma = run_simulation()
    print("Option price at S=K:", price)
    print("Exact Black-Scholes price:", exact)
    print("Absolute error:", error)
