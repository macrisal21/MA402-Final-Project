# MA402 Final Project

## Summary 
This project focuses on the numerical solution of the Black–Scholes partial differential equation (PDE) for pricing a European call option. The equation is given by

$$
\frac{\partial V}{\partial t} + \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + r S \frac{\partial V}{\partial S} - r V = 0
$$

where $V(S,t)$ is the option price, $S$ is the asset price, $t$ is time, $\sigma$ is volatility, and $r$ is the risk-free rate.

At maturity $t = T$, the payoff is

$$
V(S, T) = \max(S - K, 0)
$$

where $K$ is the strike price.

The boundary conditions are

$$
V(0, t) = 0, \quad
V(S_{\max}, t) \approx S_{\max} - K e^{-r(T - t)}.
$$

To solve this numerically, we switch to the variable $\tau = T - t$ so we can step forward in time. The spatial derivatives are approximated using central finite differences, which gives a tridiagonal system.

For time stepping, we use the backward Euler method, which leads to linear systems of the form

$$
(I - \Delta t \, L) V^{n+1} = V^n
$$

These systems are assembled as sparse matrices using `petsc4py`, and solved using PETSc’s KSP solver. In this case, LU factorization is used for a direct solve.

The solution is computed on a grid in space and time, and the final option price is compared to the exact Black–Scholes formula to check accuracy.

Overall, this project combines finite differences with PETSc to solve a common problem in quantitative finance.

## AI Translation Experience

While the AI provided a useful starting point, the script required significant debugging and modification. I had to fix syntax issues and ensure that the output was correct. The AI helped speed up the initial translation, but a lot of work went in to refining the script.
