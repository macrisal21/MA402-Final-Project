# MA402-Final-Project
This project focuses on the numerical solution of the Black–Scholes partial differential equation (PDE) for pricing a European call option. The equation is given by

$$ \frac{\partial V}{\partial t}

    \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}
    r S \frac{\partial V}{\partial S}

    r V = 0 $$

where V ( S , t ) is the option price, S is the asset price, t is time, σ is volatility, and r is the risk-free rate.

At maturity t = T , the payoff is

V ( S , T ) = max ( S − K , 0 ) ,

where K is the strike price. The boundary conditions are

V ( 0 , t ) = 0 , V ( S max , t ) ≈ S max − K e − r ( T − t ) .

To solve this numerically, we switch to the variable τ = T − t so we can step forward in time. The spatial derivatives are approximated using central finite differences, which gives a tridiagonal system.

For time stepping, we use the backward Euler method, which leads to linear systems of the form

( I − Δ t , L ) V n + 1 = V n .

These systems are assembled as sparse matrices using petsc4py, and solved using PETSc’s KSP solver. In this case, LU factorization is used for a direct solve.

The solution is computed on a grid in space and time, and the final option price is compared to the exact Black–Scholes formula to check accuracy.

Overall, this project combines finite differences with PETSc to solve a common problem in quantitative finance.
