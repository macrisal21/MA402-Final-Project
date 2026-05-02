# MA402 Final Project: Black–Scholes Option Pricing with PETSc

## Summary 
This project solves the Black–Scholes equation for a European call option:

$$
\frac{\partial V}{\partial t} + \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + r S \frac{\partial V}{\partial S} - r V = 0,
\qquad
V(S,T)=\max(S-K,0),
$$

where $V(S,t)$ is the option value, $S$ is the stock price, $K$ is the strike price, $\sigma$ is the volatility, $r$ is the risk-free interest rate, and $T$ is the maturity time. Since the exact Black–Scholes formula is known for a European call option, the numerical solver can be checked by comparing the computed option price with the analytical price.

To solve the PDE numerically, the time variable is transformed using

$$
\tau = T-t,
$$

so that the solver starts from the known payoff at maturity and steps forward in $\tau$. The corresponding asset-price domain is

$$
0 \leq S \leq S_{\max},
$$

with boundary conditions

$$
V(0,t)=0,
\qquad
V(S_{\max},t)\approx S_{\max}-Ke^{-r(T-t)}.
$$

The transformed PDE is discretized using finite differences in the stock-price direction and backward Euler time stepping. This produces a sparse linear system at each time step of the form

$$
A V^{n+1}=b,
$$

where $A$ comes from the finite difference approximation of the Black–Scholes operator and $b$ contains the option values from the previous time step, adjusted for boundary conditions.

The resulting systems are solved using PETSc's KSP solver through `petsc4py`. In this project, `ksp.setFromOptions()` allows solver options to be changed at runtime without editing the code. The line `ksp.getPC().setType("lu")` selects the preconditioner attached to the KSP solver and sets it to LU factorization, so the sparse linear system is solved using a direct factorization method. After the solve is complete, `ksp.getConvergedReason()` is used to check whether the PETSc solver successfully converged or whether the solve failed for some reason.

Together, these three functions help configure, execute, and verify the sparse linear solves that come from the Black–Scholes finite difference discretization.

## Connection to PETSc Functions

The key component of this project is understanding how Python code interacts with its underlying PETSc implementations. In particular, we will examine how petsc4py wraps several PETSc functions written in C and implemented through Cython.

As introduced in the above section, we take a closer look at three `petsc4py` functions used in this solver:
1) `ksp.setFromOptions()`, which wraps the PETSc C function `KSPSetFromOptions()`
2) `ksp.getPC().setType("lu")`, where `ksp.getPC()` wraps `KSPGetPC()` and `pc.setType("lu")` wraps `PCSetType()`
3) `ksp.getConvergedReason()`, which wraps the PETSc C function `KSPGetConvergedReason()`

These functions are critical for configuring and solving the linear systems derived from the finite difference discretization. By analyzing these functions, this project helps bridge the gap between Python code and its underlying PETSc infrastructure.

## AI Translation Experience

For this project, I used ChatGPT to help translate the structure of a PETSc-style C solver into a Python implementation using `petsc4py`. The AI was useful for generating an initial version of the code, especially for setting up PETSc objects such as `Mat`, `Vec`, `KSP`, and `PC`. It also helped identify where the main PETSc solver calls should appear in the workflow.

The AI translation experience was unique, because I was not translating a specific tutorial from the PETSc C tutorials line-by-line into Python. Instead, I used `ex2.c` as a reference for the PETSc KSP workflow and then built my own Black–Scholes solver around that structure. This made the translation process more challenging because the mathematical model, finite difference discretization, boundary conditions, and verification step had to be adapted to my own solver.

So, as one might expect, the AI-generated code needed a lot of debugging. I had to manually fix the Black–Scholes discretization, in particular the finite difference coefficients, the backward Euler time-stepping structure, and the boundary conditions at $S=0$ and $S=S_{\max}$. I also had to debug the PETSc setup, which included the matrix assembly, assigning the operator to the KSP solver, setting the LU preconditioner, and checking the solver convergence reason.

A major part of the experience was learning that translating the code was not just changing syntax from C to Python. The numerical method still had to be understood and verified. In this project, that meant comparing the computed option price with the exact Black–Scholes formula and making sure the PETSc functions were being used in the correct part of the solve process.

Overall, the AI translation was helpful as a starting point, but the final solver required manual debugging and verification against the analytical solution.

**PETSc tutorial reference**: https://petsc.org/main/src/ksp/ksp/tutorials/ex2.c.html

## How to Run

This project can be run either as a standalone Python script or through the Jupyter notebook.

### Option 1: Run the Python Module

1) Create and activate a virtual environment with PETSc and petsc4py installed.

Note: A standard Python virtual environment is not sufficient by itself, since `petsc4py` depends on the PETSc library, which must already be installed and configured.

This project was developed in an Ubuntu environment, where PETSc and `petsc4py` are more reliably supported. Installing `petsc4py` directly on other systems (such as Windows) may fail due to configuration issues.

The virtual environment ensures that Python, PETSc, and `petsc4py` are aligned, preventing installation conflicts.

2) Install remaining dependencies:
```markdown
pip install numpy matplotlib
```

3) Run `python tutorial_module.py`
```markdown
python tutorial_module.py
```

### Option 2: Run the Jupyter Notebook
1) Open `tutorial_presentation.ipynb`
```markdown
jupyter notebook tutorial_presentation.ipynb
```

2) Run all cells to:
- execute the solver
- visualize the numerical vs. exact solution
- explore how the method behaves
