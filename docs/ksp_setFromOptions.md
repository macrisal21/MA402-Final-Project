# `KSP.setFromOptions()`

## Docstring

```python
def setFromOptions(self):
    """
    Configure the KSP solver using PETSc runtime options.

    This method tells PETSc to check the options database and apply any
    solver settings that were provided externally. These options can include
    the solver type, preconditioner type, tolerances, and monitoring options.
    This is useful because the behavior of the solver can be changed without
    editing the Python source code.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    This method does not solve the linear system. It only configures the
    KSP solver before the call to ``ksp.solve(bvec, xvec)``.

    This method should be called after creating the KSP object and setting
    the main solver structure, but before calling ``solve()``.

    Minimal Working Example
    -----------------------
    Basic usage inside a PETSc KSP setup:

    >>> ksp = PETSc.KSP().create()
    >>> ksp.setOperators(A)
    >>> ksp.setType("preonly")
    >>> ksp.getPC().setType("lu")
    >>> ksp.setFromOptions()
    >>> ksp.solve(bvec, xvec)

    Illustrative Example
    --------------------
    In the Black-Scholes solver, ``ksp.setFromOptions()`` allows the user
    to adjust the PETSc solver from the command line without modifying
    ``tutorial_module.py``. For example, solver monitoring can be enabled
    externally when solving the sparse linear systems from the finite
    difference discretization.

    References
    ----------
    PETSc Manual:
    https://petsc.org/release/manualpages/KSP/KSPSetFromOptions/

    PETSc Source:
    https://gitlab.com/petsc/petsc/-/blob/main/src/ksp/ksp/interface/itcl.c
    """
```

## Source Mapping

The `petsc4py` method `ksp.setFromOptions()` is implemented in `petsc4py/src/PETSc/KSP.pyx` as the Python method `setFromOptions()`. In the Cython wrapper, this method calls `CHKERR(KSPSetFromOptions(self.ksp))`, which directly wraps the PETSc C function `KSPSetFromOptions(KSP ksp)`. This function is declared in the PETSc header file `petsc/include/petscksp.h` and implemented in the source file `petsc/src/ksp/ksp/interface/itcl.c`.
