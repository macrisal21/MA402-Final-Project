# `KSP.getConvergedReason()`

## Docstring

```python
def getConvergedReason(self):
    """
    Return the reason why the KSP solver terminated.

    This method is called after ``ksp.solve(bvec, xvec)`` and returns a
    convergence reason. Positive values indicate convergence, while
    negative values indicate divergence or failure.

    Parameters
    ----------
    None

    Returns
    -------
    int
        An integer indicating the convergence reason.

    Notes
    -----
    This method does not solve the linear system. It only reports what
    happened during the most recent call to ``ksp.solve()``.

    Minimal Working Example
    -----------------------
    Basic usage after a PETSc KSP solve:

    >>> ksp.solve(bvec, xvec)
    >>> reason = ksp.getConvergedReason()
    >>> if reason <= 0:
    ...     raise RuntimeError(f"KSP failed, reason = {reason}")

    Illustrative Example
    --------------------
    In the Black-Scholes solver, ``ksp.getConvergedReason()`` is used after
    each linear solve to verify that PETSc successfully solved the sparse
    system before updating the option-value vector.

    References
    ----------
    PETSc Manual:
    https://petsc.org/release/manualpages/KSP/KSPGetConvergedReason/

    PETSc Source:
    https://gitlab.com/petsc/petsc/-/blob/main/src/ksp/ksp/interface/itfunc.c
    """
```

## Source Mapping

The `petsc4py` method `ksp.getConvergedReason()` is implemented in `petsc4py/src/PETSc/KSP.pyx` as the Python method `getConvergedReason()`. In the Cython wrapper, this method calls `CHKERR(KSPGetConvergedReason(self.ksp, &reason))`, which directly wraps the PETSc C function `KSPGetConvergedReason(KSP ksp, KSPConvergedReason *reason)`. This function is declared in the PETSc header file `petsc/include/petscksp.h` and implemented in the source file `petsc/src/ksp/ksp/interface/itfunc.c`.
