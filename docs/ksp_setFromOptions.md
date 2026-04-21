# ksp.setFromOptions()

Configure a PETSc Krylov solver using runtime options.

Updates the `KSP` object using options stored in the PETSc options
database. In this project, this function is used in the Black–Scholes
solver after specifying the solver type and preconditioner so that
additional solver settings can be applied without modifying the code.

Parameters
----------
ksp : petsc4py.PETSc.KSP
    Krylov solver object whose settings are being updated.

Returns
-------
None

Notes
-----
The corresponding PETSc C function is
`KSPSetFromOptions(KSP ksp)`.

Declared in
    petsc/include/petscksp.h

Implemented in
    petsc/src/ksp/ksp/interface/itcl.c

PETSc GitHub source
    https://github.com/petsc/petsc/blob/main/src/ksp/ksp/interface/itcl.c

In the Black–Scholes solver, this function is called after setting
`"preonly"` and `"lu"` so that any external PETSc options (such as
tolerances or monitors) can still be applied before solving each
linear system.

Examples
--------
```python
>>> from petsc4py import PETSc
>>> ksp = PETSc.KSP().create()
>>> ksp.setType("preonly")
>>> ksp.getPC().setType("lu")
>>> ksp.setFromOptions()
```
