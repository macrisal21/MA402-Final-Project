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
`KSPSetFromOptions`, for which more information can be found here:
https://petsc.org/main/manualpages/KSP/KSPSetFromOptions/

Examples
--------
```python
>>> from petsc4py import PETSc
>>> ksp = PETSc.KSP().create()
>>> ksp.setType("preonly")
>>> ksp.getPC().setType("lu")
>>> ksp.setFromOptions()
```
