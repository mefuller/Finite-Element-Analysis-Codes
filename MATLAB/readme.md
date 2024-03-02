This code is comprised of the following functions:

BmatdetJ                Computes the B matrix and the determinate of the Jacobian matrix for an element.
Cmat_transvIcoV3b       Computes the C matrix for the constitutive model.
extrapnodesV1_4         Extrapolates computed stress quantities to the element nodes.
inputfile_fracture      An example of an input file that specifies the FEA problem.  A mode II crack tip analysis is shown in this case.
makeMesh                A function that creates a radial mesh with logarithmic element spacing.  Useful for studying the stress around a crack tip.
makeRectMesh            A function that creates a rectangular shaped mesh region with quadrilateral elements.
RK42_delS               Used by the solve to compute the stress increment.
sol_nlFEM_18            The solver function which computes the FEA solution.

The code is run by creating a simple script that calls the sol_nlFEM_18 function and supplies it with inputs (see the Matlab live script file "crack_stress.mlx" as a good example of this implementation).  The an input script in the form of a *.m file can be used to generate the necessary inputs.  See "inputfile_fracture* supplied as an example.
