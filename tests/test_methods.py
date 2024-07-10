#   ---------------------------------------------------------------------------------
#   Copyright (c) Microsoft Corporation. All rights reserved.
#   Licensed under the MIT License. See LICENSE in project root for information.
#   ---------------------------------------------------------------------------------
"""This is a sample python file for testing functions from the source code."""
from __future__ import annotations

from finite_element_analysis_codes.nlFEMPY import *


def test_cylindrical():
    """
    Verify a simple transformation from Cartessian to cylindrical, 2D
    """
    a,b = r_theta(1, 0)
    assert a==1.
    assert b==0.

# write a bunch more tests here...

# Test code
def imported():
    start = timeit.default_timer()

    mesher1 = Mesher()
    # mesher1.set_params([2,2], [[3,3], [3,3]], mesher1.coords_quarterCircle(1), [3], [[4,7,4,5]], surfs=[[0,2], [0,6], [6,7], [5,2]])
    mesher1.set_params([1,1], [[2],[2]], mesher1.coords_Quad(1, 1), surfs=[[0,1], [0,2], [2,3]])
    mesher1.create()

    mesher1.nodes[:,4] = np.array([[.6, .2]])

    mesh1 = Mesh()
    mesh1.make_mesh(mesher1)
    steel = Material_model([30e6, 0.30], "linear elastic, plane strain")
    mesh1.assign_material(steel)
    K = Global_K_matrix(mesh1)
    mesh1.plot()

    E = Strain(mesh1)
    dE = delta_Strain(mesh1)
    S = Stress(mesh1)
    U = Displacement(mesh1)
    T = Global_T_matrix(mesh1)
    F = Global_F_matrix(mesh1)

    topsurf1 = mesher1.surfs[2]
    bottomsurf = mesher1.surfs[0]
    sidesurf = mesher1.surfs[1]

    BC1 = Boundary_condition(K)
    BC1.apply_BC(sidesurf, np.zeros(len(sidesurf)), 'U1')
    BC1.apply_BC(bottomsurf, np.zeros(len(sidesurf)), 'U2')

    F.apply_traction(topsurf1, -10000, 'y')
                
    solution = Standard(K, T, F, BC1, S, E, U, mesh1)
    solution.start()

    stop = timeit.default_timer()
    print('Elapsed time: ', f'{stop-start} seconds.')

    # E.compute(U.return_all())
    # S.compute(U.return_all())

    plot_result(mesh1, S, 'S22', U)
    plot_result(mesh1, S, 'S11', U)
    plot_result(mesh1, S, 'S12', U)

    plot_result(mesh1, U, 'U1', U=U)
    plot_result(mesh1, U, 'U2', U=U)