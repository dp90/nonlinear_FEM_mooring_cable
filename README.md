# nonlinear_FEM_mooring_cable

Application of nonlinear FEM (rod elements) to a mooring cable (so geometrical nonlinearities) loaded by water, 
ship movements and gravity. A schematization of the system is available in the PNG-file Schematization_mooring_cable.

The Matlab files contain the following functionality:
1. Main.m: sets parameters and creates plots.
2. get_initial_conditions.m finds a first guess to the numerical solver with.
3. compute_stiffness_K.m computes the stiffness matrix K, representing the stiffness of the elements, which is updated
    at every iteration of the numerical solver to account for the 'current' geometry of the cable.
4. compute_external_forces.m computes the external forces at every iteration of the numerical solver to account for the
    'current' geometry of the cable.
5. newton_raphson.m is the numerical solver that iterates until the difference in position between two subsequent 
    iterations gets below a given tolerance. 
    


