Mixed Quantum Classical (MQC) Dynamics

	-By Ignacio Loaiza G.
	Contact: ignacio.loaiza@mail.utoronto.ca

A Julia implementation of MQC methods for molecular dynamics.
Includes Ehrenfest Mean-Field (EH), Tully's Fewest Switches Surface Hopping (FSSH)(1) (in both adiabatic and diabatic representations), Collective Modes (with and without surface hopping) (CM2, CM3)(2), Collective Mode Frictional Surface Hopping (CMFSH)(3), Surface Hopping in Excited Ehrenfest Potential (SHEEP)(4), Born-Oppenheimer (BO), and Split Operator (SO) exact dynamics for benchmarking. Now also includes friction with memory (FRIC)(2,5).


Supports parallelization in shared memory (i.e. 1 node).


	-DYNAMICS.jl is used for running dynamics and/or plotting from the resulting saved data. It must be run with an Initial_data .jl driver, the included examples reproduce the data from (3)single set of initial conditions under the MQC methods. (Still, an initial wigner distribution can be used)

	-SPLIT_OPERATOR.jl is used for the SO propagation with a driver dile in SO_Initial_data.

	-K_SIMULATIONS.jl is used for simulating under an array of initial momenta.

	-SINGLE_TEST.jl and FULL_MEMORY.jl are for doing a single trajectory integration, testing purposes.

	-For data visualization use the Data_treatment Jupyter notebooks, using any of the output .h5 files in data.


Check the folders Initial_conditions and SO_Initial_conditions for the initial conditions files. 


References.

	(1) - J. C. Tully, Molecular dynamics with electronic transitions. J. Chem. Phys. 93, 1061 (1990).

	(2) - I. G. Ryabinkin, and A. F. Izmaylov, Mixed Quantum-Classical Dynamics Using Collective Electronic Variables: A Better Alternative to Electronic Friction Theories. J. Phys. Chem. Lett., 8 (2), pp 440–444 (2017).

	(3) - I. Loaiza, and A. F. Izmaylov, Collective Mode Frictional Surface Hopping: Mixed Quantum-Classical dynamics for molecules on metallic surfaces. under preparation.

	(4) - S. A. Fischer, C. T. Chapman, and X. Li, Surface hopping with Ehrenfest excited potential. J. Chem. Phys., 135, 144102 (2011).

	(5) - M. Head-Gordon, J. C. Tully, Molecular dynamics with electronic frictions. J. Chem. Phys. 103, 23 (1990).
