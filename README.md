Mixed Quantum Classical (MQC) Dynamics

	-By Ignacio Loaiza G.
	Contact: ignacio.loaiza@mail.utoronto.ca

A Julia implementation of some MQC methods for molecular dynamics.
Includes Ehrenfest Mean-Field (EH), Tully's Fewest Switches Surface Hopping (FSSH)(1) (in both adiabatic and diabatic representations), Collective Modes (with and without surface hopping) (CM2, CM3, CM2_FSSH and CM3_FSSH)(2,3), Surface Hopping in Excited Ehrenfest Potential (SHEEP)(5), Born-Oppenheimer (BO), and Split Operator (SO) exact dynamics for benchmarking.


Supports parallelization in 1 node.


	-DYNAMICS.jl is used for a single set of initial conditions under the MQC methods. (Still, an initial wigner distribution can be used)

	-SPLIT_OPERATOR.jl is used for the SO propagation.

	-K_SIMULATIONS.jl is used for simulating under an array of initial momenta.

	-SINGLE_TEST.jl and FULL_MEMORY.jl are for doing a single trajectory integration, testing purposes.

	-For data visualization use the Data_treatment Jupyter notebooks, using any of the output .h5 files in data.


Check the folders Initial_conditions and SO_Initial_conditions for the initial conditions files. For the implementation of a new potential, it's easier to add it to the potential.jl file (or make a new .jl file, just remember to add it in the include.jl file!).

References.

(1) - Tully, J. C., Molecular dynamics with electronic transitions. J. Chem. Phys. 93, 1061 (1990)
(2) - Ryabinkin, I. G., and A. F. Izmaylov, Mixed Quantum-Classical Dynamics Using Collective Electronic Variables: A Better Alternative to Electronic Friction Theories. J. Phys. Chem. Lett., 8 (2), pp 440–444 (2017)
