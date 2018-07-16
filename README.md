Mixed Quantum Classical (MQC) Dynamics

	-By Ignacio Loaiza G.
	Contact: ignacio.loaiza@mail.utoronto.ca

A Julia implementation of some MQC methods for molecular dynamics.
Includes Ehrenfest, Tully's Fewest Switches Surface Hopping (FSSH) (in both adiabatic and diabatic representations), Collective Modes (with and without surface hopping) and Split Operator (SO) exact dynamics for debugging.


Currently supports parallelization in 1 node. 


	-DYNAMICS.jl is used for a single set of initial conditions under the MQC methods.

	-SPLIT_OPERATOR.jl is used for the SO propagation.

	-K_SIMULATIONS.jl is used for simulating under an array of initial momenta.

	-DATA_TREATMENT.jl is IDE friendly and has all the necessary commands for the treatment of all the generated data.


Check the folders Initial_conditions and SO_Initial_conditions for the initial conditions files. For the implementation of a new potential, it's easier to add it to the potential.jl file (or make a new .jl file, just remember to add it in the include.jl file!).


