#include module, very practical
using SharedArrays
using LinearAlgebra
using HDF5
using Distributions

include("general_functions.jl")
include("potentials.jl")
include("types.jl")
include("adiabatic_values.jl")
include("type_builder.jl")
include("distributions.jl")
include("diff_eq.jl")
include("runge_step.jl")
include("hop.jl")
include("integration.jl")
include("plotting.jl")
include("energy_checker.jl")
include("data_functions.jl")
