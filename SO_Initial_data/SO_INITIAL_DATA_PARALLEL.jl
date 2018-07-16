######### INITIAL DATA FOR SPLIT OPERATOR DYNAMICS
using HDF5 #package used for saving data

######### Define global variables and integration parameters
@everywhere const nsts=90
@everywhere const mass=2000
@everywhere const δt=0.1
flags=100

@everywhere include("/gpfs/fs0/scratch/i/izmaylov/nacho/Complete_code/general_functions.jl")
@everywhere include("/gpfs/fs0/scratch/i/izmaylov/nacho/Complete_code/potentials.jl")
@everywhere include("/gpfs/fs0/scratch/i/izmaylov/nacho/Complete_code/split_operator_functions_parallel.jl")

######### Define potential and initial conditions
R0=-24
p0=0
C0=zeros(nsts)
C0[1]=1
xmin=0.2;
xmax=30
tf=5000;
pow=13; #2^pow is the number of total points, 2^10=1024, 2^13=8k
@everywhere pot(R)=pot_layered(R,0.015)
potname="layered_015" #will save this in ./data folder, adding SO_ prefix, R0 and p0. This name should have information about the potential and its parameters only


########## Stuff for text output and extra needed stuff
file="./data/SO_"*potname*"_R0($R0)_p0($p0).h5"
@everywhere potential(R)=SO_pot(pot,R);

h5write(file,"R0",R0)
h5write(file,"p0",p0)
h5write(file,"C0",C0)
########### Text output
println("The timestep is")
@show δt
println("The mass is")
@show mass
println("The final time, in a.u. is")
@show tf
println("The number of grid points is")
@show 2^pow
println("The name of the file for the savename is")
@show file
println("The boxsize is")
@show boxsize
