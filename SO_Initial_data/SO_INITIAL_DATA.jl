######### INITIAL DATA FOR SPLIT OPERATOR DYNAMICS
using HDF5 #package used for saving data

######### Define global variables and integration parameters
const nsts=2
const mass=2000
const δt=1.0
flags=100

include("../potentials.jl")
include("../split_operator_functions.jl")

######### Define potential and initial conditions
R0=-10
p0=8.5
C0=zeros(nsts)
C0[1]=1
boxsize=50;
tf=3500;
pow=13; #2^pow is the number of total points, 2^10=1024, 2^13=8k
pot(R)=pot_extended_coupling(R)
potname="extended_coupling" #will save this in ./data folder, adding SO_ prefix, R0 and p0. This name should have information about the potential and its parameters only


########## Stuff for text output and extra needed stuff
file="./data/SO_"*potname*"_R0($R0)_p0($p0).h5"
potential(R)=SO_pot(pot,R);

h5write(file,"R0",R0)
h5write(file,"p0",p0)
h5write(file,"C0",C0)
