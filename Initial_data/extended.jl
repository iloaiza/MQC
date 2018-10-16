######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 2; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 1.0; #timestep in atomic units
const tf=5000; #final time for integration in a.u.
Ntrajs=2; #number of trajectories for integration
flags=100; #number of checkpoints throughout integration (each checkpoint saves the current state, i.e. snapshots)

include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_extended_coupling(R)   # w=frequency,h=height,c=coupling
potname="extended"
R0=-10
p0=8.5
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations

DYN_LIST=["EH","FSSH"]

FIRST_RUN=false
