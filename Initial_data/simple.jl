######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 2; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 1.0; #timestep in atomic units
const R_min=[-10.5,10]; #estimated travel distance in a.u., will be used to determine final time (be careful if the passage through strong NAC zones is repeated!)
Ntrajs=10; #number of trajectories for integration
flags=100; #number of checkpoints throughout integration (each checkpoint saves the current state, i.e. snapshots)

include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_simple_crossing(R)
potname="simple"
R0=-10
K=collect(4:1:30)
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations

DYN_LIST=["EH","FSSH"]

#FIRST_RUN=false
#initial_dist=wigner; potname="wigner_"*potname;
