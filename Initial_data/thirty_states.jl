######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 30; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 0.01; #timestep in atomic units
const R_min=10; #point to pass for stopping integrations
Ntrajs=100; #number of trajectories for integration
flags=100; #number of checkpoints throughout integration (each checkpoint saves the current state, i.e. snapshots)


include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_thirty_states(R)
potname="thirty_states"
R0=1
K=collect(20:5:50)
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations

DYN_LIST=["EH","FSSH","SHEEP","CM2","CM3","CM2_FSSH","CM3_FSSH"]
#DYN_LIST=["CM2","CM3"]

FIRST_RUN=true         #true by default

#initial_dist=wigner    #constant_dist by default
