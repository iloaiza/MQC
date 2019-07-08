######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 30; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 0.01; #timestep in atomic units
const R_min=8; #point to pass for stopping integrations
Ntrajs=1000; #number of trajectories for integration
flags=100; #number of checkpoints throughout integration (each checkpoint saves the current state, i.e. snapshots)
tf=10000

include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_thirty_states(R)
potname="thirty_states"
R0=1
K=collect(20:5:50)
p0=20
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations



#DYN_LIST=["CMFSH"]
DYN_LIST=["EH","FSSH","SHEEP","CMFSH"]
#DYN_LIST=["EH"]
#RelTol[3]*=1e-1
#AbsTol[3]*=1e-1
#dt_min=1e-4

high_verbose=true
time_print=true
FIRST_RUN=false         #true by default

#initial_dist=wigner    #constant_dist by default
