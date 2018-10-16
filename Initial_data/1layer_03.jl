######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 9; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 0.1; #timestep in atomic units
const tf=50; #final time for integration in a.u.
Ntrajs=2; #number of trajectories for integration
flags=100; #number of checkpoints throughout integration (each checkpoint saves the current state, i.e. snapshots)

include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_1layer(R,0.03)
potname="1layer_03"
R0=-24
p0=0
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations

DYN_LIST=["EH","FSSH","CM2","CM2_FSSH"]
#DYN_LIST=["CM2_FSSH"]

FIRST_RUN=false
