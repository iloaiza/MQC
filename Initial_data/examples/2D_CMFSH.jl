######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 20; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 0.01; #timestep in atomic units
const tf = 3000; #final time for integration in a.u.
Ntrajs = 10000; #number of trajectories for integration
const flags = 10; #number of checkpoints throughout integration (each checkpoint saves the current state, i.e. snapshots)

include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_2D_sinus(R,0.005)   # w=frequency,h=height,c=coupling
potname="2D_CMFSH"
R0=[-pi,-pi]
p0=[13.6,0]
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations

#DYN_LIST=["EH","FSSH","SHEEP","CM2","CM3","CM2_FSSH","CM3_FSSH","CMFSH","CM3_FSSH_FRIC"]
#DYN_LIST=["FSSH","SHEEP","CM2_FSSH","CM3_FSSH","CMFSH","CM3_FSSH_FRIC"]
#DYN_LIST=["CM2","CM3"]
#DYN_LIST=["EH"]
DYN_LIST=["CMFSH"]
dt_min*=1e-2
#AbsTol*=1e-1
#RelTol*=1e-1

#high_verbose=true
#time_print=true

#SHEEP_REL=[collect(1:10),collect(1:20)]

initial_dist=wigner_vel_2D

FIRST_RUN=false;
DO_DYN = false
