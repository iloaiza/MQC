######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 90; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 0.1; #timestep in atomic units
const tf=5000; #final time for integration in a.u.
#Ntrajs=2; #number of trajectories for integration
#flags=100; #number of checkpoints throughout integration (each checkpoint saves the current state, i.e. snapshots)

include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_layered(R,0.03)
potname="layered_03"
R0=-22.5
p0=0
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations

DYN_LIST=["CM2","CM3"]
#SHEEP_REL=[collect(1:10),collect(11:20),collect(21:90)]

#FIRST_RUN=false
