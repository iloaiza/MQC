######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 2; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 1.0; #timestep in atomic units
const tf=5000; #final time for integration in a.u.
Ntrajs=1000; #number of trajectories for integration
flags=100; #number of checkpoints throughout integration (each checkpoint saves the current state, i.e. snapshots)

include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_sinus(R,0.5,0.02,0.01)   # w=frequency,h=height,c=coupling
potname="sinus_01"
R0=-22
p0=12.5
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations
