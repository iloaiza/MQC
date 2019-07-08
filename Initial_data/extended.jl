######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 2; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 1.0; #timestep in atomic units
const R_min=[-10.5,10]; #final time for integration in a.u.
Ntrajs=100; #number of trajectories for integration

include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_extended_coupling(R)   # w=frequency,h=height,c=coupling
potname="extended"
R0=-10
K=collect(4:2:31)
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations

#RelTol[3]=1e-6
#AbsTol[3]=1e-7
DYN_LIST=["EH","FSSH"]
#dt_min=0.1
#sanity_checks=false
#sanity_breaks=false
#fixed_step=true
#time_print=true
#fixed_step=true


#FIRST_RUN=false
#initial_dist=wigner; potname="wigner_"*potname;
