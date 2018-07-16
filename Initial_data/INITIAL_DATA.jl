######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 9; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 0.1; #timestep in atomic units
const tf=5000; #final time for integration in a.u.
Ntrajs=400; #number of trajectories for integration
flags=100; #number of checkpoints throughout integration (each checkpoint saves the current state, i.e. snapshots)

include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_1layer(R,0.03)   # w=frequency,h=height,c=coupling
potname="1layer_03"
R0=-24
p0=0
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations

D2_0=[1.0+0im,0]
D3_0=[1.0+0im,0,0]

EH_DYN=true
FSSH_DYN=true
FSSH_DIA_DYN=true
CM2_DYN=true
CM3_DYN=true
