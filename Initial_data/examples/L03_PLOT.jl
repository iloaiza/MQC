######### INITIAL DATA FILE FOR MQC DYNAMICS
######### Define global variables and integration parameters
const nsts = 90; #number of electronic states
const mass = 2000 #mass in atomic units
const dt = 0.5; #initial timestep in atomic units
const tf=4000; #final time for integration in a.u.
Ntrajs=2000; #number of trajectories for integration
const flags = 500
include("../include.jl")

######### Define potential and initial conditions
potential(R)=pot_layered(R,0.03)
potname="L03"
R0=-22.5
p0=0
C0=zeros(Complex,nsts)
C0[1]=1
ast0=1; #initial adiabatic state for surface hopping simulations


DYN_LIST=["CMFSH","BO","FSSH","SHEEP","EH","FRIC"]

FIRST_RUN=false

DO_DYN = false
BO_COMP = true
EL_PLOTS = true
HISTO_PLOTS = false
