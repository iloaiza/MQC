abstract type CL_state end #used for classical kind of variables which do not need to carry information on electronic variables
abstract type MF_state end #these types only need to carry one integration for constant initial conditions since they will not have any stochastic elements
abstract type SH_state end

struct C_state #classical state, has classical nuclei information
    R
    p
    NDOFs :: Int #number of nuclear (i.e. classical) degrees of freedom, used for scaling towards higher dimensions
    mem #general variable, can hold anything for memory throughout the dynamics
end

struct Q_state #quantum state, has all electronic information. used even for classical methods
    C
    E
    W
    F
    Ua
    Γ
end

struct ODE_state #carries the ODE parameters for runge-kutta integration
    Rdot
    pdot
    Cdot
    memdot
end

struct CM2_extra
    z
    tnorm #normalized couplings vector
    wvec #epsilons vector
end

struct CM3_extra
    z
    zbar
    tnorm #normalized couplings vector
    tnorm2
    wvec #epsilons vector
    wvec2
end


struct BO_state <: CL_state
    cl :: C_state
    el :: Q_state
    ODE :: ODE_state
    prefix :: String
end

struct EH_state <: MF_state
    cl :: C_state
    el :: Q_state
    ODE :: ODE_state
    prefix :: String
end

struct FSSH_state <: SH_state
    cl :: C_state
    el :: Q_state
    ODE :: ODE_state
    ast :: Int
    prefix :: String
end

struct FSSH_dia_state <: SH_state
    cl :: C_state
    el :: Q_state
    ODE :: ODE_state
    ast :: Int
    prefix :: String
end

struct CM2_state <: MF_state
    cl :: C_state
    el :: Q_state
    ODE :: ODE_state
    CM2 :: CM2_extra
    prefix :: String
end

struct CM3_state <: MF_state
    cl :: C_state
    el :: Q_state
    ODE :: ODE_state
    CM3 :: CM3_extra
    prefix :: String
end

struct CM2_FSSH_state <: SH_state
    cl :: C_state
    el :: Q_state
    ODE :: ODE_state
    CM2 :: CM2_extra
    ast :: Int
    prefix :: String
end

struct CM3_FSSH_state <: SH_state
    cl :: C_state
    el :: Q_state
    ODE :: ODE_state
    CM3 :: CM3_extra
    ast :: Int
    prefix :: String
end


#A CONSTANT NAMED SHEEP_REL MUST BE DEFINED, WITH ARRAYS OF THE STATES THAT GO TOGETHER:
#DEFAULT IS GROUND-STATE - EVERYTHING ELSE MUST BE ADDED HERE
#IT CAN ALSO BE ADDED IN THE INIT FILE, MAKE SURE TO ADD AFTER INCLUDE PART
SHEEP_REL=[[1],collect(2:nsts)]

struct SHEEP_state <: SH_state
    cl :: C_state
    el :: Q_state
    ODE :: ODE_state
    ast :: Int
    prefix :: String
end

struct FRIC_state <: CL_state
    cl :: C_state
    el :: Q_state
    ODE :: ODE_state
    prefix :: String
end

###################################### TYPES META INFORMATION, UPDATE WHEN ADDING NEW METHOD! #####################
#ADD NEW METHOD, SPECIFY IF IT'S CLASSICAL, MEAN-FIELD OR SURFACE HOPPING
#ADD IN TWO POINTS: LIST, NSTS USED

CL_LIST=["BO","FRIC"]
MF_LIST=["EH","CM2","CM3"]
SH_LIST=["FSSH","FSSH_dia","CM2_FSSH","CM3_FSSH","SHEEP"]

EH_sts=nsts
CM2_sts=2
CM3_sts=3
FSSH_sts=nsts
FSSH_dia_sts=nsts
CM2_FSSH_sts=2
CM3_FSSH_sts=3
SHEEP_sts=nsts

#IMPORT METHOD LISTS INTO METHOD_LIST
METHOD_LIST=copy(CL_LIST)
append!(METHOD_LIST,MF_LIST)
append!(METHOD_LIST,SH_LIST)

#ASIGN DYN VARIABLE "DYN" VALUE, USEFUL FOR META PROGRAMING IN DYNAMICAL CODE
for (i,DYN) in enumerate(METHOD_LIST)
    dyn_string="$DYN=METHOD_LIST[$i]"
    eval(Meta.parse(dyn_string))
end

#CREATE STATE BUILDING FUNCTIONS
#####   CL  #####
CL_string="function builder_CL_state(R,p,prefix,NDOFs=length(R),mem=0);"
CL_if_strings=String[]

for DYN in CL_LIST
    push!(CL_if_strings,"""if prefix == $DYN;   S=$(DYN)_state_builder(R,p,NDOFs,mem);  end;    """)
end
CL_string=CL_string*prod(CL_if_strings)*"return S;  end"
eval(Meta.parse(CL_string))

#####   MF  #####
MF_string="function builder_MF_state(R,p,C,prefix,Uold=0,NDOFs=length(R),mem=0);"
MF_if_strings=String[]

for DYN in MF_LIST
    push!(MF_if_strings,"""if prefix == $DYN;   S=$(DYN)_state_builder(R,p,C,Uold,NDOFs,mem);  end;    """)
end
MF_string=MF_string*prod(MF_if_strings)*"return S;  end"
eval(Meta.parse(MF_string))
####    SH  #####
SH_string="function builder_SH_state(R,p,C,ast,prefix,Uold=0,NDOFs=length(R),mem=0);"
SH_if_strings=String[]

for DYN in SH_LIST
    push!(SH_if_strings,"""if prefix == $DYN;   S=$(DYN)_state_builder(R,p,C,ast,Uold,NDOFs,mem);  end;    """)
end
SH_string=SH_string*prod(SH_if_strings)*"return S;  end"
eval(Meta.parse(SH_string))


###################################################################################
#= HOW TO MAKE A NEW TYPE

1: DECIDE THE FAMILY:
    A: CLASSICAL: ONLY NUCLEAR POSITIONS AND MOMENTA ARE TRACKED
    B: MEAN-FIELD: CLASSICAL + ELECTRONIC COEFFICIENTS
    C: SURFACE HOPPING: MEAN-FIELD + CURRENT STATE

2: ADD THE TYPE, MAKING SURE TO MARK THE FAMILY IT BELONGS TO (E.G. STRUCT NEW_MF_state <:MF_state)

3: ADD THE NUMBER OF ELECTRONIC STATES IT USES NAMED AS DYN_sts IF IT'S NOT A CLASSICAL METHOD

4: ADD THE TYPE CONSTRUCTOR IN TYPE_BUILDER, MAKE SURE TO KNOW WHAT'S THE DEFAULT ENTRY FOR FAMILY. FOR THIS, YOU WILL NEED:
    -THE ODE INFORMATION (GOES IN ODE OBJECT!)

5: IF IT'S A SURFACE-HOPPING KIND OF METHOD, ADD A NEW HOP FUNCTION FOR THE NEW METHOD WITH THE NEW CRITERION

=#
