using HDF5 #package used for saving data
using Distributions #used for initial conditions distributions (such as wigner)

input_name=ARGS[1] #when running, use julia SINGLE_TEST.jl NAME_OF_INITIAL_FILE (without .jl, must be inside Initial_data directory)
#input_name="1layer_03"

BO_DYN=false
EH_DYN=false
CM2_VANILLA_DYN=false
CM3_VANILLA_DYN=false
CM2_DYN=false
CM3_DYN=false
FIRST_RUN=true

println("Including function modules")
include("Initial_data/"*input_name*".jl")

########## Stuff for text output (i.e. do a logfile!)
S_EH=EH_state(R0,p0,C0)
file="./data/SINGLE_"*potname*"_R0($R0)_p0($p0).h5"

Ep0,~,~,~,~,~=adiabatic_values(R0);
E0=real(C0'*diagm(Ep0)*C0)[1]+p0^2/2/mass;

if FIRST_RUN
    h5write(file,"R0",R0)
    h5write(file,"p0",p0)
    h5write(file,"C0_real",real.(C0))
    h5write(file,"C0_imag",imag.(C0))
    h5write(file,"E0",E0)
end

########### Text output
println("The timestep is")
@show dt
println("The mass is")
@show mass
println("The final time, in a.u. is")
@show tf
println("The name of the file for the savename is")
@show file
println("The initial conditions are")
@show R0
@show p0
@show C0
@show E0
println("The dynamics that will be done are:")
if BO_DYN
    println("Born-Oppenheimer")
    S_BO=BO_state(R0,p0)
end
if EH_DYN
    println("Ehrenfest")
end
if CM2_VANILLA_DYN
    println("Two collective modes - vanilla version")
    S_CM2_VANILLA=CM2_VANILLA_state(R0,p0,D2_0)
end
if CM3_VANILLA_DYN
    println("Three collective modes - vanilla version")
    S_CM3_VANILLA=CM3_VANILLA_state(R0,p0,D3_0)
end
if CM2_DYN
    println("Two collective modes")
    S_CM2=CM2_state(R0,p0,D2_0)
end
if CM3_DYN
    println("Three collective modes")
    S_CM3=CM3_state(R0,p0,D3_0)
end

###########    Born-Oppenheimer    ############
if BO_DYN
    println("BEGINNING BO INTEGRATION")
    te0=time()
    T_bo,R_bo,p_bo=single_integration(tf,S_BO);   #for single trajectory integration
    #EH_sanity_check(Rvec_e,pvec_e,C_e);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
    tef=time()
    println("FINISHED BO INTEGRATION FOR $(round(tef-te0,3))s, SAVING...")
    BO_save(T_bo,R_bo,p_bo,file)
end


###########    EHRENFEST    ############
if EH_DYN
    println("BEGINNING EHRENFEST INTEGRATION")
    te0=time()
    T_eh,R_eh,p_eh,C_eh=single_integration(tf,S_EH);   #for single trajectory integration
    #EH_sanity_check(Rvec_e,pvec_e,C_e);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
    tef=time()
    println("FINISHED EHRENFEST INTEGRATION FOR $(round(tef-te0,3))s, SAVING...")
    EH_save(T_eh,R_eh,p_eh,C_eh,file)
end


###########    CM2_vanilla    ############
if CM2_VANILLA_DYN
    println("BEGINNING CM2 VANILLA INTEGRATION")
    te0=time()
    T_cm2_vanilla,R_cm2_vanilla,p_cm2_vanilla,D_cm2_vanilla=single_integration(tf,S_CM2_VANILLA);   #for single trajectory integration
    #CM2_sanity_check(R_cm2,p_cm2,D_cm2);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
    tef=time()
    println("FINISHED CM2 VANILLA INTEGRATION FOR $(round(tef-te0,3))s, SAVING...")
    CM2_VANILLA_save(T_cm2_vanilla,R_cm2_vanilla,p_cm2_vanilla,D_cm2_vanilla,file)
end


###########    CM3_vanilla    ############
if CM3_VANILLA_DYN
    println("BEGINNING CM3 VANILLA INTEGRATION")
    te0=time()
    T_cm3_vanilla,R_cm3_vanilla,p_cm3_vanilla,D_cm3_vanilla=single_integration(tf,S_CM3_VANILLA);   #for single trajectory integration
    #CM3_sanity_check(R_cm2,p_cm2,D_cm2);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
    tef=time()
    println("FINISHED CM3 VANILLA INTEGRATION FOR $(round(tef-te0,3))s, SAVING...")
    CM3_VANILLA_save(T_cm3_vanilla,R_cm3_vanilla,p_cm3_vanilla,D_cm3_vanilla,file)
end

###########    CM2    ############
if CM2_DYN
    println("BEGINNING CM2 INTEGRATION")
    te0=time()
    T_cm2,R_cm2,p_cm2,D_cm2=single_integration(tf,S_CM2);   #for single trajectory integration
    #CM2_sanity_check(R_cm2,p_cm2,D_cm2);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
    tef=time()
    println("FINISHED CM2 INTEGRATION FOR $(round(tef-te0,3))s, SAVING...")
    CM2_save(T_cm2,R_cm2,p_cm2,D_cm2,file)
end


###########    CM3    ############
if CM3_DYN
    println("BEGINNING CM3 INTEGRATION")
    te0=time()
    T_cm3,R_cm3,p_cm3,D_cm3=single_integration(tf,S_CM3);   #for single trajectory integration
    #CM3_sanity_check(R_cm3,p_cm3,D_cm3);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
    tef=time()
    println("FINISHED CM3 INTEGRATION FOR $(round(tef-te0,3))s, SAVING...")
    CM3_save(T_cm3,R_cm3,p_cm3,D_cm3,file)
end


println("FINISHED SUCCESSFULLY!")
