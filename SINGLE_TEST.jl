using HDF5 #package used for saving data
using Distributions #used for initial conditions distributions (such as wigner)

input_name=ARGS[1] #when running, use julia SINGLE_TEST.jl NAME_OF_INITIAL_FILE (without .jl, must be inside Initial_data directory)
#input_name="extended"

BO_DYN=false
EH_DYN=false
CM2_VANILLA_DYN=false
CM3_VANILLA_DYN=false
CM2_DYN=false
CM3_DYN=false
FIRST_RUN=true
ast0=1 #default starting state is ground state

println("Including function modules")
include("Initial_data/"*input_name*".jl")

########## Stuff for text output (i.e. do a logfile!)
file="./data/SINGLE_"*potname*"_R0($R0)_p0($p0).h5"

S_EH=EH_state_builder(R0,p0,C0)
Ep0,~,~,~,~,~=adiabatic_values(R0);
E0=real(C0'*Diagonal(Ep0)*C0)[1]+p0^2/2/mass;

if FIRST_RUN
    h5write(file,"R0",R0)
    h5write(file,"p0",p0)
    h5write(file,"C0_real",real.(C0))
    h5write(file,"C0_imag",imag.(C0))
    h5write(file,"E0",E0)
    h5write(file,"DYN_LIST",DYN_LIST)
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
for DYN in DYN_LIST
    dyn_sts=eval(Meta.parse("$(DYN)_sts"))
    if DYN in CL_LIST
        println(DYN)
        S_string="S_$(DYN)=$(DYN)_state_builder($R0,$p0)"
        eval(Meta.parse(S_string))
    elseif DYN in MF_LIST
        println(DYN)
        if dyn_sts!=nsts
            println("Warning! Still no implementation for $C0, using [1,0,...] instead")
            D0=zeros(Complex,dyn_sts)
            D0[1]=1
        else
            D0=C0
        end
        S_string="S_$(DYN)=$(DYN)_state_builder($R0,$p0,$D0)"
        eval(Meta.parse(S_string))
    elseif DYN in SH_LIST
        println(DYN)
        if dyn_sts!=nsts
            println("Warning! Still no implementation for $C0, using [1,0,...] instead")
            D0=zeros(Complex,dyn_sts)
            D0[1]=1
        else
            D0=C0
        end
        S_string="S_$(DYN)=$(DYN)_state_builder($R0,$p0,$D0,$ast0)"
        eval(Meta.parse(S_string))
    else
        println("$DYN dynamics are not implemented in the lists, check the types.jl file, in the end, in the META section")
    end
end

for DYN in DYN_LIST
    if DYN in CL_LIST
        println("BEGINNING $DYN INTEGRATION")
        te0=time()
        integ_string="single_integration(tf,S_$(DYN))"
        T,R,P=eval(Meta.parse(integ_string))
        tef=time()
        println("FINISHED $DYN INTEGRATION FOR $(round(tef-te0,digits=3))s, SAVING...")
        CL_save(T,R,P,file,DYN)
    elseif DYN in MF_LIST
        println("BEGINNING $DYN INTEGRATION")
        te0=time()
        integ_string="single_integration(tf,S_$(DYN))"
        T,R,P,C=eval(Meta.parse(integ_string))
        tef=time()
        println("FINISHED $DYN INTEGRATION FOR $(round(tef-te0,digits=3))s, SAVING...")
        MF_save(T,R,P,C,file,DYN)
    elseif DYN in SH_LIST
        println("BEGINNING $DYN INTEGRATION")
        te0=time()
        integ_string="single_integration(tf,S_$(DYN))"
        T,R,P,C,AST=eval(Meta.parse(integ_string))
        tef=time()
        println("FINISHED $DYN INTEGRATION FOR $(round(tef-te0,digits=3))s, SAVING...")
        SH_save(T,R,P,C,AST,file,DYN)
    end
end

println("FINISHED SUCCESSFULLY!")
