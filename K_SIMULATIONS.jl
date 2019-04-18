########################################
##########     MQC FILE   ##############
########################################
#for parallel computing, run using 'julia -p N ENERGY_SIMULATIONS.jl INITIAL_DATA_FILENAME'
#otherwise, just run using 'julia ENERGY_SIMULATIONS.jl INITIAL_DATA_FILENAME'
#the INITIAL_DATA_FILENAME is a .jl file that is in the Initial_data folder
tic00=time()
GEN_METHOD = "K_SIMULATIONS"

using HDF5 #package used for saving data
using Distributions #used for initial conditions distributions (such as wigner)

@everywhere input_name=remotecall_fetch(i->ARGS[i],1,1)


initial_dist=false
FIRST_RUN=true
tmax=10000 #default maximum time, ~241 fs

println("Including function modules")
@everywhere include("Initial_data/"*input_name*".jl")
if initial_dist==false
    initial_dist=constant_dist
end

msize=length(K) #size of momentum vector
dirname="data/K_"*potname*"_R0($R0)"
if FIRST_RUN
    mkdir(dirname)
end

if initial_dist==false
    initial_dist=constant_dist
end

NDOFs=length(R0)
if NDOFs>1
    println("Warning: the R_min stopping point is not implemented for several dimensions")
end

########### Text output
println("The timestep is")
@show dt
println("The mass is")
@show mass
println("The stop point is")
@show R_min
println("The number of trajectories is")
@show Ntrajs
println("The initial conditions are")
@show R0
@show C0
println("The momentum array is")
@show K
println("The maximum time per trajectory is")
@show tmax
println("The initial distribution for the trajectories is")
@show initial_dist
println("The dynamics that will be done are:")
for DYN in DYN_LIST
    dyn_sts=eval(Meta.parse("$(DYN)_sts"))
    if DYN in CL_LIST
        R_string="FR_"*DYN*"=SharedArray{Float64}($msize,$NDOFs)"
        P_string="FP_"*DYN*"=SharedArray{Float64}($msize,$NDOFs)"
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        println("$DYN classical dynamics will be done")
    elseif DYN in MF_LIST
        R_string="FR_"*DYN*"=SharedArray{Float64}($msize,$NDOFs)"
        P_string="FP_"*DYN*"=SharedArray{Float64}($msize,$NDOFs)"
        pop_string="Fpop_"*DYN*"=SharedArray{Float64}($msize,$dyn_sts)"
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        eval(Meta.parse(pop_string))
        println("$DYN mean-field dynamics will be done")
    elseif DYN in SH_LIST
        R_string="FR_"*DYN*"=SharedArray{Float64}($Ntrajs,$msize,$NDOFs)"
        P_string="FP_"*DYN*"=SharedArray{Float64}($Ntrajs,$msize,$NDOFs)"
        pop_string="Fpop_"*DYN*"=SharedArray{Float64}($Ntrajs,$msize,$dyn_sts)"
        ast_string="Fast_"*DYN*"=SharedArray{Float64}($Ntrajs,$msize)"
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        eval(Meta.parse(pop_string))
        eval(Meta.parse(ast_string))
        println("$DYN surface hopping dynamics will be done")
    else
        println("$DYN was in the list and cannot be found in the list, add it in the types.jl file, by the end in the META part")
    end
end


for p0 in K
    tic0=time()
    println("Beginning with K=$p0")
    ########## Stuff for text output (i.e. do a logfile!)
    S_EH=EH_state_builder(R0,p0,C0)
    file=dirname*"/p0($p0).h5"

    Ep0,_,_,_,_,_=adiabatic_values(R0);
    E0=real(C0'*Diagonal(Ep0)*C0)[1]+p0^2/2/mass;

    if FIRST_RUN
        h5write(file,"R0",R0)
        h5write(file,"p0",p0)
        h5write(file,"C0_real",real.(C0))
        h5write(file,"C0_imag",imag.(C0))
        h5write(file,"E0",E0)
        h5write(file,"DYN_LIST",DYN_LIST)
    end

    for DYN in DYN_LIST
        dyn_sts=eval(Meta.parse("$(DYN)_sts"))
        println("BEGINNING $DYN INTEGRATION")
        if eval(Meta.parse("@isdefined $(DYN)_mem")) #block for filling memory variable with 0 if it is not defined
            eval(Meta.parse("mem=$(DYN)_mem"))
            println("""Using memory $(DYN)_mem = $(mem)""")
        else
            eval(Meta.parse("$(DYN)_mem=0"))
            println("Memory not defined for $(DYN), make sure it's intentional")
        end
        if DYN in CL_LIST
            te0=time()
            single_dist_string="single_distance_integration($R_min,$(DYN)_state_builder($R0,$p0,0,$NDOFs,$(DYN)_mem),$tmax)"
            T,S=eval(Meta.parse(single_dist_string))
            tef=time()
            println("FINISHED $DYN INTEGRATION FOR $(round(tef-te0,digits=3))s, SAVING...")
            CL_save(T,S.cl.R,S.cl.p,file,DYN)
        elseif DYN in MF_LIST
            te0=time()
            if dyn_sts!=nsts
                println("Warning: replacing C0 for [1,0,..] for this method since implementation is not made")
                D0=zeros(Complex,dyn_sts)
                D0[1]=1
            else
                D0=C0
            end
            single_dist_string="single_distance_integration($R_min,$(DYN)_state_builder($R0,$p0,$D0,0,$NDOFs,$(DYN)_mem),$tmax)"
            T,S=eval(Meta.parse(single_dist_string))
            tef=time()
            println("FINISHED $DYN INTEGRATION FOR $(round(tef-te0,digits=3))s, SAVING...")
            MF_save(T,S.cl.R,S.cl.p,S.el.C,file,DYN)
        elseif DYN in SH_LIST
            te0=time()
            if dyn_sts!=nsts
                println("Warning: replacing C0 for [1,0,..] for this method since implementation is not made")
                D0=zeros(Complex,dyn_sts)
                D0[1]=1
            else
                D0=C0
            end
            many_dist_string="many_distance_integration($R_min,$(DYN)_state_builder($R0,$p0,$D0,$ast0,0,$NDOFs,$(DYN)_mem),$Ntrajs,$tmax)"
            Ts,Rs,Ps,Cs,ASTs=eval(Meta.parse(many_dist_string))
            tef=time()
            println("FINISHED $DYN INTEGRATION FOR $(round(tef-te0,digits=3))s, SAVING...")
            T=[t for t in Ts]
            R=[r for r in Rs]
            P=[p for p in Ps]
            C=[c for c in Cs]
            AST=[ast for ast in ASTs]
            SH_save(T,R,P,C,AST,file,DYN)
        end
    end
end

include("K_FILEMAKER.jl")

if plot_out #option activated in code_config or in startup file after
    include("automatic_plotting.jl")
end
