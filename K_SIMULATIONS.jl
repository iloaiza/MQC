########################################
##########     MQC FILE   ##############
########################################
#for parallel computing, run using 'julia -p N ENERGY_SIMULATIONS.jl INITIAL_DATA_FILENAME'
#otherwise, just run using 'julia ENERGY_SIMULATIONS.jl INITIAL_DATA_FILENAME'
#the INITIAL_DATA_FILENAME is a .jl file that is in the Initial_data folder
tic00=time()

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

    for DYN in DYN_LIST
        dyn_sts=eval(Meta.parse("$(DYN)_sts"))
        println("BEGINNING $DYN INTEGRATION")
        if DYN in CL_LIST
            te0=time()
            single_dist_string="single_distance_integration($R_min,$(DYN)_state_builder($R0,$p0),$tmax)"
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
            single_dist_string="single_distance_integration($R_min,$(DYN)_state_builder($R0,$p0,$D0),$tmax)"
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
            many_dist_string="many_distance_integration($R_min,$(DYN)_state_builder($R0,$p0,$D0,$ast0),$Ntrajs,$tmax)"
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


println("Finished for all K array, passing to single save file...")
for (pnum,p0) in enumerate(K)
    file=dirname*"/p0($p0).h5"
    if NDOFs == 1
        for DYN in DYN_LIST
            dyn_sts=eval(Meta.parse("$(DYN)_sts"))
            if DYN in CL_LIST
                T,R,P=CL_read(file,DYN)
                R_string="FR_$(DYN)[$pnum,1]=$R[end]"
                P_string="FP_$(DYN)[$pnum,1]=$P[end]"
                eval(Meta.parse(R_string))
                eval(Meta.parse(P_string))
            elseif DYN in MF_LIST
                T,R,P,C=MF_read(file,DYN)
                R_string="FR_$(DYN)[$pnum,1]=$R[end]"
                P_string="FP_$(DYN)[$pnum,1]=$P[end]"
                for st in 1:dyn_sts
                    pop_string="Fpop_$(DYN)[$pnum,$st]=abs2.($C[$st])"
                    eval(Meta.parse(pop_string))
                end
                eval(Meta.parse(R_string))
                eval(Meta.parse(P_string))
            elseif DYN in SH_LIST
                for traj in 1:Ntrajs
                    T,R,P,C,AST=SH_read(file,DYN)
                    R_string="FR_$(DYN)[$traj,$pnum,1].=$R[end,$traj]"
                    P_string="FP_$(DYN)[$traj,$pnum,1].=$P[end,$traj]"
                    for st in 1:dyn_sts
                        pop_string="Fpop_$(DYN)[$traj,$pnum,$st]=abs2.($C[$traj,$st])"
                        eval(Meta.parse(pop_string))
                    end
                    ast_string="Fast_$(DYN)[$traj,$pnum]=$AST[$traj]"
                    eval(Meta.parse(R_string))
                    eval(Meta.parse(P_string))
                    eval(Meta.parse(ast_string))
                end #for traj in 1:Ntrajs
            end #if DYN in SH list
        end #list of DYN
    else #if NDOFs != 1
        for DYN in DYN_LIST
            dyn_sts=eval(Meta.parse("$(DYN)_sts"))
            if DYN in CL_LIST
                T,R,P=CL_read(file,DYN)
                R_string="FR_$(DYN)[$pnum,:]=$R[end,:]"
                P_string="FP_$(DYN)[$pnum,:]=$P[end,:]"
                eval(Meta.parse(R_string))
                eval(Meta.parse(P_string))
            elseif DYN in MF_LIST
                T,R,P,C=MF_read(file,DYN)
                R_string="FR_$(DYN)[$pnum,:]=$R[end,:]"
                P_string="FP_$(DYN)[$pnum,:]=$P[end,:]"
                for st in 1:dyn_sts
                    pop_string="Fpop_$(DYN)[$pnum,$st]=abs2.($C[$st])"
                    eval(Meta.parse(pop_string))
                end
                eval(Meta.parse(R_string))
                eval(Meta.parse(P_string))
            elseif DYN in SH_LIST
                for traj in 1:Ntrajs
                    T,R,P,C,AST=SH_read(file,DYN)
                    R_string="FR_$(DYN)[$traj,$pnum,:].=$R[end,$traj,:]"
                    P_string="FP_$(DYN)[$traj,$pnum,:].=$P[end,$traj,:]"
                    for st in 1:dyn_sts
                        pop_string="Fpop_$(DYN)[$traj,$pnum,$st]=abs2.($C[$traj,$st])"
                        eval(Meta.parse(pop_string))
                    end
                    ast_string="Fast_$(DYN)[$traj,$pnum]=$AST[$traj]"
                    eval(Meta.parse(R_string))
                    eval(Meta.parse(P_string))
                    eval(Meta.parse(ast_string))
                end #for traj in 1:Ntrajs
            end #if DYN in SH list
        end #list of DYN
    end #if NDOFs==1
end #for pnums


"Starting save file"
file=dirname*"/K_SIMULATION.h5"
if FIRST_RUN
    h5write(file,"K",K)
end

for DYN in DYN_LIST
    if DYN in CL_LIST
        R_string="fr=[r for r in FR_$(DYN)]"
        P_string="fp=[p for p in FP_$(DYN)]"
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        general_K_save([fr,fp],file,DYN)
    elseif DYN in MF_LIST
        R_string="fr=[r for r in FR_$(DYN)]"
        P_string="fp=[p for p in FP_$(DYN)]"
        pop_string="fpop=[pop for pop in Fpop_$(DYN)]"
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        eval(Meta.parse(pop_string))
        general_K_save([fr,fp,fpop],file,DYN)
    elseif DYN in SH_LIST
        R_string="fr=[r for r in FR_$(DYN)]"
        P_string="fp=[p for p in FP_$(DYN)]"
        pop_string="fpop=[pop for pop in Fpop_$(DYN)]"
        ast_string="fast=[ast for ast in Fast_$(DYN)]"
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        eval(Meta.parse(pop_string))
        eval(Meta.parse(ast_string))
        general_K_save([fr,fp,fpop,fast],file,DYN)
    end
end

println("FINISHED SIMULATIONS TOTAL TIME WAS!")
println("$(time()-tic00) seconds")

include("K_FILEMAKER")
