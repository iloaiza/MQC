META_DICT=Dict()

tic00=time()

INSIDE_SIMULATIONS=true
#INSIDE_SIMULATIONS=false                        #define false for running on its own, true (i.e. comment this line) for automatic running inside K_SIMULATIONS

if !INSIDE_SIMULATIONS
        println("Warning: INSIDE_SIMULATIONS marked as false inside K_FILEMAKER, make sure this is intentional")
        include("data_functions.jl")
        using SharedArrays
        using Distributed
        FROM_INITIAL_FILE=true
        FIRST_RUN=false
        if FROM_INITIAL_FILE
                INITIAL_NAME="thirty_states"
                jl_file=pwd()*"/Initial_data/"*INITIAL_NAME*".jl"
                @show jl_file
                include(jl_file)
                dir_name="data/K_"*potname*"_R0($R0)"
                NDOFs=length(R0)
                @show K
        else
                K=1
                dir_name=1
                DYN_LIST=["EH"]
        end
end

msize=length(K)


for DYN in DYN_LIST
    dyn_sts=eval(Meta.parse("$(DYN)_sts"))
    if DYN in CL_LIST
        println("$DYN classical dynamics will be saved")
        R_string="FR_"*DYN*"=SharedArray{Float64}($msize,$NDOFs)"
        P_string="FP_"*DYN*"=SharedArray{Float64}($msize,$NDOFs)"
        meta_string_R="""get!(META_DICT,"FR_$(DYN)","[Ksize,NDOFs]")"""
        meta_string_P="""get!(META_DICT,"FP_$(DYN)","[Ksize,NDOFs]")"""
        eval(Meta.parse(meta_string_R))
        eval(Meta.parse(meta_string_P))
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
    elseif DYN in MF_LIST
        println("$DYN mean-field dynamics will be saved")
        R_string="FR_$(DYN)=SharedArray{Float64}($msize,$NDOFs)"
        P_string="FP_$(DYN)=SharedArray{Float64}($msize,$NDOFs)"
        pop_string="Fpop_$(DYN)=SharedArray{Float64}($msize,$dyn_sts)"
        meta_string_R="""get!(META_DICT,"FR_$(DYN)","[Ksize,NDOFs]")"""
        meta_string_P="""get!(META_DICT,"FP_$(DYN)","[Ksize,NDOFs]")"""
        meta_string_pop="""get!(META_DICT,"Fpop_$(DYN)","[Ksize,electronic_states]")"""
        eval(Meta.parse(meta_string_pop))
        eval(Meta.parse(meta_string_R))
        eval(Meta.parse(meta_string_P))
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        eval(Meta.parse(pop_string))
    elseif DYN in SH_LIST
        println("$DYN surface hopping dynamics will be saved")
        R_string="FR_"*DYN*"=SharedArray{Float64}($Ntrajs,$msize,$NDOFs)"
        P_string="FP_"*DYN*"=SharedArray{Float64}($Ntrajs,$msize,$NDOFs)"
        pop_string="Fpop_"*DYN*"=SharedArray{Float64}($Ntrajs,$msize,$dyn_sts)"
        ast_string="Fast_"*DYN*"=SharedArray{Float64}($Ntrajs,$msize)"
        meta_string_R="""get!(META_DICT,"FR_$(DYN)","[Ntrajs,Ksize,NDOFs]")"""
        meta_string_P="""get!(META_DICT,"FP_$(DYN)","[Ntrajs,Ksize,NDOFs]")"""
        meta_string_pop="""get!(META_DICT,"Fpop_$(DYN)","[Ntrajs,Ksize,electronic_states]")"""
        meta_string_ast="""get!(META_DICT,"Fast_$(DYN)","[Ntrajs,Ksize]")"""
        eval(Meta.parse(meta_string_ast))
        eval(Meta.parse(meta_string_pop))
        eval(Meta.parse(meta_string_R))
        eval(Meta.parse(meta_string_P))
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        eval(Meta.parse(pop_string))
        eval(Meta.parse(ast_string))
    else
        println("$DYN was in DYN_LIST and cannot be found in the list, add it in the types.jl file, by the end in the META part")
    end
end



for (pnum,p0) in enumerate(K)
    file_name=dir_name*"/p0($p0).h5"
    for DYN in DYN_LIST
                dyn_sts=eval(Meta.parse("$(DYN)_sts"))
                if DYN in CL_LIST
                    T,R,P=CL_read(file_name,DYN)
                    for dof in 1:NDOFs
                            R_string="FR_$(DYN)[$pnum,$dof]=$R[$dof]"
                            P_string="FP_$(DYN)[$pnum,$dof]=$P[$dof]"
                            eval(Meta.parse(R_string))
                            eval(Meta.parse(P_string))
                    end
                elseif DYN in MF_LIST
                    T,R,P,C=MF_read(file_name,DYN)
                    for dof in 1:NDOFs
                            R_string="FR_$(DYN)[$pnum,$dof]=$R[$dof]"
                            P_string="FP_$(DYN)[$pnum,$dof]=$P[$dof]"
                            eval(Meta.parse(R_string))
                            eval(Meta.parse(P_string))
                    end
                    for st in 1:dyn_sts
                        pop_string="Fpop_$(DYN)[$pnum,$st]=abs2.($C[$st])"
                        eval(Meta.parse(pop_string))
                    end
                elseif DYN in SH_LIST
                    T,R,P,C,AST=SH_read(file_name,DYN)
                    for dof in 1:NDOFs
                                R_string="FR_$(DYN)[:,$pnum,$dof].=$R[:,end,$dof]"
                                P_string="FP_$(DYN)[:,$pnum,$dof].=$P[:,end,$dof]"
                                eval(Meta.parse(R_string))
                                eval(Meta.parse(P_string))
                    end
                    for st in 1:dyn_sts
                        pop_string="Fpop_$(DYN)[:,$pnum,$st].=abs2.($C[:,$st])"
                        eval(Meta.parse(pop_string))
                    end
                    ast_string="Fast_$(DYN)[:,$pnum].=$AST[:]"
                    eval(Meta.parse(ast_string))
                end #if DYN in SH list
                #println("FINISHED p=$p0 FOR $DYN")
    end #list of DYN
    #println("FINISHED FOR p=$p0")
end #for pnums



"Starting save file"
file_name=dir_name*"/K_SIMULATION.h5"
if FIRST_RUN
    h5write(file_name,"K",K)
    h5write(file_name,"META_DICT",string(META_DICT))
end

for DYN in DYN_LIST
    if DYN in CL_LIST
        R_string="fr=[r for r in FR_$(DYN)]"
        P_string="fp=[p for p in FP_$(DYN)]"
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        general_K_save([fr,fp],file_name,DYN)
    elseif DYN in MF_LIST
        R_string="fr=[r for r in FR_$(DYN)]"
        P_string="fp=[p for p in FP_$(DYN)]"
        pop_string="fpop=[pop for pop in Fpop_$(DYN)]"
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        eval(Meta.parse(pop_string))
        general_K_save([fr,fp,fpop],file_name,DYN)
    elseif DYN in SH_LIST
        R_string="fr=[r for r in FR_$(DYN)]"
        P_string="fp=[p for p in FP_$(DYN)]"
        pop_string="fpop=[pop for pop in Fpop_$(DYN)]"
        ast_string="fast=[ast for ast in Fast_$(DYN)]"
        eval(Meta.parse(R_string))
        eval(Meta.parse(P_string))
        eval(Meta.parse(pop_string))
        eval(Meta.parse(ast_string))
        general_K_save([fr,fp,fpop,fast],file_name,DYN)
    end
end



println("FINISHED SIMULATIONS, TOTAL SAVE TIME WAS!")
println("$(time()-tic00) seconds")
