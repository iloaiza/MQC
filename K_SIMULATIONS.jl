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
BO_DYN=false
EH_DYN=false
FSSH_DYN=false
FSSH_DIA_DYN=false
CM2_VANILLA_DYN=false
CM3_VANILLA_DYN=false
CM2_DYN=false
CM3_DYN=false
CM2_FSSH_DYN=false
CM3_FSSH_DYN=false
CM2_FSSH_VANILLA_DYN=false
CM3_FSSH_VANILLA_DYN=false
FIRST_RUN=true

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


########### Text output
println("The timestep is")
@show dt
println("The mass is")
@show mass
println("The travel distance, in a.u. is")
@show travel_distance
println("The number of trajectories is")
@show Ntrajs
println("The initial conditions are")
@show R0
@show C0
println("The momentum array is")
@show K
println("The dynamics that will be done are:")
if BO_DYN
    println("Born-Oppenheimer")
    FR_BO=SharedArray{Float64}(msize,Ntrajs)
    FM_BO=SharedArray{Float64}(msize,Ntrajs)
end
if EH_DYN
    println("Ehrenfest")
    FR_EH=SharedArray{Float64}(msize,Ntrajs)
    FM_EH=SharedArray{Float64}(msize,Ntrajs)
    Fpop_EH=SharedArray{Float64}(msize,Ntrajs,nsts)
end
if FSSH_DYN
    println("Adiabatic FSSH")
    FR_FSSH=SharedArray{Float64}(msize,Ntrajs)
    FM_FSSH=SharedArray{Float64}(msize,Ntrajs)
    Fpop_FSSH=SharedArray{Float64}(msize,Ntrajs,nsts)
    Fast_FSSH=SharedArray{Float64}(msize,Ntrajs)
end
if FSSH_DIA_DYN
    println("Diabatic FSSH")
    FR_FSSHd=SharedArray{Float64}(msize,Ntrajs)
    FM_FSSHd=SharedArray{Float64}(msize,Ntrajs)
    Fpop_FSSHd=SharedArray{Float64}(msize,Ntrajs,nsts)
    Fast_FSSHd=SharedArray{Float64}(msize,Ntrajs)
end
if CM2_VANILLA_DYN
    println("Two collective modes - vanilla version")
    FR_CM2_V=SharedArray{Float64}(msize,Ntrajs)
    FM_CM2_V=SharedArray{Float64}(msize,Ntrajs)
    Fpop_CM2_V=SharedArray{Float64}(msize,Ntrajs,2)
end
if CM3_VANILLA_DYN
    println("Three collective modes - vanilla version")
    FR_CM3_V=SharedArray{Float64}(msize,Ntrajs)
    FM_CM3_V=SharedArray{Float64}(msize,Ntrajs)
    Fpop_CM3_V=SharedArray{Float64}(msize,Ntrajs,3)
end
if CM2_DYN
    println("Two collective modes")
    FR_CM2=SharedArray{Float64}(msize,Ntrajs)
    FM_CM2=SharedArray{Float64}(msize,Ntrajs)
    Fpop_CM2=SharedArray{Float64}(msize,Ntrajs,2)
end
if CM3_DYN
    println("Three collective modes")
    FR_CM3=SharedArray{Float64}(msize,Ntrajs)
    FM_CM3=SharedArray{Float64}(msize,Ntrajs)
    Fpop_CM3=SharedArray{Float64}(msize,Ntrajs,3)
end
if CM2_FSSH_VANILLA_DYN
    println("Surface hopping in vanilla two collective modes")
    FR_CM2_F_V=SharedArray{Float64}(msize,Ntrajs)
    FM_CM2_F_V=SharedArray{Float64}(msize,Ntrajs)
    Fpop_CM2_F_V=SharedArray{Float64}(msize,Ntrajs,2)
    Fast_CM2_F_V=SharedArray{Float64}(msize,Ntrajs)
end
if CM3_FSSH_VANILLA_DYN
    println("Surface hopping in vanilla three collective modes")
    FR_CM3_F_V=SharedArray{Float64}(msize,Ntrajs)
    FM_CM3_F_V=SharedArray{Float64}(msize,Ntrajs)
    Fpop_CM3_F_V=SharedArray{Float64}(msize,Ntrajs,3)
    Fast_CM3_F_V=SharedArray{Float64}(msize,Ntrajs)
end


for p0 in K
    tic0=time()
    tf=round(mass*travel_distance/p0)
    ########## Stuff for text output (i.e. do a logfile!)
    S_EH=EH_state(R0,p0,C0)
    file=dirname*"/p0($p0).h5"

    Ep0,~,~,~,~,~=adiabatic_values(R0);
    E0=real(C0'*diagm(Ep0)*C0)[1]+p0^2/2/mass;

    if FIRST_RUN
        h5write(file,"R0",R0)
        h5write(file,"p0",p0)
        h5write(file,"C0_real",real.(C0))
        h5write(file,"C0_imag",imag.(C0))
        h5write(file,"E0",E0)
    end

    ###########    BO    ############
    if BO_DYN
        println("BEGINNING BO INTEGRATIONS")
        te0=time()
        #Tf_e,Rvec_e,pvec_e,C_e=single_integration(tf,S_EH);   #for single trajectory integration
        #EH_sanity_check(Rvec_e,pvec_e,C_e);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
        T_BO_s,R_BO_s,P_BO_s=dist_BO_integration(tf,R0,p0,Ntrajs,initial_dist);
        tef=time()
        println("FINISHED BO INTEGRATIONS FOR $(round(tef-te0,3))s, SAVING...")
        T_BO=zeros(T_BO_s)
        T_BO.=T_BO_s
        R_BO=zeros(R_BO_s)
        R_BO.=R_BO_s
        P_BO=zeros(P_BO_s)
        P_BO.=P_BO_s
        BO_save(T_BO,R_BO,P_BO,file)
    end

    ###########    EHRENFEST    ############
    if EH_DYN
        println("BEGINNING EHRENFEST INTEGRATIONS")
        te0=time()
        #Tf_e,Rvec_e,pvec_e,C_e=single_integration(tf,S_EH);   #for single trajectory integration
        #EH_sanity_check(Rvec_e,pvec_e,C_e);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
        T_EH_s,R_EH_s,P_EH_s,C_EH_s=dist_EH_integration(tf,R0,p0,C0,Ntrajs,initial_dist);
        tef=time()
        println("FINISHED EHRENFEST INTEGRATIONS FOR $(round(tef-te0,3))s, SAVING...")
        T_EH=zeros(T_EH_s)
        T_EH.=T_EH_s
        R_EH=zeros(R_EH_s)
        R_EH.=R_EH_s
        P_EH=zeros(P_EH_s)
        P_EH.=P_EH_s
        C_EH=zeros(C_EH_s)
        C_EH.=C_EH_s
        EH_save(T_EH,R_EH,P_EH,C_EH,file)
    end

    ########### SURFACE HOPPING  ###########
    if FSSH_DYN
        println("BEGINNING FSSH INTEGRATIONS")
        tf0=time()
        #Tf_f,Rvec_f,pvec_f,C_f,Ast_f=single_integration(tf,S_FSSH);   #for single trajectory integration
        #FSSH_sanity_check(Rvec_f,pvec_f,C_f,Ast_f);    #sanity check for FSSH to know if timestep is good. needs to run single trajectory first. FSSH should have at least 1 hop for debugging!
        T_FSSH_s,R_FSSH_s,P_FSSH_s,C_FSSH_s,AST_FSSH_s=dist_FSSH_integration(tf,R0,p0,C0,ast0,Ntrajs,initial_dist);
        tff=time()
        println("FINISHED FSSH INTEGRATIONS FOR $(round(tff-tf0,3))s, SAVING...")
        T_FSSH=zeros(T_FSSH_s)
        T_FSSH.=T_FSSH_s
        R_FSSH=zeros(R_FSSH_s)
        R_FSSH.=R_FSSH_s
        P_FSSH=zeros(P_FSSH_s)
        P_FSSH.=P_FSSH_s
        C_FSSH=zeros(C_FSSH_s)
        C_FSSH.=C_FSSH_s
        AST_FSSH=zeros(AST_FSSH_s)
        AST_FSSH.=AST_FSSH_s
        FSSH_save(T_FSSH,R_FSSH,P_FSSH,C_FSSH,AST_FSSH,file)
    end

    ########### DIABATIC SURFACE HOPPING  ###########
    if FSSH_DIA_DYN
        println("BEGINNING DIABATIC FSSH INTEGRATIONS")
        tf0=time()
        #Tf_fd,Rvec_fd,pvec_fd,C_fd,Dst_fd=single_integration(tf,S_FSSH_dia);   #for single trajectory integration
        #FSSH_dia_sanity_check(Rvec_fd,pvec_fd,C_fd,Dst_fd);    #sanity check for FSSH to know if timestep is good. needs to run single trajectory first. FSSH should have at least 1 hop for debugging!
        T_FSSHd_s,R_FSSHd_s,P_FSSHd_s,C_FSSHd_s,DST_FSSHd_s=dist_FSSH_dia_integration(tf,R0,p0,C0,ast0,Ntrajs,initial_dist);
        tff=time()
        println("FINISHED DIABATIC FSSH INTEGRATIONS FOR $(round(tff-tf0,3))s, SAVING...")
        T_FSSHd=zeros(T_FSSHd_s)
        T_FSSHd.=T_FSSHd_s
        R_FSSHd=zeros(R_FSSHd_s)
        R_FSSHd.=R_FSSHd_s
        P_FSSHd=zeros(P_FSSHd_s)
        P_FSSHd.=P_FSSHd_s
        C_FSSHd=zeros(C_FSSHd_s)
        C_FSSHd.=C_FSSHd_s
        DST_FSSHd=zeros(DST_FSSHd_s)
        DST_FSSHd.=DST_FSSHd_s
        FSSH_dia_save(T_FSSHd,R_FSSHd,P_FSSHd,C_FSSHd,DST_FSSHd,file)
    end

    ###########    CM2_VANILLA    ############
    if CM2_VANILLA_DYN
        println("BEGINNING CM2 VANILLA INTEGRATIONS")
        te0=time()
        #Tf_cm2,Rvec_cm2,pvec_cm2,D_cm2=single_integration(tf,S_CM2);   #for single trajectory integration
        #CM2_sanity_check(Rvec_cm2,pvec_cm2,D_cm2);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
        T_CM2_VANILLA_s,R_CM2_VANILLA_s,P_CM2_VANILLA_s,D_CM2_VANILLA_s=dist_CM2_VANILLA_integration(tf,R0,p0,D2_0,Ntrajs,initial_dist);
        tef=time()
        println("FINISHED CM2 VANILLA INTEGRATIONS FOR $(round(tef-te0,3))s, SAVING...")
        T_CM2_VANILLA=zeros(T_CM2_VANILLA_s)
        T_CM2_VANILLA.=T_CM2_VANILLA_s
        R_CM2_VANILLA=zeros(R_CM2_VANILLA_s)
        R_CM2_VANILLA.=R_CM2_VANILLA_s
        P_CM2_VANILLA=zeros(P_CM2_VANILLA_s)
        P_CM2_VANILLA.=P_CM2_VANILLA_s
        D_CM2_VANILLA=zeros(D_CM2_VANILLA_s)
        D_CM2_VANILLA.=D_CM2_VANILLA_s
        CM2_VANILLA_save(T_CM2_VANILLA,R_CM2_VANILLA,P_CM2_VANILLA,D_CM2_VANILLA,file)
    end

    ###########    CM3_VANILLA    ############
    if CM3_VANILLA_DYN
        println("BEGINNING CM3 VANILLA INTEGRATIONS")
        te0=time()
        #Tf_cm2,Rvec_cm2,pvec_cm2,D_cm2=single_integration(tf,S_CM2);   #for single trajectory integration
        #CM2_sanity_check(Rvec_cm2,pvec_cm2,D_cm2);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
        T_CM3_VANILLA_s,R_CM3_VANILLA_s,P_CM3_VANILLA_s,D_CM3_VANILLA_s=dist_CM3_VANILLA_integration(tf,R0,p0,D3_0,Ntrajs,initial_dist);
        tef=time()
        println("FINISHED CM3 VANILLA INTEGRATIONS FOR $(round(tef-te0,3))s, SAVING...")
        T_CM3_VANILLA=zeros(T_CM3_VANILLA_s)
        T_CM3_VANILLA.=T_CM3_VANILLA_s
        R_CM3_VANILLA=zeros(R_CM3_VANILLA_s)
        R_CM3_VANILLA.=R_CM3_VANILLA_s
        P_CM3_VANILLA=zeros(P_CM3_VANILLA_s)
        P_CM3_VANILLA.=P_CM3_VANILLA_s
        D_CM3_VANILLA=zeros(D_CM3_VANILLA_s)
        D_CM3_VANILLA.=D_CM3_VANILLA_s
        CM3_VANILLA_save(T_CM3_VANILLA,R_CM3_VANILLA,P_CM3_VANILLA,D_CM3_VANILLA,file)
    end


    ###########    CM2    ############
    if CM2_DYN
        println("BEGINNING CM2 INTEGRATIONS")
        te0=time()
        #Tf_cm2,Rvec_cm2,pvec_cm2,D_cm2=single_integration(tf,S_CM2);   #for single trajectory integration
        #CM2_sanity_check(Rvec_cm2,pvec_cm2,D_cm2);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
        T_CM2_s,R_CM2_s,P_CM2_s,D_CM2_s=dist_CM2_integration(tf,R0,p0,D2_0,Ntrajs,initial_dist);
        tef=time()
        println("FINISHED CM2 INTEGRATIONS FOR $(round(tef-te0,3))s, SAVING...")
        T_CM2=zeros(T_CM2_s)
        T_CM2.=T_CM2_s
        R_CM2=zeros(R_CM2_s)
        R_CM2.=R_CM2_s
        P_CM2=zeros(P_CM2_s)
        P_CM2.=P_CM2_s
        D_CM2=zeros(D_CM2_s)
        D_CM2.=D_CM2_s
        CM2_save(T_CM2,R_CM2,P_CM2,D_CM2,file)
    end

    ###########    CM3    ############
    if CM3_DYN
        println("BEGINNING CM3 INTEGRATIONS")
        te0=time()
        #Tf_cm3,Rvec_cm3,pvec_cm3,D_cm3=single_integration(tf,S_CM3);   #for single trajectory integration
        #CM3_sanity_check(Rvec_cm3,pvec_cm3,D_cm3);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
        T_CM3_s,R_CM3_s,P_CM3_s,D_CM3_s=dist_CM3_integration(tf,R0,p0,D3_0,Ntrajs,initial_dist);
        tef=time()
        println("FINISHED CM3 INTEGRATIONS FOR $(round(tef-te0,3))s, SAVING...")
        T_CM3=zeros(T_CM3_s)
        T_CM3.=T_CM3_s
        R_CM3=zeros(R_CM3_s)
        R_CM3.=R_CM3_s
        P_CM3=zeros(P_CM3_s)
        P_CM3.=P_CM3_s
        D_CM3=zeros(D_CM3_s)
        D_CM3.=D_CM3_s
        CM3_save(T_CM3,R_CM3,P_CM3,D_CM3,file)
    end


    ###########  CM2_FSSH  ############
    if CM2_FSSH_DYN
        println("BEGINNING CM2 FSSH INTEGRATIONS")
        te0=time()
        #Tf_cm2_fssh,Rvec_cm2_fssh,pvec_cm2_fssh,D_cm2_fssh=single_integration(tf,S_CM2_fssh);   #for single trajectory integration
        #CM2_sanity_check(Rvec_cm2_fssh,pvec_cm2_fssh,D_cm2_fssh);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
        T_CM2_FSSH_s,R_CM2_FSSH_s,P_CM2_FSSH_s,D_CM2_FSSH_s=dist_CM2_FSSH_integration(tf,R0,p0,D2_0,Ntrajs,initial_dist);
        tef=time()
        println("FINISHED CM2 FSSH INTEGRATIONS FOR $(round(tef-te0,3))s, SAVING...")
        T_CM2_FSSH=zeros(T_CM2_FSSH_s)
        T_CM2_FSSH.=T_CM2_FSSH_s
        R_CM2_FSSH=zeros(R_CM2_FSSH_s)
        R_CM2_FSSH.=R_CM2_FSSH_s
        P_CM2_FSSH=zeros(P_CM2_FSSH_s)
        P_CM2_FSSH.=P_CM2_FSSH_s
        D_CM2_FSSH=zeros(D_CM2_FSSH_s)
        D_CM2_FSSH.=D_CM2_FSSH_s
        CM2_FSSH_save(T_CM2_FSSH,R_CM2_FSSH,P_CM2_FSSH,D_CM2_FSSH,file)
    end

    ###########  CM3_FSSH  ############
    if CM3_FSSH_DYN
        println("BEGINNING CM3 FSSH INTEGRATIONS")
        te0=time()
        #Tf_cm2_fssh,Rvec_cm2_fssh,pvec_cm2_fssh,D_cm2_fssh=single_integration(tf,S_CM2_fssh);   #for single trajectory integration
        #CM2_sanity_check(Rvec_cm2_fssh,pvec_cm2_fssh,D_cm2_fssh);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
        T_CM3_FSSH_s,R_CM3_FSSH_s,P_CM3_FSSH_s,D_CM3_FSSH_s=dist_CM3_FSSH_integration(tf,R0,p0,D3_0,Ntrajs,initial_dist);
        tef=time()
        println("FINISHED CM3 FSSH INTEGRATIONS FOR $(round(tef-te0,3))s, SAVING...")
        T_CM3_FSSH=zeros(T_CM3_FSSH_s)
        T_CM3_FSSH.=T_CM3_FSSH_s
        R_CM3_FSSH=zeros(R_CM3_FSSH_s)
        R_CM3_FSSH.=R_CM3_FSSH_s
        P_CM3_FSSH=zeros(P_CM3_FSSH_s)
        P_CM3_FSSH.=P_CM3_FSSH_s
        D_CM3_FSSH=zeros(D_CM3_FSSH_s)
        D_CM3_FSSH.=D_CM3_FSSH_s
        CM3_FSSH_save(T_CM3_FSSH,R_CM3_FSSH,P_CM3_FSSH,D_CM3_FSSH,file)
    end

    ###########  CM2_FSSH_VANILLA  ############
    if CM2_FSSH_VANILLA_DYN
        println("BEGINNING CM2 FSSH VANILLA INTEGRATIONS")
        te0=time()
        #Tf_cm2_fssh,Rvec_cm2_fssh,pvec_cm2_fssh,D_cm2_fssh=single_integration(tf,S_CM2_fssh);   #for single trajectory integration
        #CM2_sanity_check(Rvec_cm2_fssh,pvec_cm2_fssh,D_cm2_fssh);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
        T_CM2_FSSH_VANILLA_s,R_CM2_FSSH_VANILLA_s,P_CM2_FSSH_VANILLA_s,D_CM2_FSSH_VANILLA_s=dist_CM2_FSSH_VANILLA_integration(tf,R0,p0,D2_0,Ntrajs,initial_dist);
        tef=time()
        println("FINISHED CM2 FSSH VANILLA INTEGRATIONS FOR $(round(tef-te0,3))s, SAVING...")
        T_CM2_FSSH_VANILLA=zeros(T_CM2_FSSH_VANILLA_s)
        T_CM2_FSSH_VANILLA.=T_CM2_FSSH_VANILLA_s
        R_CM2_FSSH_VANILLA=zeros(R_CM2_FSSH_VANILLA_s)
        R_CM2_FSSH_VANILLA.=R_CM2_FSSH_VANILLA_s
        P_CM2_FSSH_VANILLA=zeros(P_CM2_FSSH_VANILLA_s)
        P_CM2_FSSH_VANILLA.=P_CM2_FSSH_VANILLA_s
        D_CM2_FSSH_VANILLA=zeros(D_CM2_FSSH_VANILLA_s)
        D_CM2_FSSH_VANILLA.=D_CM2_FSSH_VANILLA_s
        CM2_FSSH_VANILLA_save(T_CM2_FSSH_VANILLA,R_CM2_FSSH_VANILLA,P_CM2_FSSH_VANILLA,D_CM2_FSSH_VANILLA,file)
    end

    ###########  CM3_FSSH_VANILLA  ############
    if CM3_FSSH_VANILLA_DYN
        println("BEGINNING CM3 FSSH VANILLA INTEGRATIONS")
        te0=time()
        #Tf_cm2_fssh,Rvec_cm2_fssh,pvec_cm2_fssh,D_cm2_fssh=single_integration(tf,S_CM2_fssh);   #for single trajectory integration
        #CM2_sanity_check(Rvec_cm2_fssh,pvec_cm2_fssh,D_cm2_fssh);    #sanity check for Ehrenfest (also works for FSSH) to know if timestep is good. needs to run single trajectory first
        T_CM3_FSSH_VANILLA_s,R_CM3_FSSH_VANILLA_s,P_CM3_FSSH_VANILLA_s,D_CM3_FSSH_VANILLA_s=dist_CM3_FSSH_VANILLA_integration(tf,R0,p0,D3_0,Ntrajs,initial_dist);
        tef=time()
        println("FINISHED CM3 FSSH VANILLA INTEGRATIONS FOR $(round(tef-te0,3))s, SAVING...")
        T_CM3_FSSH_VANILLA=zeros(T_CM3_FSSH_VANILLA_s)
        T_CM3_FSSH_VANILLA.=T_CM3_FSSH_VANILLA_s
        R_CM3_FSSH_VANILLA=zeros(R_CM3_FSSH_VANILLA_s)
        R_CM3_FSSH_VANILLA.=R_CM3_FSSH_VANILLA_s
        P_CM3_FSSH_VANILLA=zeros(P_CM3_FSSH_VANILLA_s)
        P_CM3_FSSH_VANILLA.=P_CM3_FSSH_VANILLA_s
        D_CM3_FSSH_VANILLA=zeros(D_CM3_FSSH_VANILLA_s)
        D_CM3_FSSH_VANILLA.=D_CM3_FSSH_VANILLA_s
        CM3_FSSH_VANILLA_save(T_CM3_FSSH_VANILLA,R_CM3_FSSH_VANILLA,P_CM3_FSSH_VANILLA,D_CM3_FSSH_VANILLA,file)
    end

    println("FINISHED SUCCESSFULLY FOR p0=$p0")
    tic1=time()
    println("The total time of running is, so far, $(round(tic1-tic00,4))s")
end



println("Finished for all K array, passing to single save file...")
for (pnum,p0) in enumerate(K)
    file=dirname*"/p0($p0).h5"
    if BO_DYN
        T,R,P=BO_read(file)
        FR_BO[pnum,:].=R[end,1,:]
        FM_BO[pnum,:].=P[end,1,:]
    end
    if EH_DYN
        T,R,P,C=EH_read(file)
        FR_EH[pnum,:].=R[end,1,:]
        FM_EH[pnum,:].=P[end,1,:]
        for i in 1:nsts
            Fpop_EH[pnum,:,i].=abs2.(C[end,i,:])
        end
    end
    if FSSH_DYN
        T,R,P,C,AST=FSSH_read(file)
        FR_FSSH[pnum,:].=R[end,1,:]
        FM_FSSH[pnum,:].=P[end,1,:]
        for i in 1:nsts
            Fpop_FSSH[pnum,:,i].=abs2.(C[end,i,:])
        end
        Fast_FSSH[pnum,:].=AST[end,:]
    end
    if FSSH_DIA_DYN
        T,R,P,C,AST=FSSH_dia_read(file)
        FR_FSSHd[pnum,:].=R[end,1,:]
        FM_FSSHd[pnum,:].=P[end,1,:]
        for i in 1:nsts
            Fpop_FSSHd[pnum,:,i].=abs2.(C[end,i,:])
        end
        Fast_FSSHd[pnum,:].=AST[end,:]
    end
    if CM2_DYN
        T,R,P,D=CM2_read(file)
        FR_CM2[pnum,:].=R[end,1,:]
        FM_CM2[pnum,:].=P[end,1,:]
        for i in 1:2
            Fpop_CM2[pnum,:,i].=abs2.(D[end,i,:])
        end
    end
    if CM3_DYN
        T,R,P,D=CM3_read(file)
        FR_CM3[pnum,:].=R[end,1,:]
        FM_CM3[pnum,:].=P[end,1,:]
        for i in 1:3
            Fpop_CM3[pnum,:,i].=abs2.(D[end,i,:])
        end
    end
    if CM2_VANILLA_DYN
        T,R,P,D=CM2_VANILLA_read(file)
        FR_CM2_V[pnum,:].=R[end,1,:]
        FM_CM2_V[pnum,:].=P[end,1,:]
        for i in 1:2
            Fpop_CM2_V[pnum,:,i].=abs2.(D[end,i,:])
        end
    end
    if CM3_VANILLA_DYN
        T,R,P,D=CM3_VANILLA_read(file)
        FR_CM3_V[pnum,:].=R[end,1,:]
        FM_CM3_V[pnum,:].=P[end,1,:]
        for i in 1:3
            Fpop_CM3_V[pnum,:,i].=abs2.(D[end,i,:])
        end
    end
    if CM2_FSSH_DYN
        T,R,P,D,AST=CM2_FSSH_read(file)
        FR_CM2_F[pnum,:].=R[end,1,:]
        FM_CM2_F[pnum,:].=P[end,1,:]
        for i in 1:2
            Fpop_CM2_F[pnum,:,i].=abs2.(D[end,i,:])
        end
        Fast_CM2_F[pnum,:].=AST[end,:]
    end
    if CM3_FSSH_DYN
        T,R,P,D=CM3_FSSH_read(file)
        FR_CM3_F[pnum,:].=R[end,1,:]
        FM_CM3_F[pnum,:].=P[end,1,:]
        for i in 1:3
            Fpop_CM3_F[pnum,:,i].=abs2.(D[end,i,:])
        end
        Fast_CM3_F[pnum,:].=AST[end,:]
    end
    if CM2_FSSH_VANILLA_DYN
        T,R,P,D,AST=CM2_FSSH_VANILLA_read(file)
        FR_CM2_F_V[pnum,:].=R[end,1,:]
        FM_CM2_F_V[pnum,:].=P[end,1,:]
        for i in 1:2
            Fpop_CM2_F_V[pnum,:,i].=abs2.(D[end,i,:])
        end
        Fast_CM2_F_V[pnum,:].=AST[end,:]
    end
    if CM3_FSSH_VANILLA_DYN
        T,R,P,D=CM3_FSSH_VANILLA_read(file)
        FR_CM3_F_V[pnum,:].=R[end,1,:]
        FM_CM3_F_V[pnum,:].=P[end,1,:]
        for i in 1:3
            Fpop_CM3_F_V[pnum,:,i].=abs2.(D[end,i,:])
        end
        Fast_CM3_F_V[pnum,:].=AST[end,:]
    end
end

"Starting save file"
file=dirname*"/K_SIMULATION.h5"
h5write(file,"K",K)
if BO_DYN
    println("Born-Oppenheimer")
    fr_bo=[r for r in FR_BO]
    fm_bo=[p for p in FM_BO]
    general_K_save([fr_bo,fm_bo],file,"BO")
end
if EH_DYN
    println("Ehrenfest")
    fr_eh=[r for r in FR_EH]
    fm_eh=[p for p in FM_EH]
    fpop_eh=[pop for pop in Fpop_EH]
    general_K_save([fr_eh,fm_eh,fpop_eh],file,"EH")
end
if FSSH_DYN
    println("Adiabatic FSSH")
    fr_fssh=[r for r in FR_FSSH]
    fm_fssh=[p for p in FM_FSSH]
    fpop_fssh=[pop for pop in Fpop_FSSH]
    fast_fssh=[ast for ast in Fast_FSSH]
    general_K_save([fr_fssh,fm_fssh,fpop_fssh,fast_fssh],file,"FSSH")
end
if FSSH_DIA_DYN
    println("Diabatic FSSH")
    fr_fsshd=[r for r in FR_FSSHd]
    fm_fsshd=[p for p in FM_FSSHd]
    fpop_fsshd=[pop for pop in Fpop_FSSHd]
    fast_fsshd=[ast for ast in Fast_FSSHd]
    general_K_save([fr_fsshd,fm_fsshd,fpop_fsshd,fast_fsshd],file,"FSSHd")
end
if CM2_VANILLA_DYN
    println("Two collective modes - vanilla version")
    fr_cm2_v=[r for r in FR_CM2_V]
    fm_cm2_v=[p for p in FM_CM2_V]
    fpop_cm2_v=[pop for pop in Fpop_CM2_V]
    general_K_save([fr_cm2_v,fm_cm2_v,fpop_cm2_v],file,"CM2_V")
end
if CM3_VANILLA_DYN
    println("Three collective modes - vanilla version")
    fr_cm3_v=[r for r in FR_CM3_V]
    fm_cm3_v=[p for p in FM_CM3_V]
    fpop_cm3_v=[pop for pop in Fpop_CM3_V]
    general_K_save([fr_cm3_v,fm_cm3_v,fpop_cm3_v],file,"CM3_V")
end
if CM2_DYN
    println("Two collective modes")
    fr_cm2=[r for r in FR_CM2]
    fm_cm2=[p for p in FM_CM2]
    fpop_cm2=[pop for pop in Fpop_CM2]
    general_K_save([fr_cm2,fm_cm2,fpop_cm2],file,"CM2")
end
if CM3_DYN
    println("Three collective modes")
    fr_cm3=[r for r in FR_CM3]
    fm_cm3=[p for p in FM_CM3]
    fpop_cm3=[pop for pop in Fpop_CM3]
    general_K_save([fr_cm3,fm_cm3,fpop_cm3],file,"CM3")
end
if CM2_FSSH_VANILLA_DYN
    println("Surface hopping in vanilla two collective modes")
    fr_cm2_f_v=[r for r in FR_CM2_F_V]
    fm_cm2_f_v=[p for p in FM_CM2_F_V]
    fpop_cm2_f_v=[pop for pop in Fpop_CM2_F_V]
    fast_cm2_f_v=[ast for ast in Fast_CM2_F_V]
    general_K_save([fr_cm2_f_v,fm_cm2_f_v,fpop_cm2_f_v,fast_cm2_f_v],file,"CM2_F_V")
end
if CM3_FSSH_VANILLA_DYN
    println("Surface hopping in vanilla three collective modes")
    fr_cm3_f_v=[r for r in FR_CM3_F_V]
    fm_cm3_f_v=[p for p in FM_CM3_F_V]
    fpop_cm3_f_v=[pop for pop in Fpop_CM3_F_V]
    fast_cm3_f_v=[ast for ast in Fast_CM3_F_V]
    general_K_save([fr_cm3_f_v,fm_cm3_f_v,fpop_cm3_f_v,fast_cm3_f_v],file,"CM3_F_V")
end

println("FINISHED SIMULATIONS TOTAL TIME WAS!")
println("$(time()-tic00) seconds")
