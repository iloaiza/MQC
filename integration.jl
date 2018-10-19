function single_integration(tf,S::CL_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(range(1,stop=steps,length=flags+1)))]
    Rvec=zeros(Float64,flags+1,S.cl.NDOFs)
    pvec=zeros(Float64,flags+1,S.cl.NDOFs)

    counter=2
    Rvec[1,:].=S.cl.R
    pvec[1,:].=S.cl.p
    for i in 2:steps
        S=runge_step(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.cl.R
            pvec[counter,:].=S.cl.p
            counter+=1
        end
    end
    Rvec[end,:].=S.cl.R
    pvec[end,:].=S.cl.p
    return Tf,Rvec,pvec
end

function single_integration(tf,S::MF_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(range(1,stop=steps,length=flags+1)))]
    Rvec=zeros(Float64,flags+1,S.cl.NDOFs)
    pvec=zeros(Float64,flags+1,S.cl.NDOFs)
    tr_sts=length(S.el.C) #number of electronic states that will be tracked over te dynamics
    C=zeros(Complex,flags+1,tr_sts)

    counter=2
    Rvec[1,:].=S.cl.R
    pvec[1,:].=S.cl.p
    C[1,:].=S.el.C
    #t00=time()
    for i in 2:steps
        S=runge_step(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.cl.R
            pvec[counter,:].=S.cl.p
            C[counter,:].=S.el.C
            counter+=1
            #= Uncomment comment line for printing debug mode
            @show T[i]
            @show S.cl.R
            println("Time since last check: $(round(time()-t00,digits=3))")
            t00=time()
            ## =#
        end
    end
    Rvec[end,:].=S.cl.R
    pvec[end,:].=S.cl.p
    C[end,:].=S.el.C
    return Tf,Rvec,pvec,C
end

function single_integration(tf,S::SH_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(range(1,stop=steps,length=flags+1)))]
    Rvec=zeros(Float64,flags+1,S.cl.NDOFs)
    pvec=zeros(Float64,flags+1,S.cl.NDOFs)
    tr_sts=length(S.el.C) #number of electronic states to track
    C=zeros(Complex,flags+1,tr_sts)
    Ast=zeros(Int,flags+1)

    counter=2
    Rvec[1,:].=S.cl.R
    pvec[1,:].=S.cl.p
    C[1,:].=S.el.C
    Ast[1]=S.ast
    for i in 2:steps
        S=runge_step(S)
        hop!(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.cl.R
            pvec[counter,:].=S.cl.p
            C[counter,:].=S.el.C
            Ast[counter]=S.ast
            counter+=1
        end
    end
    Rvec[end,:].=S.cl.R
    pvec[end,:].=S.cl.p
    C[end,:].=S.el.C
    Ast[end]=S.ast

    return Tf,Rvec,pvec,C,Ast
end


function wigner_CL_integration(tf,R0,p0,prefix,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[builder_CL_state(INITIAL[1,i],INITIAL[2,i],prefix) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(range(1,stop=steps,length=flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @distributed for i in 1:Ntrajs
            @time Tf,Rvec,pvec=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC
end

function wigner_MF_integration(tf,R0,p0,C0,prefix,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    C_VEC=SharedArray{Complex128}(flags+1,nsts,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[builder_MF_state(INITIAL[1,i],INITIAL[2,i],C0,prefix) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(range(1,stop=steps,length=flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @distributed for i in 1:Ntrajs
            @time Tf,Rvec,pvec,C=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            C_VEC[:,:,i].=C
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,C_VEC
end

function wigner_SH_integration(tf,R0,p0,C0,ast0,prefix,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    tr_sts=length(C0)
    C_VEC=SharedArray{Complex}(flags+1,tr_sts,Ntrajs)
    AST_VEC=SharedArray{Int64}(flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[builder_SH_state(INITIAL[1,i],INITIAL[2,i],C0,ast0,prefix) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(range(1,stop=steps,length=flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @distributed for i in 1:Ntrajs
            @time Tf,Rvec,pvec,C,Ast=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            C_VEC[:,:,i].=C
            AST_VEC[:,i].=Ast
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,C_VEC,AST_VEC
end

function dist_CL_integration(tf,R0,p0,prefix,Ntrajs,DIST,flags=100)
    if DIST==constant_dist
        println("Warning: running many classical trajectories with same initial conditions, one trajectory is enough")
    end
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[builder_CL_state(INITIAL[1,i],INITIAL[2,i],prefix) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(range(1,stop=steps,length=flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @distributed for i in 1:Ntrajs
            @time Tf,Rvec,pvec=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC
end

function dist_MF_integration(tf,R0,p0,C0,prefix,Ntrajs,DIST,flags=100)
    if DIST==constant_dist
        println("Warning: running many classical trajectories with same initial conditions, one trajectory is enough")
    end
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    tr_sts=length(C0)
    C_REAL_VEC=SharedArray{Float64}(flags+1,tr_sts,Ntrajs)
    C_IMAG_VEC=SharedArray{Float64}(flags+1,tr_sts,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[builder_MF_state(INITIAL[1,i],INITIAL[2,i],C0,prefix) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(range(1,stop=steps,length=flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @distributed for i in 1:Ntrajs
            @time Tf,Rvec,pvec,C=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            C_REAL_VEC[:,:,i].=Float64.(real.(C))
            C_IMAG_VEC[:,:,i].=Float64.(imag.(C))
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

    C_VEC=C_REAL_VEC.+1im*C_IMAG_VEC

        return TF,R_VEC,P_VEC,C_VEC
end

function dist_SH_integration(tf,R0,p0,C0,ast0,prefix,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    tr_sts=length(C0)
    C_REAL_VEC=SharedArray{Float64}(flags+1,tr_sts,Ntrajs)
    C_IMAG_VEC=SharedArray{Float64}(flags+1,tr_sts,Ntrajs)
    AST_VEC=SharedArray{Int64}(flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[builder_SH_state(INITIAL[1,i],INITIAL[2,i],C0,ast0,prefix) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(range(1,stop=steps,length=flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @distributed for i in 1:Ntrajs
            @time Tf,Rvec,pvec,C,Ast=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            C_REAL_VEC[:,:,i].=Float64.(real.(C))
            C_IMAG_VEC[:,:,i].=Float64.(imag.(C))
            AST_VEC[:,i].=Ast
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

    C_VEC=C_REAL_VEC.+1im*C_IMAG_VEC

        return TF,R_VEC,P_VEC,C_VEC,AST_VEC
end


########### TRAVEL DISTANCE INTEGRATIONS

function single_distance_integration(R_min,S::CL_state,tmax=10000)
    tf=0
    if length(R_min)==1
        while S.cl.R-R_min<0 && tf<tmax
            tf+=dt
            S=runge_step(S)
        end
    elseif length(R_min)==2
        while S.cl.R-R_min[2]<0 && S.cl.R-R_min[1]>0 && tf<tmax
            tf+=dt
            S=runge_step(S)
        end
    else
        error("R_min is neither a number nor an interval!")
    end

    return tf,S
end

function single_distance_integration(R_min,S::MF_state,tmax=10000)
    tf=0
    if length(R_min)==1
        while S.cl.R-R_min<0 && tf<tmax
            tf+=dt
            S=runge_step(S)
        end
    elseif length(R_min)==2
        while S.cl.R-R_min[2]<0 && S.cl.R-R_min[1]>0 && tf<tmax
            tf+=dt
            S=runge_step(S)
        end
    else
        error("R_min is neither a number nor an interval!")
    end

    return tf,S
end

function single_distance_integration(R_min,S::SH_state,tmax=10000)
    tf=0
    if length(R_min)==1
        while S.cl.R-R_min<0 && tf<tmax
            tf+=dt
            S=runge_step(S)
            hop!(S)
        end
    elseif length(R_min)==2
        while S.cl.R-R_min[2]<0 && S.cl.R-R_min[1]>0 && tf<tmax
            tf+=dt
            S=runge_step(S)
            hop!(S)
        end
    else
        error("R_min is neither a number nor an interval!")
    end

    if tf==tmax
        println("Warning! This trajectory was ended due to time limit")
    end
    return tf,S
end

function many_distance_integration(R_min,S::SH_state,Ntrajs,tmax=10000)
    #this kind of integration will only be useful for surface hopping schemes,
    #since mean-field schemes only need to be integrated once
    S_array=[S for k in 1:Ntrajs]
    T_array=SharedArray{Float64}(Ntrajs)
    R_array=SharedArray{Float64}(Ntrajs,S.cl.NDOFs)
    P_array=SharedArray{Float64}(Ntrajs,S.cl.NDOFs)
    tr_sts=eval(Meta.parse("$(S.prefix)_sts"))
    C_real_array=SharedArray{Float64}(Ntrajs,tr_sts)
    C_imag_array=SharedArray{Float64}(Ntrajs,tr_sts)
    AST_array=SharedArray{Int64}(Ntrajs)
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @distributed for i in 1:Ntrajs
        @time T_array[i],S=single_distance_integration(R_min,S_array[i],tmax)
        REMAINING[1]+=-1
        println("$(REMAINING[1]) trajectories remaining")
        R_array[i,:].=S.cl.R
        P_array[i,:].=S.cl.p
        C_real_array[i,:].=real.(S.el.C)
        C_imag_array[i,:].=imag.(S.el.C)
        AST_array[i,]=S.ast
    end

    C_array=C_real_array+1im.*C_imag_array
    return T_array,R_array,P_array,C_array,AST_array
end


#########################################################
function full_memory_integration(tf,S::CL_state)
    T=collect(0:dt:tf)
    steps=length(T)
    S_ARRAY=[S for i in 1:steps]

    #t00=time()
    for i in 2:steps
        S_ARRAY[i]=runge_step(S_ARRAY[i-1])
    end

    return T,S_ARRAY
end

function full_memory_integration(tf,S::MF_state)
    T=collect(0:dt:tf)
    steps=length(T)
    S_ARRAY=[S for i in 1:steps]

    #t00=time()
    for i in 2:steps
        S_ARRAY[i]=runge_step(S_ARRAY[i-1])
    end

    return T,S_ARRAY
end

function full_memory_integration(tf,S::SH_state)
    T=collect(0:dt:tf)
    steps=length(T)
    S_ARRAY=[S for i in 1:steps]

    #t00=time()
    for i in 2:steps
        S_ARRAY[i]=runge_step(S_ARRAY[i-1])
        hop!(S_ARRAY[i])
    end

    return T,S_ARRAY
end
