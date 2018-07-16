function single_integration(tf,S::BO_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    for i in 2:steps
        S=runge_step(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    return Tf,Rvec,pvec
end

function single_integration(tf,S::EH_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    C=zeros(Complex128,flags+1,nsts)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    C[1,:].=S.C
    for i in 2:steps
        S=runge_step(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            C[counter,:].=S.C
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    C[end,:].=S.C
    return Tf,Rvec,pvec,C
end

function single_integration(tf,S::FSSH_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    C=zeros(Complex128,flags+1,nsts)
    Ast=zeros(Int,flags+1)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    C[1,:].=S.C
    Ast[1]=S.ast
    for i in 2:steps
        S=runge_step(S)
        hop!(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            C[counter,:].=S.C
            Ast[counter]=S.ast
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    C[end,:].=S.C
    Ast[end]=S.ast

    return Tf,Rvec,pvec,C,Ast
end

function single_integration(tf,S::FSSH_dia_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    C=zeros(Complex128,flags+1,nsts)
    Dst=zeros(Int,flags+1)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    C[1,:].=S.C
    Dst[1]=S.dst
    for i in 2:steps
        S=runge_step(S)
        hop!(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            C[counter,:].=S.C
            Dst[counter]=S.dst
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    C[end,:].=S.C
    Dst[end]=S.dst

    return Tf,Rvec,pvec,C,Dst
end

function single_integration(tf,S::CM2_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    D=zeros(Complex128,flags+1,2)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    D[1,:].=S.D
    for i in 2:steps
        S=runge_step(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            D[counter,:].=S.D
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    D[end,:].=S.D
    return Tf,Rvec,pvec,D
end

function single_integration(tf,S::CM3_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    D=zeros(Complex128,flags+1,3)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    D[1,:].=S.D
    for i in 2:steps
        S=runge_step(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            D[counter,:].=S.D
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    D[end,:].=S.D
    return Tf,Rvec,pvec,D
end

function single_integration(tf,S::CM2_VANILLA_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    D=zeros(Complex128,flags+1,2)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    D[1,:].=S.D
    for i in 2:steps
        S=runge_step(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            D[counter,:].=S.D
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    D[end,:].=S.D
    return Tf,Rvec,pvec,D
end

function single_integration(tf,S::CM3_VANILLA_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    D=zeros(Complex128,flags+1,3)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    D[1,:].=S.D
    for i in 2:steps
        S=runge_step(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            D[counter,:].=S.D
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    D[end,:].=S.D
    return Tf,Rvec,pvec,D
end

function single_integration(tf,S::CM2_FSSH_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    D=zeros(Complex128,flags+1,2)
    AST=zeros(Int,flags+1)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    D[1,:].=S.D
    AST[1].=S.ast
    for i in 2:steps
        S=runge_step(S)
        hop!(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            D[counter,:].=S.D
            AST[counter]=S.ast
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    D[end,:].=S.D
    AST[end]=S.ast
    return Tf,Rvec,pvec,D,AST
end

function single_integration(tf,S::CM3_FSSH_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    D=zeros(Complex128,flags+1,3)
    AST=zeros(Int,flags+1)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    D[1,:].=S.D
    AST[1].=S.ast
    for i in 2:steps
        S=runge_step(S)
        hop!(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            D[counter,:].=S.D
            AST[counter]=S.ast
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    D[end,:].=S.D
    AST[end]=S.ast
    return Tf,Rvec,pvec,D,AST
end

function single_integration(tf,S::CM2_FSSH_VANILLA_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    D=zeros(Complex128,flags+1,2)
    AST=zeros(Int,flags+1)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    D[1,:].=S.D
    AST[1].=S.ast
    for i in 2:steps
        S=runge_step(S)
        hop!(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            D[counter,:].=S.D
            AST[counter]=S.ast
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    D[end,:].=S.D
    AST[end]=S.ast
    return Tf,Rvec,pvec,D,AST
end

function single_integration(tf,S::CM3_FSSH_VANILLA_state,flags=100)
    T=collect(0:dt:tf)
    steps=length(T)
    if steps<flags+1
        error("There are more steps ($steps) than flags ($flags) for saving! Check your timestep, final time and number of flags...")
    end
    Tf=T[Int.(round.(linspace(1,steps,flags+1)))]
    Rvec=zeros(Float64,flags+1,S.NDOFs)
    pvec=zeros(Rvec)
    D=zeros(Complex128,flags+1,3)
    AST=zeros(Int,flags+1)

    counter=2
    Rvec[1,:].=S.R
    pvec[1,:].=S.p
    D[1,:].=S.D
    AST[1].=S.ast
    for i in 2:steps
        S=runge_step(S)
        hop!(S)
        if Tf[counter]==T[i]
            Rvec[counter,:].=S.R
            pvec[counter,:].=S.p
            D[counter,:].=S.D
            AST[counter]=S.ast
            counter+=1
        end
    end
    Rvec[end,:].=S.R
    pvec[end,:].=S.p
    D[end,:].=S.D
    AST[end]=S.ast
    return Tf,Rvec,pvec,D,AST
end

function wigner_BO_integration(tf,R0,p0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[BO_state(INITIAL[1,i],INITIAL[2,i]) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
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

function wigner_EH_integration(tf,R0,p0,C0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    C_VEC=SharedArray{Complex128}(flags+1,nsts,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[EH_state(INITIAL[1,i],INITIAL[2,i],C0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
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

function wigner_FSSH_integration(tf,R0,p0,C0,ast0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    C_VEC=SharedArray{Complex128}(flags+1,nsts,Ntrajs)
    AST_VEC=SharedArray{Int64}(flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[FSSH_state(INITIAL[1,i],INITIAL[2,i],C0,ast0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
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

function wigner_FSSH_dia_integration(tf,R0,p0,C0,ast0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    C_VEC=SharedArray{Complex128}(flags+1,nsts,Ntrajs)
    DST_VEC=SharedArray{Int64}(flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[FSSH_dia_state(INITIAL[1,i],INITIAL[2,i],C0,ast0,NDOFs,true) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,C,Dst=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            C_VEC[:,:,i].=C
            DST_VEC[:,i].=Dst
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,C_VEC,DST_VEC
end

function wigner_CM2_VANILLA_integration(tf,R0,p0,D0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,2,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[CM2_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function wigner_CM3_VANILLA_integration(tf,R0,p0,D0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,3,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[CM3_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function wigner_CM2_integration(tf,R0,p0,D0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,2,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[CM2_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function wigner_CM3_integration(tf,R0,p0,D0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,3,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[CM3_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function wigner_CM2_FSSH_integration(tf,R0,p0,D0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,2,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[CM2_FSSH_state(INITIAL[1,i],INITIAL[2,i],D0,1) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function wigner_CM3_FSSH_integration(tf,R0,p0,D0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,3,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[CM3_FSSH_state(INITIAL[1,i],INITIAL[2,i],D0,1) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function wigner_CM2_FSSH_VANILLA_integration(tf,R0,p0,D0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,2,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[CM2_FSSH_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0,1) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function wigner_CM3_FSSH_VANILLA_integration(tf,R0,p0,D0,Ntrajs,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,3,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],wigner,Ntrajs)
    S0=[CM3_FSSH_state(INITIAL[1,i],INITIAL[2,i],D0,1) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function dist_BO_integration(tf,R0,p0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[BO_state(INITIAL[1,i],INITIAL[2,i]) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
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

function dist_EH_integration(tf,R0,p0,C0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    C_VEC=SharedArray{Complex128}(flags+1,nsts,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[EH_state(INITIAL[1,i],INITIAL[2,i],C0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
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

function dist_FSSH_integration(tf,R0,p0,C0,ast0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    C_VEC=SharedArray{Complex128}(flags+1,nsts,Ntrajs)
    AST_VEC=SharedArray{Int64}(flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    #S0=[FSSH_state(INITIAL[1,i],INITIAL[2,i],C0,ast0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,C,Ast=single_integration(tf,FSSH_state(INITIAL[1,i],INITIAL[2,i],C0,ast0),flags)
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

function dist_FSSH_dia_integration(tf,R0,p0,C0,ast0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    C_VEC=SharedArray{Complex128}(flags+1,nsts,Ntrajs)
    DST_VEC=SharedArray{Int64}(flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[FSSH_dia_state(INITIAL[1,i],INITIAL[2,i],C0,ast0,NDOFs,true) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,C,Dst=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            C_VEC[:,:,i].=C
            DST_VEC[:,i].=Dst
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,C_VEC,DST_VEC
end

function dist_CM2_VANILLA_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,2,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM2_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function dist_CM3_VANILLA_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,3,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM3_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function dist_CM2_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,2,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM2_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function dist_CM3_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,3,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM3_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function dist_CM2_FSSH_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,2,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM2_FSSH_state(INITIAL[1,i],INITIAL[2,i],D0,1) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function dist_CM3_FSSH_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,3,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM3_FSSH_state(INITIAL[1,i],INITIAL[2,i],D0,1) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function dist_CM2_FSSH_VANILLA_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,2,Ntrajs)
    AST_VEC=SharedArray{Float64}(flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM2_FSSH_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0,1) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D,Ast=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            AST_VEC[:,i].=Ast
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC,AST_VEC
end

function dist_CM3_FSSH_VANILLA_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    @everywhere NDOFs=length(R0)
    R_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    P_VEC=SharedArray{Float64}(flags+1,NDOFs,Ntrajs)
    D_VEC=SharedArray{Complex128}(flags+1,3,Ntrajs)
    AST_VEC=SharedArray{Float64}(flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM3_FSSH_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0,1) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=SharedArray(T[Int.(round.(linspace(1,steps,flags+1)))])
    REMAINING=SharedArray{Int64}(1)
    REMAINING[1]=Ntrajs

    @sync @parallel for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D,Ast=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            AST_VEC[:,i].=Ast
            REMAINING[1]+=-1
            println("$(REMAINING[1]) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC,AST_VEC
end

function series_BO_integration(tf,R0,p0,Ntrajs,DIST,flags=100)
    NDOFs=length(R0)
    R_VEC=zeros(flags+1,NDOFs,Ntrajs)
    P_VEC=zeros(flags+1,NDOFs,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[BO_state(INITIAL[1,i],INITIAL[2,i]) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=T[Int.(round.(linspace(1,steps,flags+1)))]
    REMAINING=Ntrajs

    for i in 1:Ntrajs
            @time Tf,Rvec,pvec=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            REMAINING+=-1
            println("$(REMAINING) trajectories remaining")
    end

        return TF,R_VEC,P_VEC
end

function series_EH_integration(tf,R0,p0,C0,Ntrajs,DIST,flags=100)
    NDOFs=length(R0)
    R_VEC=zeros(flags+1,NDOFs,Ntrajs)
    P_VEC=zeros(flags+1,NDOFs,Ntrajs)
    C_VEC=zeros(Complex,flags+1,nsts,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[EH_state(INITIAL[1,i],INITIAL[2,i],C0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=T[Int.(round.(linspace(1,steps,flags+1)))]
    REMAINING=Ntrajs

    for i in 1:Ntrajs
            @time Tf,Rvec,pvec,C=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            C_VEC[:,:,i].=C
            REMAINING+=-1
            println("$(REMAINING) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,C_VEC
end

function series_FSSH_integration(tf,R0,p0,C0,ast0,Ntrajs,DIST,flags=100)
    NDOFs=length(R0)
    R_VEC=zeros(flags+1,NDOFs,Ntrajs)
    P_VEC=zeros(flags+1,NDOFs,Ntrajs)
    C_VEC=zeros(Complex,flags+1,nsts,Ntrajs)
    AST_VEC=zeros(Int,flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[FSSH_state(INITIAL[1,i],INITIAL[2,i],C0,ast0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=T[Int.(round.(linspace(1,steps,flags+1)))]
    REMAINING=Ntrajs

    for i in 1:Ntrajs
            @time Tf,Rvec,pvec,C,Ast=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            C_VEC[:,:,i].=C
            AST_VEC[:,i].=Ast
            REMAINING+=-1
            println("$(REMAINING) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,C_VEC,AST_VEC
end

function series_FSSH_dia_integration(tf,R0,p0,C0,ast0,Ntrajs,DIST,flags=100)
    NDOFs=length(R0)
    R_VEC=zeros(flags+1,NDOFs,Ntrajs)
    P_VEC=zeros(flags+1,NDOFs,Ntrajs)
    C_VEC=zeros(Complex,flags+1,nsts,Ntrajs)
    DST_VEC=zeros(Int,flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[FSSH_dia_state(INITIAL[1,i],INITIAL[2,i],C0,ast0,NDOFs,true) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=T[Int.(round.(linspace(1,steps,flags+1)))]
    REMAINING=Ntrajs

    for i in 1:Ntrajs
            @time Tf,Rvec,pvec,C,Dst=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            C_VEC[:,:,i].=C
            DST_VEC[:,i].=Dst
            REMAINING+=-1
            println("$(REMAINING) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,C_VEC,DST_VEC
end

function series_CM2_VANILLA_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    NDOFs=length(R0)
    R_VEC=zeros(flags+1,NDOFs,Ntrajs)
    P_VEC=zeros(flags+1,NDOFs,Ntrajs)
    D_VEC=zeros(Complex,flags+1,2,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM2_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=T[Int.(round.(linspace(1,steps,flags+1)))]
    REMAINING=Ntrajs

    for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING+=-1
            println("$(REMAINING) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function dist_CM3_VANILLA_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    NDOFs=length(R0)
    R_VEC=zeros(flags+1,NDOFs,Ntrajs)
    P_VEC=zeros(flags+1,NDOFs,Ntrajs)
    D_VEC=zeros(Complex,flags+1,3,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM3_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=T[Int.(round.(linspace(1,steps,flags+1)))]
    REMAINING=Ntrajs

    for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING+=-1
            println("$(REMAINING) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function series_CM2_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    NDOFs=length(R0)
    R_VEC=zeros(flags+1,NDOFs,Ntrajs)
    P_VEC=zeros(flags+1,NDOFs,Ntrajs)
    D_VEC=zeros(Complex,flags+1,2,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM2_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=T[Int.(round.(linspace(1,steps,flags+1)))]
    REMAINING=Ntrajs

    for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING+=-1
            println("$(REMAINING) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function dist_CM3_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    NDOFs=length(R0)
    R_VEC=zeros(flags+1,NDOFs,Ntrajs)
    P_VEC=zeros(flags+1,NDOFs,Ntrajs)
    D_VEC=zeros(Complex,flags+1,3,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM3_state(INITIAL[1,i],INITIAL[2,i],D0) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=T[Int.(round.(linspace(1,steps,flags+1)))]
    REMAINING=Ntrajs

    for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            REMAINING+=-1
            println("$(REMAINING) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC
end

function series_CM2_FSSH_VANILLA_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    NDOFs=length(R0)
    R_VEC=zeros(flags+1,NDOFs,Ntrajs)
    P_VEC=zeros(flags+1,NDOFs,Ntrajs)
    D_VEC=zeros(Complex,flags+1,2,Ntrajs)
    AST_VEC=zeros(flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM2_FSSH_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0,1) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=T[Int.(round.(linspace(1,steps,flags+1)))]
    REMAINING=Ntrajs

    for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D,Ast=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            AST_VEC[:,i].=Ast
            REMAINING+=-1
            println("$(REMAINING) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC,AST_VEC
end

function dist_CM3_FSSH_VANILLA_integration(tf,R0,p0,D0,Ntrajs,DIST,flags=100)
    NDOFs=length(R0)
    R_VEC=zeros(flags+1,NDOFs,Ntrajs)
    P_VEC=zeros(flags+1,NDOFs,Ntrajs)
    D_VEC=zeros(Complex,flags+1,3,Ntrajs)
    AST_VEC=zeros(flags+1,Ntrajs)
    INITIAL=initial_1D_distribution([R0,p0],DIST,Ntrajs)
    S0=[CM3_FSSH_VANILLA_state(INITIAL[1,i],INITIAL[2,i],D0,1) for i in 1:Ntrajs]

    T=collect(0:dt:tf)
    steps=length(T)
    TF=T[Int.(round.(linspace(1,steps,flags+1)))]
    REMAINING=Ntrajs

    for i in 1:Ntrajs
            @time Tf,Rvec,pvec,D,Ast=single_integration(tf,S0[i],flags)
            if length(TF)<length(Tf)
                push!(TF,Tf[end])
            end
            R_VEC[:,:,i].=Rvec
            P_VEC[:,:,i].=pvec
            D_VEC[:,:,i].=D
            AST_VEC[:,i].=Ast
            REMAINING+=-1
            println("$(REMAINING) trajectories remaining")
    end

        return TF,R_VEC,P_VEC,D_VEC,AST_VEC
end
