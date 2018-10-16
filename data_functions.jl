#=
FUNCTIONS FOR DATA SAVING, READING AND DATA ANALYSIS (TULLY-LIKE PLOTS FOR FINAL POPULATIONS)
=#

function CL_save(Tf,Rvec,pvec,file,prefix)
    h5write(file,"T_"*prefix,Tf)
    h5write(file,"R_"*prefix,Rvec)
    h5write(file,"P_"*prefix,pvec)
end


function MF_save(Tf,Rvec,pvec,C,file,prefix)
    Creal=Float64.(real.(C))
    Cim=Float64.(imag.(C))
    h5write(file,"T_"*prefix,Tf)
    h5write(file,"R_"*prefix,Rvec)
    h5write(file,"P_"*prefix,pvec)
    h5write(file,"C_REAL_"*prefix,Creal)
    h5write(file,"C_IMAG_"*prefix,Cim)
end

function SH_save(Tf,Rvec,pvec,C,AST,file,prefix)
    Creal=Float64.(real.(C))
    Cim=Float64.(imag.(C))
    h5write(file,"T_"*prefix,Tf)
    h5write(file,"R_"*prefix,Rvec)
    h5write(file,"P_"*prefix,pvec)
    h5write(file,"C_REAL_"*prefix,Creal)
    h5write(file,"C_IMAG_"*prefix,Cim)
    h5write(file,"AST_"*prefix,AST)
end

function CL_read(file,prefix)
    T=h5read(file,"T_"*prefix)
    R=h5read(file,"R_"*prefix)
    P=h5read(file,"P_"*prefix)
    return T,R,P
end

function MF_read(file,prefix)
    T=h5read(file,"T_"*prefix)
    R=h5read(file,"R_"*prefix)
    P=h5read(file,"P_"*prefix)
    C_REAL=h5read(file,"C_REAL_"*prefix)
    C_IMAG=h5read(file,"C_IMAG_"*prefix)
    C=C_REAL+1im.*C_IMAG
    return T,R,P,C
end

function SH_read(file,prefix)
    T=h5read(file,"T_"*prefix)
    R=h5read(file,"R_"*prefix)
    P=h5read(file,"P_"*prefix)
    C_REAL=h5read(file,"C_REAL_"*prefix)
    C_IMAG=h5read(file,"C_IMAG_"*prefix)
    AST=h5read(file,"AST_"*prefix)
    C=C_REAL+1im.*C_IMAG
    return T,R,P,C,AST
end

function CL_K_save(A,file,prefix)
    h5write(file,"FR_"*prefix,A[1])
    h5write(file,"FP_"*prefix,A[2])
end

function MF_K_save(A,file,prefix)
    h5write(file,"FR_"*prefix,A[1])
    h5write(file,"FM_"*prefix,A[2])
    h5write(file,"Fpop_"*prefix,A[3])
end

function SH_K_save(A,file,prefix)
    h5write(file,"FR_"*prefix,A[1])
    h5write(file,"FM_"*prefix,A[2])
    h5write(file,"Fpop_"*prefix,A[3])
    h5write(file,"Fast_"*prefix,A[4])
end

function general_K_save(A,file,prefix)
    if length(A)==2
        CL_K_save(A,file,prefix)
    elseif length(A)==3
        MF_K_save(A,file,prefix)
    elseif length(A)==4
        SH_K_save(A,file,prefix)
    else
        error("Check our input, it should have 2, 3 or 4 entries!")
    end
end

function general_K_read(file,prefix)
    A=Array[]
    if prefix in METHOD_LIST
        FR=h5read(file,"FR_"*prefix)
        FM=h5read(file,"FM_"*prefix)
        push!(A,FR)
        push!(A,FM)
        if prefix in MF_LIST || prefix in SH_LIST
            Fpop=h5read(file,"Fpop_"*prefix)
            push!(A,Fpop)
            if prefix in SH_LIST
                Fast=h5read(file,"Fast_"*prefix)
                push!(A,Fast)
            end
        end
    end

    return A
end


function refl_trans_SH(K,FR,Fast,crossing_point)
    #this function is used for treating K_SIMULATIONS outputs
    #for tully-like plots. It separates transmitted, reflected ground state
    #and transmitted excited state
    GS_TRANS_FSSH=zeros(size(K))
    GS_REFL_FSSH=zeros(size(K))
    ET_TRANS_FSSH=zeros(size(K))
    for k in 1:length(K)
        for i in 1:Ntrajs
            if FR_FSSH[k,i]>crossing_point && Fast_FSSH[k,i]!=1
                ET_TRANS_FSSH[k]+=1
            elseif FR_FSSH[k,i]>crossing_point && Fast_FSSH[k,i]==1
                GS_TRANS_FSSH[k]+=1
            else
                GS_REFL_FSSH[k]+=1
            end
        end
    end
    GS_TRANS_FSSH=GS_TRANS_FSSH/Ntrajs
    GS_REFL_FSSH=GS_REFL_FSSH/Ntrajs
    ET_TRANS_FSSH=ET_TRANS_FSSH/Ntrajs

    return GS_TRANS_FSSH,GS_REFL_FSSH,ET_TRANS_FSSH
end

function GS_SH_prob(K,Fast)
    GS_PROB=zeros(size(K))
    for k in 1:length(size(K))
        for traj in 1:Ntrajs
            if Fast[traj,k]==1
                GS_PROB[k]+=1
            end
        end
    end

    return GS_PROB./Ntrajs
end

function GS_MF_prob(K,Fpop)
    GS_PROB=zeros(Float64,size(K))
    for k in 1:length(K)
        GS_PROB[k]+=Fpop[k,1]
    end

    return GS_PROB
end

is_true(R,P)=true

function multi_MF_prob(K,A,sts,condition=is_true)
    #sts is an array with all the electronic states to track
    #condition(A) is a function that returns some true value as a function of R and P
    num_sts=length(sts)
    k_size=length(K)
    MULTI_PROB=zeros(Float64,k_size,num_sts)
    FR=A[1]
    FP=A[2]
    Fpop=A[3]

    if length(A)==3
        for st in sts
            for k in 1:k_size
                if condition(FR[k,:],FP[k,:])==true
                    MULTI_PROB[k,st]+=Fpop[k,st]
                end
            end
        end
    elseif length(A)==4
        Ntrajs=length(FR[:,1,1])
        for traj in 1:Ntrajs
            for st in sts
                for k in 1:k_size
                    if condition(FR[traj,k,:],FP[traj,k,:])==true
                        MULTI_PROB[k,st]+=Fpop[traj,k,st]
                    end
                end
            end
        end
        MULTI_PROB=MULTI_PROB/Ntrajs
    end

    return MULTI_PROB

end

function multi_SH_prob(K,A,sts,condition=is_true)
    #sts is an array with all the electronic states to track
    #condition(A) is a function that returns some true value as a function of R and P
    num_sts=length(sts)
    k_size=length(K)
    MULTI_PROB=zeros(Float64,k_size,num_sts)
    FR=A[1]
    FP=A[2]
    Fpop=A[3]
    Fast=A[4]
    Ntrajs=length(FR[:,1,1])

    for traj in 1:Ntrajs
        for st in sts
            for k in 1:k_size
                if Fast[traj,k]==st
                    if condition(FR[traj,k,:],FP[traj,k,:])==true
                        MULTI_PROB[k,st]+=1
                    end
                end
            end
        end
    end

    return MULTI_PROB/Ntrajs

end
