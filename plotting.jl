function energies_4_plot(Rmin,Rmax,res,ex=nsts)
    #the resulting array will have a resolution of res^NDOFs
    global NDOFs=length(Rmin) #noted as N sometimes for short
    if length(res)!=NDOFs
        res=res*ones(Int,NDOFs)
    end
    #build NxNx...xN (N times) matrix with d dimensional vectors in each entry

    if NDOFs!=1
        BASE=[range(Rmin[k],stop=Rmax[k],length=res[k]) for k in 1:NDOFs]
        global R=NxN(BASE,ex)
        forstring=""
        endstring=""
        if NDOFs!=1
            forstring="for i$NDOFs in 1:$(res[NDOFs])   "
            endstring=endstring*"end;   "
            for i in 1:NDOFs-2
                forstring=forstring*"for i$(NDOFs-i) in 1:$(res[NDOFs-i])   "
                endstring=endstring*"end;   "
            end
        end
        forstring=forstring*"for (i1,q) in enumerate(R[:"
        endstring=endstring*"end;   "
        for i in 2:NDOFs
            forstring=forstring*",i$i"
        end
        forstring=forstring*"])  E,Γ,F,W,Ua,dHd=adiabatic_values(q[1:NDOFs],NDOFs); R[i1"
        for i in 2:NDOFs
            forstring=forstring*",i$i"
        end
        forstring=forstring*"][NDOFs+1:NDOFs+$ex].=E;   "
        forstring=forstring*endstring
        eval(Meta.parse(forstring))
    else
        R=zeros(Float64,res,1+ex)
        BASE=linspace(xmin,xmax,res)
        R[:,1].=BASE
        for i in 1:res
            E,Γ,F,W,Ua,dHd=adiabatic_values(R[i,1],1)
            R[i,2:ex+1].=E
        end
    end

    return R
end

function all_4_1D_plot(xmin,xmax,res,basest=1,targetsts=2)
    R=zeros(Float64,res,1+nsts+length(targetsts))
    BASE=range(xmin,stop=xmax,length=res)
    R[:,1].=BASE
    C=zeros(nsts)
    C[1]=1
    S=EH_state_builder(xmin,0,C)
    for i in 1:res
        S=EH_state_builder(R[i,1],0,C,S.el.Ua)
        R[i,2:nsts+1].=S.el.E
        R[i,nsts+2:end].=S.el.Γ[1][basest,targetsts]
    end

    return R
end

function couplings_4_plot(R,basest=1,targetst=2)
    ntargets=length(targetst)
    COUPLINGS=[zeros(size(R)) for i in 1:ntargets]
    C=zeros(nsts)
    C[1]=1
    S=EH_state_builder(R[1],0,C)

    for (i,q) in enumerate(R)
        S=EH_state_builder(q,0,C,S.el.Ua)
        for (ntar,tst) in enumerate(targetst)
            COUPLINGS[ntar][i]=S.el.Γ[1][basest,tst]
        end
    end

    return COUPLINGS
end

function diabatic_1D_plot(Rmin,Rmax,res,basest=1,targetsts=2)
    R=zeros(Float64,res,1+nsts+length(targetsts))
    BASE=range(xmin,stop=xmax,length=res)
    R[:,1].=BASE
    for i in 1:res
        H,dH=diabatic_values(R[i,1])
        R[i,2:nsts+1].=diag(H)
        R[i,nsts+2:end].=H[basest,targetsts]
    end

    return R
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
    for k in 1:length(K)
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
