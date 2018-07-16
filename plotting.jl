function energies_4_plot(Rmin,Rmax,res,ex=nsts)
    #the resulting array will have a resolution of res^NDOFs
    global NDOFs=length(Rmin) #noted as N sometimes for short
    if length(res)!=NDOFs
        res=res*ones(Int,NDOFs)
    end
    #build NxNx...xN (N times) matrix with d dimensional vectors in each entry

    if NDOFs!=1
        BASE=[linspace(Rmin[k],Rmax[k],res[k]) for k in 1:NDOFs]
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
        eval(parse(forstring))
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

function all_4_1D_plot(Rmin,Rmax,res,basest=1,targetsts=2)
    R=zeros(Float64,res,1+nsts+length(targetsts))
    BASE=linspace(xmin,xmax,res)
    R[:,1].=BASE
    for i in 1:res
        E,Γ,F,W,Ua,dHd=adiabatic_values(R[i,1],1)
        R[i,2:nsts+1].=E
        R[i,nsts+2:end].=Γ[1][basest,targetsts]
    end

    return R
end

function couplings_4_plot(R,basest=1,targetst=2)
    ntargets=length(targetst)
    COUPLINGS=[zeros(size(R)) for i in 1:ntargets]
    C=zeros(nsts)
    C[1]=1
    S=EH_state(R[1],0,C)

    for (i,q) in enumerate(R)
        S=EH_state(q,0,C,S.Ua)
        for (ntar,tst) in enumerate(targetst)
            COUPLINGS[ntar][i]=S.Γ[1][basest,tst]
        end
    end

    return COUPLINGS
end
