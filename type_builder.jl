function phase_tracker(Uold,Ua)
    if Uold==0     #the first run, there is a previous state
        return false,Ua
    else 
        phasechange=false
        phase=zeros(Int,nsts)
        for st in 1:nsts
            if Uold[1,st]!=0
                phase[st]=round(Int,Ua[1,st]/Uold[1,st])
            end
            if phase[st] != 1 && phase[st] !=-1
                #println("First warning! Using full phase tracking with round")
                phase[st]=round.(Int,Diagonal((Ua')*Uold)[st,st])
                if phase[st] != 1 && phase[st] !=-1
                    #println("Second warning! Phase not 1,-1, using sign phase")
                    phase[st]=Int(sign(Diagonal((Ua')*Uold)[st,st]))
                    if phase[st] != 1 && phase[st] !=-1
                        #println("Third warning! Using symmetrized double phase")
                        phase[st]=Int(sign(Diagonal((Ua')*Uold)[st,st]+Diagonal((Uold')*Ua)[st,st]))
                        if phase[st] != 1 && phase[st] !=-1
                            println("Warning! Phase not tracked for st=$st, resetting its value...")
                            phase[st]=1
                        end
                    end
                end
            end
            Ua[:,st] *= phase[st]
            if phase[st]==-1 && !phasechange
                phasechange=true
            end
        end
        

        return phasechange,Ua
    end
end

function C_state_builder(R,p,NDOFs=length(R),mem=0)
    return C_state(R,p,NDOFs,mem)
end

function Q_state_builder(C,R,Uold=0,NDOFs=length(R))
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)

    phasechange,Ua=phase_tracker(Uold,Ua)
    if phasechange
        for i in 1:NDOFs
            F[i]=(Ua')*dHd[i]*Ua
            Γ[i]=-F[i]./(W+Diagonal(ones(nsts)))
            Γ[i]=Γ[i]-Diagonal(Γ[i])
        end
    end

    return Q_state(C,E,W,F,Ua,Γ)
end

function BO_state_builder(R,p,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false)
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(0,R,0,NDOFs)

    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-el.F[k][1,1]
    end
    ODE=ODE_state(p/mass,pdot,0,0)

    return BO_state(cl,el,ODE,"BO",extra)
end

function EH_state_builder(R,p,C,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)

    Rdot=p/mass
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-real(C'*el.F[k]*C)[1]
    end
    NACs=zeros(nsts,nsts)
    for dof in 1:NDOFs
        NACs += Rdot[dof]*el.Γ[dof]
    end
    Cdot=-(1im*Diagonal(el.E)+NACs)
    #@show R,p,sum(abs2.(C))
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return EH_state(cl,el,ODE,"EH",extra)
end

function FSSH_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)

    Rdot=p/mass
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-el.F[k][ast,ast]
    end
    Cdot=-(1im*Diagonal(el.E)+sum(Rdot.*el.Γ))
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return FSSH_state(cl,el,ODE,ast,"FSSH",extra)
end

function FSSH_dia_state_builder(R,p,C,ast=1,NDOFs=length(R),mem=0,extra=Any[],first=false)
    if first==true
        _,_,_,_,Ua,_=adiabatic_values(R,NDOFs)
        C=Ua'*C
        Cprob=abs2.(C)
        for i in 2:nsts
            Cprob[i]+=Cprob[i-1]
        end
        chi=rand()
        for i in 1:nsts
            if chi<Cprob[i]
                ast=i
                break
            end
        end
    end

    cl=C_state_builder(R,p,NDOFs,mem)
    V,dV=diabatic_values(R)
    el=Q_state_builder(C,V,0,dV,0,0)

    Rdot=p/mass
    Cdot=[-1im*sum([V[k,j]*C[j] for j in 1:nsts]) for k in 1:nsts]
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-dV[k][ast,ast]
    end
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return FSSH_dia_state(cl,el,ODE,ast,"FSSH_dia",extra)
end

function CM2_extra(p,Γ,W,NDOFs)
    z,tnorm,wvec=CM2_additional_values(p,Γ,W,NDOFs)
    return CM2_extra(z,tnorm,wvec)
end

function CM3_extra(p,Γ,W,NDOFs)
    z,tnorm,wvec=CM2_additional_values(p,Γ,W,NDOFs)
    zbar,wvec2,tnorm2=CM3_additional_values(wvec,tnorm,z,NDOFs)

    return CM3_extra(z,zbar,tnorm,tnorm2,wvec,wvec2)
end

function CM2_state_builder(R,p,C,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)
    CM2=CM2_extra(p,el.Γ,el.W,NDOFs)

    F2=[zeros(2,2) for k in 1:NDOFs]

    for k in 1:NDOFs
        F2[k][1,2]=sum(el.F[k][2:end,1].*CM2.tnorm)
        F2[k][2,1]=F2[k][1,2]
        F2[k][1,1]=el.F[k][1,1]
        F2[k][2,2]=el.F[k][1,1]
        #F2[k][2,2]=sum([el.F[k][st,st].*(CM2.tnorm[st-1]^2) for st in 2:nsts]) + el.F[k][1,1]
    end

    Weff=sum(CM2.wvec.*(CM2.tnorm.^2))
    Cdot=[0 -CM2.z;CM2.z 1im*Weff]

    Rdot=p/mass
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-real(el.C'*F2[k]*el.C)[1]
    end
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return CM2_state(cl,el,ODE,CM2,"CM2",extra)
end

function CM3_state_builder(R,p,C,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)
    CM3=CM3_extra(p,el.Γ,el.W,NDOFs)

    barW=sum(CM3.wvec.*(CM3.tnorm.^2))
    barW2=sum(CM3.wvec2.*(CM3.tnorm2.^2))
    Cdot=[0 -CM3.z 0;CM3.z 1im*barW -1im*CM3.zbar;0 -1im*CM3.zbar 1im*barW2]

    F3=[zeros(3,3) for k in 1:NDOFs]

    for k in 1:NDOFs
        F3[k][1,2]=sum(el.F[k][2:end,1].*CM3.tnorm)
        F3[k][1,3]=sum((el.F[k][2,1].*CM3.tnorm[2:end]-CM3.tnorm[1].*el.F[k][3:end,1]).*CM3.tnorm2)
        F3[k][2,1]=F3[k][1,2]
        F3[k][3,1]=F3[k][1,3]
        F3[k][1,1]=el.F[k][1,1]
        F3[k][2,2]=el.F[k][1,1]
        F3[k][3,3]=el.F[k][1,1]
    end

    Rdot=p/mass
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-real(C'*F3[k]*C)[1]
    end
    ODE=ODE_state(Rdot,pdot,Cdot,0)


    return CM3_state(cl,el,ODE,CM3,"CM3",extra)
end

function CM2_FSSH_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)
    CM2=CM2_extra(p,el.Γ,el.W,NDOFs)

    Weff=sum(CM2.wvec.*(CM2.tnorm.^2))
    Cdot=[0 -CM2.z;CM2.z 1im*Weff]

    if length(extra)==0 #first run
        extra=[sum(abs2.(p))/2/mass + el.E[1],zeros(nsts-1)] #first entry holds initial energy, second entry holds state for jumping memory
    end

    Rdot=p/mass
    pdot=zeros(NDOFs)

    for k in 1:NDOFs
        if ast==1
            pdot[k]=-el.F[k][1,1]
        else
            for st in 2:nsts
                #pdot[k]+=-el.F[k][st,st]*(extra[2][st-1]^2)
                for st2 in 2:nsts
                    pdot[k]+=-el.F[k][st,st2]*extra[2][st-1]*extra[2][st2-1]
                end
                #pdot[k]+=-el.F[k][st,st]*(CM2.tnorm[st-1]^2)
            end
        end
    end
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return CM2_FSSH_state(cl,el,ODE,CM2,ast,"CM2_FSSH",extra)
end

function CM3_FSSH_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)
    CM3=CM3_extra(p,el.Γ,el.W,NDOFs)

    barW=sum(CM3.wvec.*(CM3.tnorm.^2))
    barW2=sum(CM3.wvec2.*(CM3.tnorm2.^2))
    Cdot=[0 -CM3.z 0;CM3.z 1im*barW -1im*CM3.zbar;0 -1im*CM3.zbar 1im*barW2]

    Rdot=p/mass
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-el.F[k][1]
    end
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return CM3_FSSH_state(cl,el,ODE,CM3,ast,"CM3_FSSH",extra)
end

function CMFSH_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)
    CM2=CM2_extra(p,el.Γ,el.W,NDOFs)

    Weff=sum(CM2.wvec.*(CM2.tnorm.^2))
    Cdot=[0 -CM2.z;CM2.z 1im*Weff]

    Rdot=p/mass
    pdot=zeros(NDOFs)

    E_k=zeros(nsts-2)
    F_k=zeros(nsts-2)

    if mem==0
        mem=zeros(2*(nsts-2)+1) # keeps 2*(nsts-2) for fric memory, end for energy lost due to friction
        cl=C_state_builder(R,p,NDOFs,mem)
        if ast>2
            ast=2
        end
    else
        for i in 1:nsts-2
            E_k[i]=mem[i]
            F_k[i]=mem[i+nsts-2]
        end
    end

    U=zeros(nsts,nsts)
    U[1,1]=1
    for i in 1:nsts-1
        U[2,i+1]=CM2.tnorm[i]
        U[i+1,2]=U[2,i+1]
    end
    for i in 3:nsts
        U[i,i]=-CM2.tnorm[1]
    end

    G=zeros(nsts,nsts,NDOFs)
    for dof in 1:NDOFs
        G[:,:,dof]=(U')*el.Γ[dof]*U
    end

    TAUtilde=zeros(nsts-2)
    for dof in 1:NDOFs
        TAUtilde .+= Rdot[dof]*G[2,3:end,dof]
    end
    Etilde=(U')*Diagonal(el.E)*U
    Wtilde=zeros(nsts-1)
    for st in 2:nsts
        Wtilde[st-1]= Etilde[st,st]-Etilde[2,2]
    end

    E_kdot=zeros(nsts-2)
    F_kdot=zeros(nsts-2)
    memdot=zeros(2*(nsts-2)+1) # add end for tracking disipated energy due to friction
    for i in 1:nsts-2
        E_kdot[i]=TAUtilde[i]*Wtilde[i+1]*abs2(C[2])-Wtilde[i+1]*F_k[i]
        F_kdot[i]=Wtilde[i+1]*E_k[i]
        memdot[i]=E_kdot[i]
        memdot[i+nsts-2]=F_kdot[i]
    end

    pdot_fric=zeros(NDOFs)
    for k in 1:NDOFs
        if ast==1
            pdot[k]=-el.F[k][1,1]
        else
            for st in 2:nsts
                pdot[k]+=-el.F[k][st,st]*extra[st-1]
            end
        end
        for st in 1:nsts-2
            pdot_fric[k]+=-2*G[st+2,1,k]*E_k[st]*abs(C[2])
        end
    end
    pdot+=pdot_fric

    memdot[end]=sum(p .*pdot_fric)/mass #track energy lost due to friction

    ODE=ODE_state(Rdot,pdot,Cdot,memdot)

    return CMFSH_state(cl,el,ODE,CM2,ast,"CMFSH",extra)
end

function CM3_FSSH_FRIC_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)
    CM2=CM2_extra(p,el.Γ,el.W,NDOFs)

    Weff=sum(CM2.wvec.*(CM2.tnorm.^2))
    Cdot=[0 -CM2.z;CM2.z 1im*Weff]

    Rdot=p/mass
    pdot=zeros(NDOFs)

    E_k=zeros(nsts-2)
    F_k=zeros(nsts-2)
    if mem==0 #first run, initialize ek's and fk's
        mem=zeros(2*(nsts-2)+1) # keeps 2*(nsts-2) for fric memory, end-1 for energy lost due to friction and end for total energy (invariant)
        cl=C_state_builder(R,p,NDOFs,mem)
    else
        for i in 1:nsts-2 #could be defined more simply just as zeros, but construction is shown for future reference
            E_k[i]=mem[i]
            F_k[i]=mem[i+nsts-2]
        end
    end

    if length(extra)==0 #first run
        extra=zeros(nsts-1)
    end

    U=zeros(nsts,nsts)
    U[1,1]=1
    for i in 1:nsts-1
        U[2,i+1]=CM2.tnorm[i]
        U[i+1,2]=U[2,i+1]
    end
    for i in 3:nsts
        U[i,i]=-CM2.tnorm[1]
    end

    G=zeros(nsts,nsts,NDOFs)
    for dof in 1:NDOFs
        G[:,:,dof]=(U')*el.Γ[dof]*U
    end

    #TAU=CM2.tnorm .*CM2.z
    TAUtilde=zeros(nsts-2)
    for dof in 1:NDOFs
        TAUtilde .+= Rdot[dof]*G[2,3:end,dof]
    end
    Etilde=(U')*Diagonal(el.E)*U
    Wtilde=zeros(nsts-1)
    for st in 2:nsts
        Wtilde[st-1]= Etilde[st,st]-Etilde[2,2]
    end

    E_kdot=zeros(nsts-2)
    F_kdot=zeros(nsts-2)
    memdot=zeros(2*(nsts-2)+1) # add end for tracking disipated energy due to friction
    for i in 1:nsts-2
        #E_kdot[i]=TAU[i+1]*CM2.wvec[i+1]*abs2(C[2])-CM2.wvec[i+1]*F_k[i]
        E_kdot[i]=TAUtilde[i]*Wtilde[i+1]*abs2(C[2])-Wtilde[i+1]*F_k[i]
        #F_kdot[i]=CM2.wvec[i+1]*E_k[i]
        F_kdot[i]=Wtilde[i+1]*E_k[i]
        memdot[i]=E_kdot[i]
        memdot[i+nsts-2]=F_kdot[i]
    end

    pdot_fric=zeros(NDOFs)
    for k in 1:NDOFs
        if ast==1
            pdot[k]=-el.F[k][1,1]
        else
            for st in 2:nsts
                pdot[k]+=-el.F[k][st,st]*extra[st-1]
            end
        end
        for st in 1:nsts-2
            #pdot_fric[k]+=-2*G[st+2,1,k]*E_k[st]*abs(C[2])
            for dof2 in 1:NDOFs
                pdot_fric[k]+=-2*G[st+2,1,k]*G[st+2,1,dof2]*abs(Wtilde[st+1])*Rdot[dof2]
            end
        end
    end
    pdot+=pdot_fric

    memdot[end]=sum(p .*pdot_fric)/mass #track energy lost due to friction

    ODE=ODE_state(Rdot,pdot,Cdot,memdot)

    return CM3_FSSH_FRIC_state(cl,el,ODE,CM2,ast,"CM3_FSSH_FRIC",extra)
end

function SHEEP_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)

    F=zeros(NDOFs)
    pdot=zeros(NDOFs)

    if ast==1
        for dof in 1:NDOFs
            F[dof]=-el.F[dof][1,1]
        end
    elseif ast>=2
        ast=2
        rho00=abs2(C[1])
        drho=zeros(NDOFs)
        Ees=sum(el.E[2:end].*abs2.(C[2:end]))/(1-rho00)
        for dof in 1:NDOFs
            drho[dof]=2*real(conj(C[1])*sum(el.Γ[dof][1,:].*C))
            Feh=-real(C'*el.F[dof]*C)[1] #Ehrenfest force
            #Fcross=C[1]*sum(conj.(C).*el.F[dof][:,1])+conj(C[1])*sum(C.*el.F[dof][1,:])
            Fgs=-el.F[dof][1,1]
            #F[dof]=real(Feh-abs2(C[1])*Fgs)/(1-abs2(C[1]))
            F[dof]=(Feh-rho00*Fgs+drho[dof]*(Ees-el.E[1]))/(1-rho00)
        end
    end

    for dof in 1:NDOFs
        pdot[dof]=F[dof]
    end

    Rdot=p/mass
    Cdot=-(1im*Diagonal(el.E)+sum(Rdot.*el.Γ))
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return SHEEP_state(cl,el,ODE,ast,"SHEEP",extra)
end

function FRIC_state_builder(R,p,Uold=0,NDOFs=length(R),mem=0,extra=Any[],first=false)    #this is markovian friction with no memory!!
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(0,R,Uold,NDOFs)

    pdot=zeros(NDOFs)

    E_k=zeros(nsts-1)
    F_k=zeros(nsts-1)
    if mem==0 #first run, initialize ek's and fk's
        mem=zeros(2*(nsts-1)+1) #last memory entrance will track lost energy due to friction
    else
        for i in 1:nsts-1 #could be defined more simply just as zeros, but construction is shown for future reference
            E_k[i]=mem[i]
            F_k[i]=mem[i+nsts-1]
        end
    end

    Rdot=p/mass
    E_kdot=zeros(nsts-1)
    F_kdot=zeros(nsts-1)
    memdot=zeros(2*(nsts-1)+1)
    for i in 1:nsts-1
        E_kdot[i]=sum(Rdot.*el.Γ)[i+1,1]*el.W[i+1,1]-el.W[i+1,1]*F_k[i]
        F_kdot[i]=el.W[i+1,1]*E_k[i]
        memdot[i]=E_kdot[i]
        memdot[i+nsts-1]=F_kdot[i]
    end

    pdot_fric=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-el.F[k][1,1]
        for st in 1:nsts-1
            fric=-2*el.Γ[k][st+1,1]*E_k[st]
            pdot_fric[k]+=fric
            pdot[k]+=fric
        end
    end
    memdot[end] = sum(p.*pdot_fric)/mass
    ODE=ODE_state(Rdot,pdot,0,memdot)

    return FRIC_state(cl,el,ODE,"FRIC",extra)
end

################################################################################

function CM2_initial_builder(C)
    return [C[1],sqrt(1-abs2(C[1]))],[]
end
function CMFSH_initial_builder(C)
    D0=[C[1],sqrt(1-abs2(C[1]))]
    return D0,abs2.(C[2:end])
end


function init_builder(DYN,C)
    if eval(Meta.parse("@isdefined $(DYN)_initial_builder"))
        return eval(Meta.parse("$(DYN)_initial_builder($C)"))
    else
        println("Warning, haven't defined transformation for initial electonic condition of $DYN, using [1,0,...] instead")
        D0=zeros(eval(Meta.parse("$(DYN)_sts")))
        D0[1]=1
        return D0,[]
    end
end