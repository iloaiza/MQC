function phase_tracker(Uold,Ua)
    if Uold==0     #the first run, there is a previous state
        return false,Ua
    else
        phase=round.(Int,Diagonal(Ua'*Uold)[diagind(Ua)])
            if  prod(abs.(phase))!=1
            println("Warning, not al phases are +-1")# at position")
            #@show R
            println("Using sign phases instead of round. Consider changing timestep if warning is repeated many times...")
            phase=sign.(Diagonal(Ua'*Uold)[diagind(Ua)])
            if prod(phase)==0
                println("Second warning: the phase sign has a zero, resetting phase on those values")
                for (i,p) in enumerate(phase)
                    if p==0
                        phase[i]=1
                    end
                end
            end
        end

        phasechange=0
        for i in phase
            if i==-1
                phasechange=1
                break
            end
        end

        if phasechange==0
            return false,Ua
        else
            Ua=Ua*Diagonal(phase)
            return true,Ua
        end
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
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+Diagonal(ones(nsts)))
            Γ[i]=Γ[i]-Diagonal(Γ[i])
        end
    end

    return Q_state(C,E,W,F,Ua,Γ)
end

function BO_state_builder(R,p,Uold=0,NDOFs=length(R),mem=0)
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(0,R,0,NDOFs)

    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-el.F[k][1]
    end
    ODE=ODE_state(p/mass,pdot,0,0)

    return BO_state(cl,el,ODE,"BO")
end

function EH_state_builder(R,p,C,Uold=0,NDOFs=length(R),mem=0) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)

    Rdot=p/mass
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-real(el.C'*el.F[k]*el.C)[1]
    end
    Cdot=-(1im*Diagonal(el.E)+sum(Rdot.*el.Γ))
    ODE=ODE_state(p/mass,pdot,Cdot,0)

    return EH_state(cl,el,ODE,"EH")
end

function FSSH_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)

    Rdot=p/mass
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-el.F[k][ast,ast]
    end
    Cdot=-(1im*Diagonal(el.E)+sum(Rdot.*el.Γ))
    ODE=ODE_state(p/mass,pdot,Cdot,0)

    return FSSH_state(cl,el,ODE,ast,"FSSH")
end

function FSSH_dia_state_builder(R,p,C,ast=1,NDOFs=length(R),mem=0,first=false)
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
    Cdot=[-1im*sum(V[k,j]*C[j] for j in 1:nsts) for k in 1:nsts]
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-dV[k][ast,ast]
    end
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return FSSH_dia_state(cl,el,ODE,ast,"FSSH_dia")
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

function CM2_state_builder(R,p,C,Uold=0,NDOFs=length(R),mem=0) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)
    CM2=CM2_extra(p,el.Γ,el.W,NDOFs)

    F2=[zeros(2,2) for k in 1:NDOFs]

    for k in 1:NDOFs
        F2[k][1,2]=sum(el.F[k][2:end,1].*CM2.tnorm)
        F2[k][2,1]=F2[k][1,2]
        F2[k][1,1]=el.F[k][1,1]
        F2[k][2,2]=el.F[k][1,1]
    end

    Weff=sum(CM2.wvec.*(CM2.tnorm.^2))
    Cdot=[0 -CM2.z;CM2.z 1im*Weff]

    Rdot=p/mass
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-real(el.C'*F2[k]*el.C)[1]
    end
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return CM2_state(cl,el,ODE,CM2,"CM2")
end

function CM3_state_builder(R,p,C,Uold=0,NDOFs=length(R),mem=0) #Uold reresents the Ua value for the previous state, used for sign consistency
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


    return CM3_state(cl,el,ODE,CM3,"CM3")
end

function CM2_FSSH_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)
    CM2=CM2_extra(p,el.Γ,el.W,NDOFs)

    Weff=sum(CM2.wvec.*(CM2.tnorm.^2))
    Cdot=[0 -CM2.z;CM2.z 1im*Weff]

    Rdot=p/mass
    pdot=zeros(NDOFs)
    for k in 1:NDOFs
        pdot[k]=-el.F[k][1]
    end
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return CM2_FSSH_state(cl,el,ODE,CM2,ast,"CM2_FSSH")
end

function CM3_FSSH_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0) #Uold reresents the Ua value for the previous state, used for sign consistency
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

    return CM3_FSSH_state(cl,el,ODE,CM3,ast,"CM3_FSSH")
end

function CM2_FSSH_FRIC_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0) #Uold reresents the Ua value for the previous state, used for sign consistency
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
        mem=zeros(2*(nsts-2))
        cl=C_state_builder(R,p,NDOFs,mem)
    else
        for i in 1:nsts-2 #could be defined more simply just as zeros, but construction is shown for future reference
            E_k[i]=mem[i]
            F_k[i]=mem[i+nsts-2]
        end
    end

    TAU=CM2.tnorm.*CM2.z
    E_kdot=zeros(nsts-2)
    F_kdot=zeros(nsts-2)
    memdot=zeros(2*(nsts-2))
    for i in 1:nsts-2
        E_kdot[i]=TAU[i+1]*CM2.wvec[i+1]*abs2(C[2])-CM2.wvec[i+1]*F_k[i]
        F_kdot[i]=CM2.wvec[i+1]*E_k[i]
        memdot[i]=E_kdot[i]
        memdot[i+nsts-2]=F_kdot[i]
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
    Udag=U'

    G=zeros(nsts,nsts,NDOFs)
    for dof in 1:NDOFs
        G[:,:,NDOFs]=Udag*el.Γ[dof]*U
    end


    for k in 1:NDOFs
        pdot[k]=-el.F[k][1]
        for st in 1:nsts-2
            pdot[k]+=-2*G[st+2,1,k]*E_k[st]*abs(C[2])
        end
    end
    ODE=ODE_state(Rdot,pdot,Cdot,memdot)

    return CM2_FSSH_FRIC_state(cl,el,ODE,CM2,ast,"CM2_FSSH_FRIC")
end

function CM3_FSSH_FRIC_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)
    CM3=CM3_extra(p,el.Γ,el.W,NDOFs)

    barW=sum(CM3.wvec.*(CM3.tnorm.^2))
    barW2=sum(CM3.wvec2.*(CM3.tnorm2.^2))
    Cdot=[0 -CM3.z 0;CM3.z 1im*barW -1im*CM3.zbar;0 -1im*CM3.zbar 1im*barW2]

    Rdot=p/mass
    pdot=zeros(NDOFs)

    E_k=zeros(nsts-3)
    F_k=zeros(nsts-3)
    if mem==0 #first run, initialize ek's and fk's
        mem=zeros(2*(nsts-3))
        cl=C_state_builder(R,p,NDOFs,mem)
    else
        for i in 1:nsts-3 #could be defined more simply just as zeros, but construction is shown for future reference
            E_k[i]=mem[i]
            F_k[i]=mem[i+nsts-3]
        end
    end

    TAU=CM3.tnorm2.*CM3.zbar
    E_kdot=zeros(nsts-3)
    F_kdot=zeros(nsts-3)
    memdot=zeros(2*(nsts-3))
    for i in 1:nsts-3
        E_kdot[i]=TAU[i+1]*CM3.wvec2[i+1]*abs2(C[3])-CM3.wvec[i+1]*F_k[i]
        F_kdot[i]=CM3.wvec2[i+1]*E_k[i]
        memdot[i]=E_kdot[i]
        memdot[i+nsts-3]=F_kdot[i]
    end

    U2=zeros(nsts,nsts)
    U3=zeros(nsts,nsts)
    U2[1,1]=1
    for i in 1:nsts-1
        U2[2,i+1]=CM3.tnorm[i]
        U2[i+1,2]=U2[2,i+1]
    end
    for i in 3:nsts
        U2[i,i]=-CM3.tnorm[1]
    end
    U3[1,1]=CM3.zbar
    U3[2,2]=1
    for i in 1:nsts-2
        U3[3,i+2]=CM3.tnorm2[i]
        U3[i+2,3]=U3[3,i+2]
    end
    for i in 4:nsts
        U3[i,i]=-CM3.tnorm2[1]
    end
    U=U2*U3
    Udag=U'

    G=zeros(nsts,nsts,NDOFs)
    for dof in 1:NDOFs
        G[:,:,NDOFs]=Udag*el.Γ[dof]*U
    end

    for k in 1:NDOFs
        pdot[k]=-el.F[k][1]
        for st in 1:nsts-3
            pdot[k]+=-2*G[st+3,2,k]*E_k[st]*abs(C[3])
        end
    end

    ODE=ODE_state(Rdot,pdot,Cdot,memdot)

    return CM3_FSSH_FRIC_state(cl,el,ODE,CM3,ast,"CM3_FSSH_FRIC")
end

function SHEEP_state_builder(R,p,C,ast,Uold=0,NDOFs=length(R),mem=0) #Uold reresents the Ua value for the previous state, used for sign consistency
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(C,R,Uold,NDOFs)

    F=zeros(NDOFs)
    pdot=zeros(NDOFs)
    for dof in 1:NDOFs
        avg_pop=0
        for avg_st in SHEEP_REL[ast]
            c2=abs2(C[avg_st])
            avg_pop+=c2
            F[dof]+=c2*el.F[dof][avg_st,avg_st]
        end
        F[dof]=F[dof]/avg_pop
        pdot[dof]=-F[dof]
    end

    Rdot=p/mass
    Cdot=-(1im*Diagonal(el.E)+sum(Rdot.*el.Γ))
    ODE=ODE_state(Rdot,pdot,Cdot,0)

    return SHEEP_state(cl,el,ODE,ast,"SHEEP")
end

#IMPLEMENTATION UNDER WAY!!
function FRIC_state_builder(R,p,Uold=0,NDOFs=length(R),mem=0)    #this is markovian friction with no memory!!
    cl=C_state_builder(R,p,NDOFs,mem)
    el=Q_state_builder(0,R,Uold,NDOFs)

    pdot=zeros(NDOFs)

    E_k=zeros(nsts-1)
    F_k=zeros(nsts-1)
    if mem==0 #first run, initialize ek's and fk's
        mem=zeros(2*(nsts-1))
    else
        for i in 1:nsts-1 #could be defined more simply just as zeros, but construction is shown for future reference
            E_k[i]=mem[i]
            F_k[i]=mem[i+nsts-1]
        end
    end

    Rdot=p/mass
    E_kdot=zeros(nsts-1)
    F_kdot=zeros(nsts-1)
    memdot=zeros(2*(nsts-1))
    for i in 1:nsts-1
        E_kdot[i]=sum(Rdot.*el.Γ)[i+1,1]*el.W[i+1,1]-el.W[i+1,1]*F_k[i]
        F_kdot[i]=el.W[i+1,1]*E_k[i]
        memdot[i]=E_kdot[i]
        memdot[i+nsts-1]=F_kdot[i]
    end


    for k in 1:NDOFs
        pdot[k]=-el.F[k][1]
        for st in 1:nsts-1
            pdot[k]+=-2*el.Γ[k][st+1,1]*E_k[st]
        end
    end
    ODE=ODE_state(Rdot,pdot,0,memdot)

    return FRIC_state(cl,el,ODE,"FRIC")
end

################################################################################

function CM2_initial_builder(R,p,C)
    D0=[C[1],sqrt(1-abs2(C[1]))]
end

function CM3_initial_builder(R,p,C)
    D0=[C[1],sqrt(1-abs2(C[1]))]
end
