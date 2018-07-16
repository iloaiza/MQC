function phase_tracker(Uold,Ua)
    if Uold==0     #the first run, there is a previous state
        return false,Ua
    else
        phase=round.(diag(Ua'*Uold))
            if  prod(abs.(phase))!=1
            println("Warning, not al phases are +-1")# at position")
            #@show R
            println("Using sign phases instead of round")
            phase=sign.(diag(Ua'*Uold))
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
            Ua=Ua*diagm(phase)
            return true,Ua
        end
    end
end

function BO_state(R,p,NDOFs=length(R))
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)

    return BO_state(R,p,E[1],F[1],NDOFs)
end

function EH_state(R,p,C,Uold=0,NDOFs=length(R)) #Uold reresents the Ua value for the previous state, used for sign consistency
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)
    phasechange,Ua=phase_tracker(Uold,Ua)

    if phasechange
        for i in 1:NDOFs
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+eye(nsts))
            Γ[i]=Γ[i]-diagm(diag(Γ[i]))
        end
    end

    return EH_state(R,p,C,E,Γ,F,Ua,NDOFs)
end

function FSSH_state(R,p,C,ast,Uold=0,NDOFs=length(R)) #Uold reresents the Ua value for the previous state, used for sign consistency
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)
    phasechange,Ua=phase_tracker(Uold,Ua)

    if phasechange
        for i in 1:NDOFs
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+eye(nsts))
            Γ[i]=Γ[i]-diagm(diag(Γ[i]))
        end
    end

    return FSSH_state(R,p,C,E,Γ,F,Ua,ast,NDOFs)
end

function FSSH_dia_state(R,p,C,dst=1,NDOFs=length(R),first=false)
    if first==true
        ~,~,~,~,Ua,~=adiabatic_values(R,NDOFs)
        C=Ua'*C
        Cprob=abs2.(C)
        for i in 2:nsts
            Cprob[i]+=Cprob[i-1]
        end
        chi=rand()
        for i in 1:nsts
            if chi<Cprob[i]
                dst=i
                break
            end
        end
    end

    Hd,dHd=diabatic_values(R)

    return FSSH_dia_state(R,p,C,Hd,dHd,dst,NDOFs)
end

function CM2_VANILLA_state(R,p,D,Uold=0,NDOFs=length(R)) #Uold reresents the Ua value for the previous state, used for sign consistency
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)
    phasechange,Ua=phase_tracker(Uold,Ua)

    if phasechange
        for i in 1:NDOFs
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+eye(nsts))
            Γ[i]=Γ[i]-diagm(diag(Γ[i]))
        end
    end
    z,tnorm,wvec=CM2_additional_values(p,Γ,W,NDOFs)

    F2=[zeros(2,2) for k in 1:NDOFs]

    for k in 1:NDOFs
        F2[k][1,2]=sum(F[k][2:end,1].*tnorm)
        F2[k][2,1]=F2[k][1,2]
        F2[k][1,1]=F[k][1,1]
        F2[k][2,2]=F[k][1,1]
    end

    Weff=sum(wvec.*(tnorm.^2))
    Dmat=[0 -z;z 1im*Weff]


    return CM2_VANILLA_state(R,p,D,E,Ua,F2,z,tnorm,wvec,Dmat,NDOFs)
end

function CM3_VANILLA_state(R,p,D,Uold=0,NDOFs=length(R)) #Uold reresents the Ua value for the previous state, used for sign consistency
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)
    phasechange,Ua=phase_tracker(Uold,Ua)

    if phasechange
        for i in 1:NDOFs
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+eye(nsts))
            Γ[i]=Γ[i]-diagm(diag(Γ[i]))
        end
    end
    z,tnorm,wvec=CM2_additional_values(p,Γ,W,NDOFs)
    zbar,wvec2,tnorm2=CM3_additional_values(wvec,tnorm,z,NDOFs)
    barW=sum(wvec.*(tnorm.^2))
    barW2=sum(wvec2.*(tnorm2.^2))
    Dmat=[0 -z 0;z 1im*barW -1im*zbar;0 -1im*zbar 1im*barW2]

    F3=[zeros(3,3) for k in 1:NDOFs]

    for k in 1:NDOFs
        F3[k][1,2]=sum(F[k][2:end,1].*tnorm)
        F3[k][1,3]=sum((F[k][2,1].*tnorm[2:end]-tnorm[1].*F[k][3:end,1]).*tnorm2)
        #= F3[k][1,2]+=F[k][2,1]*tnorm[1]
        for j in 3:nsts
            #F3[k][1,2]+=F[k][j,1]*tnorm[j-1]
            #F3[k][1,3]+=(F[k][2,1]*tnorm[j-1]-F[k][j,1]*tnorm[1])*tnorm2[j-2]
        end  #   =#
        F3[k][2,1]=F3[k][1,2]
        F3[k][3,1]=F3[k][1,3]
        F3[k][1,1]=F[k][1,1]
        F3[k][2,2]=F[k][1,1]
        F3[k][3,3]=F[k][1,1]
    end

    return CM3_VANILLA_state(R,p,D,E,Ua,F3,z,zbar,tnorm,tnorm2,wvec,wvec2,Dmat,NDOFs)
end

function CM2_FSSH_VANILLA_state(R,p,D,ast,Uold=0,NDOFs=length(R)) #Uold reresents the Ua value for the previous state, used for sign consistency
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)
    phasechange,Ua=phase_tracker(Uold,Ua)

    if phasechange
        for i in 1:NDOFs
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+eye(nsts))
            Γ[i]=Γ[i]-diagm(diag(Γ[i]))
        end
    end
    z,tnorm,wvec=CM2_additional_values(p,Γ,W,NDOFs)

    Fgs=[F[k][1,1] for k in 1:NDOFs]

    Weff=sum(wvec.*(tnorm.^2))
    Dmat=[0 -z;z 1im*Weff]


    return CM2_FSSH_VANILLA_state(R,p,D,E,Ua,Fgs,z,tnorm,Weff,Dmat,ast,NDOFs)
end

function CM3_FSSH_VANILLA_state(R,p,D,ast,Uold=0,NDOFs=length(R)) #Uold reresents the Ua value for the previous state, used for sign consistency
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)
    phasechange,Ua=phase_tracker(Uold,Ua)

    if phasechange
        for i in 1:NDOFs
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+eye(nsts))
            Γ[i]=Γ[i]-diagm(diag(Γ[i]))
        end
    end
    z,tnorm,wvec=CM2_additional_values(p,Γ,W,NDOFs)
    zbar,wvec2,tnorm2=CM3_additional_values(wvec,tnorm,z,NDOFs)
    barW=sum(wvec.*(tnorm.^2))
    barW2=sum(wvec2.*(tnorm2.^2))
    Dmat=[0 -z 0;z 1im*barW -1im*zbar;0 -1im*zbar 1im*barW2]

    Fgs=[F[k][1,1] for k in 1:NDOFs]

    return CM3_FSSH_VANILLA_state(R,p,D,E,Ua,Fgs,z,zbar,tnorm,tnorm2,barW,barW2,Dmat,ast,NDOFs)
end

function CM2_state(R,p,D,Uold=0,NDOFs=length(R)) #Uold reresents the Ua value for the previous state, used for sign consistency
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)
    phasechange,Ua=phase_tracker(Uold,Ua)

    if phasechange
        for i in 1:NDOFs
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+eye(nsts))
            Γ[i]=Γ[i]-diagm(diag(Γ[i]))
        end
    end
    z,tnorm,wvec=CM2_additional_values(p,Γ,W,NDOFs)
    F2=reduced_force(F,tnorm,z,NDOFs,nsts)
    Weff=sum(wvec.*(tnorm.^2))
    Dmat=[0 -z;z 1im*Weff]

    #F_off_diag=z*mass/p.*[0 1;1 0]*((E[2:end]'*tnorm.^2)[1]-E[1])
    #F2+=-F_off_diag

    return CM2_state(R,p,D,E,Ua,F2,z,tnorm,wvec,Dmat,NDOFs)
end

function CM3_state(R,p,D,Uold=0,NDOFs=length(R)) #Uold reresents the Ua value for the previous state, used for sign consistency
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)
    phasechange,Ua=phase_tracker(Uold,Ua)

    if phasechange
        for i in 1:NDOFs
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+eye(nsts))
            Γ[i]=Γ[i]-diagm(diag(Γ[i]))
        end
    end
    z,tnorm,wvec=CM2_additional_values(p,Γ,W,NDOFs)
    zbar,wvec2,tnorm2=CM3_additional_values(wvec,tnorm,z,NDOFs)

    F2=collective_force(F,tnorm,z,NDOFs,nsts)
    F23=reduced_force(gs_remove(F2,NDOFs,nsts),tnorm2,zbar,NDOFs,nsts-1)

    F3=[zeros(3,3) for k in 1:NDOFs]
    for k in 1:NDOFs
        F3[k][1,1:2].=F2[k][1,1:2]
        F3[k][2,1]=F3[k][1,2]
        F3[k][2:3,2:3].=F23[k][1:2,1:2]
        for j in 3:nsts
            F3[k][1,3]+=-tnorm2[j-2]*F2[k][1,j]
        end
        F3[k][3,1]=F3[k][1,3]
    end
    barW=sum(wvec.*(tnorm.^2))
    barW2=sum(wvec2.*(tnorm2.^2))
    Dmat=[0 -z 0;z 1im*barW -1im*zbar;0 -1im*zbar 1im*barW2]


    return CM3_state(R,p,D,E,Ua,F3,z,zbar,tnorm,tnorm2,wvec,wvec2,Dmat,NDOFs)
end

function CM2_FSSH_state(R,p,D,ast,Uold=0,NDOFs=length(R)) #Uold reresents the Ua value for the previous state, used for sign consistency
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)
    phasechange,Ua=phase_tracker(Uold,Ua)

    if phasechange
        for i in 1:NDOFs
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+eye(nsts))
            Γ[i]=Γ[i]-diagm(diag(Γ[i]))
        end
    end
    z,tnorm,wvec=CM2_additional_values(p,Γ,W,NDOFs)
    F2=diagonal_force(F,tnorm,z,NDOFs,nsts)
    Weff=sum(wvec.*(tnorm.^2))
    Dmat=[0 -z;z 1im*Weff]

    return CM2_FSSH_state(R,p,D,E,Ua,F2,z,tnorm,wvec,Dmat,NDOFs,ast)
end

function CM3_FSSH_state(R,p,D,ast,Uold=0,NDOFs=length(R)) #Uold reresents the Ua value for the previous state, used for sign consistency
    E,Γ,F,W,Ua,dHd=adiabatic_values(R,NDOFs)
    phasechange,Ua=phase_tracker(Uold,Ua)

    if phasechange
        for i in 1:NDOFs
            F[i]=Ua'*dHd[i]*Ua
            Γ[i]=-F[i]./(W+eye(nsts))
            Γ[i]=Γ[i]-diagm(diag(Γ[i]))
        end
    end
    z,tnorm,wvec=CM2_additional_values(p,Γ,W,NDOFs)
    zbar,wvec2,tnorm2=CM3_additional_values(wvec,tnorm,z,NDOFs)

    F2=collective_force(F,tnorm,z,NDOFs,nsts)
    F23=diagonal_force(gs_remove(F2,NDOFs,nsts),tnorm2,zbar,NDOFs,nsts-1)

    F3=[zeros(3) for k in 1:NDOFs]
    for k in 1:NDOFs
        F3[k][1]=F2[k][1,1]
        F3[k][2:3].=F23[k][1:2]
    end
    barW=sum(wvec.*(tnorm.^2))
    barW2=sum(wvec2.*(tnorm2.^2))
    Dmat=[0 -z 0;z 1im*barW -1im*zbar;0 -1im*zbar 1im*barW2]


    return CM3_FSSH_state(R,p,D,E,Ua,F3,z,zbar,tnorm,tnorm2,wvec,wvec2,Dmat,NDOFs,ast)
end
