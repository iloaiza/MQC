function energy(S) #default returns false (for methods that don't track energy)
    return false
end

function energy(S::EH_state)
    Epot=sum(abs2.(S.el.C).*S.el.E)
    Ekin=sum(abs2.(S.cl.p))/2/mass

    return Epot+Ekin
end

function energy(S::FSSH_state)
    Epot=S.el.E[S.ast]
    Ekin=sum(abs2.(S.cl.p))/2/mass

    return Epot+Ekin
end

function energy(S::FSSH_dia_state)
    Epot=S.el.E[S.ast,S.ast]
    Ekin=sum(abs2.(S.cl.p))/2/mass

    return Epot+Ekin
end


function energy(S::BO_state)
    Epot=S.el.E[1]
    Ekin=sum(abs2.(S.cl.p))/2/mass

    return Epot+Ekin
end

function energy(S::FRIC_state)
    Epot=S.el.E[1]
    Ekin=sum(abs2.(S.cl.p))/2/mass
    Elost=S.cl.mem[end]

    return Epot+Ekin-Elost
end

function energy(S::SHEEP_state)
    ast_array=SHEEP_REL[S.ast]
    rho_kk=sum(abs2.(S.el.C[ast_array]))
    Epot=[0.0]
    for ast in ast_array
        c2=abs2(S.el.C[ast])
        Epot[1]+=c2*S.el.E[ast]
    end
    Epot[1]=Epot[1]/rho_kk
    Ekin=sum(abs2.(S.cl.p))/2/mass

    return Epot[1]+Ekin
end

function energy(S::CM2_state)
    E1=S.el.E[1]
    E2=E1+sum(S.CM2.wvec.*(S.CM2.tnorm.^2))
    Epot=sum(abs2.(S.el.C).*[E1,E2])
    Ekin=sum(abs2.(S.cl.p))/2/mass
    @show Epot, Ekin, Epot+Ekin
    return false
end

function energy(S::CM3_state)
    #CM3 doesn't conserve energy!
    return false
end

function energy(S::CM2_FSSH_state)
    if S.ast == 1
        Epot=S.el.E[1]
    else
        Epot=S.el.E[1]+sum([S.el.W[k,1]*(S.CM2.tnorm[k-1]^2) for k in 2:nsts])
    end
    Ekin=sum(abs2.(S.cl.p))/2/mass
    @show Epot, Ekin, Epot+Ekin
    return false
end

function energy(S::CM2_FSSH_FRIC_state)
    Ekin=sum(abs2.(S.cl.p))/2/mass
    Elost=S.cl.mem[end-1]
    if S.ast == 1
        Epot=S.el.E[1]
    else
        Epot=S.cl.mem[end]-Ekin+Elost
    end
    #@show Epot, Ekin, Elost, Epot+Ekin-Elost-S.cl.mem[end]
    return Epot+Ekin-Elost
end

#sanity check subroutine
function health_check(S,E0)
    normbreak=false
    if E0 == false
        dE = 0
    else
        E = energy(S)
        dE = abs(E-E0)/abs(E0)
    end
    if S.prefix in MF_LIST || S.prefix in SH_LIST
        Cnorm=sum(abs2.(S.el.C))
        dnorm = abs(1-Cnorm)
        if dnorm > tol
            normbreak=true
            @show Cnorm
            println("Warning: norm breaking!")
        end
    end
    if dE>tol || normbreak
        @show S.cl
        @show E
        @show E0
        @show dE
        message="Warning! "
        if dE>tol
            message=message*"Energy conservation being violated beyond tolerance $tol. "
        end
        if normbreak
            message=message*"Norm breaking being violated beyond tolerance $tol."
        end
        if sanity_breaks
            error(message)
        else
            println(message)
        end
    end
end
