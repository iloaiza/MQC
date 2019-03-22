function energy(S) #default returns false (for methods that don't track energy)
    return false
end

function energy(S::EH_state)
    E_pot=sum(abs2.(S.el.C).*S.el.E)
    E_kin=sum(abs2.(S.cl.p))/2/mass

    return E_pot+E_kin
end

function energy(S::FSSH_state)
    E_pot=S.el.E[S.ast]
    E_kin=sum(abs2.(S.cl.p))/2/mass

    return E_pot+E_kin
end

function energy(S::FSSH_dia_state)
    E_pot=S.el.E[S.ast,S.ast]
    E_kin=sum(abs2.(S.cl.p))/2/mass

    return E_pot+E_kin
end


function energy(S::BO_state)
    E_pot=S.el.E[1]
    E_kin=sum(abs2.(S.cl.p))/2/mass

    return E_pot+E_kin
end

function energy(S::FRIC_state)
    E_pot=S.el.E[1]
    E_kin=sum(abs2.(S.cl.p))/2/mass
    E_lost=S.cl.mem[end]

    return E_pot+E_kin-E_lost
end

function energy(S::SHEEP_state)
    ast_array=SHEEP_REL[S.ast]
    rho_kk=sum(abs2.(S.el.C[ast_array]))
    Epot=0.0
    for ast in ast_array
        c2=abs2(S.el.C[ast])
        Epot+=c2*S.el.E[ast]
    end
    Epot=Epot/rho_kk
    E_kin=sum(abs2.(S.cl.p))/2/mass

    return E_pot+E_kin-E_lost
end

function energy(S::CM2_state)
    #CM2 doesn't conserve energy
    return false
end

function energy(S::CM3_state)
    #CM3 doesn't conserve energy!
    return false
end

#sanity check subroutine
function health_check(S,E0)
    E = energy(S)
    if E==false #method cannot track energy
        dE=0
    else
        dE = abs(E-E0)/abs(E0)
    end
    if S.prefix in MF_LIST || S.prefix in SH_LIST
        Cnorm=sum(abs2.(S.el.C))
        dnorm = abs(1-Cnorm)
        if dnorm > tol
            dE=tol+1
            @show Cnorm
        end
    end
    if dE>tol
        @show S.cl
        @show E
        @show E0
        @show dE
        error("Warning: energy conservation and/or norm breaking being violated beyond tolerance $tol")
    end
end
