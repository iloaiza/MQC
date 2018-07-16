function EH_energy_check(R,P,C)
    steps=length(R)
    Ep=zeros(steps,nsts) #PES value at each step
    Eeff=zeros(steps) #effective potential energy at each step
    K=P.^2/2/mass #kinetic energy at each step

    for (i,r) in enumerate(R)
        energy,~,~,~,~,~=adiabatic_values(r)
        Ep[i,:].=energy[:]
        Eeff[i].=real.(C[i,:]'*diagm(energy)*C[i,:])
    end

    return Ep,Eeff,K
end

function FSSH_energy_check(R,P,AST)
    steps=length(R)
    Ep=zeros(steps,nsts) #PES value at each step
    Eeff=zeros(steps) #effective potential energy at each step
    K=P.^2/2/mass #kinetic energy at each step

    for (i,r) in enumerate(R)
        energy,~,~,~,~,~=adiabatic_values(r)
        Ep[i,:].=energy[:]
        Eeff[i]=energy[AST[i]]
    end

    return Ep,Eeff,K
end

function EH_sanity_check(R,P,C)
    ~,Eeff,K=EH_energy_check(R,P,C)
    Eini=Eeff[1]+K[1]
    steps=length(R)
    DE=[Eeff[i]+K[i]-Eini for i in 1:steps]
    E_error=maximum(abs.(DE))
    NORM=[sum(abs2.(C[i,:])) for i in 1:steps]
    de_neg=minimum(DE)
    de_pos=maximum(DE)
    if abs(de_neg)>abs(de_pos)
        e_string="The loss of energy is "
    else
        e_string="The gain of energy is "
    end
    E_perc=round(E_error/Eini,4)*100
    norm_max=maximum(NORM)
    norm_min=minimum(NORM)

    println(e_string*"ΔE=$(round(E_error,6)), which represents $E_perc% of the initial energy E0=$(round(Eini,6))")
    println("           (this is a fraction of $(E_perc/100))")
    println("The maximum norm is Nmax=$(round(norm_max,6)), while the minimum is Nmin=$(round(norm_min,6))")
end

function FSSH_sanity_check(R,P,C,AST)
    ~,Eeff,K=FSSH_energy_check(R,P,AST)
    Eini=Eeff[1]+K[1]
    steps=length(R)
    DE=[Eeff[i]+K[i]-Eini for i in 1:steps]
    E_error=maximum(abs.(DE))
    NORM=[sum(abs2.(C[i,:])) for i in 1:steps]
    de_neg=minimum(DE)
    de_pos=maximum(DE)
    if abs(de_neg)>abs(de_pos)
        e_string="The loss of energy is "
    else
        e_string="The gain of energy is "
    end
    E_perc=round(E_error/Eini,4)*100
    norm_max=maximum(NORM)
    norm_min=minimum(NORM)

    HOP_COUNT=0;
    for i in 2:steps
        if AST[i]!=AST[i-1]
            HOP_COUNT+=1
        end
    end

    if HOP_COUNT==0
        println("Warning! There were no hops in this trajectory")
    else
        println("There were $HOP_COUNT hops in this trajectory")
    end

    println(e_string*"ΔE=$(round(E_error,6)), which represents $E_perc% of the initial energy E0=$(round(Eini,6))")
    println("           (this is a fraction of $(E_perc/100))")
    println("The maximum norm is Nmax=$(round(norm_max,6)), while the minimum is Nmin=$(round(norm_min,6))")
end
