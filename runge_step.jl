function runge_step(S::BO_state,δt=dt,NDOFs=S.NDOFs)
    R0=S.R
    p0=S.p

    dR1,dp1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1

    R2=R0+kR1/2
    p2=p0+kp1/2
    S=BO_state(R2,p2,NDOFs)
    dR2,dp2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2

    R3=R0+kR2/2
    p3=p0+kp2/2
    S=BO_state(R3,p3,NDOFs)
    dR3,dp3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3

    R4=R0+kR3
    p4=p0+kp3
    S=BO_state(R4,p4,NDOFs)
    dR4,dp4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6

    return BO_state(Rnew,pnew,NDOFs)
end

function runge_step(S::EH_state,δt=dt,NDOFs=S.NDOFs)
    R0=S.R
    p0=S.p
    C0=S.C

    dR1,dp1,dC1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kC1=δt*dC1

    R2=R0+kR1/2
    p2=p0+kp1/2
    C2=C0+kC1/2
    S=EH_state(R2,p2,C2,S.Ua,NDOFs)
    dR2,dp2,dC2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kC2=δt*dC2

    R3=R0+kR2/2
    p3=p0+kp2/2
    C3=C0+kC2/2
    S=EH_state(R3,p3,C3,S.Ua,NDOFs)
    dR3,dp3,dC3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kC3=δt*dC3

    R4=R0+kR3
    p4=p0+kp3
    C4=C0+kC3
    S=EH_state(R4,p4,C4,S.Ua,NDOFs)
    dR4,dp4,dC4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kC4=δt*dC4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Cnew=C0+kC1/6+kC2/3+kC3/3+kC4/6

    return EH_state(Rnew,pnew,Cnew,S.Ua,NDOFs)
end

function runge_step(S::FSSH_state,δt=dt,NDOFs=S.NDOFs)
    R0=S.R
    p0=S.p
    C0=S.C

    dR1,dp1,dC1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kC1=δt*dC1

    R2=R0+kR1/2
    p2=p0+kp1/2
    C2=C0+kC1/2
    S=FSSH_state(R2,p2,C2,S.ast,S.Ua,NDOFs)
    dR2,dp2,dC2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kC2=δt*dC2

    R3=R0+kR2/2
    p3=p0+kp2/2
    C3=C0+kC2/2
    S=FSSH_state(R3,p3,C3,S.ast,S.Ua,NDOFs)
    dR3,dp3,dC3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kC3=δt*dC3

    R4=R0+kR3
    p4=p0+kp3
    C4=C0+kC3
    S=FSSH_state(R4,p4,C4,S.ast,S.Ua,NDOFs)
    dR4,dp4,dC4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kC4=δt*dC4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Cnew=C0+kC1/6+kC2/3+kC3/3+kC4/6

    return FSSH_state(Rnew,pnew,Cnew,S.ast,S.Ua)
end

function runge_step(S::FSSH_dia_state,δt=dt,NDOFs=S.NDOFs)
    R0=S.R
    p0=S.p
    C0=S.C

    dR1,dp1,dC1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kC1=δt*dC1

    R2=R0+kR1/2
    p2=p0+kp1/2
    C2=C0+kC1/2
    S=FSSH_dia_state(R2,p2,C2,S.dst,NDOFs)
    dR2,dp2,dC2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kC2=δt*dC2

    R3=R0+kR2/2
    p3=p0+kp2/2
    C3=C0+kC2/2
    S=FSSH_dia_state(R3,p3,C3,S.dst,NDOFs)
    dR3,dp3,dC3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kC3=δt*dC3

    R4=R0+kR3
    p4=p0+kp3
    C4=C0+kC3
    S=FSSH_dia_state(R4,p4,C4,S.dst,NDOFs)
    dR4,dp4,dC4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kC4=δt*dC4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Cnew=C0+kC1/6+kC2/3+kC3/3+kC4/6

    return FSSH_dia_state(Rnew,pnew,Cnew,S.dst)
end

function runge_step(S::CM2_VANILLA_state,δt=dt)
    R0=S.R
    p0=S.p
    D0=S.D

    dR1,dp1,dD1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kD1=δt*dD1

    R2=R0+kR1/2
    p2=p0+kp1/2
    D2=D0+kD1/2
    S=CM2_VANILLA_state(R2,p2,D2,S.Ua,S.NDOFs)
    dR2,dp2,dD2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kD2=δt*dD2

    R3=R0+kR2/2
    p3=p0+kp2/2
    D3=D0+kD2/2
    S=CM2_VANILLA_state(R3,p3,D3,S.Ua,S.NDOFs)
    dR3,dp3,dD3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kD3=δt*dD3

    R4=R0+kR3
    p4=p0+kp3
    D4=D0+kD3
    S=CM2_VANILLA_state(R4,p4,D4,S.Ua,S.NDOFs)
    dR4,dp4,dD4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kD4=δt*dD4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Dnew=D0+kD1/6+kD2/3+kD3/3+kD4/6

    return CM2_VANILLA_state(Rnew,pnew,Dnew,S.Ua,S.NDOFs)
end

function runge_step(S::CM3_VANILLA_state,δt=dt)
    R0=S.R
    p0=S.p
    D0=S.D

    dR1,dp1,dD1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kD1=δt*dD1

    R2=R0+kR1/2
    p2=p0+kp1/2
    D2=D0+kD1/2
    S=CM3_VANILLA_state(R2,p2,D2,S.Ua,S.NDOFs)
    dR2,dp2,dD2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kD2=δt*dD2

    R3=R0+kR2/2
    p3=p0+kp2/2
    D3=D0+kD2/2
    S=CM3_VANILLA_state(R3,p3,D3,S.Ua,S.NDOFs)
    dR3,dp3,dD3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kD3=δt*dD3

    R4=R0+kR3
    p4=p0+kp3
    D4=D0+kD3
    S=CM3_VANILLA_state(R4,p4,D4,S.Ua,S.NDOFs)
    dR4,dp4,dD4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kD4=δt*dD4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Dnew=D0+kD1/6+kD2/3+kD3/3+kD4/6

    return CM3_VANILLA_state(Rnew,pnew,Dnew,S.Ua,S.NDOFs)
end

function runge_step(S::CM2_state,δt=dt)
    R0=S.R
    p0=S.p
    D0=S.D

    dR1,dp1,dD1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kD1=δt*dD1

    R2=R0+kR1/2
    p2=p0+kp1/2
    D2=D0+kD1/2
    S=CM2_state(R2,p2,D2,S.Ua,S.NDOFs)
    dR2,dp2,dD2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kD2=δt*dD2

    R3=R0+kR2/2
    p3=p0+kp2/2
    D3=D0+kD2/2
    S=CM2_state(R3,p3,D3,S.Ua,S.NDOFs)
    dR3,dp3,dD3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kD3=δt*dD3

    R4=R0+kR3
    p4=p0+kp3
    D4=D0+kD3
    S=CM2_state(R4,p4,D4,S.Ua,S.NDOFs)
    dR4,dp4,dD4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kD4=δt*dD4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Dnew=D0+kD1/6+kD2/3+kD3/3+kD4/6

    return CM2_state(Rnew,pnew,Dnew,S.Ua,S.NDOFs)
end

function runge_step(S::CM3_state,δt=dt)
    R0=S.R
    p0=S.p
    D0=S.D

    dR1,dp1,dD1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kD1=δt*dD1

    R2=R0+kR1/2
    p2=p0+kp1/2
    D2=D0+kD1/2
    S=CM3_state(R2,p2,D2,S.Ua,S.NDOFs)
    dR2,dp2,dD2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kD2=δt*dD2

    R3=R0+kR2/2
    p3=p0+kp2/2
    D3=D0+kD2/2
    S=CM3_state(R3,p3,D3,S.Ua,S.NDOFs)
    dR3,dp3,dD3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kD3=δt*dD3

    R4=R0+kR3
    p4=p0+kp3
    D4=D0+kD3
    S=CM3_state(R4,p4,D4,S.Ua,S.NDOFs)
    dR4,dp4,dD4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kD4=δt*dD4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Dnew=D0+kD1/6+kD2/3+kD3/3+kD4/6

    return CM3_state(Rnew,pnew,Dnew,S.Ua,S.NDOFs)
end

function runge_step(S::CM2_FSSH_state,δt=dt)
    R0=S.R
    p0=S.p
    D0=S.D

    dR1,dp1,dD1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kD1=δt*dD1

    R2=R0+kR1/2
    p2=p0+kp1/2
    D2=D0+kD1/2
    S=CM2_FSSH_state(R2,p2,D2,S.ast,S.Ua,S.NDOFs)
    dR2,dp2,dD2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kD2=δt*dD2

    R3=R0+kR2/2
    p3=p0+kp2/2
    D3=D0+kD2/2
    S=CM2_FSSH_state(R3,p3,D3,S.ast,S.Ua,S.NDOFs)
    dR3,dp3,dD3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kD3=δt*dD3

    R4=R0+kR3
    p4=p0+kp3
    D4=D0+kD3
    S=CM2_FSSH_state(R4,p4,D4,S.ast,S.Ua,S.NDOFs)
    dR4,dp4,dD4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kD4=δt*dD4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Dnew=D0+kD1/6+kD2/3+kD3/3+kD4/6

    return CM2_FSSH_state(Rnew,pnew,Dnew,S.ast,S.Ua,S.NDOFs)
end

function runge_step(S::CM3_FSSH_state,δt=dt)
    R0=S.R
    p0=S.p
    D0=S.D

    dR1,dp1,dD1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kD1=δt*dD1

    R2=R0+kR1/2
    p2=p0+kp1/2
    D2=D0+kD1/2
    S=CM3_FSSH_state(R2,p2,D2,S.ast,S.Ua,S.NDOFs)
    dR2,dp2,dD2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kD2=δt*dD2

    R3=R0+kR2/2
    p3=p0+kp2/2
    D3=D0+kD2/2
    S=CM3_FSSH_state(R3,p3,D3,S.ast,S.Ua,S.NDOFs)
    dR3,dp3,dD3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kD3=δt*dD3

    R4=R0+kR3
    p4=p0+kp3
    D4=D0+kD3
    S=CM3_FSSH_state(R4,p4,D4,S.ast,S.Ua,S.NDOFs)
    dR4,dp4,dD4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kD4=δt*dD4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Dnew=D0+kD1/6+kD2/3+kD3/3+kD4/6

    return CM3_FSSH_state(Rnew,pnew,Dnew,S.ast,S.Ua,S.NDOFs)
end

function runge_step(S::CM2_FSSH_VANILLA_state,δt=dt)
    R0=S.R
    p0=S.p
    D0=S.D

    dR1,dp1,dD1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kD1=δt*dD1

    R2=R0+kR1/2
    p2=p0+kp1/2
    D2=D0+kD1/2
    S=CM2_FSSH_VANILLA_state(R2,p2,D2,S.ast,S.Ua,S.NDOFs)
    dR2,dp2,dD2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kD2=δt*dD2

    R3=R0+kR2/2
    p3=p0+kp2/2
    D3=D0+kD2/2
    S=CM2_FSSH_VANILLA_state(R3,p3,D3,S.ast,S.Ua,S.NDOFs)
    dR3,dp3,dD3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kD3=δt*dD3

    R4=R0+kR3
    p4=p0+kp3
    D4=D0+kD3
    S=CM2_FSSH_VANILLA_state(R4,p4,D4,S.ast,S.Ua,S.NDOFs)
    dR4,dp4,dD4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kD4=δt*dD4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Dnew=D0+kD1/6+kD2/3+kD3/3+kD4/6

    return CM2_FSSH_VANILLA_state(Rnew,pnew,Dnew,S.ast,S.Ua,S.NDOFs)
end

function runge_step(S::CM3_FSSH_VANILLA_state,δt=dt)
    R0=S.R
    p0=S.p
    D0=S.D

    dR1,dp1,dD1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kD1=δt*dD1

    R2=R0+kR1/2
    p2=p0+kp1/2
    D2=D0+kD1/2
    S=CM3_FSSH_VANILLA_state(R2,p2,D2,S.ast,S.Ua,S.NDOFs)
    dR2,dp2,dD2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kD2=δt*dD2

    R3=R0+kR2/2
    p3=p0+kp2/2
    D3=D0+kD2/2
    S=CM3_FSSH_VANILLA_state(R3,p3,D3,S.ast,S.Ua,S.NDOFs)
    dR3,dp3,dD3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kD3=δt*dD3

    R4=R0+kR3
    p4=p0+kp3
    D4=D0+kD3
    S=CM3_FSSH_VANILLA_state(R4,p4,D4,S.ast,S.Ua,S.NDOFs)
    dR4,dp4,dD4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kD4=δt*dD4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Dnew=D0+kD1/6+kD2/3+kD3/3+kD4/6

    return CM3_FSSH_VANILLA_state(Rnew,pnew,Dnew,S.ast,S.Ua,S.NDOFs)
end
