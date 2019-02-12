function runge_step(S::CL_state,δt=dt,NDOFs=S.cl.NDOFs)
    R0=S.cl.R
    p0=S.cl.p
    mem0=S.cl.mem

    dR1,dp1,dmem1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kmem1=δt*dmem1

    R2=R0+kR1/2
    p2=p0+kp1/2
    mem2=mem0+kmem1/2
    S=builder_CL_state(R2,p2,S.prefix,NDOFs,mem2)
    dR2,dp2,dmem2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kmem2=δt*dmem2

    R3=R0+kR2/2
    p3=p0+kp2/2
    mem3=mem0+kmem2/2
    S=builder_CL_state(R3,p3,S.prefix,NDOFs,mem3)
    dR3,dp3,dmem3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kmem3=δt*dmem3

    R4=R0+kR3
    p4=p0+kp3
    mem4=mem0+kmem3
    S=builder_CL_state(R4,p4,S.prefix,NDOFs,mem4)
    dR4,dp4,dmem4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kmem4=δt*dmem4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    memnew=mem0+kmem1/6+kmem2/3+kmem3/3+kmem4/6


    return builder_CL_state(Rnew,pnew,S.prefix,NDOFs,memnew)
end

function runge_step(S::MF_state,δt=dt,NDOFs=S.cl.NDOFs)
    R0=S.cl.R
    p0=S.cl.p
    C0=S.el.C
    mem0=S.cl.mem

    dR1,dp1,dC1,dmem1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kC1=δt*dC1
    kmem1=δt*dmem1

    R2=R0+kR1/2
    p2=p0+kp1/2
    C2=C0+kC1/2
    mem2=mem0+kmem1/2
    S=builder_MF_state(R2,p2,C2,S.prefix,S.el.Ua,NDOFs,mem2)
    dR2,dp2,dC2,dmem2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kC2=δt*dC2
    kmem2=δt*dmem2

    R3=R0+kR2/2
    p3=p0+kp2/2
    C3=C0+kC2/2
    mem3=mem0+kmem2/2
    S=builder_MF_state(R3,p3,C3,S.prefix,S.el.Ua,NDOFs,mem3)
    dR3,dp3,dC3,dmem3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kC3=δt*dC3
    kmem3=δt*dmem3

    R4=R0+kR3
    p4=p0+kp3
    C4=C0+kC3
    mem4=mem0+kmem3
    S=builder_MF_state(R4,p4,C4,S.prefix,S.el.Ua,NDOFs,mem4)
    dR4,dp4,dC4,dmem4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kC4=δt*dC4
    kmem4=δt*dmem4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Cnew=C0+kC1/6+kC2/3+kC3/3+kC4/6
    memnew=mem0+kmem1/6+kmem2/3+kmem3/3+kmem4/6

    return builder_MF_state(Rnew,pnew,Cnew,S.prefix,S.el.Ua,NDOFs,memnew)
end

function runge_step(S::SH_state,δt=dt,NDOFs=S.cl.NDOFs)
    R0=S.cl.R
    p0=S.cl.p
    C0=S.el.C
    mem0=S.cl.mem

    dR1,dp1,dC1,dmem1=diff_eq(S)
    kR1=δt*dR1
    kp1=δt*dp1
    kC1=δt*dC1
    kmem1=δt*dmem1

    R2=R0+kR1/2
    p2=p0+kp1/2
    C2=C0+kC1/2
    mem2=mem0+kmem1/2
    S=builder_SH_state(R2,p2,C2,S.ast,S.prefix,S.el.Ua,NDOFs,mem2)
    dR2,dp2,dC2,dmem2=diff_eq(S)
    kR2=δt*dR2
    kp2=δt*dp2
    kC2=δt*dC2
    kmem2=δt*dmem2

    R3=R0+kR2/2
    p3=p0+kp2/2
    C3=C0+kC2/2
    mem3=mem0+kmem2/2
    S=builder_SH_state(R3,p3,C3,S.ast,S.prefix,S.el.Ua,NDOFs,mem3)
    dR3,dp3,dC3,dmem3=diff_eq(S)
    kR3=δt*dR3
    kp3=δt*dp3
    kC3=δt*dC3
    kmem3=δt*dmem3

    R4=R0+kR3
    p4=p0+kp3
    C4=C0+kC3
    mem4=mem0+kmem3
    S=builder_SH_state(R4,p4,C4,S.ast,S.prefix,S.el.Ua,NDOFs,mem4)
    dR4,dp4,dC4,dmem4=diff_eq(S)
    kR4=δt*dR4
    kp4=δt*dp4
    kC4=δt*dC4
    kmem4=δt*dmem4

    Rnew=R0+kR1/6+kR2/3+kR3/3+kR4/6
    pnew=p0+kp1/6+kp2/3+kp3/3+kp4/6
    Cnew=C0+kC1/6+kC2/3+kC3/3+kC4/6
    memnew=mem0+kmem1/6+kmem2/3+kmem3/3+kmem4/6

    return builder_SH_state(Rnew,pnew,Cnew,S.ast,S.prefix,S.el.Ua,NDOFs,memnew)
end
