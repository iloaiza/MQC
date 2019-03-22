function runge_state_builder(S,y)
    R,p,C,mem=y[1],y[2],y[3],y[4]
    if S.prefix in CL_LIST
        S=builder_CL_state(R,p,S.prefix,S.el.Ua,S.cl.NDOFs,mem)
    elseif S.prefix in MF_LIST
        S=builder_MF_state(R,p,C,S.prefix,S.el.Ua,S.cl.NDOFs,mem)
    else
        S=builder_SH_state(R,p,C,S.ast,S.prefix,S.el.Ua,S.cl.NDOFs,mem)
    end

    return S
end

function runge_ks(S0,tstep)
    if S0.cl.NDOFs==1
        y0=[S0.cl.R[1],S0.cl.p[1],S0.el.C,S0.cl.mem]
    else
        y0=[S0.cl.R,S0.cl.p,S0.el.C,S0.cl.mem]
    end
    S=runge_state_builder(S0,y0) #make a copy, won't change S0 until necessary

    k1=tstep*diff_eq(S)
    y1=y0.+k1/4
    S=runge_state_builder(S,y1)

    k2=tstep*diff_eq(S)
    y2=y0 .+ 3k1/32 .+ 9k2/32
    S=runge_state_builder(S,y2)

    k3=tstep*diff_eq(S)
    y3=y0 .+ 1932k1/2197 .- 7200k2/2197 .+ 7296k3/2197
    S=runge_state_builder(S,y3)

    k4=tstep*diff_eq(S)
    y4=y0 .+ 439k1/216 .- 8k2 .+ 3680k3/513 .- 845k4/4104
    S=runge_state_builder(S,y4)

    k5=tstep*diff_eq(S)
    y5=y0 .- 8k1/27 .+ 2k2 .- 3544k3/2565 .+ 1859k4/4104 .- 11k5/40
    S=runge_state_builder(S,y5)

    k6=tstep*diff_eq(S)

    return k1,k2,k3,k4,k5,k6
end

function runge45_step(S0,tstep,dt_min,force=false,eps=rk_tol)
    if tstep>dt_max
        tstep=dt_max
    end
    k1,k2,k3,k4,k5,k6 = runge_ks(S0,tstep)
    if S0.cl.NDOFs==1
        y0=[S0.cl.R[1],S0.cl.p[1],S0.el.C,S0.cl.mem]
    else
        y0=[S0.cl.R,S0.cl.p,S0.el.C,S0.cl.mem]
    end

    ynew=y0 .+ 25k1/216 .+ 1408k3/2565 .+ 2197k4/4101 .- k5/5
    znew=y0 .+ 16k1/135 .+ 6656k3/12825 .+ 28561k4/56430 .- 9k5/50 .+ 2k6/55
    R=sum([sum(abs.(znew[i]-ynew[i])) for i in 1:4])
    ds=tstep*sqrt(sqrt(rk_tol/(2*R)))
    if R<=rk_tol || force
        S0=runge_state_builder(S0,ynew)
        if S0.prefix in SH_LIST
            S0=hop!(S0,tstep)
        end
        return tstep,S0,ds
    else
        if ds<dt_min
            println("Runge-Kutta warning: convergence not achieved for minimum timestep")
            println("Continuing, be wary of results (specially if warning repeats!)...")
            @show ds
            @show dt_min
            runge45_step(S0,dt_min,dt_min,true,eps)
        else
            runge45_step(S0,ds,dt_min,false,eps)
        end

    end
end

function rk45_bigstep(S0,t0,tf,dt_ini,dt_min,eps=rk_tol)
    while t0<tf
        if tf-t0<dt_ini
            dt_ini=tf-t0
        end
        tstep,S0,dt_ini=runge45_step(S0,dt_ini,dt_min,false,eps)
        t0+=tstep
    end
    return S0,tf
end
