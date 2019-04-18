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

function runge_ks(S,tstep)
    if S.cl.NDOFs==1
        y0=[S.cl.R[1],S.cl.p[1],S.el.C,S.cl.mem]
    else
        y0=[S.cl.R,S.cl.p,S.el.C,S.cl.mem]
    end

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

function runge45_step(S0,tstep,force=false)
    if fixed_step
        k1,_,k3,k4,k5,k6 = runge_ks(S0,tstep)
        if S0.cl.NDOFs==1
            y0=[S0.cl.R[1],S0.cl.p[1],S0.el.C,S0.cl.mem]
        else
            y0=[S0.cl.R,S0.cl.p,S0.el.C,S0.cl.mem]
        end

        ynew=y0 .+ 16k1/135 .+ 6656k3/12825 .+ 28561k4/56430 .- 9k5/50 .+ 2k6/55
        S0=runge_state_builder(S0,ynew)
        if S0.prefix in SH_LIST
            S0=hop!(S0,tstep)
        end
        return tstep,S0,tstep
    else
        if tstep>dt_max
            tstep=dt_max
        end
        k1,_,k3,k4,k5,k6 = runge_ks(S0,tstep)
        if S0.cl.NDOFs==1
            y0=[S0.cl.R[1],S0.cl.p[1],S0.el.C,S0.cl.mem]
        else
            y0=[S0.cl.R,S0.cl.p,S0.el.C,S0.cl.mem]
        end

        ynew=y0 .+ 25k1/216 .+ 1408k3/2565 .+ 2197k4/4101 .- k5/5
        znew=y0 .+ 16k1/135 .+ 6656k3/12825 .+ 28561k4/56430 .- 9k5/50 .+ 2k6/55
        Rnorm2=sum(abs2.(ynew[1]-znew[1]))
        pnorm2=sum(abs2.(ynew[2]-znew[2]))
        Cnorm2=sum(abs2.(ynew[3]))+sum(abs2.(znew[3]))-2*sum(real.(ynew[3].*conj.(znew[3])))
        Cnorm2=maximum([Cnorm2,0])
        memnorm2=sum(abs2.(ynew[4]-znew[4]))
        R_err=sqrt(Rnorm2)/(RelTol[1]*sqrt(sum(abs2.(y0[1])))+AbsTol[1])
        p_err=sqrt(pnorm2)/(RelTol[2]*sqrt(sum(abs2.(y0[2])))+AbsTol[2])
        C_err=sqrt(Cnorm2)/(RelTol[3]*sqrt(sum(abs2.(y0[3])))+AbsTol[3])
        mem_err=sqrt(memnorm2)/(RelTol[4]*sqrt(sum(abs2.(y0[4])))+AbsTol[4])
        R=sqrt(abs2(R_err)+abs2(p_err)+abs2(C_err)+abs2(mem_err))
        ds=tstep*0.84*sqrt(sqrt(1/R))
        #@show R_err,p_err,C_err,mem_err,R,ds
        if R<=1 || force
            S0=runge_state_builder(S0,ynew)
            if S0.prefix in SH_LIST
                S0=hop!(S0,tstep)
            end
            return tstep,S0,ds
        else
            if ds<dt_min #dt_min defined in code_config
                println("Runge-Kutta warning: convergence not achieved for minimum timestep")
                println("Continuing, be wary of results (specially if warning repeats!)...")
                @show ds
                @show dt_min
                runge45_step(S0,dt_min,true)
            else
                runge45_step(S0,ds,false)
            end

        end
    end
end

function rk45_bigstep(S0,t0,tf,dt_ini)
    tstep=dt_ini
    dti=dt_ini
    if time_print
        T=Float64[]
    end
    while t0<tf
        tstep,S0,dti=runge45_step(S0,dti,false)
        if time_print
            push!(T,tstep)
        end
        t0+=tstep
    end

    if time_print
        @show (maximum(T),minimum(T),mean(T))
    end
    return S0,t0,tstep
end
