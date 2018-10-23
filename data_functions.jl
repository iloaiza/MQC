#=
FUNCTIONS FOR DATA SAVING, READING AND DATA ANALYSIS (TULLY-LIKE PLOTS FOR FINAL POPULATIONS)
=#

function CL_save(Tf,Rvec,pvec,file,prefix)
    h5write(file,"T_"*prefix,Tf)
    h5write(file,"R_"*prefix,Rvec)
    h5write(file,"P_"*prefix,pvec)
end


function MF_save(Tf,Rvec,pvec,C,file,prefix)
    Creal=Float64.(real.(C))
    Cim=Float64.(imag.(C))
    h5write(file,"T_"*prefix,Tf)
    h5write(file,"R_"*prefix,Rvec)
    h5write(file,"P_"*prefix,pvec)
    h5write(file,"C_REAL_"*prefix,Creal)
    h5write(file,"C_IMAG_"*prefix,Cim)
end

function SH_save(Tf,Rvec,pvec,C,AST,file,prefix)
    Creal=Float64.(real.(C))
    Cim=Float64.(imag.(C))
    h5write(file,"T_"*prefix,Tf)
    h5write(file,"R_"*prefix,Rvec)
    h5write(file,"P_"*prefix,pvec)
    h5write(file,"C_REAL_"*prefix,Creal)
    h5write(file,"C_IMAG_"*prefix,Cim)
    h5write(file,"AST_"*prefix,AST)
end

function CL_read(file,prefix)
    T=h5read(file,"T_"*prefix)
    R=h5read(file,"R_"*prefix)
    P=h5read(file,"P_"*prefix)
    return T,R,P
end

function MF_read(file,prefix)
    T=h5read(file,"T_"*prefix)
    R=h5read(file,"R_"*prefix)
    P=h5read(file,"P_"*prefix)
    C_REAL=h5read(file,"C_REAL_"*prefix)
    C_IMAG=h5read(file,"C_IMAG_"*prefix)
    C=C_REAL+1im.*C_IMAG
    return T,R,P,C
end

function SH_read(file,prefix)
    T=h5read(file,"T_"*prefix)
    R=h5read(file,"R_"*prefix)
    P=h5read(file,"P_"*prefix)
    C_REAL=h5read(file,"C_REAL_"*prefix)
    C_IMAG=h5read(file,"C_IMAG_"*prefix)
    AST=h5read(file,"AST_"*prefix)
    C=C_REAL+1im.*C_IMAG
    return T,R,P,C,AST
end

function CL_K_save(A,file,prefix)
    h5write(file,"FR_"*prefix,A[1])
    h5write(file,"FP_"*prefix,A[2])
end

function MF_K_save(A,file,prefix)
    h5write(file,"FR_"*prefix,A[1])
    h5write(file,"FM_"*prefix,A[2])
    h5write(file,"Fpop_"*prefix,A[3])
end

function SH_K_save(A,file,prefix)
    h5write(file,"FR_"*prefix,A[1])
    h5write(file,"FM_"*prefix,A[2])
    h5write(file,"Fpop_"*prefix,A[3])
    h5write(file,"Fast_"*prefix,A[4])
end

function general_K_save(A,file,prefix)
    if length(A)==2
        CL_K_save(A,file,prefix)
    elseif length(A)==3
        MF_K_save(A,file,prefix)
    elseif length(A)==4
        SH_K_save(A,file,prefix)
    else
        error("Check our input, it should have 2, 3 or 4 entries!")
    end
end

function general_K_read(file,prefix)
    A=Array[]
    if prefix in METHOD_LIST
        FR=h5read(file,"FR_"*prefix)
        FM=h5read(file,"FM_"*prefix)
        push!(A,FR)
        push!(A,FM)
        if prefix in MF_LIST || prefix in SH_LIST
            Fpop=h5read(file,"Fpop_"*prefix)
            push!(A,Fpop)
            if prefix in SH_LIST
                Fast=h5read(file,"Fast_"*prefix)
                push!(A,Fast)
            end
        end
    end

    return A
end

function full_memory_save(T,S_ARRAY,file)
    prefix=S_ARRAY[1].prefix
    steps=length(S_ARRAY)
    NDOFs=S_ARRAY[1].cl.NDOFs
    if NDOFs==1
        R=zeros(steps)
        p=zeros(steps)
    else
        R=zeros(steps,NDOFs)
        p=zeros(steps,NDOFs)
    end

    if !(prefix in CL_LIST) #not classical, need electronic coefficients tracking
        eff_sts=length(S_ARRAY[1].el.C)
        C=zeros(Complex,steps,eff_sts)
        U=zeros(steps,nsts,nsts)
    end

    for k in 1:steps
        S=S_ARRAY[k]
        if NDOFs==1
            R[k]=S.cl.R
            p[k]=S.cl.p
        else
            R[k,:].=S.cl.R
            p[k,:].=S.cl.p
        end
        if !(prefix in CL_LIST) #not classical, need electronic coefficients tracking
            C[k,:].=S.el.C
            U[k,:,:].=S.el.Ua
        end
    end
    Creal=zeros(size(C))
    Cimag=zeros(size(C))
    Creal.=real.(C)
    Cimag.=imag.(C)

    h5write(file,"R_"*prefix,R)
    h5write(file,"p_"*prefix,p)
    if !(prefix in CL_LIST)
        h5write(file,"Creal_"*prefix,Creal)
        h5write(file,"Cimag_"*prefix,Cimag)
        h5write(file,"U_"*prefix,U)
    end
    if prefix in SH_LIST
        AST=[S_ARRAY[k].ast for k in 1:steps]
        h5write(file,"AST_"*prefix,AST)
    end
    h5write(file,"T_"*prefix,T)
end

function full_memory_load(prefix,file,flags=100)
    R=h5read(file,"R_"*prefix)
    p=h5read(file,"p_"*prefix)
    T=h5read(file,"T_"*prefix)
    steps=length(T)
    NDOFs=length(R[1,:])

    kmap=Int.(round.(collect(range(1,stop=steps,length=flags))))

    if !(prefix in CL_LIST)
        C=h5read(file,"Creal_"*prefix)+1im.*h5read(file,"Cimag_"*prefix)
        U=h5read(file,"U_"*prefix)
        if prefix in SH_LIST
            AST=h5read(file,"AST_"*prefix)
        end
    end

    if NDOFs==1
        if prefix in CL_LIST
            S_ARRAY=[builder_CL_state(R[k],p[k],prefix,NDOFs) for k in kmap]
        elseif prefix in MF_LIST
            S_ARRAY=[builder_MF_state(R[k],p[k],C[k,:],prefix,U[k,:,:],NDOFs) for k in kmap]
        else
            S_ARRAY=[builder_SH_state(R[k],p[k],C[k,:],AST[k],prefix,U[k,:,:],NDOFs) for k in kmap]
        end
    else
        if prefix in CL_LIST
            S_ARRAY=[builder_CL_state(R[k,:],p[k,:],prefix,NDOFs) for k in kmap]
        elseif prefix in MF_LIST
            S_ARRAY=[builder_MF_state(R[k,:],p[k,:],C[k,:],prefix,U[k,:,:],NDOFs) for k in kmap]
        else
            S_ARRAY=[builder_SH_state(R[k,:],p[k,:],C[k,:],AST[k],prefix,U[k,:,:],NDOFs) for k in kmap]
        end
    end

    return T[kmap],S_ARRAY
end
