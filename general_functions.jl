function super_zeros(n,dims,le_type="Float64")
    zero_string="Z=[zeros($le_type,$n) for i1 in 1:$(dims[1])"
    for i in 2:length(dims)
        zero_string=zero_string*", i$i in 1:$(dims[2])"
    end
    zero_string=zero_string*"]"
    eval(parse(zero_string))

    return Z
end

function NxN(BASE,ex,N=length(BASE),le_type=Float64)
    #build (res)x(res)x...x(res) (N times) matrix with N dimensional vectors in each entry
    #BASE is an array of linspaces over which the array is built
    #ex is the number of extra dimensions to add over N
    #le_type admits string entry, change to "Complex" for storing complex stuff
    #or maybe even touples ;)
    global A_NxN=super_zeros(N+ex,length.(BASE),le_type)
    global B_NxN=BASE
    #@show size(A)
    forstring="for i1 in 1:$(length(BASE[1]))"
    for i in 2:N
        forstring=forstring*", i$i in 1:$(length(BASE[i]))"
    end
    forstring=forstring*"    A_NxN[i1"
    for i in 2:N
        forstring=forstring*",i$i"
    end
    forstring=forstring*"]=[B_NxN[1][i1]"
    for i in 2:N
        forstring=forstring*",B_NxN[$i][i$i]"
    end
    for i in 1:ex
        forstring=forstring*",0"
    end
    forstring=forstring*"];  end"
    eval(parse(forstring))

    return A
end

function BO_save(Tf,Rvec,pvec,file)
    h5write(file,"T_BO",Tf)
    h5write(file,"R_BO",Rvec)
    h5write(file,"P_BO",pvec)
end


function EH_save(Tf,Rvec,pvec,C,file)
    Creal=real.(C)
    Cim=imag.(C)
    h5write(file,"T_EH",Tf)
    h5write(file,"R_EH",Rvec)
    h5write(file,"P_EH",pvec)
    h5write(file,"C_REAL_EH",Creal)
    h5write(file,"C_IMAG_EH",Cim)
end

function FSSH_save(Tf,Rvec,pvec,C,AST,file)
    Creal=real.(C)
    Cim=imag.(C)
    h5write(file,"T_FSSH",Tf)
    h5write(file,"R_FSSH",Rvec)
    h5write(file,"P_FSSH",pvec)
    h5write(file,"C_REAL_FSSH",Creal)
    h5write(file,"C_IMAG_FSSH",Cim)
    h5write(file,"AST_FSSH",AST)
end

function FSSH_dia_save(Tf,Rvec,pvec,C,DST,file)
    Creal=real.(C)
    Cim=imag.(C)
    h5write(file,"T_FSSHd",Tf)
    h5write(file,"R_FSSHd",Rvec)
    h5write(file,"P_FSSHd",pvec)
    h5write(file,"C_REAL_FSSHd",Creal)
    h5write(file,"C_IMAG_FSSHd",Cim)
    h5write(file,"DST_FSSHd",DST)
end

function CM2_VANILLA_save(Tf,Rvec,pvec,D,file)
    Dreal=real.(D)
    Dim=imag.(D)
    h5write(file,"T_CM2_VANILLA",Tf)
    h5write(file,"R_CM2_VANILLA",Rvec)
    h5write(file,"P_CM2_VANILLA",pvec)
    h5write(file,"D_REAL_CM2_VANILLA",Dreal)
    h5write(file,"D_IMAG_CM2_VANILLA",Dim)
end

function CM3_VANILLA_save(Tf,Rvec,pvec,D,file)
    Dreal=real.(D)
    Dim=imag.(D)
    h5write(file,"T_CM3_VANILLA",Tf)
    h5write(file,"R_CM3_VANILLA",Rvec)
    h5write(file,"P_CM3_VANILLA",pvec)
    h5write(file,"D_REAL_CM3_VANILLA",Dreal)
    h5write(file,"D_IMAG_CM3_VANILLA",Dim)
end

function CM2_save(Tf,Rvec,pvec,D,file)
    Dreal=real.(D)
    Dim=imag.(D)
    h5write(file,"T_CM2",Tf)
    h5write(file,"R_CM2",Rvec)
    h5write(file,"P_CM2",pvec)
    h5write(file,"D_REAL_CM2",Dreal)
    h5write(file,"D_IMAG_CM2",Dim)
end

function CM3_save(Tf,Rvec,pvec,D,file)
    Dreal=real.(D)
    Dim=imag.(D)
    h5write(file,"T_CM3",Tf)
    h5write(file,"R_CM3",Rvec)
    h5write(file,"P_CM3",pvec)
    h5write(file,"D_REAL_CM3",Dreal)
    h5write(file,"D_IMAG_CM3",Dim)
end

function CM2_FSSH_save(Tf,Rvec,pvec,D,file)
    Dreal=real.(D)
    Dim=imag.(D)
    h5write(file,"T_CM2_FSSH",Tf)
    h5write(file,"R_CM2_FSSH",Rvec)
    h5write(file,"P_CM2_FSSH",pvec)
    h5write(file,"D_REAL_CM2_FSSH",Dreal)
    h5write(file,"D_IMAG_CM2_FSSH",Dim)
end

function CM3_FSSH_save(Tf,Rvec,pvec,D,file)
    Dreal=real.(D)
    Dim=imag.(D)
    h5write(file,"T_CM3_FSSH",Tf)
    h5write(file,"R_CM3_FSSH",Rvec)
    h5write(file,"P_CM3_FSSH",pvec)
    h5write(file,"D_REAL_CM3_FSSH",Dreal)
    h5write(file,"D_IMAG_CM3_FSSH",Dim)
end

function CM2_FSSH_VANILLA_save(Tf,Rvec,pvec,D,file)
    Dreal=real.(D)
    Dim=imag.(D)
    h5write(file,"T_CM2_FSSH_VANILLA",Tf)
    h5write(file,"R_CM2_FSSH_VANILLA",Rvec)
    h5write(file,"P_CM2_FSSH_VANILLA",pvec)
    h5write(file,"D_REAL_CM2_FSSH_VANILLA",Dreal)
    h5write(file,"D_IMAG_CM2_FSSH_VANILLA",Dim)
end

function CM3_FSSH_VANILLA_save(Tf,Rvec,pvec,D,file)
    Dreal=real.(D)
    Dim=imag.(D)
    h5write(file,"T_CM3_FSSH_VANILLA",Tf)
    h5write(file,"R_CM3_FSSH_VANILLA",Rvec)
    h5write(file,"P_CM3_FSSH_VANILLA",pvec)
    h5write(file,"D_REAL_CM3_FSSH_VANILLA",Dreal)
    h5write(file,"D_IMAG_CM3_FSSH_VANILLA",Dim)
end

function BO_read(file)
    T=h5read(file,"T_BO")
    R=h5read(file,"R_BO")
    P=h5read(file,"P_BO")
    return T,R,P
end

function EH_read(file)
    T=h5read(file,"T_EH")
    R=h5read(file,"R_EH")
    P=h5read(file,"P_EH")
    C_REAL=h5read(file,"C_REAL_EH")
    C_IMAG=h5read(file,"C_IMAG_EH")
    C=C_REAL+1im.*C_IMAG
    return T,R,P,C
end

function FSSH_read(file)
    T=h5read(file,"T_FSSH")
    R=h5read(file,"R_FSSH")
    P=h5read(file,"P_FSSH")
    C_REAL=h5read(file,"C_REAL_FSSH")
    C_IMAG=h5read(file,"C_IMAG_FSSH")
    AST=h5read(file,"AST_FSSH")
    C=C_REAL+1im.*C_IMAG
    return T,R,P,C,AST
end

function FSSH_dia_read(file)
    T=h5read(file,"T_FSSHd")
    R=h5read(file,"R_FSSHd")
    P=h5read(file,"P_FSSHd")
    C_REAL=h5read(file,"C_REAL_FSSHd")
    C_IMAG=h5read(file,"C_IMAG_FSSHd")
    DST=h5read(file,"DST_FSSHd")
    C=C_REAL+1im.*C_IMAG
    return T,R,P,C,DST
end

function CM2_VANILLA_read(file)
    T=h5read(file,"T_CM2_VANILLA")
    R=h5read(file,"R_CM2_VANILLA")
    P=h5read(file,"P_CM2_VANILLA")
    D_REAL=h5read(file,"D_REAL_CM2_VANILLA")
    D_IMAG=h5read(file,"D_IMAG_CM2_VANILLA")
    D=D_REAL+1im.*D_IMAG
    return T,R,P,D
end

function CM3_VANILLA_read(file)
    T=h5read(file,"T_CM3_VANILLA")
    R=h5read(file,"R_CM3_VANILLA")
    P=h5read(file,"P_CM3_VANILLA")
    D_REAL=h5read(file,"D_REAL_CM3_VANILLA")
    D_IMAG=h5read(file,"D_IMAG_CM3_VANILLA")
    D=D_REAL+1im.*D_IMAG
    return T,R,P,D
end

function CM2_read(file)
    T=h5read(file,"T_CM2")
    R=h5read(file,"R_CM2")
    P=h5read(file,"P_CM2")
    D_REAL=h5read(file,"D_REAL_CM2")
    D_IMAG=h5read(file,"D_IMAG_CM2")
    D=D_REAL+1im.*D_IMAG
    return T,R,P,D
end

function CM3_read(file)
    T=h5read(file,"T_CM3")
    R=h5read(file,"R_CM3")
    P=h5read(file,"P_CM3")
    D_REAL=h5read(file,"D_REAL_CM3")
    D_IMAG=h5read(file,"D_IMAG_CM3")
    D=D_REAL+1im.*D_IMAG
    return T,R,P,D
end

function CM2_FSSH_read(file)
    T=h5read(file,"T_CM2_FSSH")
    R=h5read(file,"R_CM2_FSSH")
    P=h5read(file,"P_CM2_FSSH")
    D_REAL=h5read(file,"D_REAL_CM2_FSSH")
    D_IMAG=h5read(file,"D_IMAG_CM2_FSSH")
    D=D_REAL+1im.*D_IMAG
    return T,R,P,D
end

function CM3_FSSH_read(file)
    T=h5read(file,"T_CM3_FSSH")
    R=h5read(file,"R_CM3_FSSH")
    P=h5read(file,"P_CM3_FSSH")
    D_REAL=h5read(file,"D_REAL_CM3_FSSH")
    D_IMAG=h5read(file,"D_IMAG_CM3_FSSH")
    D=D_REAL+1im.*D_IMAG
    return T,R,P,D
end

function CM2_FSSH_VANILLA_read(file)
    T=h5read(file,"T_CM2_FSSH_VANILLA")
    R=h5read(file,"R_CM2_FSSH_VANILLA")
    P=h5read(file,"P_CM2_FSSH_VANILLA")
    D_REAL=h5read(file,"D_REAL_CM2_FSSH_VANILLA")
    D_IMAG=h5read(file,"D_IMAG_CM2_FSSH_VANILLA")
    D=D_REAL+1im.*D_IMAG
    return T,R,P,D
end

function CM3_FSSH_VANILLA_read(file)
    T=h5read(file,"T_CM3_FSSH_VANILLA")
    R=h5read(file,"R_CM3_FSSH_VANILLA")
    P=h5read(file,"P_CM3_FSSH_VANILLA")
    D_REAL=h5read(file,"D_REAL_CM3_FSSH_VANILLA")
    D_IMAG=h5read(file,"D_IMAG_CM3_FSSH_VANILLA")
    D=D_REAL+1im.*D_IMAG
    return T,R,P,D
end

function general_K_save(A,file,suffix)
    #A array is made of arrays to save. A[1] is final positions, A[2] final momentum
    #A[3] are electronic populations and A[4] the effective states for SH
    h5write(file,"FR_"*suffix,A[1])
    h5write(file,"FM_"*suffix,A[2])
    if length(A)>=3
        h5write(file,"Fpop_"*suffix,A[3])
        if length(A)>=4
            h5write(file,"Fast_"*suffix,A[4])
        end
    end
end

function general_K_read(file,suffix,n)
    A=Array[]
    if n>=2
        FR=h5read(file,"FR_"*suffix)
        FM=h5read(file,"FM_"*suffix)
        push!(A,FR)
        push!(A,FM)
        if n>=3
            Fpop=h5read(file,"Fpop_"*suffix)
            push!(A,Fpop)
            if n==4
                Fast=h5read(file,"Fast_"*suffix)
                push!(A,Fast)
            end
        end
    end

    return A
end


function my_histogram(X,n)
    Xrange=linspace(minimum(X),maximum(X),n+1)
    plotting_X=linspace(Xrange[1]+step/2,Xrange[end]-step/2,n)
    H=zeros(n)

    for i in 1:n-1
        H[i]=count(x->(Xrange[i]<=x<Xrange[i+1]),X)
    end
    H[n]=count(x->(Xrange[n]<=x<=Xrange[n+1]),X)

    return plotting_X,H
end

function my_histogram(X,xmin,xmax,n)
    Xrange=linspace(xmin,xmax,n+1)
    step=Xrange[2]-Xrange[1]
    plotting_X=linspace(Xrange[1]+step/2,Xrange[end]-step/2,n)
    H=zeros(n)

    for i in 1:n-1
        H[i]=count(x->(Xrange[i]<=x<Xrange[i+1]),X)
    end
    H[n]=count(x->(Xrange[n]<=x<=Xrange[n+1]),X)

    return plotting_X,H
end

function super_histo(R,n)
    xmin=minimum(R)
    xmax=maximum(R)
    flags=length(R[:,1])
    plotting_R=zeros(n,flags)
    histo_R=zeros(n,flags)
    for i in 1:flags
        plotting_R[:,i],histo_R[:,i]=my_histogram(R[i,:],xmin,xmax,n)
    end

    return plotting_R,histo_R
end

function super_histo(R,xmin,xmax,n)
    flags=length(R[:,1])
    plotting_R=zeros(n,flags)
    histo_R=zeros(n,flags)
    for i in 1:flags
        plotting_R[:,i],histo_R[:,i]=my_histogram(R[i,:],xmin,xmax,n)
    end

    return plotting_R,histo_R
end

function integral(f,a,b,dx=1e-5)
    X=collect(a:dx:b)
    SUM=0
    for x in X
        SUM+=f(x)
    end
    SUM=SUM*dx

    return SUM
end
