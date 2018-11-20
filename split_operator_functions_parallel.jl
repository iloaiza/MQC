function SO_pot(potential_function,R)
    H,_=potential_function(R)
    return H
end

function basis_builder(xmin,xmax,npts)
    X=linspace(xmin,xmax,npts)
    dx=X[2]-X[1]

    p=linspace(-pi,pi,npts+1)
    p=p[1:end-1]/dx
    P=ifftshift(p)

    return X,dx,P
end

function wigner(R,R0,p0)
    #return exp(-(R-R0)^2*mass*ω/2)*sqrt(sqrt(mass*ω/pi))*exp(1im*R*p0)
    exp(-(R-R0)^2/2)*sqrt(sqrt(2/pi))*exp(1im*R*p0)
end

function initial_wavefunction(R0,p0,C0,X,npts=length(X))
    Ψ=zeros(Complex{Float64},npts,nsts)
    for i in 1:npts
        Ψ[i,:].=wigner(X[i],R0,p0).*C0
    end
    return SharedArray(Ψ)
end

function matrix_builder(npts,X,P,pot)
    Vmat=SharedArray(zeros(nsts,nsts,npts))
    hatW=SharedArray(zeros(nsts,npts))
    U=SharedArray(zeros(nsts,nsts,npts))
    @sync @distributed for ipt in 1:npts
        V=pot(X[ipt])
        for i in 1:nsts
            for j in 1:nsts
                Vmat[i,j,ipt]=V[i,j]
            end
        end
        E,Ua=eig(V)
        ind=sortperm(E)
        sort!(E)
        Ua=Ua[:,ind]
        hatW[:,ipt].=E[:]
        U[:,:,ipt]=Ua
    end


    hatT=0.5/mass.*P.^2
    expT=exp.(-1im*hatT*δt/2);

    expV=SharedArray{Complex128}(zeros(Complex,nsts,nsts,npts))
    @sync @distributed for ipt in 1:npts
        expV[:,:,ipt]=U[:,:,ipt]*diagm(exp.(-1im*hatW[:,ipt]*δt))*(U[:,:,ipt]')
    end

    return Vmat,hatW,U,hatT,expT,expV
end

function ad_2_dia(Ψad,U,npts=length(Ψad[:,1]))
    Ψd=zeros(Ψad)
    for i in 1:npts
        Ψd[i,:]=U[:,:,i]*Ψad[i,:]
    end
    return SharedArray(Ψd)
end

function dia_2_ad(Ψd,U,npts=length(Ψd[:,1]))
    Ψad=zeros(Ψd)
    for i in 1:npts
        Ψad[i,:]=U[:,:,i]'*Ψd[i,:]
    end
    return SharedArray(Ψad)
end

function next_step!(Ψd,expT,expV,npts)
    @sync @distributed for ist in 1:nsts
        Ψd[:,ist]=ifft(expT.*fft(Ψd[:,ist]))
    end

    @sync @distributed for ipt in 1:npts
        Ψd[ipt,:]=expV[:,:,ipt]*Ψd[ipt,:]
    end

    @sync @distributed for ist in 1:nsts
        Ψd[:,ist]=ifft(expT.*fft(Ψd[:,ist]))
    end
end

function Prob_X(Ψ,npts=length(Ψ[:,1]))
    PX=[sum(abs2.(Ψ[i,:])) for i in 1:npts]
    return PX
end

function SO_propagation(tf,pow,xmin,xmax,R0,p0,C0,potential,flags)
    t00=time()
    npts=2^pow
    X,dx,P=basis_builder(xmin,xmax,npts)
    Ψ0=initial_wavefunction(R0,p0,C0,X,npts);
    PSI_AD=[Ψ0]
    sizehint!(PSI_AD,flags+1)

    Vmat,hatW,U,hatT,expT,expV=matrix_builder(npts,X,P,potential)

    Ψd=ad_2_dia(Ψ0,U)

    T=collect(0:δt:tf)
    nsteps=length(T)
    Tf=T[1:Int(round(nsteps/flags)):nsteps]
    if Tf[end]!=tf
        push!(Tf,tf)
    end

    t0=time()
    println("Finished setup during $(round(t0-t00,4))s, starting wavefunction propagation")

    told=copy(t0)
    counter=2
    for tstep in 1:nsteps
        if Tf[counter]==T[tstep]
            tnew=time()
            push!(PSI_AD,dia_2_ad(Ψd,U))
            println("Currently at $(round(tstep/nsteps,4)*100)%")
            println("$(round(tnew-told,4))s elapsed since last % show")
            told=time()
            counter+=1
        end
        next_step!(Ψd,expT,expV,npts)
    end
    push!(PSI_AD,dia_2_ad(Ψd,U))
    tf=time()
    println("Finished dynamics. Total propagation time was $(round(tf-t0,4)/3600) hours")

    return Tf,X,P,PSI_AD
end

function SO_save(file,Tf,X,P,PSI_AD)
    h5write(file,"Tf",Tf)
    h5write(file,"X",collect(X))
    h5write(file,"P",collect(P))
    for i in 1:length(PSI_AD)
        h5write(file,"PSI_AD_REAL_$i",real.(PSI_AD[i]))
        h5write(file,"PSI_AD_IMAG_$i",imag.(PSI_AD[i]))
    end
end

function SO_read(file)
    Tf=h5read(file,"Tf")
    X=h5read(file,"X")
    P=h5read(file,"P")
    PSI_AD=[h5read(file,"PSI_AD_REAL_$i")+1im*h5read(file,"PSI_AD_IMAG_$i") for i in 1:flags+1]
    PX=Prob_X.(PSI_AD)

    return Tf,X,P,PSI_AD,PX
end

function SO_histo_builder(X,PX,res,xmin=X[1],xmax=X[end])
    Y=range(xmin,stop=xmax,length=res)
    H=zeros(size(Y))
    dY=Y[2]-Y[1]
    dX=X[2]-X[1]
    size_factor=dX/dY

    count=1
    while Y[1]>X[count]
        count+=1
    end

    for (i,y) in enumerate(Y)
        while y>=X[count]
            H[i]+=PX[count]
            count+=1
        end
    end

    return Y,H.*size_factor
end
