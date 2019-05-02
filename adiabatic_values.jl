function adiabatic_values(R,NDOFs=length(R))
    Hd,dHd=potential(R) #dHd must be an array with NDOFs nsts*nsts arrays
    E,Ua=eigen(Hd)
    ind=sortperm(E)
    sort!(E)
    Ua=Ua[:,ind]
    W=Diagonal(E)*ones(nsts,nsts)-ones(nsts,nsts)*Diagonal(E)
    F=[zeros(nsts,nsts) for k in 1:NDOFs]
    Γ=[zeros(nsts,nsts) for k in 1:NDOFs]
    for k in 1:NDOFs
        F[k]=(Ua')*dHd[k]*Ua
        Γ[k]=-F[k]./(W+Diagonal(ones(nsts)))
        Γ[k]=Γ[k]-Diagonal(Γ[k])
    end

    return E,Γ,F,W,Ua,dHd
end

function diabatic_values(R)
    Hd,dHd=potential(R)

    return Hd,dHd
end

function CM2_additional_values(p,Γ,W,NDOFs)
    NACs=sum(p/mass.*Γ)
    tvec=-NACs[2:end,1]

    z=norm(tvec)*sign(sum(tvec))
    if z==0
        tnorm=zeros(size(tvec))
    else
        tnorm=tvec./abs(z)
    end
    wvec=W[2:end,1]

    return z,tnorm,wvec
end

function CM3_additional_values(wvec,tnorm,z,NDOFs)
    if z==0
        kvec=zeros(size(tnorm))
        K=1
    else
        kvec=tnorm./tnorm[1]
        K=1/tnorm[1]
    end
    w_2=wvec[2:end]
    k_2=kvec[2:end]

    wvec2=w_2.-(k_2.^2/K^2).*(w_2.-wvec[1])
    tvec2=k_2.*(w_2.-wvec[1])./K^2

    tmax=maximum(tvec2)
    tmin=minimum(tvec2)
    if abs(tmin)>abs(tmax)
        zbarsign=sign(tmin)
    else
        zbarsign=sign(tmax)
    end
    zbar=norm(tvec2)*zbarsign

    if zbar==0
        tnorm2=zeros(size(tvec2))
    else
        tnorm2=tvec2./zbar  #carries zbar sign! needed for equation consistency
    end

    return zbar,wvec2,tnorm2
end


function collective_force(F,tnorm,z,NDOFs,sts=length(tnorm)+1)
    #this function does the full transformation of the symmetric matrix F into collective space under the unitary U assumption
    F2=[zeros(sts,sts) for k in 1:NDOFs]
    for k in 1:NDOFs
        F2[k][1,1]=F[k][1,1]                                                                                                                    #f00
        for j in 3:sts
            F2[k][1,j]=F[k][1,2]*tnorm[j-1]-F[k][1,j]*tnorm[1]                                                                                  #f0j
            F2[k][j,1]=F2[k][1,j]
            for m in 2:sts
                F2[k][2,j]+=tnorm[j-1]*F[k][2,m]*tnorm[m-1]-tnorm[1]*F[k][m,j]*tnorm[m-1]                                                      #f1j
            end                                                       #
            for m in 3:sts
                F2[k][j,m]=tnorm[m-1]*(F[k][2,2]*tnorm[j-1]-F[k][2,j]*tnorm[1])-tnorm[1]*(F[k][2,m]*tnorm[j-1]-F[k][j,m]*tnorm[1])             #fjm
            end
        end
        for j in 2:sts
            F2[k][j,2]=F2[k][2,j]
            F2[k][1,2]+=F[k][1,j]*tnorm[j-1]                                                                                                   #f01
            for m in 2:sts
                F2[k][2,2]+=F[k][j,m]*tnorm[j-1]*tnorm[m-1]                                                                                    #f11
            end
        end
        F2[k][2,1]=F2[k][1,2]
    end

    return F2
end

function reduced_force(F,tnorm,z,NDOFs,sts=length(tnorm)+1)
    #this function only calculates the 2x2 upper left block of the collective variable transformed force
    F2=[zeros(2,2) for k in 1:NDOFs]

    for k in 1:NDOFs
        F2[k][1,1]=F[k][1,1]
        for j in 2:sts
            F2[k][1,2]+=F[k][1,j]*tnorm[j-1]
            for m in 2:sts
                F2[k][2,2]+=F[k][j,m]*tnorm[j-1]*tnorm[m-1]
            end
        end
        F2[k][2,1]=F2[k][1,2]
    end

    return F2
end

function diagonal_force(F,tnorm,z,NDOFs,sts=length(tnorm)+1)
    #this function only calculates the first 2 diagonal elements of the collective force matrix
    F2=[zeros(2) for k in 1:NDOFs]

    for k in 1:NDOFs
        F2[k][1]=F[k][1,1]
        for j in 2:sts
            for m in 2:sts
                F2[k][2]+=F[k][j,m]*tnorm[j-1]*tnorm[m-1]
            end
        end
    end

    return F2
end

function gs_remove(F,NDOFs=length(F),sts=length(F[1][:,1]))
    F2=[zeros(sts-1,sts-1) for k in 1:NDOFs]

    for k in 1:NDOFs
        F2[k][:,:].=F[k][2:end,2:end]
    end

    return F2
end
