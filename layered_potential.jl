const Omega_layered=0.004/sqrt(2)
const thickness_layered=9

dist=0.16/2
offset_delta=0.0001
offset=offset_delta.*kron(ones(9,1),0:thickness_layered)

alpha=0.004


pos=dist.*(-4:1:4)
D=kron(pos,ones(thickness_layered+1,1))+offset
ALPHA_MAT=zeros(9*(thickness_layered+1),9*(thickness_layered+1))
for i in 2:9*(thickness_layered+1)
    ALPHA_MAT[i,i-1]=alpha
    ALPHA_MAT[i-1,i]=alpha
end
for i in 1:thickness_layered-1
    ALPHA_MAT[10i+1,10i]=0
    ALPHA_MAT[10i,10i+1]=0
end
const Dreal_layered=D./(mass*Omega_layered^2)
const A_layered=ALPHA_MAT

function Vd(R,D)
    return 0.5*mass*Omega_layered^2*(R-D)^2
end

function Vtot(R,beta=0.03,gamma=0.02,delta=0.01)
    Bmat=beta.*eye(thickness_layered+1)
    Cmat=gamma.*eye(thickness_layered+1)
    Dmat=delta.*eye(thickness_layered+1)

    V=zeros(9*(thickness_layered+1),9*(thickness_layered+1))
    kron_mat=zeros(thickness_layered,thickness_layered)

    for i in 1:9*(thickness_layered+1)
        V[i,i]=Vd(R,Dreal_layered[i])
    end

    for i in 2:9
        kron_mat[i-1,i]=1
        kron_mat[i,i-1]=1
    end
    V+=kron(kron_mat,Bmat)

    kron_mat=zeros(9,9)
    for i in 3:9
        kron_mat[i-2,i]=1
        kron_mat[i,i-2]=1
    end
    V+=kron(kron_mat,Cmat)

    kron_mat=ones(9,9)
    for j in 1:3
        for i in j:9
            kron_mat[i-(j-1),i]=0
            kron_mat[i,i-(j-1)]=0
        end
    end
    V+=kron(kron_mat,Dmat)

    return V+A_layered
end

function dVd(R,D)
    return mass*Omega_layered^2*(R-D)
end

function dVtot(R)

    dV=zeros(9*(thickness_layered+1),9*(thickness_layered+1))

    for i in 1:9*(thickness_layered+1)
        dV[i,i]=dVd(R,Dreal_layered[i])
    end

    return dV
end

function pot_layered(R,beta=0.03,gamma=0.02,delta=0.01)
    return Vtot(R,beta,gamma,delta),[dVtot(R)]
end
