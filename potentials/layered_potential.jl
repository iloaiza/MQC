function pot_1layer(R,β=0.03,N=4,M=2000,Ω=0.0028,D=5,γ=0.02,δ=0.01)
    V=δ*ones(2N+1,2N+1)
    dV=zeros(size(V))

    icount=0
    for i in -N:N
        icount+=1
        V[icount,icount]=M*Ω^2*(R+i*D)^2/2
        dV[icount,icount]=M*Ω^2*(R+i*D)
    end

    for i in 1:2N
        V[i,i+1]=β
        V[i+1,i]=β
    end

    for i in 1:2N-1
        V[i,i+2]=γ
        V[i+2,i]=γ
    end

    return V,[dV]
end


function pot_sinus(R,w,h,c)#width, height and coupling
    H=zeros(2,2)
    dH=zeros(2,2)

    H[1,1]=h*sin(w*R)
    H[2,2]=-h*sin(w*R)
    H[1,2]=c
    H[2,1]=c

    dH[1,1]=h*w*cos(w*R)
    dH[2,2]=-h*w*cos(w*R)
    dH[1,2]=0
    dH[2,1]=0

    return H,[dH]
end

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


function pot_2D_sinus(R,c,l=10,d=1e-4,h=0.02,w=0.5,RL=[1e-2,1e-2])#w=width, h=height, c=normal coupling, d=intra_layer coupling, l=layers, RL=layer-layer dist
    H=zeros(2*l,2*l)
    dHx=zeros(2*l,2*l)
    dHy=zeros(2*l,2*l)

    for layer in 1:l
        R_layer=R-(layer-1)*RL
        H[layer,layer]=h*sin(w*R_layer[1])*sin(w*R_layer[2])
        H[l+layer,l+layer]=-h*sin(w*R_layer[1])*sin(w*R_layer[2])
        H[layer,l+layer]=c
        H[l+layer,layer]=c

        dHx[layer,layer]=h*w*cos(w*R_layer[1])*sin(w*R_layer[2])
        dHy[layer,layer]=h*w*cos(w*R_layer[2])*sin(w*R_layer[1])
        dHx[layer+l,layer+l]=-h*w*cos(w*R_layer[1])*sin(w*R_layer[2])
        dHy[layer+l,layer+l]=-h*w*cos(w*R_layer[2])*sin(w*R_layer[1])

        if layer<l
            H[layer,layer+1]=d
            H[layer+1,layer]=d
            H[l+layer,l+layer+1]=d
            H[l+layer+1,l+layer]=d
        end
    end


    return H,[dHx,dHy]
end
