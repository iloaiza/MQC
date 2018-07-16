function pot_LVC_2D(R,W=[1,1],x0=1,gamma=1)
    #wx=W[1],wy=W[2]
    H=zeros(Float64,2,2)
    dHx=zeros(H)
    dHy=zeros(H)
    c=gamma*x0*W[1]^2

    H[1,1]=0.5*(W[1]^2*(R[1]+x0)^2+W[2]^2*R[2]^2)
    H[2,2]=0.5*(W[1]^2*(R[1]-x0)^2+W[2]^2*R[2]^2)
    H[1,2]=c*R[2]
    H[2,1]=H[1,2]

    dHx[1,1]=W[1]^2*(R[1]+x0)
    dHx[2,2]=W[1]^2*(R[1]-x0)

    dHy[1,1]=W[2]^2*R[2]
    dHy[2,2]=dHy[1,1]
    dHy[1,2]=c
    dHy[2,1]=c

    return H,[dHx,dHy]
end

function pot_fric_2D_model_C(R,w=0.0002,wy=0.001,g=20.6097,ya=5,DGR0=-0.0038,DGL0=4DGR0,E=1e-6,G0=2e-5,G1=0.0001,K=5,n=0)
    H=zeros(Float64,2,2)
    dHx=zeros(H)
    dHy=zeros(H)

    x=R[1]
    y=R[2]
    exp_term=exp(-K*(y-n))
    da_sqrt=sqrt((mass*wy^2*y*ya+0.5*(DGL0-DGR0))^2+E^2)

    H[1,1]=0.5*mass*w^2*(x^2+y^2)
    H[2,2]=0.5*mass*(w^2*(x-g)^2+wy^2*(y^2+ya^2))+0.5*(DGL0+DGR0)-da_sqrt
    H[1,2]=G0+(G1-G0)/(1+exp_term)
    H[2,1]=H[1,2]

    dHx[1,1]=mass*w^2*x
    dHx[2,2]=mass*w^2*(x-g)

    dHy[1,1]=mass*w^2*y
    dHy[2,2]=mass*wy^2*y-mass*wy^2*ya*(mass*wy^2*y*ya+0.5*(DGL0-DGR0))/da_sqrt
    dHy[1,2]=(G1-G0)*K*y*exp_term/(1+exp_term)^2
    dHy[2,1]=dHy[1,2]

    return H,[dHx,dHy]
end

function pot_spin_boson(R,C,D,w) #C is the coupling strength, D is the distance between parabola centers and w is the harmonic force strength
    H=zeros(Float64,2,2)
    dH=zeros(Float64,2,2)

    H[1,1]=w^2/2*(R+D/2)^2
    H[2,2]=w^2/2*(R-D/2)^2
    H[1,2]=C
    H[2,1]=C

    dH[1,1]=w^2*(R+D/2)
    dH[2,2]=w^2*(R-D/2)

    return H,[dH]
end

function pot_simple_crossing(R,A=0.01,B=1.6,C=0.0045,D=1.0)
    H=zeros(Float64,2,2)
    dH=zeros(Float64,2,2)

    H[1,1]=sign(R)*A*(1-exp(-B*abs(R)))
    H[2,2]=-H[1,1]
    H[1,2]=C*exp(-D*R^2)
    H[2,1]=H[1,2]

    dH[1,1]=A*B*exp(-B*abs(R))
    dH[2,2]=-dH[1,1]
    dH[1,2]=-C*D*2R*exp(-D*R^2)
    dH[2,1]=dH[1,2]

    return H,[dH]
end

function pot_dual_crossing(R,A=0.1,B=0.28,C=0.015,D=0.06,E0=0.05)
    H=zeros(Float64,2,2)
    dH=zeros(Float64,2,2)

    #H[1,1]=0
    H[2,2]=-A*exp(-B*R^2)+E0
    H[1,2]=C*exp(-D*R^2)
    H[2,1]=H[1,2]

    #dH[1,1]=0
    dH[2,2]=A*B*2R*exp(-B*R^2)
    dH[1,2]=-C*D*2R*exp(-D*R^2)
    dH[2,1]=dH[1,2]

    return H,[dH]
end

function pot_extended_coupling(R,A=6e-4,B=0.1,C=0.9)
    H=zeros(Float64,2,2)
    dH=zeros(Float64,2,2)

    H[1,1]=A
    H[2,2]=-A
    if R<=0
        H[1,2]=B*exp(C*R)
        dH[1,2]=C*B*exp(C*R)
    else
        H[1,2]=B*(2-exp(-C*R))
        dH[1,2]=B*C*exp(-C*R)
    end
    H[2,1]=H[1,2]

    #dH[1,1]=0
    #dH[2,2]=0
    dH[2,1]=dH[1,2]

    return H,[dH]
end

function pot_1layer(R,β=0.03,N=4,M=2000,Ω=0.0028,D=5,γ=0.02,δ=0.01)
    V=δ*ones(2N+1,2N+1)
    dV=zeros(V)

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

include("layered_potential.jl")

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

Xn = [0.0, 1.46799, 1.85059, 1.89366, 1.93795, 1.98358, 2.03068, 2.07939, 2.12987, 2.1823, 2.23689, 2.29388, 2.35354, 2.41619, 2.48221, 2.55206, 2.62627, 2.70552, 2.79061, 2.8826, 2.98281, 3.093, 3.21554, 3.35376, 3.51253, 3.69941, 3.92705, 4.21916, 4.629, 5.32647]
#values of intersections for thirty states potential function
function pot_thirty_states(R)
    H=zeros(30,30)
    dH=zeros(30,30)

    H[1,1]=0.38*(1-exp(1-R))^2
    dH[1,1]=2*0.38*(1-exp(1-R))*exp(1-R)
    H[2,2]=exp(-2R)
    dH[2,2]=-2exp(-2R)
    H[1,2]=0.05*exp(-(R-Xn[2])^2)
    H[2,1]=H[1,2]
    dH[1,2]=-0.1*(Xn[2]-R)*exp(-(R-Xn[2])^2)
    dH[2,1]=dH[1,2]
    for n in 3:30
        H[n,n]=exp(-2R)+(n+7)/100
        dH[n,n]=-2exp(-2R)
        H[1,n]=0.05*exp(-(R-Xn[n])^2)
        H[n,1]=H[1,n]
        dH[1,n]=-0.1*(Xn[n]-R)*exp(-(R-Xn[n])^2)
        dH[n,1]=dH[1,n]
    end

    return H,dH
end

NaCl_coupling(y)=0.0387./((y-8.222).^2+0.0778^2)
function pot_NaCl(R)
    x =R*0.52;
    alpha = -integral(NaCl_coupling,1,R)*0.52;
    U=zeros(2,2)
    U[1,1] = cos(alpha); U[1,2] = sin(alpha); U[2,1] = -U[1,2] ; U[2,2] = cos(alpha);
    H_ad=zeros(2,2)
    H_diab=zeros(2,2)

    De = 4.82; re = 2; be  =.8522; ce = 1.182;
    Vg = De*(1-ce*exp(-be*(x-re))).^2/27.2;
    a=-0.386 ;b =      0.0294  ;c =       86.76 ;d =       1.468  ;r =       6.01  ;
    Ve = (a*x+b*x.^2+c*exp(-d*x.^2)+ r)/27.2;
    H_ad[1,1] = Vg; H_ad[2,2] = Ve;
    H_diab = U'*H_ad*U;

    return H_diab,zeros(2)
end
