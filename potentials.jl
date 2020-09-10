#=
POTENTIALS, SOME GENERAL ONES PLUS ADDS ALL INSIDE POTENTIALS FOLDER
REMEMBER THE POTENTIAL RETURN MUST BE H,[dH], WITH H THE DIABATIC HAMILTONIAN AND [dH] AN ARRAY WITH THE PARTIAL DERIVATIVES
[dH] SIZE WILL THEN BE GIVEN BY THE NUCLEAR DEGREES OF FREEDOM (NDOFs)
=#

function pot_LVC_2D(R,W=[1,1],x0=1,gamma=1)
    #wx=W[1],wy=W[2]
    H=zeros(Float64,2,2)
    dHx=zeros(2,2)
    dHy=zeros(2,2)
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

function pot_LZ(R,d=1,c=0.5)
    H=[d*R c;c -d*R]
    dH=[d 0;0 -d]

    return H,[dH]
end


for jl_POT in readdir("./potentials/")
    include("./potentials/"*jl_POT)
end
