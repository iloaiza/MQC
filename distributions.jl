using Distributions
#srand(123)  ####uncomment for reproducible set of random numbers

function constant_dist(A)
    return A
end

function abs_2_2D(v) #turns a velocity magnitude into a random 2D vector
    theta=2pi*rand()
    vx=v*cos(theta)
    vy=v*sin(theta)
    return vx,vy
end

function wigner(A) #A contains the parameters for Wigner distribution. returns a position and a momentum under the Wigner distribution. 1D function
    xm=A[1] #mean position
    pm=A[2] #mean momentum
    σx=1/sqrt(2)
    σp=1/sqrt(2)
    X=rand(Normal(xm,σx))
    P=rand(Normal(pm,σp))
    return [X,P]
end

function wigner_std(A,σx=1/sqrt(2),σp=1/sqrt(2)) #A=[xm,pm] contains the parameters for Wigner distribution, set σx and σp here for variable standard deviation. returns a position and a momentum under the Wigner distribution
    xm=A[1] #mean position
    pm=A[2] #mean momentum
    X=zeros(Float64,size(xm))
    P=zeros(Float64,size(pm))
    for i in 1:length(xm)
        X[i]=rand(Normal(xm[i],σx))
        P[i]=rand(Normal(pm[i],σp))
    end
    return [X,P]
end

function wigner_vel_2D(A) #A is an array with two entries, A=[Rm,Pm] where Rm=[xm,ym], Pm=[pxm,pym]
    #this will not take into account the direction of the momentum vector, just the magnitude and take a random direction
    Rm=A[1] #mean position
    Pm=sqrt(sum(abs2.(A[2]))) #mean momentum
    xm=Rm[1]
    ym=Rm[2]
    σr=1/sqrt(2)
    σp=1/sqrt(2)
    X=rand(Normal(xm,σr))
    Y=rand(Normal(ym,σr))
    P=rand(Normal(Pm,σp))
    px,py=abs_2_2D(P)

    return [[X,Y],[px,py]]
end

function wigner_dir_2D(A) #same as wigner_vel_2D, but takes gaussian distribution along each direction of the momentum vector
    Rm=A[1] #mean position
    Pm=A[2] #mean momentum
    xm=Rm[1]
    ym=Rm[2]
    pxm=Pm[1]
    pxy=Pm[2]
    σr=1/sqrt(2)
    σp=1/sqrt(2)
    X=rand(Normal(xm,σr))
    Y=rand(Normal(ym,σr))
    Px=rand(Normal(pxm,σp))
    Py=rand(Normal(pxy,σp))

    return [[X,Y],[Px,Py]]
end

function boltzmann_velocity(A) #returns a velocity under a boltzmann distribution with A=[mass,temperature]
    mass=A[1] #mass
    T=A[2] #temperature
    V=zeros(2)
    v=rand(Chi(3))
    v=v*sqrt(2*T/mass)
    return v
end

function uniform_velocity(A) #returns a velocity under a uniform distribution between A=[vmin,vmax]
    vmin=A[1]
    vmax=A[2]
    v=rand(Uniform(vmin,vmax))
    return v
end

function initial_1D_distribution(A,dist_method,N) #A has the important parameters for dist_method. N is the number of initial parameters to calculate
    if length(A[1])==1
        I1=dist_method(A)
        out_size=length(I1)
        I=zeros(out_size,N)
        I[:,1].=I1
        for i in 2:N
            I[:,i].=dist_method(A)
        end

        return I
    else #multiple NDOFs
        NDOFs=length(A[1])
        I=[dist_method(A) for n in 1:N]
        I_fin=[zeros(NDOFs) for i in 1:2, n in 1:N]
        for n in 1:N
            I_fin[1,n]=I[n][1]
            I_fin[2,n]=I[n][2]
        end

        return I_fin
    end
end
