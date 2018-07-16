using Distributions
#srand(123)  ####uncomment for reproducible set of random numbers

function constant_dist(A)
    return A
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

function wigner_std(A) #A=[xm,pm,sigma_x,sigma_y] contains the parameters for Wigner distribution with variable standard deviation. returns a position and a momentum under the Wigner distribution
    xm=A[1] #mean position
    pm=A[2] #mean momentum
    σx=A[3] #x standard deviation
    σp=A[4] #p standard deviation
    X=zeros(Float64,size(xm))
    P=zeros(Float64,size(pm))
    for i in 1:length(xm)
        X[i]=rand(Normal(xm[i],σx))
        P[i]=rand(Normal(pm[i],σp))
    end
    return [X,P]
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

function abs_2_2D(v) #turns a velocity magnitude into a random 2D vector
    theta=2pi*rand()
    vx=v*cos(theta)
    vy=v*sin(theta)
    return [vx,vy]
end

function initial_1D_distribution(A,dist_method,N) #A has the important parameters for dist_method. N is the number of initial parameters to calculate
    I1=dist_method(A)
    out_size=length(I1)
    I=zeros(out_size,N)
    I[:,1].=I1
    for i in 2:N
        I[:,i].=dist_method(A)
    end

    return I
end
