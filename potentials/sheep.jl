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
    dH[1,2]=0.1*(Xn[2]-R)*exp(-(R-Xn[2])^2)
    dH[2,1]=dH[1,2]
    for n in 3:30
        H[n,n]=exp(-2R)+(n+7)/100
        dH[n,n]=-2exp(-2R)
        H[1,n]=0.05*exp(-(R-Xn[n])^2)
        H[n,1]=H[1,n]
        dH[1,n]=0.1*(Xn[n]-R)*exp(-(R-Xn[n])^2)
        dH[n,1]=dH[1,n]
    end

    return H,[dH]
end

function pot_multiple_crossings(R)
    H=zeros(3,3)
    dH=zeros(3,3)

    H[2,2]=0.045
    H[3,3]=-0.1*exp(-0.28*R^2)+0.05
    H[1,3]=0.015*exp(-0.06*R^2)
    H[3,1]=H[1,3]
    H[2,3]=0.035*exp(-0.06*R^2)
    H[3,2]=H[2,3]

    dH[3,3]=0.056*R*exp(-0.28*R^2)
    dH[1,3]=-0.0018*R*exp(-0.06*R^2)
    dH[3,1]=dH[1,3]
    dH[2,3]=-0.0042*R*exp(-0.06*R^2)
    dH[3,2]=dH[2,3]

    return H,[dH]
end
