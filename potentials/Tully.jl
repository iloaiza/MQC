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
