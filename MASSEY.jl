const mass=20000
include("include.jl")

function LZ_estimator(R,R0,p0,P,chi=1,inist=1,finst=2)
    C0=zeros(nsts)
    C0[inist]=1
    E0,Γ0,F0,W0,Ua0,dHd0=adiabatic_values(R0)
    E,Γ,F,W,Ua,dHd=adiabatic_values(R)
    Hd,dHd=diabatic_values(R)
    dHd=dHd[1]
    DF=dHd[finst,finst]-dHd[inist,inist]
    TE0=real(C0'*E0*C0)[inist]+p0^2/2/mass
    w21=E[finst]-E[inist]
    Rdot=sqrt(2/mass*(TE0+P*w21-E[finst]))
    c=sqrt(abs(DF)*Rdot*log(1/(1-P))/2/pi)
    cexp=Hd[inist,finst]
    Pexp=1-exp(-2pi*cexp^2/Rdot/abs(DF))
    println("The estimated coupling is $c, while the real coupling is $cexp")
    println("The Landau Zener probability is $Pexp")
    Fa=dHd[inist,inist]
    Fb=dHd[finst,finst]
    Popt=1-exp(-chi)
    println("The optimal probability ratio is $Popt")
    Rdotopt=sqrt(2/mass*(TE0+Popt*w21-E[finst]))
    copt=sqrt(abs(DF)*Rdotopt/2/pi)
    println("The optimal coupling parameter would then be $copt")
    dia_M=2pi*cexp^2/Rdot/abs(DF)
    ad_M=1/abs(Rdot*Γ[1][1,2]/w21)
    @show dia_M
    @show ad_M
    @show(ad_M-dia_M)
    @show(pi*ad_M/dia_M)
    @show Γ[1][inist,finst]
    @show abs(DF)/(4*cexp)
end

function adiabatic_LZ_estimator(R,R0,p0,P=1,inist=1,finst=2)
    C0=zeros(nsts)
    C0[inist]=1
    E0,Γ0,F0,W0,Ua0,dHd0=adiabatic_values(R0)
    E,Γ,F,W,Ua,dHd=adiabatic_values(R)
    Eini=E0[1]+p0^2/2/mass
    p=sqrt(2*mass*(Eini-P*E[1]-(1-P)*E[2]))
    Rdot=p/mass
    xi=pi*(E[2]-E[1])/4/Rdot/abs(Γ[1][1,2])
    @show p
    @show Γ[1][1,2]
    PROB=1-exp(-xi)
    return PROB
end

function quiet_LZ_estimator(R,R0,p0,P,chi=1,inist=1,finst=2)
    C0=zeros(nsts)
    C0[inist]=1
    E0,Γ0,F0,W0,Ua0,dHd0=adiabatic_values(R0)
    E,Γ,F,W,Ua,dHd=adiabatic_values(R)
    Hd,dHd=diabatic_values(R)
    dHd=dHd[1]
    DF=dHd[finst,finst]-dHd[inist,inist]
    TE0=real(C0'*E0*C0)[inist]+p0^2/2/mass
    w21=E[finst]-E[inist]
    Rdot=sqrt(2/mass*(TE0+P*w21-E[finst]))
    c=sqrt(abs(DF)*Rdot*log(1/(1-P))/2/pi)
    cexp=Hd[inist,finst]
    Pexp=1-exp(-2pi*cexp^2/Rdot/abs(DF))
    Fa=dHd[inist,inist]
    Fb=dHd[finst,finst]
    Popt=1-exp(-chi)
    Rdotopt=sqrt(2/mass*(TE0+Popt*w21-E[finst]))
    copt=sqrt(abs(DF)*Rdotopt/2/pi)
    dia_M=2pi*cexp^2/Rdot/abs(DF)
    ad_M=1/abs(Rdot*Γ[1][inist,finst]/w21)
    return pi*ad_M/dia_M
end


const nsts=2
potential(R)=pot_NaCl(R)
adiabatic_LZ_estimator(11,1,0,1)
LZ_estimator(15,1,0,0.5,100)


const nsts=2
c=0.0038
potential(R)=pot_sinus(R,0.5,0.02,c)
LZ_estimator(-18.85,-22,11.5,0.5,100)

const nsts=9
beta=0.015
potential(R)=pot_1layer(R,beta)
LZ_estimator(-17.5,-24,0,0.5)

const nsts=2
potential(R)=pot_simple_crossing(R)
LZ_estimator(0,-10,11.5,0.5)

const nsts=90
beta=0.0126
potential(R)=pot_layered(R,beta)
LZ_estimator(-17.5,-24,0,0.5,1,1,11)


using Plots
plotlyjs()
const nsts=2
C=logspace(-2,-5,1000)
F=zeros(C)
for (i,c) in enumerate(C)
    global potential(R)=pot_sinus(R,0.5,0.02,c)
    F[i]=quiet_LZ_estimator(-18.85,-22,11.5,0.5,1,1,2)
end
plot(C,F,xscale=:log)


const nsts=9
C=linspace(0.005,0.5,1000)
F=zeros(C)
for (i,c) in enumerate(C)
    global potential(R)=pot_1layer(R,c)
    F[i]=quiet_LZ_estimator(-17.5,-24,0,0.5,1,1,2)
end
plot(C,F)


const nsts=90
C=linspace(0.002,0.05,1000)
F=zeros(C)
for (i,c) in enumerate(C)
    global potential(R)=pot_layered(R,c)
    F[i]=quiet_LZ_estimator(-17.5,-24,0,0.5,1,1,11)
end
plot(C,F,xscale=:log)
