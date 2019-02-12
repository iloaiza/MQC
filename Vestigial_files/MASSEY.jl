const mass=2000
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
    Rdot_adia=sqrt(2/mass*(TE0-E[inist]))
    c=sqrt(abs(DF)*Rdot*log(1/(1-P))/2/pi)
    c_adia=sqrt(abs(DF)*Rdot_adia*log(1/(1-P))/2/pi)
    cexp=Hd[inist,finst]
    Pexp=1-exp(-2pi*cexp^2/Rdot/abs(DF))
    println("The estimated coupling is $c, while the real coupling is $cexp")
    println("From the GS energy, the estimated coupling is $c_adia")
    println("The Landau Zener probability is $Pexp")
    Fa=dHd[inist,inist]
    Fb=dHd[finst,finst]
    Popt=1-exp(-chi)
    println("The optimal probability ratio is $Popt")
    Rdotopt=sqrt(2/mass*(TE0+Popt*w21-E[finst]))
    Rdotopt_adia=sqrt(2/mass*(TE0-E[inist]))
    copt=sqrt(abs(DF)*Rdotopt/2/pi)
    copt_adia=sqrt(abs(DF)*Rdotopt_adia/2/pi)
    println("The optimal coupling parameter would then be $copt")
    println("The optimal coupling parameter from GS energy would be $copt_adia")
    dia_M=2pi*cexp^2/Rdot/abs(DF)
    ad_M=pi*w21/abs(4*Rdot*Γ[1][inist,finst])
    opt_d12=pi*w21/4/Rdotopt_adia
    @show dia_M
    @show ad_M
    @show(ad_M-dia_M)
    @show(pi*ad_M/dia_M)
    @show Γ[1][inist,finst]
    @show opt_d12
    @show abs(DF)/(4*cexp)
end

function Massey_calculator(R,R0,p0,inist=1,finst=2)
    C0=zeros(nsts)
    C0[inist]=1
    E0,Γ0,F0,W0,Ua0,dHd0=adiabatic_values(R0)
    E,Γ,F,W,Ua,dHd=adiabatic_values(R)
    Hd,dHd=diabatic_values(R)
    dHd=dHd[1]
    DF=dHd[finst,finst]-dHd[inist,inist]
    TE0=real(C0'*E0*C0)[inist]+p0^2/2/mass
    w21=E[finst]-E[inist]
    Rdot_adia=sqrt(2/mass*(TE0-E[inist]))
    Fa=dHd[inist,inist]
    Fb=dHd[finst,finst]
    Rdot_adia=sqrt(2/mass*(TE0-E[inist]))
    ad_M=pi*w21/abs(4*Rdot_adia*Γ[1][inist,finst])
    @show ad_M
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
c=0.0044
potential(R)=pot_sinus(R,0.5,0.02,c)
Massey_calculator(-18.85,-22,15)


const nsts=9
beta=0.01078
potential(R)=pot_1layer(R,beta)
LZ_estimator(-17.5,-24,0,0.5)

const nsts=2
potential(R)=pot_simple_crossing(R)
LZ_estimator(0,-10,11.5,0.5)

const nsts=90
beta=0.003
potential(R)=pot_layered(R,beta)
Massey_calculator(-17.5,-22.5,0,1,11)
LZ_estimator(-17.5,-22.5,0,1,1,1,11)


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
