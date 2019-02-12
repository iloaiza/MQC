############################# USING STUFF
using HDF5
using Plots
plotlyjs()

############################# INPUT STUFF
input_name="24layered_015"

#filename=""
############################# LOAD SAVE FILES
include("Initial_data/"*input_name*".jl")
file="MQC/data/"*potname*"_R0($R0)_p0($p0).h5"
#file="MQC/data/"*filename*".h5"
SOfile="MQC/data/SO_"*potname*"_R0($R0)_p0($p0).h5"
SOfile="/home/nacho/Desktop/Desktop/HEAVY_DATA/SO_"*potname*"_R0($R0)_p0($p0).h5"

############################# HISTOGRAM STUFF
HISTO_RES=150;
xmin=-30;
xmax=32;
###SO

include("split_operator_functions_parallel.jl")
T_SO,X_SO,P_SO,PSI_AD_SO,PX_SO=SO_read(SOfile);           #LOAD SPLIT OPERATOR
Y,PX0=SO_histo_builder(X_SO,PX_SO[1],HISTO_RES,xmin,xmax);
Y,PXf=SO_histo_builder(X_SO,PX_SO[end],HISTO_RES,xmin,xmax);
PXf=PXf./maximum(PX0)
PX0=PX0./maximum(PX0)
###

E0=h5read(file,"E0");                                     #LOAD INITIAL ENERGY
T_EH,R_EH,P_EH,C_EH=EH_read(file);                        #LOAD EHRENFEST
T_FSSH,R_FSSH,P_FSSH,C_FSSH,AST_FSSH=FSSH_read(file);     #LOAD SURFACE HOPPING
T_CM2_V,R_CM2_V,P_CM2_V,D_CM2_V=CM2_VANILLA_read(file);
T_CM3_V,R_CM3_V,P_CM3_V,D_CM3_V=CM3_VANILLA_read(file);
T_CM2,R_CM2,P_CM2,D_CM2=CM2_read(file);
T_CM3,R_CM3,P_CM3,D_CM3=CM3_read(file);
T_CM2_F,R_CM2_F,P_CM2_F,D_CM2_F=CM2_FSSH_read(file);
T_CM3_F,R_CM3_F,P_CM3_F,D_CM3_F=CM3_FSSH_read(file);
T_CM2_F_V,R_CM2_F_V,P_CM2_F_V,D_CM2_F_V=CM2_FSSH_VANILLA_read(file);
T_CM3_F_V,R_CM3_F_V,P_CM3_F_V,D_CM3_F_V=CM3_FSSH_VANILLA_read(file);



R_eh,HR_eh=super_histo(R_EH[:,1,:],xmin,xmax,HISTO_RES);
R_fssh,HR_fssh=super_histo(R_FSSH[:,1,:],xmin,xmax,HISTO_RES);
R_cm2,HR_cm2=super_histo(R_CM2[:,1,:],xmin,xmax,HISTO_RES);
R_cm3,HR_cm3=super_histo(R_CM3[:,1,:],xmin,xmax,HISTO_RES);
R_cm2_v,HR_cm2_v=super_histo(R_CM2_V[:,1,:],xmin,xmax,HISTO_RES);
R_cm3_v,HR_cm3_v=super_histo(R_CM3_V[:,1,:],xmin,xmax,HISTO_RES);
R_cm2_f,HR_cm2_f=super_histo(R_CM2_F[:,1,:],xmin,xmax,HISTO_RES);
R_cm3_f,HR_cm3_f=super_histo(R_CM3_F[:,1,:],xmin,xmax,HISTO_RES);
R_cm2_f_v,HR_cm2_f_v=super_histo(R_CM2_F_V[:,1,:],xmin,xmax,HISTO_RES);
R_cm3_f_v,HR_cm3_f_v=super_histo(R_CM3_F_V[:,1,:],xmin,xmax,HISTO_RES);
                                        #NORMALIZE THE HISTOGRAMS
HR_eh=HR_eh./maximum(HR_eh[:,1])
HR_fssh=HR_fssh./maximum(HR_fssh[:,1])
HR_cm2_v=HR_cm2_v./maximum(HR_cm2_v[:,1])
HR_cm3_v=HR_cm3_v./maximum(HR_cm3_v[:,1])
HR_cm2=HR_cm2./maximum(HR_cm2[:,1])
HR_cm3=HR_cm3./maximum(HR_cm3[:,1])
HR_cm2_f=HR_cm2_f./maximum(HR_cm2_f[:,1])
HR_cm3_f=HR_cm3_f./maximum(HR_cm3_f[:,1])
HR_cm2_f_v=HR_cm2_f_v./maximum(HR_cm2_f_v[:,1])
HR_cm3_f_v=HR_cm3_f_v./maximum(HR_cm3_f_v[:,1]);


############### HISTOGRAM PLOTTING
GUIDEFONT=font(24,"Helvetica");
TICKFONT=font(24,"Helvetica");
P=plot(Y,PX0,line=(2.0,:dash),color=:grey);
plot!(Y,PXf,line=(2.5,:dot),color=:black);
#P=plot(R_eh[:,1],HR_eh[:,1],label="EH_i",line=(1,:dash),color=:black);
plot!(R_eh[:,end],-HR_eh[:,end],label="EH_f",line=(2,:solid),color=:blue);
plot!(R_fssh[:,end],HR_fssh[:,end],label="FSSH_f",line=(2.5,:solid),color=:red);
#plot!(R_fssh_d[:,end],HR_fssh_d[:,end],label="FSSH_d_f",line=(1.5,:solid),color=:green);
plot!(R_cm2[:,end],HR_cm2[:,end],label="CM2",line=(2,:solid))
plot!(R_cm3[:,end],HR_cm3[:,end],label="CM3",line=(2,:solid))
plot!(R_cm2_v[:,end],HR_cm2_v[:,end],label="CM2_V",line=(2,:solid))
plot!(R_cm3_v[:,end],HR_cm3_v[:,end],label="CM3_V",line=(2,:solid))
plot!(R_cm2_f[:,end],HR_cm2_f[:,end],label="CM2_F",line=(2,:solid))
plot!(R_cm3_f[:,end],HR_cm3_f[:,end],label="CM3_F",line=(2,:solid))
plot!(R_cm2_f_v[:,end],HR_cm2_f_v[:,end],label="CM2_F_V",line=(2,:solid))
plot!(R_cm3_f_v[:,end],-HR_cm3_f_v[:,end],label="CM3_F_V",line=(2,:solid))
plot!(xlabel="Position (a.u.)",ylabel="Nuclear distribution",xlims=(xmin,xmax));
plot!(xguidefont = GUIDEFONT,xtickfont=TICKFONT,yguidefont = GUIDEFONT,ytickfont=TICKFONT);
plot!(legend=false);
display(P)
ylims!(-0.25,0.25)
#plot!(legend=true)


using Interact
@manipulate for i in 1:101
    P=plot(R_eh[:,1],HR_eh[:,1],label="EH_i",line=(1,:dash),color=:gray);
    plot!(R_eh[:,i],HR_eh[:,i],label="EH_f",line=(2,:solid),color=:blue);
    plot!(R_fssh[:,i],HR_fssh[:,i],label="FSSH_f",line=(2.5,:solid),color=:red);
    #plot!(R_fssh_d[:,end],HR_fssh_d[:,end],label="FSSH_d_f",line=(1.5,:solid),color=:green);
    plot!(R_cm2[:,i],HR_cm2[:,i],label="CM2",line=(2,:solid))
    plot!(R_cm3[:,i],HR_cm3[:,i],label="CM3",line=(2,:solid))
    plot!(R_cm2_v[:,i],HR_cm2_v[:,i],label="CM2_V",line=(2,:solid))
    plot!(R_cm3_v[:,i],HR_cm3_v[:,i],label="CM3_V",line=(2,:solid))
    plot!(R_cm2_f[:,i],HR_cm2_f[:,i],label="CM2_F",line=(2,:solid))
    plot!(R_cm3_f[:,i],HR_cm3_f[:,i],label="CM3_F",line=(2,:solid))
    plot!(R_cm2_f_v[:,i],HR_cm2_f_v[:,i],label="CM2_F_V",line=(2,:solid))
    plot!(R_cm3_f_v[:,i],HR_cm3_f_v[:,i],label="CM3_F_V",line=(2,:solid))
    plot!(xlabel="Position (a.u.)",ylabel="Nuclear distribution",xlims=(xmin,xmax));
    plot!(xguidefont = GUIDEFONT,xtickfont=TICKFONT,yguidefont = GUIDEFONT,ytickfont=TICKFONT);
    display(P)
end
################################ ENERGY STUFF
xmin=-24;
xmax=24;
res=1000;
for x in linspace(xmin,xmax,100)
    #@show x
    H,dH=pot_NaCl(x)
end
R=all_4_1D_plot(xmin,xmax,res);
Xbase=R[:,1];
E=R[:,2:1+nsts];
G=R[:,2+nsts:end];
num_couplings=length(G[1,:])
S_EH=EH_state(R0,p0,C0)
Ep0,~,~,~,~,~=adiabatic_values(R0);
E0=real(C0'*diagm(Ep0)*C0)[1]+p0^2/2/mass;
############### ENERGY PLOTS
COLOURS=[:blue,:red,:green,:cyan,:purple];
GUIDEFONT=font(24,"Helvetica");
TICKFONT=font(24,"Helvetica");
P=plot(Xbase,E[:,1],line=(1.5,:solid),color=:blue)
for i in 2:nsts
    plot!(Xbase,E[:,i],line=(1.5,:solid),color=COLOURS[mod1(i,length(COLOURS))])
end
NACfact=400
for i in 1:num_couplings
    plot!(Xbase,G[:,i]/NACfact,line=(2.0,:dash),color=:green)
end
plot!(x->E0,line=(2.0,:dash),colour=:black);
plot!(xlabel="Position (a.u.)",ylabel="Energy (a.u.)",xlims=(xmin,xmax));
plot!(xguidefont = GUIDEFONT,xtickfont=TICKFONT,yguidefont = GUIDEFONT,ytickfont=TICKFONT);
plot!(legend=false);
title!("NAC/$NACfact",titlefont = GUIDEFONT)
ylims!(-0.011,0.1)
display(P)
#yaxis!()

################ EH and SH PES
GUIDEFONT=font(24,"Helvetica");
TICKFONT=font(24,"Helvetica");
Efact=1000
P=plot(Xbase,E[:,1].*Efact,line=(1.5,:solid),color=:blue);
plot!(Xbase,E[:,2].*Efact,line=(1.5,:solid),color=:blue);

ktraj1=1;
ktraj2=3;
Reh=R_EH[:,1,ktraj1]
Peh=P_EH[:,1,ktraj1]
Ceh=C_EH[:,:,ktraj1]
Rf1=R_FSSH[:,1,ktraj1]
Pf1=P_FSSH[:,1,ktraj1]
Af1=AST_FSSH[:,ktraj1]
Rf2=R_FSSH[:,1,ktraj2]
Pf2=P_FSSH[:,1,ktraj2]
Af2=AST_FSSH[:,ktraj2]
Af2[35:end].+=1;
#~,Eeh,~=EH_energy_check(R_EH[:,1,ktraj1],P_EH[:,1,ktraj1],C_EH[[:,:,ktraj1]]);
~,Eeh,~=EH_energy_check(Reh,Peh,Ceh);
#~,Efssh1,~=FSSH_energy_check(R_FSSH[:,1,ktraj1],P_FSSH[:,ktraj1],AST_FSSH[[:,ktraj1]]);
#~,Efssh2,~=FSSH_energy_check(R_FSSH[:,1,ktraj1],P_FSSH[:,ktraj1],AST_FSSH[[:,ktraj2]]);
~,Efssh1,~=FSSH_energy_check(Rf1,Pf1,Af1);
~,Efssh2,~=FSSH_energy_check(Rf2,Pf2,Af2);


plot!(R_EH[:,ktraj1],Eeh.*Efact,line=(3,:dot),color=:black);
plot!(R_FSSH[:,ktraj1],Efssh1.*Efact,line=(3,:dot),color=:red);
plot!(R_FSSH[:,ktraj2],Efssh2.*Efact,line=(3,:dot),color=:green);
plot!(legend=false,xlabel="Position (a.u.)",ylabel="Energy (a.u.x10⁻³)");
plot!(xguidefont = GUIDEFONT,xtickfont=TICKFONT,yguidefont = GUIDEFONT,ytickfont=TICKFONT);
display(P)
xlims!(-10,11);


############## 2D ENERGY PLOTS
Xbase=[R[i][1] for i in 1:res];
E=[R[i][k] for i in 1:res, k in 2:nsts+1]; #E[:,i] will be the i-th adiabatic PES

X=linspace(-20,25,225);
Y=linspace(-20,25,225);
E1=[R[i,j][3] for i in 1:length(X), j in 1:length(Y)];
E2=[R[i,j][4] for i in 1:length(X), j in 1:length(Y)];

surface(X,Y,E1)
surface!(X,Y,E2)
# =
R=energies_4_plot([-20,-20],[25,25],[225,225]);
GAMMA=zeros(Float64,size(R))
for (i,r) in enumerate(R)
    ~,G,~,~,~,~=adiabatic_values(r[1:NDOFs],NDOFs)
    GAMMA[i]=G[2][1,2]
end
=#


############## SINGLE TRAJECTORIES
input_name="1layer_01056"
include("Initial_data/"*input_name*".jl")
file="MQC/data/SINGLE_"*potname*"_R0($R0)_p0($p0).h5"
T_eh,R_eh,p_eh,C_eh=EH_read(file)
T_cm2,R_cm2,p_cm2,D_cm2=CM2_read(file)
T_cm3,R_cm3,p_cm3,D_cm3=CM3_read(file)
T_cm2_fssh,R_cm2_fssh,p_cm2_fssh,



############## K_ENERGY_SIMULATIONS
input_name="simple"
#filename=""
##
include("Initial_data/"*input_name*".jl")

dirname="MQC/data/K_"*potname*"_R0($R0)"
file=dirname*"/K_SIMULATION.h5"

A_EH=general_K_read(file,"EH",3)
FR_EH=A_EH[1]
FM_EH=A_EH[2]
Fpop_EH=A_EH[3]

A_FSSH=general_K_read(file,"FSSH",4)
FR_FSSH=A_FSSH[1]
FM_FSSH=A_FSSH[2]
Fpop_FSSH=A_FSSH[3]
Fast_FSSH=A_FSSH[4]


histogram(FR_FSSH[7,:])
GS_FSSH=GS_MF_prob(K,Fpop_FSSH)
GS_FSSH=GS_SH_prob(K,Fast_FSSH)
GS_EH=GS_MF_prob(K,Fpop_EH)


GS_TRANS_FSSH,GS_REFL_FSSH,ET_TRANS_FSSH=refl_trans_SH(K,FR_FSSH,Fast_FSSH,2)
plot(K,ET_TRANS_FSSH)
plot!(K,GS_REFL_FSSH)
plot!(K,GS_TRANS_FSSH)

plot(K,GS_FSSH)
plot!(K,GS_EH)
ylims!(0.999,1.001)

################## CHANGE SAVE FILES
using HDF5
input_name="24layered_015"
include("Initial_data/"*input_name*".jl")
newfile="MQC/data/"*potname*"_R0($R0)_p0($p0).h5"
NEWNAME="layered_015_24"
file="MQC/data/"*NEWNAME*".h5"
E0=h5read(file,"E0");                                     #LOAD INITIAL ENERGY
T_EH,R_EH,P_EH,C_EH=EH_read(file);                        #LOAD EHRENFEST
T_FSSH,R_FSSH,P_FSSH,C_FSSH,AST_FSSH=FSSH_read(file);     #LOAD SURFACE HOPPING
T_CM2_V,R_CM2_V,P_CM2_V,D_CM2_V=CM2_VANILLA_read(file);
T_CM3_V,R_CM3_V,P_CM3_V,D_CM3_V=CM3_VANILLA_read(file);
T_CM2,R_CM2,P_CM2,D_CM2=CM2_read(file);
T_CM3,R_CM3,P_CM3,D_CM3=CM3_read(file);
T_CM2_F,R_CM2_F,P_CM2_F,D_CM2_F=CM2_FSSH_read(file);
T_CM3_F,R_CM3_F,P_CM3_F,D_CM3_F=CM3_FSSH_read(file);
T_CM2_F_V,R_CM2_F_V,P_CM2_F_V,D_CM2_F_V=CM2_FSSH_VANILLA_read(file);
T_CM3_F_V,R_CM3_F_V,P_CM3_F_V,D_CM3_F_V=CM3_FSSH_VANILLA_read(file);

h5write(newfile,"E0",E0)
EH_save(T_EH,R_EH,P_EH,C_EH,newfile)
FSSH_save(T_FSSH,R_FSSH,P_FSSH,C_FSSH,AST_FSSH,newfile)
CM2_save(T_CM2,R_CM2,P_CM2,D_CM2,newfile)
CM3_save(T_CM3,R_CM3,P_CM3,D_CM3,newfile)
CM2_VANILLA_save(T_CM2_V,R_CM2_V,P_CM2_V,D_CM2_V,newfile)
CM3_VANILLA_save(T_CM3_V,R_CM3_V,P_CM3_V,D_CM3_V,newfile)
CM2_FSSH_save(T_CM2_F,R_CM2_F,P_CM2_F,D_CM2_F,newfile)
CM3_FSSH_save(T_CM3_F,R_CM3_F,P_CM3_F,D_CM3_F,newfile)
CM2_FSSH_VANILLA_save(T_CM2_F_V,R_CM2_F_V,P_CM2_F_V,D_CM2_F_V,newfile)
CM3_FSSH_VANILLA_save(T_CM3_F_V,R_CM3_F_V,P_CM3_F_V,D_CM3_F_V,newfile)
