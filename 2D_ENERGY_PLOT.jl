using Plots
using HDF5
using Distributed
using LinearAlgebra
#using Interact
pyplot()

if length(ARGS) == 0
    Warning("No input name, use Initial data file for loading potential")
    println("2D_CMFSH default")
    input_name = "2D_CMFSH"
else
    input_name = ARGS[1]
end
    
include("./Initial_data/"*input_name*".jl")
include("include.jl")

using LaTeXStrings
using Plots.PlotMeasures
using ColorSchemes

println("Loaded packages, starting plotting routine...")
#Config input for plot
FONT = font(60)
FONT2 = font(40)
L_MARG=[0mm 0mm]
B_MARG=[0mm 0mm]
SIZE=[1880,1220]

# More parameters
res=500
yrange=[-25,18]
xrange=[-25,18]
STATES = 1:20 #STATES[1] must always be 1, always plot GS
cScheme = :thermal

println("Starting energies calculation...")
@time R=energies_4_plot([xrange[1],yrange[1]],[xrange[2],yrange[2]],[res,res]);
@time Xbase=[R[i][1] for i in 1:res];
@time E=[R[i][k] for i in 1:res, k in 2:nsts+1]; #E[:,i] will be the i-th adiabatic PES


X=range(xrange[1],stop=xrange[2],length=res);
Y=range(yrange[1],stop=yrange[2],length=res);

for st in STATES
    E_string="E$(st)=[R[i,j][2+$(st)] for i in 1:length(X), j in 1:length(Y)]"
    eval(Meta.parse(E_string))
end


println("Starting plot...")
P3 = surface(Y,X,E1,color=cScheme)
for n in STATES[2:end]
    surface!(Y,X,eval(Meta.parse("E$n")),color=cScheme)
end
surface!(xtickfont=FONT2, ytickfont=FONT2, zticks=false,size=[1880,1280],
    legendfont=FONT,xlabel="\n \n x",ylabel="\n \n y",zlabel="\n Energy (a.u.)",xguidefont=FONT,
    yguidefont=FONT,zguidefont=FONT,camera=(-22,23),left_margin=L_MARG,bottom_margin=B_MARG)
savefig(P3,"2D_POT.png")

println("Finished 2D plotting!")
