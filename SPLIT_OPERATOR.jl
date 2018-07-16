########################################
##########     SO FILE    ##############
########################################

input_name=ARGS[1]

println("Including function modules")
include("SO_Initial_data/"*input_name*".jl")
########SETUP DONE##################
println("Starting SO dynamics...")
Tf,X,P,PSI_AD=SO_propagation(tf,pow,xmin,xmax,R0,p0,C0,potential,flags);

println("Saving...")
SO_save(file,Tf,X,P,PSI_AD)

println("DONE!")


#=          ###########  Animation of the wavepacket!
using Plots
plotlyjs()
PX=Prob_X.(PSI_AD)
anim = @animate for i in 1:length(PX)
    plot(X,PX[i])
end
gif(anim,file*".gif",fps=15)
=#
