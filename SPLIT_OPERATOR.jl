########################################
##########     SO FILE    ##############
########################################

input_name=ARGS[1]

println("Including function modules")

input_name=split(input_name,".")[1] #INPUT NAME SHOULD NEVER HAVE DOTS!
include("SO_Initial_data/"*input_name*".jl")

########### Text output
println("The timestep is")
@show δt
println("The mass is")
@show mass
println("The final time, in a.u. is")
@show tf
println("The number of grid points is")
@show 2^pow
println("The name of the file for the savename is")
@show file
println("The boxsize is")
@show xmax
@show xmin
@show xmax-xmin

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
