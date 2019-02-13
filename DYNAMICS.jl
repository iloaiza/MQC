########################################
##########     MQC FILE   ##############
########################################
#for parallel computing, run using 'julia -p N DYNAMICS.jl INITIAL_DATA_FILENAME'
#otherwise, just run using 'julia DYNAMICS.jl INITIAL_DATA_FILENAME'
#the INITIAL_DATA_FILENAME is a .jl file that is in the Initial_data folder
tic0=time()

using Distributed

#push!(ARGS,"24layered_03")  #UNCOMMENT FOR RUNNING INSIDE IDE
#cd("MQC")
@everywhere input_name=remotecall_fetch(i->ARGS[i],1,1)

FIRST_RUN=true
initial_dist=false

println("Including function modules")
if split(input_name,".")[end]=="jl"
    @everywhere include("Initial_data/"*input_name)
else
    @everywhere include("Initial_data/"*input_name*".jl")
end

if initial_dist==false
    initial_dist=wigner
end

########## Stuff for text output (i.e. do a logfile!)
S_EH=EH_state_builder(R0,p0,C0)
filename="./data/"*potname*"_R0($R0)_p0($p0).h5"

Ep0,_,_,_,_,_=adiabatic_values(R0);
E0=real(C0'*Diagonal(Ep0)*C0)[1]+sum(abs2.(p0))/2/mass;


if FIRST_RUN
    h5write(filename,"R0",R0)
    h5write(filename,"p0",p0)
    h5write(filename,"C0_real",real.(C0))
    h5write(filename,"C0_imag",imag.(C0))
    h5write(filename,"E0",E0)
    h5write(filename,"DYN_LIST",DYN_LIST)
end

########### Text output
println("The timestep is")
@show dt
println("The mass is")
@show mass
println("The final time, in a.u. is")
@show tf
println("The number of trajectories is")
@show Ntrajs
println("The name of the file for the savename is")
@show filename
println("The initial conditions are")
@show R0
@show p0
@show C0
@show E0

NDOFs=length(R0)
println("The dynamics that will be done are:")
for DYN in DYN_LIST
    dyn_sts=eval(Meta.parse("$(DYN)_sts"))
    if eval(Meta.parse("@isdefined $(DYN)_mem")) #block for filling memory variable with 0 if it is not defined
        eval(Meta.parse("mem=$(DYN)_mem"))
        println("""Using memory $(DYN)_mem = $(mem)""")
    else
        eval(Meta.parse("$(DYN)_mem=0"))
        println("Memory not defined for $(DYN), make sure it's intentional")
    end
    if DYN in CL_LIST
        println("$DYN CLASSICAL INTEGRATION")
        S_string="S_$(DYN)=builder_CL_state($R0,$p0,$DYN,0,$NDOFs,$(DYN)_mem)"
        eval(Meta.parse(S_string))
    elseif DYN in MF_LIST
        println("$DYN MEAN-FIELD INTEGRATION")
        if dyn_sts!=nsts
            println("Warning! Still no implementation for $C0 initial condition for $DYN, using [1,0,...] instead")
            D0=zeros(Complex,dyn_sts)
            D0[1]=1
        else
            D0=C0
        end
        S_string="S_$(DYN)=builder_MF_state($R0,$p0,$D0,$DYN,0,$NDOFs,$(DYN)_mem)"
        eval(Meta.parse(S_string))
    elseif DYN in SH_LIST
        println("$DYN SURFACE HOPPING INTEGRATION")
        if dyn_sts!=nsts
            println("Warning! Still no implementation for $C0 initial condition for $DYN, using [1,0,...] instead")
            D0=zeros(Complex,dyn_sts)
            D0[1]=1
        else
            D0=C0
        end
        S_string="S_$(DYN)=builder_SH_state($R0,$p0,$D0,$ast0,$DYN,0,$NDOFs,$(DYN)_mem)"
        eval(Meta.parse(S_string))
    else
        println("$DYN has no implementation, check the list in types.jl in the end section of META")
    end
end

for DYN in DYN_LIST
    dyn_sts=eval(Meta.parse("$(DYN)_sts"))
    if DYN in CL_LIST
        println("BEGINNING $DYN CLASSICAL INTEGRATION")
        te0=time()
        if initial_dist==constant_dist
            println("Constant distribution: only one propagation needed for classical methods")
            integ_string="single_integration(tf,S_$(DYN))"
            T,R,P=eval(Meta.parse(integ_string))
        else
            integ_string="dist_CL_integration($tf,$R0,$p0,$(DYN)_mem,$DYN,$Ntrajs,$initial_dist)"
            Ts,Rs,Ps=eval(Meta.parse(integ_string))
            T=zeros(Float64,size(Ts))
            R=zeros(Float64,size(Rs))
            P=zeros(Float64,size(Ps))
            T.=Ts
            R.=Rs
            P.=Ps
        end
        tef=time()
        println("FINISHED $DYN INTEGRATION FOR $(round(tef-te0,digits=3))s, SAVING...")
        CL_save(T,R,P,filename,DYN)
    elseif DYN in MF_LIST
        println("BEGINNING $DYN MEAN-FIELD INTEGRATION")
        te0=time()
        if dyn_sts!=nsts
            println("Warning! Still no implementation for $C0 initial condition for $DYN, using [1,0,...] instead")
            D0=zeros(Complex,dyn_sts)
            D0[1]=1
        else
            D0=C0
        end
        if initial_dist==constant_dist
            println("Constant distribution: only one propagation needed for mean-field methods")
            integ_string="single_integration(tf,S_$(DYN))"
            T,R,P,C=eval(Meta.parse(integ_string))
        else
            integ_string="dist_MF_integration($tf,$R0,$p0,$D0,$(DYN)_mem,$DYN,$Ntrajs,$initial_dist)"
            Ts,Rs,Ps,Cs=eval(Meta.parse(integ_string))
            T=zeros(Float64,size(Ts))
            R=zeros(Float64,size(Rs))
            P=zeros(Float64,size(Ps))
            C=zeros(Complex,size(Cs))
            T.=Ts
            R.=Rs
            P.=Ps
            C.=Cs
        end
        tef=time()
        println("FINISHED $DYN INTEGRATION FOR $(round(tef-te0,digits=3))s, SAVING...")
        MF_save(T,R,P,C,filename,DYN)
    elseif DYN in SH_LIST
        println("BEGINNING $DYN SURFACE-HOPPING INTEGRATION")
        te0=time()
        if dyn_sts!=nsts
            println("Warning! Still no implementation for $C0, using [1,0,...] instead")
            D0=zeros(Complex,dyn_sts)
            D0[1]=1
        else
            D0=C0
        end
        integ_string="dist_SH_integration($tf,$R0,$p0,$D0,$ast0,$(DYN)_mem,$DYN,$Ntrajs,$initial_dist)"
        Ts,Rs,Ps,Cs,ASTs=eval(Meta.parse(integ_string))
        T=zeros(Float64,size(Ts))
        R=zeros(Float64,size(Rs))
        P=zeros(Float64,size(Ps))
        C=zeros(Complex,size(Cs))
        AST=zeros(Int,size(ASTs))
        T.=Ts
        R.=Rs
        P.=Ps
        C.=Cs
        AST.=ASTs
        tef=time()
        println("FINISHED $DYN INTEGRATION FOR $(round(tef-te0,digits=3))s, SAVING...")
        SH_save(T,R,P,C,AST,filename,DYN)
    end
end


println("FINISHED SUCCESSFULLY!")
tic1=time()
println("The total time of running was $(round(tic1-tic0,digits=4))s")
