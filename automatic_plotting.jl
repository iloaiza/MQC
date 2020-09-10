println("Loading plotting packages...")
using Plots
using Plots.PlotMeasures
using ColorSchemes
#eval(Meta.parse(plot_method))
pyplot()
println("Finished loading!")
FONT=font(40)
SIZE=[1780,880]
L_MARG=[15mm 0mm]
B_MARG=[10mm 0mm]
SIZE2D=[1880,1280]
L_MARG_2D=[10mm 0mm]
B_MARG_2D=[8mm 0mm]
plotname="./plots/"
@show plotname
@show filename
COLOR_LIST = [:black,:red,:blue,:green,:purple]
LINE_LIST = [:solid,:solid,:dash,:dash,:dot]
MARK_LIST = [:circle,:xcross,:square,:dtriangle,:utriangle]
MARKER_SIZE = 20
FLAG_SKIP = 3 #1 for not skipping flags, >1 for skipping flags in BO deviation plot
SCATTER_SKIP = 40 #for only adding scatter markers every SCATTER_SKIP steps
@show FLAG_SKIP
@show SCATTER_SKIP

if GEN_METHOD == "DYNAMICS"
    if length(R0) == 1 #1 nuclear dimension
        for dyn in DYN_LIST
            println("Loading $dyn...")
            if dyn in CL_LIST
                str="""T_$(dyn),R_$(dyn),P_$(dyn)=CL_read("$filename","$dyn")"""
            elseif dyn in MF_LIST
                str="""T_$(dyn),R_$(dyn),P_$(dyn),C_$(dyn)=MF_read("$filename","$dyn")"""
            elseif dyn in SH_LIST
                str="""T_$(dyn),R_$dyn,P_$(dyn),C_$(dyn),AST_$(dyn)=SH_read("$filename","$dyn")"""
            end
            eval(Meta.parse(str))
        end

        if HISTO_PLOTS
            XLIMS=[R0+10.0,R0-10.0]
            for dyn in DYN_LIST
                min_str = "minimum(R_$(dyn))"
                max_str = "maximum(R_$(dyn))"
                mini = eval(Meta.parse(min_str))
                maxi = eval(Meta.parse(max_str))
                if XLIMS[1] > mini
                    XLIMS[1] = mini
                end
                if XLIMS[2] < maxi
                    XLIMS[2] = maxi
                end
            end

            ΔX=XLIMS[2]-XLIMS[1]
            XLIMS[1]=XLIMS[1]-0.05*ΔX
            XLIMS[2]=XLIMS[2]+0.05*ΔX

            HLIMS=[1.0,-1,1] #third one holds mirror for plot clarity
            for dyn in DYN_LIST
                println("Building histograms for $dyn")
                str="Rbase_$dyn,HR_$dyn=super_histo(R_$dyn[:,1,:],$(XLIMS[1]),$(XLIMS[2]),$HISTO_RES)"
                eval(Meta.parse(str))
                str="HR_$dyn=$(HLIMS[3]) .* HR_$dyn ./ maximum(HR_$dyn[:,1])"
                eval(Meta.parse(str))
                HLIMS[3] = HLIMS[3]*-1

                min_str = "minimum(HR_$(dyn)[:,end])"
                max_str = "maximum(HR_$(dyn)[:,end])"
                mini = eval(Meta.parse(min_str))
                maxi = eval(Meta.parse(max_str))
                if HLIMS[1] > mini
                    HLIMS[1] = mini
                end
                if HLIMS[2] < maxi
                    HLIMS[2] = maxi
                end
            end

            dyn_1=DYN_LIST[1]
            if plot_ini
                plot_str = """plot(Rbase_$(dyn_1)[:,1],HR_$(dyn_1)[:,1],label="$(dyn_1)_i",line=(2,:dash),color=:black)"""
            else
                plot_str = """plot()"""
            end
            flagTot = length(eval(Meta.parse("HR_$dyn_1[1,:]")))
            for flag in 1:flagTot
                P=eval(Meta.parse(plot_str))
                dynCount = [0]
                for dyn in DYN_LIST
                    dynCount[1] += 1
                    plot_str = """plot!(Rbase_$(dyn)[:,$flag],HR_$(dyn)[:,$flag],label="$(dyn)",line=(3.5,LINE_LIST[$(dynCount[1])],COLOR_LIST[$(dynCount[1])]))"""
                    eval(Meta.parse(plot_str))
                end
                plot!(xlabel="Position (a.u.)",xlims=XLIMS,ylabel="Nuclear distribution");
                plot!(xtickfont = FONT,xguidefont=FONT,ytickfont = FONT,yguidefont=FONT,size=SIZE,legendfont=FONT,
                    left_margin=L_MARG,bottom_margin=B_MARG);
                if flag != 1
                    plot!(legend=false)
                end
                #ylims!(hmin-0.1,hmax+0.1)
                savefig(P,plotname*potname*"R0($R0)_p0($p0).$flag.png")
                println("figure for flag $flag saved!")
            end
        end
    elseif length(R0) == 2 #2 nuclear dimensions
        HISTO_RES=100;
        mins=[-25;-25]
        maxs=[19;19]
        dynNum = length(DYN_LIST)
        for dyn in DYN_LIST
            println("Loading $dyn...")
            if dyn in CL_LIST
                str="""T_$(dyn),R_$(dyn),P_$(dyn)=CL_read("$filename","$dyn")"""
            elseif dyn in MF_LIST
                str="""T_$(dyn),R_$(dyn),P_$(dyn),C_$(dyn)=MF_read("$filename","$dyn")"""
            elseif dyn in SH_LIST
                str="""T_$(dyn),R_$dyn,P_$(dyn),C_$(dyn),AST_$(dyn)=SH_read("$filename","$dyn")"""
            end
            eval(Meta.parse(str))
        end

        println("Finished loading, starting histograms...")

        for dyn in DYN_LIST
            tic=time()
            println("$dyn histogram...")
            #eval(Meta.parse("@show R_$dyn"))
            str="Rbase_$dyn,HR_$dyn=multi_d_histo(R_$dyn,$mins,$maxs,$HISTO_RES)"
            eval(Meta.parse(str))
            #eval(Meta.parse("@show HR_$dyn"))
            #str="HR_$dyn=HR_$dyn./maximum(HR_$dyn[:,:,1])"
            #eval(Meta.parse("@show HR_$dyn"))
            #eval(Meta.parse(str))
            println("Finished $dyn histogram after $(time()-tic) seconds!")
        end


        flagTot = length(eval(Meta.parse("HR_$(DYN_LIST[1])[1,1,:]")))
        for dyn in DYN_LIST
            for flag in 1:flagTot
                println("Starting contour plot for $dyn at flag $flag, plotname is $plotname$potname.$dyn.$flag.png")
                P = eval(Meta.parse("contour(Rbase_$(dyn)[1][:,$flag],Rbase_$(dyn)[2][:,$flag],HR_$(dyn)[:,:,$flag],linealpha=0,fill=true,line=false);"))#scale=:exp,c=cgrad([:darkblue, :yellow, :orange, :red], [0.01, 0.1, 0.4, 0.9]));"))
                #P = eval(Meta.parse("heatmap(HR_$(dyn)[:,:,$flag],cgrad([:orange, :blue], [0.1, 0.3, 0.8]),size=SIZE2D);"))
                println("Plotted contour, changing style...")
                #plot!(xtickfont=FONT,ytickfont=FONT,yguidefont=FONT,xguidefont=FONT,legendfont=FONT,size=SIZE2D,left_margin=L_MARG_2D,bottom_margin=B_MARG_2D)
                plot!(xticks=false,yticks=false,yguidefont=FONT,xguidefont=FONT,legendfont=FONT,size=SIZE2D,left_margin=L_MARG_2D,bottom_margin=B_MARG_2D)
                #xticks!([(4k+1) for k in -6:6])
                #yticks!([(4k+1) for k in -6:6])
                if flag != 1
                    plot!(legend=false)
                end
                savefig(P,"$plotname$potname.$dyn.$flag.png")
                println("Plot saved!")
            end
        end

        #= Surface plot of histogram instead of contour
        for dyn in DYN_LIST
            for flag in 1:flagTot
                println("Starting surface plot for $dyn at flag $flag, plotname is $plotname$potname.$dyn.surface.$flag.png")
                P = eval(Meta.parse("surface(Rbase_$(dyn)[1][:,$flag],Rbase_$(dyn)[2][:,$flag],HR_$(dyn)[:,:,$flag]);"))#scale=:exp,c=cgrad([:darkblue, :yellow, :orange, :red], [0.01, 0.1, 0.4, 0.9]));"))
                #P = eval(Meta.parse("heatmap(HR_$(dyn)[:,:,$flag],cgrad([:orange, :blue], [0.1, 0.3, 0.8]),size=SIZE2D);"))
                println("Plotted contour, changing style...")
                plot!(xtickfont=FONT,ytickfont=FONT,zguidefont=FONT,yguidefont=FONT,xguidefont=FONT,size=SIZE2D,c=:thermal,
                    ztickfont=FONT,xlabel="\n \n x",ylabel="\n \n y",zlabel="\n \n Number of trajectories")
                savefig(P,"$plotname$potname.$dyn.surface.$flag.png")
                println("Plot saved!")
            end
        end
        # =#
    else #NDOFs != 1,2
        println("Still no plot implementation for NDOFs⩾3 (number of nuclear degrees of freedom)")
    end
end

if BO_COMP
    Ndyns = length(DYN_LIST) - 1
    dynCount = [0]
    Rdiffs = zeros(checkpoints+1,Ndyns)
    for dyn in DYN_LIST
        if dyn != "BO"
            dynCount[1] += 1
            #r_string="rm_$(dyn)=zeros($(checkpoints + 1))"
            #eval(Meta.parse(r_string))
            for i in 1:checkpoints+1
                r_string = "mean(R_$(dyn)[$i,1,:])"
                Rdiffs[i,dynCount] .= eval(Meta.parse(r_string)) - R_BO[i]
            end
        end
    end
    P = plot()
    dynCount[1] = 0
    for dyn in DYN_LIST
        if dyn != "BO"
            dynCount[1] += 1
            plot!(eval(Meta.parse("T_$(dyn)[1:$FLAG_SKIP:end]")),Rdiffs[1:FLAG_SKIP:end,dynCount], label=dyn, line=(3.5,COLOR_LIST[dynCount[1]]))
            scatter!(eval(Meta.parse("T_$(dyn)[1:$SCATTER_SKIP:end]")),Rdiffs[1:SCATTER_SKIP:end,dynCount], label="",marker=(MARKER_SIZE,COLOR_LIST[dynCount[1]],MARK_LIST[dynCount[1]]))
        end
    end

    plot!(legend=:bottomleft,xtickfont = FONT,xguidefont=FONT,ytickfont = FONT,yguidefont=FONT,size=SIZE,legendfont=FONT,
        xlabel="Time (a.u.)",ylabel="Nonadiabatic deviation",left_margin=L_MARG,bottom_margin=B_MARG)
    savefig(P,"$plotname$potname.VS_BO.png")
    println("Born-Oppenheimer difference plotted and saved!")
    P = plot()
    for k in 1:length(MARK_LIST)
        plot!([0,1],[k,k], line=(3.5,COLOR_LIST[k]))
        scatter!([0.5],[k],marker=(MARKER_SIZE,MARK_LIST[k],COLOR_LIST[k]))
    end
    plot!(ylims=(-0.5,length(MARK_LIST)+0.5))
    savefig(P,"$plotname$potname.VS_BO.LEGEND.png")
    println("Legend plot for Born-Oppenheimer saved!")
end

if EL_PLOTS == true
    println("Starting electronic dynamics")
    Ndyns = length(DYN_LIST)
    elDynList = String[] 
    for dyn in DYN_LIST
        if !(dyn in CL_LIST)
            push!(elDynList,dyn)
            println(dyn * " included for plotting")
        end
    end
    dynTot = length(elDynList)

    if !(@isdefined(Ntrajs))
        const Ntrajs = 1
    else
        @show Ntrajs
    end
    C = zeros(Complex,checkpoints+1,Ntrajs)
    T = zeros(checkpoints+1,dynTot)
    Cmean = zeros(Complex,checkpoints+1,nsts,dynTot)
    dynNum = [0]
    println("Starting Cmean calculations...")
    for dyn in elDynList
        println("Loading $(dyn)")
        dynNum[1] += 1
        T[:,dynNum[1]] = eval(Meta.parse("T_$(dyn)"))
        if EL_ES_PLOTS
            stTot = eval(Meta.parse(dyn*"_sts"))
        else
            stTot = 1
        end
        for st in 1:stTot
            println("Running over state $st of $stTot")
            if dyn in MF_LIST
                #C[:,st,1,dynNum[1]] = eval(Meta.parse("C_$(dyn)[:,$st]"))
                #Cmean[:,st,dynNum[1]] = C[:,st,1,dynNum[1]]
                Cmean[:,st,dynNum[1]] = eval(Meta.parse("C_$(dyn)[:,$st]"))
            else
                C .= eval(Meta.parse("C_$(dyn)[:,$st,:]"))
                for i in 1:checkpoints+1
                    Cmean[i,st,dynNum[1]] = mean(C[i,:])
                end
            end
        end
        
    end

    P = plot()
    dynNum[1] = 0
    for dyn in elDynList
        dynNum[1] += 1
        plot!(T[:,dynNum[1]],abs2.(Cmean[:,1,dynNum[1]]), label=dyn, line=(3.5,COLOR_LIST[dynNum[1]]))
        scatter!(T[1:SCATTER_SKIP:end,dynNum[1]],abs2.(Cmean[1:SCATTER_SKIP:end,1,dynNum[1]]), label="",marker=(MARKER_SIZE,COLOR_LIST[dynNum[1]],MARK_LIST[dynNum[1]]))
    end

    plot!(legend=false,xtickfont = FONT,xguidefont=FONT,ytickfont = FONT,yguidefont=FONT,
        size=SIZE,legendfont=FONT,xlabel="Time (a.u.) ",ylabel="Ground state population",left_margin=L_MARG,bottom_margin=B_MARG)
    savefig(P,"$plotname$potname.electronic_GS.png")

    if EL_ES_PLOTS
        P = plot()
        dynNum[1] = 0
        for dyn in elDynList
            dynNum[1] += 1
            for st in 1:eval(Meta.parse(dyn*"_sts"))
                plot!(T[:,dynNum[1]],abs2.(Cmean[:,st,dynNum[1]]), label=dyn*"_$st", line=2.5)
            end
        end

        plot!(legend=:bottomleft,xtickfont = FONT,xguidefont=FONT,ytickfont = FONT,yguidefont=FONT,size=SIZE,legendfont=FONT,xlabel="Time (a.u.)",ylabel="|Electronic populations|^2")
        savefig(P,"$plotname$potname.electronic_ALL.png")
        println("Electronic populations plotted and saved!")
    end
end

function isFricGood(T,R,P,C,dynName,tol)
    flagNum = length(T)
    stNum = length(C[1,:])
    Fk_IND = zeros(flagNum,stNum)
    for tnum in 1:flagNum
        S = builder_MF_state(R[tnum],P[tnum],C[tnum,:],dynName)
        for k in 1:stNum
            fdiff = abs(S.el.F[1][1,1]-S.el.F[1][k,k])
            if fdiff > tol
                Fk_IND[tnum,k] = 1
            end
        end
    end

    P = plot()
    for st in 1:20
        plot!(R,Fk_IND[:,st])
    end

    return P
end

if GEN_METHOD == "SINGLE" #running single trajectory
    println("Starting single trajectory plotting routine")
    if NDOFs != 1
        println("No single trajectory implementation for more than 1 nuclear degrees of freedom")
    else
        for dyn in DYN_LIST #LOAD
            println("Loading $dyn...")
            if dyn in CL_LIST
                str="""T_$(dyn),R_$(dyn),P_$(dyn)=CL_read("$filename","$dyn")"""
            elseif dyn in MF_LIST
                str="""T_$(dyn),R_$(dyn),P_$(dyn),C_$(dyn)=MF_read("$filename","$dyn")"""
            elseif dyn in SH_LIST
                str="""T_$(dyn),R_$dyn,P_$(dyn),C_$(dyn),AST_$(dyn)=SH_read("$filename","$dyn")"""
            end
            eval(Meta.parse(str))
        end

        println("Starting time-dependent plots")
        P = plot(xtickfont=FONT,ytickfont=FONT,yguidefont=FONT,xguidefont=FONT,legendfont=FONT,size=SIZE)
        plot!(xlabel="Time (a.u.)",ylabel="R (a.u.)")
        for dyn in DYN_LIST
            eval(Meta.parse("""plot!(T_$(dyn),R_$(dyn),label="$(dyn)")"""))
        end
        savefig(P,"$(plotname)TD-R.png")
        println("R plot saved as $plotname.R.png")

        P = plot(xtickfont=FONT,ytickfont=FONT,yguidefont=FONT,xguidefont=FONT,legendfont=FONT,size=SIZE)
        plot!(xlabel="Time (a.u.)",ylabel="p (a.u.)")
        for dyn in DYN_LIST
            eval(Meta.parse("""plot!(T_$(dyn),P_$(dyn),label="$(dyn)")"""))
        end
        savefig(P,"$(plotname)TD-p.png")
        println("p plot saved as $plotname.p.png")

        P = plot(xtickfont=FONT,ytickfont=FONT,yguidefont=FONT,xguidefont=FONT,legendfont=FONT,size=SIZE)
        plot!(xlabel="Time (a.u.)",ylabel="|C0|²")
        for dyn in DYN_LIST
            eval(Meta.parse("""plot!(T_$(dyn),abs2.(C_$(dyn)[:,1]),label="$(dyn)")"""))
        end
        savefig(P,"$(plotname)TD-C0.png")
        println("c plot saved!")


        println("Starting position-dependent plots")
        P = plot(xtickfont=FONT,ytickfont=FONT,yguidefont=FONT,xguidefont=FONT,legendfont=FONT,size=SIZE)
        plot!(xlabel="Position R (a.u.)",ylabel="Momentum p (a.u.)")
        for dyn in DYN_LIST
            eval(Meta.parse("""plot!(R_$(dyn),P_$(dyn),label="$(dyn)")"""))
        end
        savefig(P,"$(plotname)RD-P.png")
        println("R plot saved as $plotname.R.png")

        P = plot(xtickfont=FONT,ytickfont=FONT,yguidefont=FONT,xguidefont=FONT,legendfont=FONT,size=SIZE)
        plot!(xlabel="Position (a.u.)",ylabel="|C0|²")
        for dyn in DYN_LIST
            eval(Meta.parse("""plot!(R_$(dyn),abs2.(C_$(dyn)[:,1]),label="$(dyn)")"""))
        end
        savefig(P,"$(plotname)RD-C0.png")
        println("c plot saved!")


        for dyn in DYN_LIST
            println("Starting with $dyn")
            if dyn == "SHEEP"
                eval(Meta.parse("Cgs_SHEEP = abs2.(C_SHEEP)"))
            else
                L = eval(Meta.parse("length(T_$dyn)"))
                eval(Meta.parse("Cgs_$(dyn) = zeros($L,2)"))
                if eval(Meta.parse("$(dyn)_sts")) == 2
                    fricSts = 1
                else
                    fricSts = 10
                end
                #eval(Meta.parse("count = 0"))
                for i in 1:L
                    #eval(Meta.parse("count++"))
                    eval(Meta.parse("Cgs_$(dyn)[$i,1] = sum(abs2.(C_$(dyn)[$i,1:$fricSts]))"))
                end
            end
            println("Finished $dyn")
        end

        println("Finished treatment, starting plot...")
        P = plot(xtickfont=FONT,ytickfont=FONT,yguidefont=FONT,xguidefont=FONT,legendfont=FONT,size=SIZE)
        plot!(xlabel="Position (a.u.)",ylabel="|C0|2")
        for dyn in DYN_LIST
            eval(Meta.parse("""plot!(R_$(dyn),Cgs_$(dyn)[:,1],label="$(dyn)")"""))
        end
        savefig(P,"$(plotname)RD-Cgs.png")
        println("cgs plot saved!")

        #=
        println("Starting friction analysis on EH trajectory...")
        tol = 0.1
        while tol != 0
            println("Enter tolerance. Set tol=0 to exit...")
            tol = Meta.parse(readline())
            @show tol
            P = isFricGood(T_EH,R_EH,P_EH,C_EH,"EH",tol)
            savefig(P,"$(plotname)Fric-$(tol).png")
            display(P)
        end
        # =#


    end #NDOFs == 1
end

function GS_SH_prob(K,Fast)
    GS_PROB=zeros(size(K))
    for k in 1:length(K)
        for traj in 1:Ntrajs
            if Fast[traj,k]==1
                GS_PROB[k]+=1
            end
        end
    end

return GS_PROB./Ntrajs
end


if GEN_METHOD == "K_SIMULATIONS"
    if length(R0) == 1 #one nuclear dimension
        ######## LOAD SAVE FILES
        for DYN in DYN_LIST
            A_string="A_$(DYN)=general_K_read(filename,$(DYN))"
            R_string="FR_$(DYN)=A_$(DYN)[1]"
            P_string="FP_$(DYN)=A_$(DYN)[2]"
            eval(Meta.parse(A_string))
            eval(Meta.parse(R_string))
            eval(Meta.parse(P_string))
            if DYN in MF_LIST || DYN in SH_LIST
                pop_string="Fpop_$(DYN)=A_$(DYN)[3]"
                eval(Meta.parse(pop_string))
                if DYN in SH_LIST
                    ast_string="Fast_$(DYN)=A_$(DYN)[4]"
                    eval(Meta.parse(ast_string))
                end
            end
        end

            ########### GROUND STATE PLOTTER. PLOTS FINAL GS POPULATIONS FOR DYNAMICS IN DYN_LIST
        GS_PLOT=plot(size=SIZE)
        for DYN in DYN_LIST
            if DYN in MF_LIST
                GS_string="GS_$(DYN)=GS_MF_prob(K,Fpop_$(DYN))"
                eval(Meta.parse(GS_string))
                plt_string="plot!(K,GS_$(DYN),label=$(DYN))"
                eval(Meta.parse(plt_string))
            elseif DYN in SH_LIST
                GS_string="GS_$(DYN)=GS_SH_prob(K,Fast_$(DYN))"
                eval(Meta.parse(GS_string))
                plt_string="plot!(K,GS_$(DYN),label=$(DYN))"
                eval(Meta.parse(plt_string))
            else
                println("$DYN doesn't have GS info, won't plot")
            end
        end

        savefig(GS_PLOT,"plots/"*potname*"_R0($R0)_GS_PLOT.png")

            ####### STATE TRANSMISSION PLOTS
        #cond(R,P) is a function that gives true if R and P are within some condition of interest for the plot, otherwise false
        cond(R,P)=R[1]>0

        P=plot(size=SIZE)
        for DYN in DYN_LIST
            dyn_sts=eval(Meta.parse("$(DYN)_sts"))
            if DYN in MF_LIST
                M_string="MULTI_$(DYN)=multi_MF_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                eval(Meta.parse(M_string))
                for st in 1:dyn_sts
                    plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                    eval(Meta.parse(plt_string))
                end
            elseif DYN in SH_LIST
                if SH_eval==true
                    M_string="MULTI_$(DYN)=multi_SH_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                    eval(Meta.parse(M_string))
                    for st in 1:dyn_sts
                        plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                        eval(Meta.parse(plt_string))
                    end
                else
                    M_string="MULTI_$(DYN)=multi_MF_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                    eval(Meta.parse(M_string))
                    for st in 1:dyn_sts
                        plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                        eval(Meta.parse(plt_string))
                    end
                end
            else
                println("$DYN doesn't have electronic states info, won't plot")
            end
        end

        savefig(P,"plots/"*potname*"_R0($R0)_STATE_TRANS.png")

            ####### STATE REFLECTION PLOTS
        #cond(R,P) is a function that gives true if R and P are within some condition of interest for the plot, otherwise false
        cond(R,P)=R[1]<0

        P=plot(size=SIZE)
        for DYN in DYN_LIST
            dyn_sts=eval(Meta.parse("$(DYN)_sts"))
            if DYN in MF_LIST
                M_string="MULTI_$(DYN)=multi_MF_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                eval(Meta.parse(M_string))
                for st in 1:dyn_sts
                    plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                    eval(Meta.parse(plt_string))
                end
            elseif DYN in SH_LIST
                if SH_eval==true
                    M_string="MULTI_$(DYN)=multi_SH_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                    eval(Meta.parse(M_string))
                    for st in 1:dyn_sts
                        plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                        eval(Meta.parse(plt_string))
                    end
                else
                    M_string="MULTI_$(DYN)=multi_MF_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                    eval(Meta.parse(M_string))
                    for st in 1:dyn_sts
                        plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                        eval(Meta.parse(plt_string))
                    end
                end
            else
                println("$DYN doesn't have electronic states info, won't plot")
            end
        end

        savefig(P,"plots/"*potname*"_R0($R0)_STATE_REFL.png")

    else
        println("Still no automatic plotting implementation for multidimensional K_SIMULATIONS")
    end

end

println("Finished plotting!")
