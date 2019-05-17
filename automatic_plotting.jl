using Plots
using Plots.PlotMeasures
eval(Meta.parse(plot_method))
GUIDEFONT=font(24,"Helvetica")
TICKFONT=font(24,"Helvetica")
L_MARG=[20mm 0mm]
B_MARG=[15mm 0mm]
SIZE=[1000,700]

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
            plot_str = """plot(Rbase_$(dyn_1)[:,1],HR_$(dyn_1)[:,1],label="$(dyn_1)_i",line=(1,:dash),color=:black,left_margin=L_MARG,bottom_margin=B_MARG,size=SIZE)"""
        else
            plot_str = """plot(left_margin=L_MARG,bottom_margin=B_MARG,size=SIZE)"""
        end
        P=eval(Meta.parse(plot_str))
        for dyn in DYN_LIST
            plot_str = """plot!(Rbase_$(dyn)[:,end],HR_$(dyn)[:,end],label="$(dyn)",line=(2,:solid))"""
            eval(Meta.parse(plot_str))
        end
        plot!(xlabel="Position (a.u.)",xlims=XLIMS)#,ylabel="Nuclear distribution");
        plot!(xguidefont = GUIDEFONT,xtickfont=TICKFONT,yguidefont = GUIDEFONT,ytickfont=TICKFONT,legend=:best);
        #ylims!(hmin-0.1,hmax+0.1)
        savefig(P,"plots/"*potname*"R0($R0)_p0($p0).png")
        println("figure saved!")
    elseif NDOFs==2
        HISTO_RES=100;
        mins=[-50;-50]
        maxs=[50;50]
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
            str="Rbase_$dyn,HR_$dyn=multi_d_histo(R_$dyn,$mins,$maxs,$HISTO_RES)"
            eval(Meta.parse(str))
            str="HR_$dyn=HR_$dyn./maximum(HR_$dyn[:,:,1])"
            eval(Meta.parse(str))
            println("Finished $dyn histogram after $(time()-tic) seconds!")
        end

        for dyn in DYN_LIST
            plt_str="contour(Rbase_$(dyn)[1][:,end],Rbase_$(dyn)[2][:,end],HR_$(dyn)[:,:,end])"
            eval(Meta.parse(plt_str))
            xticks!([(2k+1)*pi for k in -7:7])
            yticks!([(2k+1)*pi for k in -7:7])
            savefig("$(filename)_$(dyn)")
        end
    else #NDOFs != 1
        println("Still no plot implementation for multi_dimensional case")
    end
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
        GS_PLOT=plot(left_margin=L_MARG,bottom_margin=B_MARG,size=SIZE)
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

        P=plot(left_margin=L_MARG,bottom_margin=B_MARG,size=SIZE)
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

        P=plot(left_margin=L_MARG,bottom_margin=B_MARG,size=SIZE)
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
