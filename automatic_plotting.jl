using Plots
eval(Meta.parse(plot_method))
GUIDEFONT=font(24,"Helvetica");
TICKFONT=font(24,"Helvetica");

if GEN_METHOD == "DYNAMICS"
    if length(R0) ==1 #1 nuclear dimension
        E0=h5read(file,"E0");                                     #LOAD INITIAL ENERGY

        for dyn in DYN_LIST
            println("Loading $dyn...")
            if dyn in CL_LIST
                str="""T_$(dyn),R_$(dyn),P_$(dyn)=CL_read("$file","$dyn")"""
            elseif dyn in MF_LIST
                str="""T_$(dyn),R_$(dyn),P_$(dyn),C_$(dyn)=MF_read("$file","$dyn")"""
            elseif dyn in SH_LIST
                str="""T_$(dyn),R_$dyn,P_$(dyn),C_$(dyn),AST_$(dyn)=SH_read("$file","$dyn")"""
            end
            eval(Meta.parse(str))
        end

        Rmin = R0 + 10
        Rmax = R0 - 10
        for dyn in DYN_LIST
            min_str = "minimum(R_$(dyn))"
            max_str = "maximum(R_$(dyn))"
            min = eval(Meta.parse(min_str))
            max = eval(Meta.parse(max_str))
            if Rmin > min
                Rmin = min
            end
            if Rmax < max
                Rmax = max
            end
        end

        xmin=Rmin-0.05*abs(Rmin)
        xmax=Rmax+0.05*abs(Rmax)

        hmax = 1
        hmin = -1
        inv = 1
        for dyn in DYN_LIST
            println("Building histograms for $dyn")
            str="Rbase_$dyn,HR_$dyn=super_histo(R_$dyn[:,1,:],$xmin,$xmax,$HISTO_RES)"
            eval(Meta.parse(str))
            str="HR_$dyn=$inv.*HR_$dyn./maximum(HR_$dyn[:,1])"
            eval(Meta.parse(str))
            inv = inv*-1

            min_str = "minimum(HR_$(dyn))"
            max_str = "maximum(HR_$(dyn))"
            min = eval(Meta.parse(min_str))
            max = eval(Meta.parse(max_str))
            if hmin > min
                hmin = min
            end
            if hmax < max
                hmax = max
            end
        end

        plot_str = """plot(Rbase_$(dyn)[:,1],HR_$(dyn)[:,1],label="$(dyn)_i",line=(1,:dash),color=:black)"""
        P=eval(Meta.parse(plot_str))
        for dyn in DYN_LIST
            plot_str = """plot!(Rbase_$(dyn)[:,end],HR_$(dyn)[:,end],label="$(dyn)",line=(2,:solid))"""
            eval(Meta.parse(plot_str))
        end
        plot!(xlabel="Position (a.u.)",ylabel="Nuclear distribution",xlims=(xmin,xmax));
        plot!(xguidefont = GUIDEFONT,xtickfont=TICKFONT,yguidefont = GUIDEFONT,ytickfont=TICKFONT);
        ylims!(hmin-0.1,hmax+0.1)
        savefig(P,filename*".png")
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
        GS_PLOT=plot()
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

        savefig(GS_PLOT,filename*"_GS_PLOT.png")

            ####### STATE TRANSMISSION PLOTS
        #cond(R,P) is a function that gives true if R and P are within some condition of interest for the plot, otherwise false
        cond(R,P)=R[1]>0

        P=plot()
        for DYN in DYN_LIST
            dyn_sts=eval(Meta.parse("$(DYN)_sts"))
            if DYN in MF_LIST
                M_string="MULTI_$(DYN)=multi_MF_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                eval(Meta.parse(M_string))
                for st in 1:nsts
                    plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                    eval(Meta.parse(plt_string))
                end
            elseif DYN in SH_LIST
                if SH_eval==true
                    M_string="MULTI_$(DYN)=multi_SH_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                    eval(Meta.parse(M_string))
                    for st in 1:nsts
                        plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                        eval(Meta.parse(plt_string))
                    end
                else
                    M_string="MULTI_$(DYN)=multi_MF_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                    eval(Meta.parse(M_string))
                    for st in 1:nsts
                        plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                        eval(Meta.parse(plt_string))
                    end
                end
            else
                println("$DYN doesn't have electronic states info, won't plot")
            end
        end

        savefig(P,filename*"_STATE_TRANS.png")

            ####### STATE REFLECTION PLOTS
        #cond(R,P) is a function that gives true if R and P are within some condition of interest for the plot, otherwise false
        cond(R,P)=R[1]<0

        P=plot()
        for DYN in DYN_LIST
            dyn_sts=eval(Meta.parse("$(DYN)_sts"))
            if DYN in MF_LIST
                M_string="MULTI_$(DYN)=multi_MF_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                eval(Meta.parse(M_string))
                for st in 1:nsts
                    plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                    eval(Meta.parse(plt_string))
                end
            elseif DYN in SH_LIST
                if SH_eval==true
                    M_string="MULTI_$(DYN)=multi_SH_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                    eval(Meta.parse(M_string))
                    for st in 1:nsts
                        plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                        eval(Meta.parse(plt_string))
                    end
                else
                    M_string="MULTI_$(DYN)=multi_MF_prob(K,A_$(DYN),1:$dyn_sts,cond)"
                    eval(Meta.parse(M_string))
                    for st in 1:nsts
                        plt_string="""plot!(K,MULTI_$(DYN)[:,$st],label="$(DYN)_st=$st")"""
                        eval(Meta.parse(plt_string))
                    end
                end
            else
                println("$DYN doesn't have electronic states info, won't plot")
            end
        end

        savefig(P,filename*"_STATE_REFL.png")

    else
        println("Still no automatic plotting implementation for multidimensional K_SIMULATIONS")
    end

end

println("Finished plotting!")
