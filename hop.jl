function hop!(S::FSSH_state,tstep)
  probs=zeros(nsts)
  ast=S.ast
  akk=abs2(S.el.C[ast])
  for i in 1:nsts
      if i != ast
          #prob=2*real(sum([abs(S.cl.p[a])/mass*S.el.Γ[a][ast,i] for a in 1:S.cl.NDOFs])*conj(S.el.C[ast])*S.el.C[i])*tstep/akk
          prob=2*real(S.ODE.Cdot[i,ast]*conj(S.el.C[ast])*S.el.C[i])*tstep/akk
          probs[i]=maximum([0,prob])
      end
  end
  for i in 2:nsts
      probs[i]+=probs[i-1]
  end
  ξ=rand()
  for i in 1:nsts
      if ξ<probs[i] #hop successful, pass to energy check
          dij=[S.el.Γ[k][i,ast] for k in 1:S.cl.NDOFs] #NAC vector where it will hop
          if S.cl.NDOFs==1
              Nij=1
          else
              Nij=abs.(dij)/norm(dij) #normalized NAC vector
          end
          p_nac=S.cl.p.*Nij #momentum alongside NAC direction
          p_ad=S.cl.p-p_nac #adiabatic momentum
          ΔE=norm2(p_nac)/2/mass+S.el.E[ast]-S.el.E[i] #check energy alongside NAC direction
          #Eini=abs2(S.p)/2/mass+S.E[ast] #for debugging purposes (1)
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              p_nac_new=sign.(p_nac).*Nij*sqrt(2*mass*ΔE)
              pnew=p_ad+p_nac_new
              S=FSSH_state_builder(S.cl.R,pnew,S.el.C,i,S.el.Ua,S.cl.NDOFs)
              return S
              #Efin=abs2(S.p)/2/mass+S.E[i] #for debugging (1)
              #@show Eini-Efin #for debugging (1)
          end
          break
      end
  end

  return S
end


function hop!(S::FSSH_dia_state,tstep)
  probs=zeros(nsts)
  ast=S.ast
  akk=abs2(S.el.C[ast])
  for i in 1:nsts
      if i != ast
          prob=2*imag(S.el.C[i]*conj(S.el.C[ast])*S.el.E[i,ast])*tstep/akk
          probs[i]=maximum([0,prob])
      end
  end
  for i in 2:nsts
      probs[i]+=probs[i-1]
      #global PPP=probs[end]
  end
  ξ=rand()
  for i in 1:nsts
      if ξ<probs[i] #hop successful, pass to energy check
          if S.cl.NDOFs!=1
              error("Haven't implemented multiple nuclear degrees of freedom for diabatic FSSH")
          end
          ΔE=S.p^2/2/mass+S.el.E[ast,ast]-S.el.E[i,i]
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i")
              pnew=sign(S.cl.p)*sqrt(2*mass*ΔE)
              S=FSSH_dia_state_builder(S.cl.R,pnew,S.el.C,i,S.cl.NDOFs)
              return S
          end
          break
      end
  end
  return S
end


function hop!(S::CM2_FSSH_state,tstep)
  probs=zeros(2)
  ast=S.ast
  akk=abs2(S.el.C[ast])

  for i in 1:2
      if i != ast
          prob=-2*real(S.ODE.Cdot[ast,i]*conj(S.el.C[ast])*S.el.C[i])*tstep/akk
          probs[i]=maximum([0,prob])
      end
  end
  probs[2]+=probs[1]

  ξ=rand()
  for i in 1:2
      if ξ<probs[i] #hop successful, pass to energy check
          if S.cl.NDOFs!=1
              U=zeros(nsts,nsts)

              U[1,1]=1
              for i in 1:nsts-1
                  U[2,i+1]=S.CM2.tnorm[i]
                  U[i+1,2]=U[2,i+1]
              end
              for i in 3:nsts
                  U[i,i]=-S.CM2.tnorm[1]
              end

              Udag=U'

              dij=[(Udag*S.el.Γ[dof]*U)[i,ast] for dof in 1:S.cl.NDOFs]
              Nij=abs.(dij)/norm(dij)
              p_nac=S.cl.p.*Nij
          else
              p_nac=S.cl.p
              Nij=1
          end
          p_ad=S.cl.p-p_nac
          #E=[S.el.E[1],sum(S.el.E[2:end].*(S.CM2.tnorm.^2))]
          E=[S.el.E[1],sum(S.el.W[2:end,1].*(S.CM2.tnorm.^2))+S.el.E[1]]
          ΔE=norm2(p_nac)/2/mass+E[ast]-E[i] #check energy alongside NAC direction
          #Eini=abs2(S.p)/2/mass+S.E[ast] #for debugging purposes (1)
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              p_nac_new=sign.(p_nac).*Nij*sqrt(2*mass*ΔE)
              pnew=p_nac_new+p_ad
              S=CM2_FSSH_state_builder(S.cl.R,pnew,S.el.C,i,S.el.Ua,S.cl.NDOFs)
              return S
              #Efin=abs2(S.p)/2/mass+S.E[i] #for debugging (1)
              #@show Eini-Efin #for debugging (1)
          end
          break
      end
  end

  return S
end


function hop!(S::CM3_FSSH_state,tstep)
  probs=zeros(3)
  ast=S.ast
  akk=abs2(S.el.C[ast])

  for i in 1:3
      if i != ast
          prob=-2*real(S.ODE.Cdot[ast,i]*conj(S.el.C[ast])*S.el.C[i])*tstep/akk
          probs[i]=maximum([0,prob])
      end
  end
  for i in 2:3
      probs[i]+=probs[i-1]
  end
  ξ=rand()
  for i in 1:3
      if ξ<probs[i] #hop successful, pass to energy check
          if S.cl.NDOFs!=1
              U2=zeros(nsts,nsts)
              U3=zeros(nsts,nsts)

              U2[1,1]=1
              for i in 1:nsts-1
                  U2[2,i+1]=S.CM3.tnorm[i]
                  U2[i+1,2]=U2[2,i+1]
              end
              for i in 3:nsts
                  U2[i,i]=-S.CM3.tnorm[1]
              end

              U3[1,1]=S.CM3.zbar
              U3[2,2]=1
              for i in 1:nsts-2
                  U3[3,i+2]=S.CM3.tnorm2[i]
                  U3[i+2,3]=U3[3,i+2]
              end
              for i in 4:nsts
                  U3[i,i]=-S.CM3.tnorm2[1]
              end

              U=U2*U3
              Udag=U'

              dij=[(Udag*S.el.Γ[dof]*U)[i,ast] for dof in 1:S.cl.NDOFs]
              Nij=abs.(dij)/norm(dij)
              p_nac=S.cl.p.*Nij
          else
              p_nac=S.cl.p
              Nij=1
          end
          p_ad=S.cl.p-p_nac
          Ebar=[S.el.E[2]*S.CM3.tnorm[k-1]^2+S.el.E[k]*S.CM3.tnorm[1]^2 for k in 2:nsts]
          E=[S.el.E[1],sum(S.el.E[2:end].*(S.CM3.tnorm.^2)),sum(Ebar[2:end].*(S.CM3.tnorm2.^2))]
          ΔE=norm2(p_nac)/2/mass+E[ast]-E[i] #check energy alongside NAC direction
          #Eini=abs2(S.p)/2/mass+S.E[ast] #for debugging purposes (1)
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              p_nac_new=sign.(p_nac).*Nij*sqrt(2*mass*ΔE)
              pnew=p_nac_new+p_ad
              S=CM3_FSSH_state_builder(S.cl.R,pnew,S.el.C,i,S.el.Ua,S.cl.NDOFs)
              return S
              #Efin=abs2(S.p)/2/mass+S.E[i] #for debugging (1)
              #@show Eini-Efin #for debugging (1)
          #else
            #  println("frustrated hop from $ast to $i")
          end
          break
      end
  end

  return S
end

function hop!(S::CM2_FSSH_FRIC_state,tstep)
  probs=zeros(2)
  ast=S.ast
  akk=abs2(S.el.C[ast])

  for i in 1:2
      if i != ast
          prob=-2*real(S.ODE.Cdot[ast,i]*conj(S.el.C[ast])*S.el.C[i])*tstep/akk
          probs[i]=maximum([0,prob])
      end
  end
  probs[2]+=probs[1]

  ξ=rand()
  for i in 1:2
      if ξ<probs[i] #hop successful, pass to energy check
          if S.cl.NDOFs!=1
              U=zeros(nsts,nsts)

              U[1,1]=1
              for i in 1:nsts-1
                  U[2,i+1]=S.CM2.tnorm[i]
                  U[i+1,2]=U[2,i+1]
              end
              for i in 3:nsts
                  U[i,i]=-S.CM2.tnorm[1]
              end

              dij=[((U')*S.el.Γ[dof]*U)[i,ast] for dof in 1:S.cl.NDOFs]
              Nij=abs.(dij)/norm(dij)
              p_nac=S.cl.p.*Nij
          else
              p_nac=S.cl.p
              Nij=1
          end
          p_ad=S.cl.p-p_nac
          if ast == 1
            ΔE=real(norm2(p_nac)/2/mass-sum(S.el.E[2:end] .*(S.CM2.tnorm.^2)))
          else
            ΔE=real(S.extra[2]-S.el.E[1]-S.cl.mem[end])
          end
          #@show ΔE
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              if ast == 1
                p_nac_new=sign.(p_nac) .*Nij*sqrt(2*mass*ΔE)
                pnew=p_nac_new+p_ad
                stMax=2
                tauMax=abs(S.CM2.tnorm[1])
                for st in 2:nsts-1
                  tauTry=abs(S.CM2.tnorm[st])
                  if tauTry>tauMax
                    tauMax=tauTry
                    stMax=i+1
                  end
                end
                extra=Any[stMax,S.extra[2]]
                energy(S)
                S=CM2_FSSH_FRIC_state_builder(S.cl.R,pnew,S.el.C,i,S.el.Ua,S.cl.NDOFs,S.cl.mem,extra)
                energy(S)
                return S
              else
                pnew=p_ad+sign.(S.cl.p) .* sqrt(2*mass*ΔE)
                energy(S)
                extra=Any[1,S.extra[2]]
                S=CM2_FSSH_FRIC_state_builder(S.cl.R,pnew,S.el.C,i,S.el.Ua,S.cl.NDOFs,S.cl.mem,extra)
                energy(S)
                return S
              end
          end
          break
      end
  end

  return S
end


function hop!(S::CM3_FSSH_FRIC_state,tstep)
  probs=zeros(3)
  ast=S.ast
  akk=abs2(S.el.C[ast])

  for i in 1:3
      if i != ast
          prob=-2*real(S.ODE.Cdot[ast,i]*conj(S.el.C[ast])*S.el.C[i])*tstep/akk
          probs[i]=maximum([0,prob])
      end
  end
  for i in 2:3
      probs[i]+=probs[i-1]
  end
  ξ=rand()
  for i in 1:3
      if ξ<probs[i] #hop successful, pass to energy check
          if S.cl.NDOFs!=1
              U2=zeros(nsts,nsts)
              U3=zeros(nsts,nsts)

              U2[1,1]=1
              for i in 1:nsts-1
                  U2[2,i+1]=S.CM3.tnorm[i]
                  U2[i+1,2]=U2[2,i+1]
              end
              for i in 3:nsts
                  U2[i,i]=-S.CM3.tnorm[1]
              end

              U3[1,1]=S.CM3.zbar
              U3[2,2]=1
              for i in 1:nsts-2
                  U3[3,i+2]=S.CM3.tnorm2[i]
                  U3[i+2,3]=U3[3,i+2]
              end
              for i in 4:nsts
                  U3[i,i]=-S.CM3.tnorm2[1]
              end

              U=U2*U3
              Udag=U'

              dij=[(Udag*S.el.Γ[dof]*U)[i,ast] for dof in 1:S.cl.NDOFs]
              Nij=abs.(dij)/norm(dij)
              p_nac=S.cl.p.*Nij
          else
              p_nac=S.cl.p
              Nij=1
          end
          p_ad=S.cl.p-p_nac
          Ebar=[S.el.E[2]*S.CM3.tnorm[k-1]^2+S.el.E[k]*S.CM3.tnorm[1]^2 for k in 2:nsts]
          E=[S.el.E[1],sum(S.el.E[2:end].*(S.CM3.tnorm.^2)),sum(Ebar[2:end].*(S.CM3.tnorm2.^2))]
          ΔE=norm2(p_nac)/2/mass+E[ast]-E[i] #check energy alongside NAC direction
          #Eini=abs2(S.p)/2/mass+S.E[ast] #for debugging purposes (1)
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              p_nac_new=sign.(p_nac).*Nij*sqrt(2*mass*ΔE)
              pnew=p_nac_new+p_ad
              S=CM3_FSSH_state_builder(S.cl.R,pnew,S.el.C,i,S.el.Ua,S.cl.NDOFs)
              return S
              #Efin=abs2(S.p)/2/mass+S.E[i] #for debugging (1)
              #@show Eini-Efin #for debugging (1)
          #else
            #  println("frustrated hop from $ast to $i")
          end
          break
      end
  end

  return S
end

function hop!(S::SHEEP_state,tstep)
  num_tr_sts=length(SHEEP_REL)
  probs=zeros(nsts)
  ast_array=SHEEP_REL[S.ast]
  rho_kk=sum(abs2.(S.el.C[ast_array]))
  for tr_st in 1:num_tr_sts
      for ast in ast_array
          if tr_st != ast
              #prob=2*real(sum([abs(S.cl.p[a])/mass*S.el.Γ[a][ast,tr_st] for a in 1:S.cl.NDOFs])*conj(S.el.C[ast])*S.el.C[tr_st])*tstep/rho_kk
              prob=2*real(S.ODE.Cdot[tr_st,ast]*conj(S.el.C[ast])*S.el.C[tr_st])*tstep/rho_kk
              probs[tr_st]+=maximum([0,prob])
          end
      end
  end
  for tr_st in 2:num_tr_sts
      probs[tr_st]+=probs[tr_st-1]
  end
  ξ=rand()
  for tr_st in 1:num_tr_sts
      if ξ<probs[tr_st] #hop successful, pass to energy check
          dij=zeros(S.cl.NDOFs)
          for dof in 1:S.cl.NDOFs
              dij[dof]=sum(S.el.Γ[dof][tr_st,S.ast])
          end
          if S.cl.NDOFs==1
              Nij=1
          else
              Nij=abs.(dij)/norm(dij) #normalized NAC vector
          end
          p_nac=S.cl.p.*Nij #momentum alongside NAC direction
          p_ad=S.cl.p-p_nac #adiabatic momentum
          East=0.0
          Eobj=0.0
          rho_obj=0.0
          for st in SHEEP_REL[tr_st]
              c2=abs2(S.el.C[st])
              rho_obj+=c2
              Eobj+=S.el.E[st]*c2
          end
          Eobj=Eobj/rho_obj
          for ast in ast_array
              c2=abs2(S.el.C[ast])
              East+=c2*S.el.E[ast]
          end
          East=East/rho_kk
          ΔE=norm2(p_nac)/2/mass+East-Eobj #check energy alongside NAC direction
          #Eini=abs2(S.p)/2/mass+S.E[ast] #for debugging purposes (1)
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              p_nac_new=sign.(p_nac).*Nij*sqrt(2*mass*ΔE)
              pnew=p_ad+p_nac_new
              S=SHEEP_state_builder(S.cl.R,pnew,S.el.C,tr_st,S.el.Ua,S.cl.NDOFs)
              return S
              #Efin=abs2(S.p)/2/mass+S.E[i] #for debugging (1)
              #@show Eini-Efin #for debugging (1)
          end
          break
      end #if prob is enough for hop
  end

  return S
end
