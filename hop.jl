function hop!(S::FSSH_state)
  probs=zeros(nsts)
  ast=S.ast
  akk=abs2(S.C[ast])
  for i in 1:nsts
      if i != ast
          prob=2*real(sum([abs(S.p[a])/mass*S.Γ[a][ast,i] for a in 1:S.NDOFs])*conj(S.C[ast])*S.C[i])*dt/akk
          probs[i]=maximum([0,prob])
      end
  end
  for i in 2:nsts
      probs[i]+=probs[i-1]
  end
  ξ=rand()
  for i in 1:nsts
      if ξ<probs[i] #hop successful, pass to energy check
          dij=[S.Γ[k][i,ast] for k in 1:S.NDOFs] #NAC vector where it will hop
          if S.NDOFs==1
              Nij=1
          else
              Nij=dij/abs(dij) #normalized NAC vector
          end
          p_nac=S.p.*Nij #momentum alongside NAC direction
          p_ad=S.p-p_nac #adiabatic momentum
          ΔE=abs2(p_nac)/2/mass+S.E[ast]-S.E[i] #check energy alongside NAC direction
          #Eini=abs2(S.p)/2/mass+S.E[ast] #for debugging purposes (1)
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              p_nac_new=sign.(p_nac).*Nij*sqrt(2*mass*ΔE)
              S.p=p_ad+p_nac_new
              S.ast=i
              #Efin=abs2(S.p)/2/mass+S.E[i] #for debugging (1)
              #@show Eini-Efin #for debugging (1)
          end
          break
      end
  end

end


function hop!(S::FSSH_dia_state)
  probs=zeros(nsts)
  ast=S.dst
  akk=abs2(S.C[ast])
  for i in 1:nsts
      if i != ast
          prob=2*imag(S.C[i]*conj(S.C[ast])*S.V[i,ast])*dt/akk
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
          ΔE=S.p^2/2/mass+S.V[ast,ast]-S.V[i,i]
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i")
              S.p=sign(S.p)*sqrt(2*mass*ΔE)
              S.dst=i
          end
          break
      end
  end

end


function hop!(S::CM2_FSSH_state)
  probs=zeros(2)
  ast=S.ast
  akk=abs2(S.D[ast])

  for i in 1:2
      if i != ast
          prob=-2*real(S.Dmat[ast,i]*conj(S.D[ast])*S.D[i])*dt/akk
          probs[i]=maximum([0,prob])
      end
  end
  for i in 2:2
      probs[i]+=probs[i-1]
  end
  ξ=rand()
  for i in 1:2
      if ξ<probs[i] #hop successful, pass to energy check
          if S.NDOFs!=1
              error("Still haven't implemented for many NDOFs")
          end
          E=[S.E[1],sum(S.E[2:end].*(S.tnorm.^2))]
          ΔE=S.p^2/2/mass+E[ast]-E[i] #check energy alongside NAC direction
          #Eini=abs2(S.p)/2/mass+S.E[ast] #for debugging purposes (1)
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              p_new=sign(S.p)*sqrt(2*mass*ΔE)
              S.p=p_new
              S.ast=i
              #Efin=abs2(S.p)/2/mass+S.E[i] #for debugging (1)
              #@show Eini-Efin #for debugging (1)
          end
          break
      end
  end

end


function hop!(S::CM3_FSSH_state)
  probs=zeros(3)
  ast=S.ast
  akk=abs2(S.D[ast])

  for i in 1:3
      if i != ast
          prob=-2*real(S.Dmat[ast,i]*conj(S.D[ast])*S.D[i])*dt/akk
          probs[i]=maximum([0,prob])
      end
  end
  for i in 2:3
      probs[i]+=probs[i-1]
  end
  ξ=rand()
  for i in 1:3
      if ξ<probs[i] #hop successful, pass to energy check
          if S.NDOFs!=1
              error("Still haven't implemented for many NDOFs")
          end
          Ebar=[S.E[2]*S.tnorm[k-1]^2+S.E[k]*S.tnorm[1]^2 for k in 2:nsts]
          E=[S.E[1],sum(S.E[2:end].*(S.tnorm.^2)),sum(Ebar[2:end].*(S.tnorm2.^2))]
          ΔE=S.p^2/2/mass+E[ast]-E[i] #check energy alongside NAC direction
          #Eini=abs2(S.p)/2/mass+S.E[ast] #for debugging purposes (1)
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              p_new=sign(S.p)*sqrt(2*mass*ΔE)
              S.p=p_new
              S.ast=i
              #Efin=abs2(S.p)/2/mass+S.E[i] #for debugging (1)
              #@show Eini-Efin #for debugging (1)
          end
          break
      end
  end

end

function hop!(S::CM2_FSSH_VANILLA_state)
  probs=zeros(2)
  ast=S.ast
  akk=abs2(S.D[ast])

  for i in 1:2
      if i != ast
          prob=-2*real(S.Dmat[ast,i]*conj(S.D[ast])*S.D[i])*dt/akk
          probs[i]=maximum([0,prob])
      end
  end
  for i in 2:2
      probs[i]+=probs[i-1]
  end
  ξ=rand()
  for i in 1:2
      if ξ<probs[i] #hop successful, pass to energy check
          if S.NDOFs!=1
              error("Still haven't implemented for many NDOFs")
          end
          #E=[S.E[1],sum(S.E[2:end].*(S.tnorm.^2))]
          E=[S.E[1],S.E[1]+S.w11]
          ΔE=S.p^2/2/mass+E[ast]-E[i] #check energy alongside NAC direction
           #Eini=abs2(S.p)/2/mass+S.E[ast] #for debugging purposes (1)
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              p_new=sign(S.p)*sqrt(2*mass*ΔE)
              S.p=p_new
              S.ast=i
              #Efin=abs2(S.p)/2/mass+S.E[i] #for debugging (1)
              #@show Eini-Efin #for debugging (1)
          end
          break
      end
  end

end


function hop!(S::CM3_FSSH_VANILLA_state)
  probs=zeros(3)
  ast=S.ast
  akk=abs2(S.D[ast])

  for i in 1:3
      if i != ast
          prob=-2*real(S.Dmat[ast,i]*conj(S.D[ast])*S.D[i])*dt/akk
          probs[i]=maximum([0,prob])
      end
  end
  for i in 2:3
      probs[i]+=probs[i-1]
  end
  ξ=rand()
  for i in 1:3
      if ξ<probs[i] #hop successful, pass to energy check
          if S.NDOFs!=1
              error("Still haven't implemented for many NDOFs")
          end
          #Ebar=[S.E[2]*S.tnorm[k-1]^2+S.E[k]*S.tnorm[1]^2 for k in 2:nsts]
          #E=[S.E[1],sum(S.E[2:end].*(S.tnorm.^2)),sum(Ebar[2:end].*(S.tnorm2.^2))]
          E=[S.E[1],S.E[1]+S.w11,S.E[1]+S.w22]
          ΔE=S.p^2/2/mass+E[ast]-E[i] #check energy alongside NAC direction
          #Eini=abs2(S.p)/2/mass+S.E[ast] #for debugging purposes (1)
          if ΔE>=0 #energy suficient
              #println("hop from $ast to $i") #for debugging
              p_new=sign(S.p)*sqrt(2*mass*ΔE)
              S.p=p_new
              S.ast=i
              #Efin=abs2(S.p)/2/mass+S.E[i] #for debugging (1)
              #@show Eini-Efin #for debugging (1)
          end
          break
      end
  end

end
