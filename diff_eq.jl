function diff_eq(S::BO_state)
    if S.NDOFs==1
        pdot=-S.F[1]
    end

    Rdot=S.p/mass

    return Rdot,pdot
end

function diff_eq(S::EH_state)
  if S.NDOFs==1
      pdot=-real(S.C'*S.F[1]*S.C)[1]
  else
      pdot=zeros(Float64,S.NDOFs)
      for i in 1:S.NDOFs
          pdot[i]=-real(S.C'*S.F[i]*S.C)[1]
      end
  end
  Rdot=S.p/mass

  Cdot=-(1im*diagm(S.E)+sum(Rdot.*S.Γ))*S.C

  return Rdot,pdot,Cdot
end

function diff_eq(S::FSSH_state)
  if S.NDOFs==1
        pdot=-S.F[1][S.ast,S.ast]
  else
      pdot=zeros(Float64,S.NDOFs)
      for i in 1:S.NDOFs
          pdot[i]=-S.F[i][S.ast,S.ast]
      end
  end
  Rdot=S.p/mass
  Cdot=-(1im*diagm(S.E)+sum(Rdot.*S.Γ))*S.C

  return Rdot,pdot,Cdot
end

function diff_eq(S::FSSH_dia_state)
  Rdot=S.p/mass
  Cdot=[-1im*sum(S.V[k,j]*S.C[j] for j in 1:nsts) for k in 1:nsts]

  if S.NDOFs==1
      pdot=-S.dV[1][S.dst,S.dst]
  else
      pdot=[-S.dV[i][S.dst,S.dst] for i in 1:S.NDOFs]
  end

  return Rdot,pdot,Cdot
end

#=
function diff_eq(S::CM3_state)
    Rdot=S.p/mass
    if S.z==0
        tnorm=zeros(S.tvec)
    else
        tnorm=S.tvec./S.z
    end
    if S.zbar==0
        tnorm2=zeros(S.tvec2)
    else
        tnorm2=S.tvec2./S.zbar
    end

    barW=sum(S.wvec.*(tnorm.^2))
    barW2=sum(S.wvec2.*(tnorm2.^2))

    if S.NDOFs==1
        F12=sum(S.fvec[1].*tnorm)
        F13=sum(S.fvec2[1].*tnorm2)
        force_matrix=-Float64[S.F[1][1,1] F12 F13;F12 S.F[1][1,1] 0;F13 0 S.F[1][1,1]]
        pdot=real(S.D'*force_matrix*S.D)[1]
    else
        error("Still haven't implemented many NDOFs for this!!")
    end

    Dmat=-Complex[0 -S.z 0;S.z 1im*barW -1im*S.zbar;0 -1im*S.zbar 1im*barW2]
    Ddot=Dmat*S.D

    return Rdot,pdot,Ddot
end
=#

function diff_eq(S::CM2_VANILLA_state)
    Rdot=S.p/mass

    if S.NDOFs==1
        pdot=-real(S.D'*S.F2[1]*S.D)[1]
    else
        for k in 1:S.NDOFs
            pdot[k]=-real(S.D'*S.F2[k]*S.D)[1]
        end
    end
    Ddot=S.Dmat*S.D

    return Rdot,pdot,Ddot
end

function diff_eq(S::CM3_VANILLA_state)
    Rdot=S.p/mass

    if S.NDOFs==1
        pdot=-real(S.D'*S.F3[1]*S.D)[1]
    else
        for k in 1:S.NDOFs
            pdot[k]=-real(S.D'*S.F3[k]*S.D)[1]
        end
    end
    Ddot=S.Dmat*S.D

    return Rdot,pdot,Ddot
end

function diff_eq(S::CM2_FSSH_VANILLA_state)
    Rdot=S.p/mass

    if S.NDOFs==1
        pdot=-S.Fgs[1]
    else
        for k in 1:S.NDOFs
            pdot[k]=-S.Fgs[k]
        end
    end
    Ddot=S.Dmat*S.D

    return Rdot,pdot,Ddot
end

function diff_eq(S::CM3_FSSH_VANILLA_state)
    Rdot=S.p/mass

    if S.NDOFs==1
        pdot=-S.Fgs[1]
    else
        for k in 1:S.NDOFs
            pdot[k]=-S.Fgs[k]
        end
    end
    Ddot=S.Dmat*S.D

    return Rdot,pdot,Ddot
end

function diff_eq(S::CM2_state)
    Rdot=S.p/mass

    if S.NDOFs==1
        pdot=-real(S.D'*S.F2[1]*S.D)[1]
    else
        for k in 1:S.NDOFs
            pdot[k]=-real(S.D'*S.F2[k]*S.D)[1]
        end
    end
    Ddot=S.Dmat*S.D

    return Rdot,pdot,Ddot
end

function diff_eq(S::CM3_state)
    Rdot=S.p/mass

    if S.NDOFs==1
        pdot=-real(S.D'*S.F3[1]*S.D)[1]
    else
        for k in 1:S.NDOFs
            pdot[k]=-real(S.D'*S.F3[k]*S.D)[1]
        end
    end
    Ddot=S.Dmat*S.D

    return Rdot,pdot,Ddot
end

function diff_eq(S::CM2_FSSH_state)
    Rdot=S.p/mass

    if S.NDOFs==1
        pdot=-S.F2[1][S.ast]
    else
        for k in 1:S.NDOFs
            pdot[k]=-S.F2[k][ast]
        end
    end

    Ddot=S.Dmat*S.D

    return Rdot,pdot,Ddot
end

function diff_eq(S::CM3_FSSH_state)
    Rdot=S.p/mass

    if S.NDOFs==1
        pdot=-S.F3[1][S.ast]
    else
        for k in 1:S.NDOFs
            pdot[k]=-S.F3[k][S.ast]
        end
    end

    Ddot=S.Dmat*S.D

    return Rdot,pdot,Ddot
end
