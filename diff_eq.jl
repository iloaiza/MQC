function diff_eq(S::CL_state)
  if S.cl.NDOFs==1
    return S.ODE.Rdot[1],S.ODE.pdot[1]
  else
    return S.ODE.Rdot,S.ODE.pdot
  end
end

function diff_eq(S::MF_state)
  if S.cl.NDOFs==1
    return S.ODE.Rdot[1],S.ODE.pdot[1],S.ODE.Cdot
  else
    return S.ODE.Rdot,S.ODE.pdot,S.ODE.Cdot
  end
end

function diff_eq(S::SH_state)
  if S.cl.NDOFs==1
    return S.ODE.Rdot[1],S.ODE.pdot[1],S.ODE.Cdot
  else
    return S.ODE.Rdot,S.ODE.pdot,S.ODE.Cdot
  end
end
