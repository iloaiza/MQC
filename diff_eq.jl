function diff_eq(S::CL_state)
  if S.cl.NDOFs==1
    return [S.ODE.Rdot[1],S.ODE.pdot[1],0,S.ODE.memdot]
  else
    return [S.ODE.Rdot,S.ODE.pdot,0,S.ODE.memdot]
  end
end

function diff_eq(S::MF_state)
  if S.cl.NDOFs==1
    return [S.ODE.Rdot[1],S.ODE.pdot[1],S.ODE.Cdot*S.el.C,S.ODE.memdot]
  else
    return [S.ODE.Rdot,S.ODE.pdot,S.ODE.Cdot*S.el.C,S.ODE.memdot]
  end
end

function diff_eq(S::SH_state)
  if S.cl.NDOFs==1
    return [S.ODE.Rdot[1],S.ODE.pdot[1],S.ODE.Cdot*S.el.C,S.ODE.memdot]
  else
    return [S.ODE.Rdot,S.ODE.pdot,S.ODE.Cdot*S.el.C,S.ODE.memdot]
  end
end
