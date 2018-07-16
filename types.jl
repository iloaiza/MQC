struct BO_state
    R
    p
    E
    F
    NDOFs
end

struct EH_state
    R #nuclear coordinates
    p #nuclear momenta
    C #electronic coefficients
    E #adiabatic energies
    Γ #NACs. Array containing arrays of Γa, Γa being NAC along a nuclear coord.
    F #Forces. Array containing arrays of Fa, Fa being force grad. alongside a
    Ua #unitary transformation from dia to ad basis
    NDOFs #Number of nuclear degrees of freedom
end

type FSSH_state
  R
  p
  C
  E
  Γ
  F
  Ua
  ast #current adiabatic state
  NDOFs
end

type FSSH_dia_state
  R
  p
  C
  V #potential matrix (has off diagonal terms)
  dV
  dst
  NDOFs
end

struct CM2_VANILLA_state
    R
    p
    D
    E
    Ua
    F2 #transformed force in collective space 2x2 matrix
    z
    tnorm #normalized couplings vector
    wvec #epsilons vector
    Dmat #2x2 matrix for electronic evolution, ddot=Dmat*d
    NDOFs
end

struct CM3_VANILLA_state
    R
    p
    D
    E
    Ua
    F3 #transformed force in collective space 3x3 matrix
    z
    zbar
    tnorm #normalized couplings vector
    tnorm2
    wvec #epsilons vector
    wvec2
    Dmat #3x3 matrix for electronic evolution, ddot=Ddot*d
    NDOFs
end

type CM2_FSSH_VANILLA_state
    R
    p
    D
    E
    Ua
    Fgs #ground state force
    z
    tnorm #normalized couplings vector
    w11 #epsilons vector
    Dmat #2x2 matrix for electronic evolution, ddot=Dmat*d
    ast
    NDOFs
end

type CM3_FSSH_VANILLA_state
    R
    p
    D
    E
    Ua
    Fgs #ground state force
    z
    zbar
    tnorm #normalized couplings vector
    tnorm2
    w11 #epsilons vector
    w22
    Dmat #3x3 matrix for electronic evolution, ddot=Ddot*d
    ast
    NDOFs
end

struct CM2_state
    R
    p
    D
    E
    Ua
    F2 #transformed force in collective space 2x2 matrix
    z
    tnorm #normalized couplings vector
    wvec #epsilons vector
    Dmat #2x2 matrix for electronic evolution, ddot=Dmat*d
    NDOFs
end

struct CM3_state
    R
    p
    D
    E
    Ua
    F3 #transformed force in collective space 3x3 matrix
    z
    zbar
    tnorm #normalized couplings vector
    tnorm2
    wvec #epsilons vector
    wvec2
    Dmat #3x3 matrix for electronic evolution, ddot=Ddot*d
    NDOFs
end

type CM2_FSSH_state
    R
    p
    D
    E
    Ua
    F2 #transformed force in collective space 2x2 matrix
    z
    tnorm #normalized couplings vector
    wvec #epsilons vector
    Dmat #2x2 matrix for electronic evolution, ddot=Dmat*d
    NDOFs
    ast
end

type CM3_FSSH_state
    R
    p
    D
    E
    Ua
    F3 #transformed force in collective space 3x3 matrix
    z
    zbar
    tnorm #normalized couplings vector
    tnorm2
    wvec #epsilons vector
    wvec2
    Dmat #3x3 matrix for electronic evolution, ddot=Ddot*d
    NDOFs
    ast
end
