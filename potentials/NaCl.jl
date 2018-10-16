NaCl_coupling(y)=0.0387./((y-8.222).^2+0.0778^2)
function pot_NaCl(R)
    x =R*0.52;
    alpha = -integral(NaCl_coupling,1,R)*0.52;
    U=zeros(2,2)
    U[1,1] = cos(alpha); U[1,2] = sin(alpha); U[2,1] = -U[1,2] ; U[2,2] = cos(alpha);
    H_ad=zeros(2,2)
    H_diab=zeros(2,2)

    De = 4.82; re = 2; be  =.8522; ce = 1.182;
    Vg = De*(1-ce*exp(-be*(x-re))).^2/27.2;
    a=-0.386 ;b =      0.0294  ;c =       86.76 ;d =       1.468  ;r =       6.01  ;
    Ve = (a*x+b*x.^2+c*exp(-d*x.^2)+ r)/27.2;
    H_ad[1,1] = Vg; H_ad[2,2] = Ve;
    H_diab = U'*H_ad*U;

    return H_diab,[zeros(2)]
end
