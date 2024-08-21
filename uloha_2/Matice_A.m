function [A] = Matice_A(a,e,w,i,W,E,M,r,n,ss)
RzW= [cos(-W)   sin(-W)     0
     -sin(-W)   cos(-W)     0
        0           0       1];

Rzw= [cos(-w)   sin(-w)     0
     -sin(-w)   cos(-w)     0
        0           0       1];

Rxi=[1       0          0
     0  cos(-i)   sin(-i)
     0  -sin(-i)  cos(-i)];

dRzW=(-1)*[-sin(-W)   cos(-W)     0
           -cos(-W)  -sin(-W)     0
               0         0      0];

dRzw=(-1)*[-sin(-w)   cos(-w)     0
           -cos(-w)  -sin(-w)     0
                0        0      0];

dRxi=(-1)*[0       0         0
           0  -sin(-i)    cos(-i)
           0  -cos(-i)   -sin(-i)];

R=RzW*Rxi*Rzw;

da=R*[cos(E)-e+(3/2)*(a/r)*M*sin(E);
    sqrt(1-e^2)*(sin(E)-(3/2)*(a/r)*M*cos(E));
    0];
de=R*[-a*(1+a/r*(sin(E)^2));
    (((a^2)*sin(E))/(r*sqrt(1-e^2)))*(cos(E)-e);
    0];
dt0=R*[n*((a^2)/r)*sin(E);
    (-n)*((a^2)/r)*sqrt(1-e^2)*cos(E);
    0];

A=[da,    de,   dt0,    RzW*Rxi*dRzw*ss',    RzW*dRxi*Rzw*ss',   dRzW*Rxi*Rzw*ss'];
end