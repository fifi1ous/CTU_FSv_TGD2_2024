function [XYZ] = rusena_druzice(Nav_rin,GM,we)
%%
%input = [t0,t,Crs,dn,M0,Cuc,e,Cus,a,Toe,Cic,Wo,Cis,i0,Crc,w0,Wt,it0];
%%
t0= Nav_rin(1);
t=  Nav_rin(2);
Crs=Nav_rin(3);
dn= Nav_rin(4);
M0= Nav_rin(5);
Cuc=Nav_rin(6);
e=  Nav_rin(7);
Cus=Nav_rin(8);
a=  Nav_rin(9);
Toe=Nav_rin(10);
Cic=Nav_rin(11);
W0= Nav_rin(12);
Cis=Nav_rin(13);
i0= Nav_rin(14);
Crc=Nav_rin(15);
w0= Nav_rin(16);
Wt= Nav_rin(17);
it0=Nav_rin(18);


n=sqrt(GM/(a^3));
M=M0+(n+dn)*(t-t0);
E = excentricka_anomalie(e, M);
v=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
u0=v+w0;
w=w0+Cuc*cos(2*u0)+Cus*sin(2*u0);
u=v+w;
r0=a*(1-e*cos(E));
r=r0+Crc*cos(2*u)+Crs*sin(2*u);
i=i0+Cic*cos(2*u)+Cis*sin(2*u)+it0*(t-t0);
l=W0+(Wt-we)*(t-t0)-we*Toe;

r_xy=[r*cos(u);r*sin(u);0];

XYZ=[cos(-l)   sin(-l)     0
    -sin(-l)   cos(-l)     0
          0         0      1]*[1        0        0
                               0   cos(-i)  sin(-i)
                               0  -sin(-i)  cos(-i)]*r_xy;
XYZ=XYZ';

function E = excentricka_anomalie(e, M, mez)
    % ## Autor Michal Kovář ##
    % Tato funkce počítá excentrickou anomálii pomocí Newtonovy iterační metody
    % e - excentricita eliptické dráhy
    % M - střední anomálie
    % mez - volitelný parametr přesnosti výsledku... defaultně 1e-12

    % keplerova rovnice M = E - e * sin(E)
    % vyjádření exentrické anomálie E = M + e * sin(E)
    % f(E) = E - e * sin(E) - M
    % f'(E) = 1 - e * cos(E)
    % první odhad E_0 = M + (1 + e * cos(M)) * sin(M)
    % iterace = E_i = E_(i-1) - ((E_(i-1)-e*sin(M))/(1-e*cos(E_0))

    % Výpočet excentrické anomálie
    E_0 = M + (1 + e .* cos(M)) .* sin(M); % první odhad excentrické anomálie
    if nargin < 3
        mez = 1e-12; % přesnost pro ukončení cyklu
    end
    rozidl = inf(size(M));

    while any(abs(rozidl(:)) > mez)
        E_1 = E_0 - ((E_0 - e .* sin(E_0) - M) ./ (1 - e .* cos(E_0)));
        rozidl = E_0 - E_1;
        E_0 = E_1;
    end

    E = E_0;
end
end