function [ss,r_xy,E,M,r,n] = ss_ICRS_1(a,e,t0,w,i,W,t,GM)
    n=sqrt(GM/(a^3));
    M=n*(t-t0);
    E = excentricka_anomalie(e, M);
    r=a*(1-e*cos(E));
    v=atan2((sqrt(1-e^2)*sin(E)),(cos(E)-e));
    r_xy=[r*cos(v);r*sin(v);0];
    ss=[cos(-W)   sin(-W)     0
       -sin(-W)   cos(-W)     0
             0       0        1]*[1       0          0
                                  0   cos(-i)   sin(-i)
                                  0  -sin(-i)   cos(-i)]*[cos(-w)   sin(-w)     0
                                                         -sin(-w)   cos(-w)     0
                                                             0          0       1]*r_xy;
    ss=ss';
    r_xy=r_xy';
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