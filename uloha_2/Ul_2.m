clc; clear; format long G; close all
% TG2 ul2 - zadani c.21
% ---------------------

% SP3 soubor
sp3 = 'IGS0OPSFIN_20240430000_01D_15M_ORB.SP3';
erp=  'IGS0OPSFIN_20240420000_07D_01D_ERP.ERP';

% PRN druzice
prn = 'G16';

[Mat, MJD,ERP] = reed_SP3(sp3,prn,erp);
xp=ERP(2)*10^(-6);
yp=ERP(3)*10^(-6);
DUT1=ERP(4)*10^(-7);
leaps=37;
data = sp32inerc(sp3,prn,xp,yp,DUT1,leaps);
GM=398600.5;

KEP = Approximate_keppler_elements(data(2:3,:),GM);

cas=data(2:end,3);
SS_ICRS=data(2:end,4:6);
ROZ=1;
kepler=zeros(size(data,1)-2,6);

for i=1:size(kepler,1)
    while any(abs(ROZ)>0.000001)
        [ss,r_xy,E,M,r,n] = ss_ICRS_1(KEP(1),KEP(2),KEP(3),KEP(4),KEP(5),KEP(6),cas(i),GM);
        A1=Matice_A(KEP(1),KEP(2),KEP(4),KEP(5),KEP(6),E,M,r,n,r_xy);
        
        [ss1,r_xy,E,M,r,n] = ss_ICRS_1(KEP(1),KEP(2),KEP(3),KEP(4),KEP(5),KEP(6),cas(i+1),GM);
        A2=Matice_A(KEP(1),KEP(2),KEP(4),KEP(5),KEP(6),E,M,r,n,r_xy);
        
        A=[A1;A2];
        ROZ=[(SS_ICRS(i,:)-ss),(SS_ICRS(i+1,:)-ss1)]';
        dte=A^-1*ROZ;
        KEP=KEP+dte';
    end
    ROZ=1;
    if KEP(4)<0
        KEP(4)=KEP(4)+2*pi;
    end
    if KEP(5)<0
        KEP(5)=KEP(5)+2*pi;
    end
    if KEP(6)<0
        KEP(6)=KEP(6)+2*pi;
    end
    kepler(i,:)=KEP;
end
cas=cas(1:end-1);

t_hod=cas/3600;
figure(1)
subplot(2,3,1);
plot(t_hod,kepler(:,1))
title('Oskulační elementy: a')
xlabel('t [h]'); ylabel('a [km]')
grid on

% figure(2)
subplot(2,3,2);
plot(t_hod,kepler(:,2))
title('Oskulační elementy: e')
xlabel('t [h]'); ylabel('e')
grid on

% figure(3)
subplot(2,3,3);
plot(t_hod,kepler(:,3))
title('Oskulační elementy: t_0')
xlabel('t [h]'); ylabel('t_0 [s]')
grid on

% figure(4)
subplot(2,3,4);
plot(t_hod,kepler(:,4)/pi*180)
title('Oskulační elementy: \omega')
xlabel('t [h]'); ylabel('\omega [°]')
grid on

% figure(5)
subplot(2,3,5);
plot(t_hod,kepler(:,5)/pi*180)
title('Oskulační elementy: i')
xlabel('t [h]'); ylabel('i [°]')
grid on

% figure(6)
subplot(2,3,6);
plot(t_hod,kepler(:,6)/pi*180)
title('Oskulační elementy: \Omega')
xlabel('t [h]'); ylabel('\Omega [°]')
grid on