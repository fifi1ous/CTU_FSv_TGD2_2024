clc; clear; format long G
% TG2 ul3 - zadani c.21
% ---------------------

% observacni RINEX
rinexo = 'gope051a00.24o';
rinexn = 'gope051a00.24n';

% identifikator kodoveho mereni
code = 'C5X';

% vypocetni epochy
t1 = [12 46  0]; % [h min s]
t2 = [12 46 15]; % [h min s]
t3 = [12 46 30]; % [h min s]

DR="G";
DATE="2024 02 20";

GM=398600.5;
we=7292115.1467*10^-11;
c=299792458.0;

T=[t1;t2;t3];
[Vysledek,XYZ] = reed_o_rinnex(rinexo,code,DR,DATE,T);

T1=[];
for n=1:size(Vysledek,1)
    T1=[T1;Vysledek{n,1}];
end

%%merni=[t-obsrvation,t-satelite,PRN,date,T-broadcast,[d1,d2,d3],IODE,efemeries,GPSweak,distance]
[mereni,DSK,U] = get_messurment(T1,DATE,rinexn,Vysledek,c);

udaje=[];
l=[];
for n=1:size(mereni,1)
    Nav_rin = export_cell(mereni(n,:));
    XYZ_d = rusena_druzice(Nav_rin,GM,we);
    udaje=[udaje;XYZ_d*1000,DSK(n)];
    l=[l;mereni{n,end}];
end
%         udaje=[ X Y Z Dt]
%                 1 2 3  4
Nez=XYZ;

for n=1:size(T1,1)
    Nez=[Nez,0];
end

%% Vyrovnání


nuls=[];
for n=1:size(T1,1)
    u=zeros(size(Vysledek{n,2},1),size(T1,1));
    u(:,n)=1;
    nuls=[nuls;u];
end
NEZ=0;

zz=1;
while max(abs(Nez(1:3)-NEZ))>0.001
    A=[(Nez(1)-udaje(:,1))./sqrt((udaje(:,1)-Nez(1)).^2+(udaje(:,2)-Nez(2)).^2+(udaje(:,3)-Nez(3)).^2),(Nez(2)-udaje(:,2))./sqrt((udaje(:,1)-Nez(1)).^2+(udaje(:,2)-Nez(2)).^2+(udaje(:,3)-Nez(3)).^2),(Nez(3)-udaje(:,3))./sqrt((udaje(:,1)-Nez(1)).^2+(udaje(:,2)-Nez(2)).^2+(udaje(:,3)-Nez(3)).^2)];
    A=[A,nuls];
    l0=sqrt((udaje(:,1)-Nez(1)).^2+(udaje(:,2)-Nez(2)).^2+(udaje(:,3)-Nez(3)).^2)+(Nez(4:end)*nuls')'-c*udaje(:,4);
    lc=l0-l;
    
    dh=-inv(A'*A)*A'*lc;
    NEZ=Nez(1:3);
    Nez=Nez+dh';
    
    v1=A*dh+lc;
    v2=sqrt((udaje(:,1)-Nez(1)).^2+(udaje(:,2)-Nez(2)).^2+(udaje(:,3)-Nez(3)).^2)+(Nez(4:end)*nuls')'-c*udaje(:,4)-l;
    zz=zz+1;
end
Nez(4:end)=Nez(4:end)/c;

s0=sqrt((v1'*v1)/(size(A,1)-length(Nez)));
Qx=inv(A'*A);
sm=s0*sqrt(diag(Qx));
sm(4:end)=sm(4:end)/c;

