function data = sp32inerc(file,prn,xp,yp,DUT1,leaps)
% --------------------------------------------------
% pro dany sp3 soubor a druzici prepocte presne efemeridy
% druzice (souradnice X,Y,Z) v terestrickem systemu do
% inercialniho nebeskeho systemu ICRS
% 
% IN:  file ... nazev sp3 souboru, verze c (string)
%       prn ... identifikator druzice (string) (pr.: 'G01')
%    xp, yp ... souradnice pohybu polu ['']
%      DUT1 ... DUT1 ( = UT1-UTC) [s]
%     leaps ... leap seconds ( = TAI-UTC) [s]
%               hodnota je v Bulletinu C mezinar.sluzby IERS
% OUT: data ... pro kazdou epochu pro danou druzici obsahuje:
%               data = [i MJD t X Y Z]
%                 i ... poradove cislo epochy
%               MJD ... modifikovane julianske datum dne [den]
%                 t ... epocha v case UTC [s od pocatku dne]
%                 X ... X_ICRS [km]
%                 Y ... Y_ICRS [km]
%                 Z ... Z_ICRS [km]
% --------------------------------------------------
% functions used:
% Rx, Ry, Rz
% --------------------------------------------------

if ischar(file) ~= 1
    error('sp32inerc: Jmeno souboru v promenne file musi byt retezec.');
end
if ischar(prn) ~= 1
    error('sp32inerc: PRN druzice v promenne prn musi byt retezec.');
end

% -------------------
% nacteni sp3 souboru
% -------------------
fid = fopen(file);
fgets(fid);
radek = fgets(fid);
% modifikovane julianske datum dne sp3 souboru
MJD = str2num(radek(40:44));
while radek(1) ~= '*'
    radek = fgets(fid);
end

% cteni dat
t = [];
data = [];
i = 0;
sat = 0;
eof = 0; % pro konec souboru
while eof == 0
    if radek(1) == '*'  % nova epocha [h min s]
        i=i+1;
        sat=0;
        h = str2num(radek(15:16));
        min = str2num(radek(18:19));
        s = str2num(radek(21:31));
        t(i) = h*3600 + min*60 + s; % [s od pocatku dne]
        radek = fgets(fid);
    end
    while radek(1) ~= '*' & strcmp(radek(1:3),'EOF') ~= 1
        if strcmp(radek(2:4),prn) == 1
            sat = 1;
            % poloha druzice
            X = str2num(radek(5:18));  % [km]
            Y = str2num(radek(19:32)); % [km]
            Z = str2num(radek(33:46)); % [km]
            data(i,:) = [i MJD t(i) X Y Z];
        end
        radek = fgets(fid);
    end
    if sat == 0 % druzice v epose t neni
        i = i-1;
        fprintf('sp32inerc: V epose%3dh%3dmin%5.1fs druzice %s nenalezena.\n',h,min,s,prn);
    end
    eof = feof(fid);
end
fclose(fid);
fprintf('sp32inerc: Nacteno %d epoch druzice %s.\n\n',i,prn);

% -------------------------------------------
% prevod z terestricke do inercialni soustavy
% -------------------------------------------
% polohy z sp3 v terestricke soustave
XYZter = data(:,4:6)'; % [X;Y;Z] [km]
% epochy z sp3 v case GPS
tgps = data(:,3);      % [s od pocatku dne] % vektor

% GPS -> TAI
TAI = tgps + 19;       % [s od pocatku dne] % vektor
% TAI -> UTC
UTC = TAI - leaps;     % [s od pocatku dne] % vektor

% nutace (dpsi a epsA ovlivnuje i S)
% ----------------------------------
T = (MJD-51544.5)/36525;
% sklon ekliptiky IAU 1976 (neovlivneny nutaci)
epsA = 84381.448 - 46.8150*T - 0.00059*T^2 + 0.001813*T^3; % ['']
epsA = epsA/3600*pi/180; % [rad]

% Delaunayovy promenne IAU 1980
% Omega = stredni delka vystupneho uzlu Mesice
Om = 125.04455501*3600 - 6962890.2665*T + 7.4722*T^2 - 0.007702*T^3 - 0.00005939*T^4;     % ['']
Om = Om/3600*pi/180;     % [rad]
Om = mod(Om,2*pi);
%    D  = stredni elongace Mesice od Slunce
D  = 297.85019547*3600 + 1602961601.2090*T - 6.3706*T^2 + 0.006593*T^3 - 0.00003169*T^4;  % ['']
D  = D/3600*pi/180;      % [rad]
D  = mod(D,2*pi);
%    F  = stredni delka Mesice - Omega
F  = 93.27209062*3600 + 1739527262.8478*T - 12.7512*T^2 - 0.001037*T^3 + 0.00000417*T^4;  % ['']
F  = F/3600*pi/180;      % [rad]
F  = mod(F,2*pi);
%    L' = stredni anomalie Slunce
Lp = 357.52910918*3600 + 129596581.0481*T - 0.5532*T^2 - 0.000136*T^3 - 0.00001149*T^4;   % ['']
Lp = Lp/3600*pi/180;     % [rad]
Lp = mod(Lp,2*pi);
%    L  = stredni anomalie Mesice
L  = 134.96340251*3600 + 1717915923.2178*T + 31.8792*T^2 + 0.051635*T^3 - 0.00024470*T^4; % ['']
L  = L/3600*pi/180;      % [rad]
L  = mod(L,2*pi);

% nutacni uhly IAU 1980
deps = ( 9.2025 + 0.00089*T)*cos(Om) ...
     + ( 0.5736 - 0.00031*T)*cos(2*F-2*D+2*Om) ...
     + ( 0.0977 - 0.00005*T)*cos(2*F+2*Om) ...
     + (-0.0895 + 0.00005*T)*cos(2*Om) ...
     + ( 0.0054 - 0.00001*T)*cos(-Lp) ...
     + (-0.0007 + 0.00000*T)*cos(L) ...
     + ( 0.0224 - 0.00006*T)*cos(Lp+2*(F-D+Om)) ...
     + ( 0.0200            )*cos(2*F+Om) ...
     + ( 0.0129 - 0.00001*T)*cos(L+2*(F+Om)) ...
     + (-0.0095 + 0.00003*T)*cos(-Lp+2*(F-D+Om)) ...
     + (-0.0001            )*cos(-L+2*D) ...
     + (-0.0070            )*cos(2*F-2*D+Om) ...
     + (-0.0053            )*cos(-L+2*(F+Om));     % ['']
deps = deps/3600*pi/180; % [rad]
   
dpsi =(-17.1996 - 0.01742*T)*sin(Om) ...
     + (-1.3187 - 0.00016*T)*sin(2*F-2*D+2*Om) ...
     + (-0.2274 - 0.00002*T)*sin(2*F+2*Om) ...
     + ( 0.2062 + 0.00002*T)*sin(2*Om) ...
     + (-0.1426 + 0.00034*T)*sin(-Lp) ...
     + ( 0.0712 + 0.00001*T)*sin(L) ...
     + (-0.0517 + 0.00012*T)*sin(Lp+2*(F-D+Om)) ...
     + (-0.0386 - 0.00004*T)*sin(2*F+Om) ...
     + (-0.0301            )*sin(L+2*(F+Om)) ...
     + ( 0.0217 - 0.00005*T)*sin(-Lp+2*(F-D+Om)) ...
     + ( 0.0158            )*sin(-L+2*D) ...
     + ( 0.0129 + 0.00001*T)*sin(2*F-2*D+Om) ...
     + ( 0.0123            )*sin(-L+2*(F+Om));     % ['']
dpsi = dpsi/3600*pi/180; % [rad]

% nutacni matice IAU 1980
N = Rx(-epsA)*Rz(dpsi)*Rx(epsA+deps); 

% svetovy hvezdny cas S
% ---------------------
du = floor(MJD) - 51544.5; % MJD = 51544.5 pro 1.1.2000 12h UT1
tu = du / 36525;
% stredni svetovy hvezdny cas S0 pro svetovou pulnoc (pro 0h UT1)
S0h = (6*3600+41*60+50.54841) + 8640184.812866*tu + 0.093104*tu^2 - 6.2e-6*tu^3; % [h]
S0 = S0h * pi/12;  % [rad]
S0 = mod(S0,2*pi); % <0,2*pi>
mi1 = 1.002737909350795 + 5.9006e-11*tu - 5.9e-15*tu^2;
% stredni svetovy hvezdny cas
Sstredni = S0 + mi1*(UTC+DUT1)/3600*pi/12; % [rad]  % vektor
% vliv nutace
dR = 0.00264*sin(Om) + 0.000063*sin(2*Om); % ['']
R = dpsi*cos(epsA) + dR/3600*pi/180;       % [rad]
% pravy svetovy hvezdny cas
S = Sstredni + R;                          % [rad]  % vektor

% precese
% -------
T = (MJD-51544.5)/36525;
% precesni uhly IAU 1976
dzeta = 2306.2181*T + 0.30188*T^2 + 0.017998*T^3; % ['']
theta = 2004.3109*T - 0.42665*T^2 - 0.041833*T^3; % ['']
z     = 2306.2181*T + 1.09468*T^2 + 0.018203*T^3; % ['']
dzeta = dzeta/3600*pi/180; % [rad]
theta = theta/3600*pi/180; % [rad]
z     = z/3600*pi/180;     % [rad]
% precesni matice
P = Rz(dzeta)*Ry(-theta)*Rz(z);

% pohyb polu
% ----------
xp = xp/3600*(pi/180); % [rad]
yp = yp/3600*(pi/180); % [rad]

% souradnice v ITRS
% -----------------
PN = P*N;
POLE = Rx(yp)*Ry(xp)*XYZter;
for i=1:length(S)
    XYZinerc(:,i) = PN*Rz(-S(i))*POLE(:,i);
end

% vysledky
data(:,3) = UTC;
for i=1:length(UTC)
    if UTC(i)<0
        data(i,2) = data(i,2) - 1; % MJD = MJD - 1
        data(i,3) = data(i,3) + 3600*24;
    elseif UTC(i)>=3600*24
        data(i,2) = data(i,2) + 1; % MJD = MJD + 1
        data(i,3) = data(i,3) - 3600*24;
    end
end
data(:,4:6) = XYZinerc';
