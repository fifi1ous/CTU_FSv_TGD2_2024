function [mereni,DSK,U] = get_messurment(T,DATE,rinexn,Vysledek,c)
%% get_messurment
%   Input           T       -   matrix of times,    row-[hh mm ss]
%                   GM      -   constant
%                   we      -   constant
%                   DATE    -   date name in format "2024 02 20"
%                   rinexn  -   name of the navigation rinnex
%                   Vylsedek-   cell of results int this structure
%                               [time,{satelite ID, distance}]
%
%   Output          mereni  -   cell of navigation and observation rinnex
%                               [t-obsrvation,t-satelite,PRN,date,T-broadcast,[d1,d2,d3],IODE,efemeries,GPSweak,distance]
%                   DSK     -   correction of cloks
%%
mereni={};
U=[];
DSK=[];
for nn=1:size(T,1)
    jj=size(Vysledek{nn,2},1);
    for ii=1:jj
        EFE = reed_n_rinnex(rinexn,Vysledek{nn,1},Vysledek{nn,2}{ii,1},DATE);
        % U=[U;EFE];

        EFE = [EFE{1}(1:3),EFE{1}(4)*3600+EFE{1}(5)*60+EFE{1}(6),EFE(2:end)];

        t0_1=Vysledek{nn,1}(1)*3600+Vysledek{nn,1}(2)*60+Vysledek{nn,1}(3);
        %%
        Tk=t0_1-Vysledek{nn,2}{ii,2}/c;
        DSK1=EFE{3}(1)+EFE{3}(2)*(Tk-EFE{2})+EFE{3}(3)*(Tk-EFE{2})^2;
        DSK=[DSK;DSK1];
        TK=t0_1-Vysledek{nn,2}{ii,2}/c-DSK1;
        %%
        mereni=[mereni;t0_1,TK,Vysledek{nn,2}{ii,1},EFE,Vysledek{nn,2}{ii,2}];
    end    
end

function [EFE] = reed_n_rinnex(rinnexn,t,PRN,DATE)
%% reed navigation rinnex
% 
% Input:        rinnexn:    name of the navigation rinnex
%               t      :    time of the observation         example:[02 00 00]
%               PRN    :    Name of the sattelite           example:"G01"
%               DATE   :    Date of the observation         example:"2024 02 20"
%
% Output:       EFE    :    navigation efemeries in cell {[date, time], [time derivates], IODE, [efemeries], GPSWEAK}  
%%
konec=0;
fid=fopen(rinnexn,'r');

j=0;
roz2=0;

pom=char(DATE);
time=t(1)*3600+t(2)*60+t(3);
while konec==0
        radek = fgets(fid);
        if PRN==radek(1:3)
            if strcmp(radek(5:14),DATE)
                T1=str2double(radek(16:17))*3600+str2double(radek(19:20))*60+str2double(radek(22:23));
            else
                if str2double(radek(13:14))<str2double(pom(9:10))
                    T1=-(24*3600-(str2double(radek(16:17))*3600+str2double(radek(19:20))*60+str2double(radek(22:23))));
                else
                    T1=str2double(radek(16:17))*3600+str2double(radek(19:20))*60+str2double(radek(22:23))+24*3600;
                end
            end

            roz=time-T1;
            if roz<=0 && roz2>0
               if roz==0
                   u=j+1;
               end
               break
            else
                u=j+1;
            end
            roz2=roz;
        end
        j=j+1;
        konec = feof(fid);
end
fclose(fid);

fid=fopen(rinnexn,'r');
konec=0;
j=0;
konst=6;
EFE=[];
while konec==0
    j=j+1;
    radek = fgets(fid);
    if j==u
        radek=strrep(radek(5:end), 'D', 'e');
        radek=strrep(radek, 'E', 'e');
        DTstamp=str2double(strsplit(radek(1:19)));
        del=[str2double(radek(20:38)),str2double(radek(39:57)),str2double(radek(58:78))];
        for n=1:konst
            radek = fgets(fid);
            radek=strrep(radek, 'D', 'e');
            radek=strrep(radek, 'E', 'e');
            POM=[str2double(radek(4:23)),str2double(radek(24:42)),str2double(radek(43:61)),str2double(radek(62:80))];
            EFE=[EFE;POM];
        end
        break
    end
    konec = feof(fid);
end
fclose(fid);
IODE=EFE(1,1);
efe=[EFE(1,2:4),EFE(2,:),EFE(3,:),EFE(4,:),EFE(5,1)];
GPSWEAK=EFE(5,3);
EFE={DTstamp,del,IODE,efe,GPSWEAK};
end

end