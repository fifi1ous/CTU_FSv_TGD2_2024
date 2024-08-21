function [Vysledek,XYZ] = reed_o_rinnex(rinexo,code,DR,DATE,T)
%%  function which will read observe rinnex
%   Input       rinexo  -   name of the rinex file
%               code    -   name of the measured code
%               DR      -   type of satellite (G,R,E,...)   example:"G"
%               DATE    -   date of the rinnex              example:"2024 02 20"
%               T       -   matrix of measured time [hh,mm,ss], each epoch
%                           in the row

%   Output      vysledky-   cell of results int this structure
%                           [time,{satelite ID, distance}]
%               XYZ     -   Coordinates of station
%%
    t=T(:,1)*3600+T(:,2)*60+T(:,3);
    t=sort(t);
    
    T1=datestr(seconds(t),'HH:MM:SS');
    
    for n=1:size(T,1)
        T(n,:)=[str2double(T1(n,1:2)),str2double(T1(n,4:5)),str2double(T1(n,7:8))];
    end
    
    fid=fopen(rinexo,'r');
    konec = 0;
    j=1;
    dr=DR+"   ";
    pom = 0;
    pocet_mereni=0;
    K=0;
    MAX=size(T,1);
    cas_t=1;
    prefix="> ";
    sufix=" ";
    date=prefix+DATE+sufix;
    Vysledek={};
    Mereni={};
    
    if T(1,1)==0
        TIME="00 ";
    else
        if T(1,1)<10
            TIME="0"+string(T(1,1))+" ";
        else
            TIME=string(T(1,1))+" ";
        end
    end
    
    if T(1,2)==0
        TIME=TIME+"00 ";
    else
        if T(1,2)<10
            TIME=TIME+"0"+string(T(1,2))+" ";
        else
            TIME=TIME+string(T(1,2))+" ";
        end
    end
    
    if T(1,3)==0
        TIME=TIME+" 0";
    else
        TIME=TIME+string(T(1,3));
    end
    
    druzice=[];
    
    while konec==0
        radek = fgets(fid);
        if j==14
            XYZ=str2double(strsplit(radek));
            XYZ=XYZ(2:4);
        end
    
        if strfind(radek,dr)
            KON=1;
        else
            KON=0;
        end
    
        if (KON==1 || pom==1)&& K==0
            if strfind(radek,code)
                pozice=strsplit(radek);
                if pom==0
                    pozice=pozice(3:end);
                end
                pocet_mereni=find((strcmp(pozice, code)));
                pom=0;
                K=1;
            else
                pom = 1;
                pocet_mereni=pocet_mereni+13;
            end
        end
    
        if cas_t<=MAX && K==1
            if strfind(radek,date+TIME)
                radek = fgets(fid);
                while radek(1)~=">"
                    pom=0;
                    while radek(1)==DR
                        text=radek;
              
                        T_1={};
                        a=length(text);
                        T_text={text(1:3)};
                        text=text(4:end);
                        for n=1:(length(text)/16)-1
                            T_0=text(1:16);
                            T_1=[T_1,str2double(T_0(1:14))];
                            text=text(17:end);
                        end
                        T_0=text(1:16);
                        T_1=[T_1,str2double(T_0(1:14))];
                        P=NaN;
                        if length(T_1)>=pocet_mereni
                            if isnan (cell2mat(T_1(pocet_mereni)))
                            
                            else
                            Mereni=[Mereni;T_text,T_1(pocet_mereni)];
                            end
                        end
                        radek = fgets(fid);
                        pom=1;
                    end
    
                    if pom==1
                        cas_mereni={T(cas_t,:),Mereni};
                        Mereni={};
                        Vysledek=[Vysledek;cas_mereni];
                    end
    
                    if pom==1 && cas_t<MAX
                        cas_t=cas_t+1;
                        if T(cas_t,1)==0
                            TIME="00 ";
                        else
                            if T(cas_t,1)<10
                                TIME="0"+string(T(cas_t,1))+" ";
                            else
                                TIME=string(T(cas_t,1))+" ";
                            end
                        end
                        
                        if T(cas_t,2)==0
                            TIME=TIME+"00 ";
                        else
                            if T(cas_t,2)<10
                                TIME=TIME+"0"+string(T(cas_t,2))+" ";
                            else
                                TIME=TIME+string(T(cas_t,2))+" ";
                            end
                        end
                        
                        if T(cas_t,3)==0
                            TIME=TIME+" 0";
                        else
                            TIME=TIME+string(T(cas_t,3));
                        end
                    else
                        radek = fgets(fid);
                        if radek(1)~=">"
                            d1=length(radek);
                        end
                    end
                end
                d2=length(radek);
                fseek(fid,-(d2+d1),'cof');
            end
        end
    
        konec = feof(fid);
        j=j+1;
    end
    fclose(fid);
end