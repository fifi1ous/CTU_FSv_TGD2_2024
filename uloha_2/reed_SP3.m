function [M,MJD,ERP] = reed_SP3(sp3,prn,erp)
%% Function that will read and export from spr file a time series of specific satellite
% INPUT     spr3    :   Name of the sp3 file
%           prn     :   Name of the Satellite
%           erp     :   Name of the erp file
% OUTPUT    M       :   Matrix of the time series
%           MJD     :   Julian day
%           ERP     :   Parameters of ERP

%%
if length(prn)~=4
    prn="P"+ prn;
end

fid=fopen(sp3,'r');
konec = 0;
n=1;
M=[];

while konec==0
    radek = fgets(fid);
        if strfind(radek,prn)
            VEC=str2double(strsplit(radek(5:end)));
            if length(VEC)~=10
                 VEC=[VEC(2:end-1),0];
            else
                 VEC=[VEC(2:end-1)];
            end
            M=[M;VEC];
        end
        if n==2
            MJD=radek;
        end
    n=n+1;
    konec = feof(fid); 
end
fclose(fid);

MJD=str2double(strsplit(MJD));
MJD=MJD(5)+0.5;

if nargin>2
    text=num2str(MJD)+"0";
    konec = 0;
    fid=fopen(erp,'r');
    while konec==0
        radek = fgets(fid);
            if strfind(radek,text)
                ERP=str2double(strsplit(radek));
            end
        konec = feof(fid); 
    end
    ERP=ERP(1:end-1);
end
end