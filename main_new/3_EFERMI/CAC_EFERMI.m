addpath('../')
addpath('../0_Public')
parameter;
if ~exist('NCL','var')
    NCL=1;
end

EFERMI_Av=0;

EFERMI_t=zeros(length((1+NSTART+NSKIP):NSKIP:NCONF),NBEAD);
NC_loc=0;
for NC=(1+NSTART+NSKIP):NSKIP:NCONF
    NC_loc=NC_loc+1;
    for NB=1:NBEAD
    FID=fopen([ROOT_DIR,'/BEAD_',int2str(NC),'/HEAD_',int2str(NB)]);
    A=fread(FID,[3,3],'double');
    EFERMI=fread(FID,1,'double');
    TOTEN=fread(FID,1,'double');
    NKPTS=fread(FID,1,'int');ISPIN=fread(FID,1,'int');
    NTYP=fread(FID,1,'int');NITYP=zeros(1,NTYP);
    for NT=1:NTYP
        NITYP(NT)=fread(FID,1,'int');
    end
    fclose(FID); 
    EFERMI_t(NC_loc,NB)=EFERMI;
    EFERMI_Av=EFERMI_Av+EFERMI;
    end
end
EFERMI_Av=EFERMI_Av/(NBEAD*size((1+NSTART+NSKIP):NSKIP:NCONF,2))
if NCL==1
    if I_Ef==1
        save('../EFERMI.mat','EFERMI_Av');
    else
        save('../EFERMI_old.mat','EFERMI_Av');
    end
end
