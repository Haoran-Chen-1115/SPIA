%% Load results obtained before
tic;
addpath('../')
addpath('../0_Public')
parameter
L_FIXCOM = false; % whether remove the center of mass displacement
%POMASS=[1.00 1.00]; % need ion mass in the case

FID=fopen([ROOT_DIR,'/BEAD_1/HEAD_1']);
A=fread(FID,[3,3],'double');
EFERMI=fread(FID,1,'double');
TOTEN=fread(FID,1,'double');
NKPTS=fread(FID,1,'int');ISPIN=fread(FID,1,'int');
NTYP=fread(FID,1,'int');NITYP=zeros(1,NTYP);
for NT=1:NTYP
    NITYP(NT)=fread(FID,1,'int');
end
fclose(FID);
NIONS=sum(NITYP);
B=2*pi*inv(A);
OMEGA=det(A);

% open position files
MSD = zeros(length(1+NSTART:NSKIP:NCONF),NTYP);
MSD_i = zeros(length(1+NSTART:NSKIP:NCONF),sum(NITYP));
POS_t = zeros(3,NIONS,NBEAD,length(1+NSTART:NSKIP:NCONF));
DP_c_t = zeros(3,length(1+NSTART:NSKIP:NCONF));

POSION_INI = zeros(3,NIONS);
NIS = 1;
for NT=1:NTYP
    POSS=fopen([ROOT_DIR,'/BEAD_',int2str(1+NSTART),'/POS_1_',int2str(NT)]);
    POSION_INI(:,NIS:NIS+NITYP(NT)-1)=...
        fread(POSS,[3,NITYP(NT)],'double');
    NIS=NIS+NITYP(NT);
    fclose(POSS);
end

NC_loc=0;
for NC=1+NSTART:NSKIP:NCONF
    NC_loc=NC_loc+1;
    DX=zeros(1,sum(NITYP));DY=zeros(1,sum(NITYP));DZ=zeros(1,sum(NITYP));
    for NB=1:NBEAD
        POSION=zeros(3,NIONS);
        NIS=1;
        for NT=1:NTYP
            POSS=fopen([ROOT_DIR,'/BEAD_',int2str(NC),'/POS_',int2str(NB),'_',int2str(NT)]);
            POSION(:,NIS:NIS+NITYP(NT)-1)=...
                fread(POSS,[3,NITYP(NT)],'double');
            NIS=NIS+NITYP(NT);
            fclose(POSS);
        end
        
        DX_t=POSION(1,:)-POSION_INI(1,:);
        DY_t=POSION(2,:)-POSION_INI(2,:);
        DZ_t=POSION(3,:)-POSION_INI(3,:);
        
        %
        DX_t(DX_t>1/2)=-1+DX_t(DX_t>1/2);
        DY_t(DY_t>1/2)=-1+DY_t(DY_t>1/2);
        DZ_t(DZ_t>1/2)=-1+DZ_t(DZ_t>1/2);
        DX_t(DX_t<-1/2)=1+DX_t(DX_t<-1/2);
        DY_t(DY_t<-1/2)=1+DY_t(DY_t<-1/2);
        DZ_t(DZ_t<-1/2)=1+DZ_t(DZ_t<-1/2);
        
        if L_FIXCOM
            % Remove center of mass displacement
            NIS=0; DX_c=0; DY_c=0; DZ_c=0;
            for NT=1:NTYP
                DX_c = DX_c+sum(POMASS(NT)*DX_2(NIS+1:NIS+NITYP(NT)));
                DY_c = DY_c+sum(POMASS(NT)*DY_2(NIS+1:NIS+NITYP(NT)));
                DZ_c = DZ_c+sum(POMASS(NT)*DZ_2(NIS+1:NIS+NITYP(NT)));
                NIS=NIS+NITYP(NT);
            end
            DX_c=DX_c/sum(POMASS.*NITYP);
            DY_c=DY_c/sum(POMASS.*NITYP);
            DZ_c=DZ_c/sum(POMASS.*NITYP);
            DP_c_t(:,NC_loc)=[DX_c;DY_c;DZ_c];
        else
            DX_c=0; DY_c=0; DZ_c=0;
        end
        
        DX_t=DX_t-DX_c;
        DY_t=DY_t-DY_c;
        DZ_t=DZ_t-DZ_c;
        
        POS_t(:,:,NB,NC_loc) = POSION_INI + [DX_t;DY_t;DZ_t];
        
        DX=DX+DX_t;
        DY=DY+DY_t;
        DZ=DZ+DZ_t;
    end
    DX=DX/NBEAD; DY=DY/NBEAD; DZ=DZ/NBEAD;
    
    DX_2=A(1,1)*DX+A(1,2)*DY+A(1,3)*DZ;
    DY_2=A(2,1)*DX+A(2,2)*DY+A(2,3)*DZ;
    DZ_2=A(3,1)*DX+A(3,2)*DY+A(3,3)*DZ;
    
    DX=DX_2;DY=DY_2;DZ=DZ_2;
    
    DR=DX.^2+DY.^2+DZ.^2;
    
    NIS=0;
    for NT=1:NTYP
        MSD(NC_loc,NT) = sum(DR(NIS+1:NIS+NITYP(NT)))/NITYP(NT);
        MSD_i(NC_loc,NIS+1:NIS+NITYP(NT))=DR(NIS+1:NIS+NITYP(NT));
        NIS=NIS+NITYP(NT);
    end
end

%%
POS_Av = mean(POS_t,4);
DP = sum(abs(POS_Av-POSION_INI).^2,1);
