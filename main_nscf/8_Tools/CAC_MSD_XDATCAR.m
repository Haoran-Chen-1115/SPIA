addpath('../')
addpath('../0_Public')
parameter;
NSTART=0;NSKIP=1;NBEAD=16;
L_FIXCOM = false; % whether remove the center of mass displacement
%POMASS=[1.00 1.00]; % need ion mass in the case

FID=fopen('../../XDATCAR_PI');

% Element type and numbers of ions and electrons
SYS_NAME=textscan(FID, '%[^\n]', 1);
SYS_NAME=char(SYS_NAME{1});
SCALE=cell2mat(textscan(FID, '%d', 1)); % Scaling factor of lattice parameter

% Lattice parameter in A
A=cell2mat(textscan(FID, '%f %f %f', 3));
A=A';
OMEGA=det(A);
B=2*pi*inv(A); % inverse lattice B/2pi

% Elements type
ELE_TYP=textscan(FID, '%[^\n]', 1);
ELE_TYP=split(char(ELE_TYP{1}));
NTYP=length(ELE_TYP);
NITYP=cell2mat(textscan(FID,'%f',NTYP))';
NTYP=length(NITYP);
NIONS=sum(NITYP);

%% Number of configurations
lines = 0;
while ~feof(FID)
    fgetl(FID);
    lines = lines +1;
end
frewind(FID);
NCONF = floor((lines)/((NIONS+1)*NBEAD+1));
%% Read in configurations
POSION=zeros(3,NIONS*NCONF*NBEAD);

NIS=0;
textscan(FID, '%f %f %f', 0, 'headerLines', 7);
for NC=1:NCONF
    %% Reading in a ion configuration
    NCC=textscan(FID, 'Direct configuration=%d', 1);
    data = textscan(FID, '%f %f %f', NIONS*NBEAD,'HeaderLines',0);
    
    %% Setting up ion positions
    % ion positions in fractional coordinates
    coefs = [data{1}  data{2}  data{3}];
    % turn to Cartesian coordinates
    POSION(:,NIS+1:NIS+NBEAD*NIONS) = coefs';
    NIS=NIS+NBEAD*NIONS;
    
    if feof(FID)
        break;
    end
end
fclose(FID);

POSION_t=reshape(POSION,3,NIONS,NBEAD,NCONF);

%%
MSD = zeros(length(1+NSTART:NSKIP:NCONF),NTYP,1);
MSD_i = zeros(length(1+NSTART:NSKIP:NCONF),sum(NITYP),1);
POS_t = zeros(3,NIONS,length(1+NSTART:NSKIP:NCONF),1);
DP_c_t = zeros(3,length(1+NSTART:NSKIP:NCONF));

POSION_INI = POSION_t(:,:,1,1+NSTART);

% open position files
NC_loc=0;

for NC=1+NSTART:NSKIP:NCONF
    NC_loc=NC_loc+1;
    
    DX=0;DY=0;DZ=0;
    for NB=1:NBEAD
        POSION=POSION_t(:,:,NB,NC);
        
        DX_2=POSION(1,:)-POSION_INI(1,:);
        DY_2=POSION(2,:)-POSION_INI(2,:);
        DZ_2=POSION(3,:)-POSION_INI(3,:);
        
        %
        DX_2(DX_2>1/2)=-1+DX_2(DX_2>1/2);
        DY_2(DY_2>1/2)=-1+DY_2(DY_2>1/2);
        DZ_2(DZ_2>1/2)=-1+DZ_2(DZ_2>1/2);
        DX_2(DX_2<-1/2)=1+DX_2(DX_2<-1/2);
        DY_2(DY_2<-1/2)=1+DY_2(DY_2<-1/2);
        DZ_2(DZ_2<-1/2)=1+DZ_2(DZ_2<-1/2);
        
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
     
        DX_2=DX_2-DX_c;
        DY_2=DY_2-DY_c;
        DZ_2=DZ_2-DZ_c;
        
        DX=DX+DX_2;DY=DY+DY_2;DZ=DZ+DZ_2;
    end
    DX=DX/NBEAD;DY=DY/NBEAD;DZ=DZ/NBEAD;
    POS_t(:,:,NC_loc) = POSION_INI + [DX;DY;DZ];
    
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

POS_t=reshape(A*reshape(POS_t,3,[]),3,NIONS,[]);
%%
POS_Av = mean(POS_t,3);
DP = sum(abs(POS_Av-POSION_INI).^2,1);