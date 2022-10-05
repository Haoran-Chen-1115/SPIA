%% Load results obtained before
addpath('../')
addpath('../0_Public')
tic;
gpuDevice(1);
parameter;
%load('FERMI_Surface.mat','kF','idF','NGF')
%load('INDEX.mat','G','G_INDEX')

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


%% Calculate Density Correlation Function, using double precision
dr_RDF=0.01;
gr = cell(NTYP,NTYP);
rmax_RDF=A(1,1); % Maximun distance
r_RDF=0:dr_RDF:rmax_RDF;
for NT1=1:NTYP
    for NT2=NT1:NTYP
        gr{NT1,NT2}=0; % initial value for radial distribution function
    end
end

NSKIP=5;NSTART=0;NCONF=3000;NB=1;
% open position files
for NC=1+NSKIP+NSTART:NSKIP:NCONF
    %for NC=1
    NB=1;
    POSION=cell(1,NTYP);
    for NT=1:NTYP
        POSS=fopen(['../../BEAD_',int2str(NC),'/POS_',int2str(NB),'_',int2str(NT)]);
        POSION{NT}=gpuArray(...
            fread(POSS,[3,NITYP(NT)],'double'));
        fclose(POSS);
    end
    POSION_INI=POSION;
    POSION_Av=POSION;

    for NB=2:NBEAD
        POSION=cell(1,NTYP);
        for NT=1:NTYP
            POSS=fopen(['../../BEAD_',int2str(NC),'/POS_',int2str(NB),'_',int2str(NT)]);
            POSION{NT}=gpuArray(...
                fread(POSS,[3,NITYP(NT)],'double'));
            fclose(POSS);

            DX=POSION{NT}(1,:)-POSION_INI{NT}(1,:);
            DY=POSION{NT}(2,:)-POSION_INI{NT}(2,:);
            DZ=POSION{NT}(3,:)-POSION_INI{NT}(3,:);

            DX(DX>1/2)=-1+DX(DX>1/2);
            DY(DY>1/2)=-1+DY(DY>1/2);
            DZ(DZ>1/2)=-1+DZ(DZ>1/2);
            DX(DX<-1/2)=1+DX(DX<-1/2);
            DY(DY<-1/2)=1+DY(DY<-1/2);
            DZ(DZ<-1/2)=1+DZ(DZ<-1/2);

            POSION{NT}(1,:)=POSION_INI{NT}(1,:)+DX;
            POSION{NT}(2,:)=POSION_INI{NT}(2,:)+DY;
            POSION{NT}(3,:)=POSION_INI{NT}(3,:)+DZ;

            POSION_Av{NT}(1,:)=POSION_Av{NT}(1,:)+POSION{NT}(1,:);
            POSION_Av{NT}(2,:)=POSION_Av{NT}(2,:)+POSION{NT}(2,:);
            POSION_Av{NT}(3,:)=POSION_Av{NT}(3,:)+POSION{NT}(3,:);
        end
    end

    for NT=1:NTYP
         POSION{NT}= POSION_Av{NT}/NBEAD;
    end

    for NT1=1:NTYP
        for NT2=NT1:NTYP
            MAXC1=1;MAXC2=1;MAXC3=1;
            NX=3;NY=3;NZ=3;

            PX=reshape(repmat(POSION{NT2}(1,:).',1,NX)+[0:MAXC1,-MAXC1:-1],1,[]);
            PX=reshape(repmat(...
                reshape(PX...
                ,1,[],NX),1,1,NY*NZ),1,[]);

            PY=reshape(repmat(POSION{NT2}(2,:).',1,NY)+[0:MAXC2,-MAXC2:-1],1,[]);
            PY=reshape(repmat(...
                reshape(PY...
                ,1,[],NY),1,NX,NZ),1,[]);

            PZ=reshape(repmat(POSION{NT2}(3,:).',1,NZ)+[0:MAXC3,-MAXC3:-1],1,[]);
            PZ=reshape(repmat(...
                reshape(PZ...
                ,1,[],NZ),1,NX*NY,1),1,[]);

            DX=POSION{NT1}(1,:)-PX.';
            DY=POSION{NT1}(2,:)-PY.';
            DZ=POSION{NT1}(3,:)-PZ.';

            DX_2=A(1,1)*DX+A(1,2)*DY+A(1,3)*DZ;
            DY_2=A(2,1)*DX+A(2,2)*DY+A(2,3)*DZ;
            DZ_2=A(3,1)*DX+A(3,2)*DY+A(3,3)*DZ;

            DX=DX_2;DY=DY_2;DZ=DZ_2;

            DR=sqrt(DX.^2+DY.^2+DZ.^2);
            DR(abs(DR)<1E-8)=1E-9;

            DR_idx=DR(DR<rmax_RDF+dr_RDF/2 & abs(DR)>1E-5);

            if NT1==NT2
                NP=length(POSION{NT1})*(length(POSION{NT1})-1);
                rho=NP/OMEGA;
            else
                NP=length(POSION{NT1})*length(POSION{NT2});
                rho=NP/OMEGA;
            end

            grx=hist(DR_idx,r_RDF);
            grx(2:end)=grx(2:end)./((r_RDF(2:end)+dr_RDF/2).^3-(r_RDF(2:end)-dr_RDF/2).^3);
            grx=grx/rho/(4*pi/3);
            gr{NT1,NT2}=gr{NT1,NT2}+grx;

        end
    end

end

%%
for NT1=1:NTYP
    for NT2=NT1:NTYP
        gr{NT1,NT2} = gr{NT1,NT2} /length(1+NSKIP+NSTART:NSKIP:NCONF);
    end
end
