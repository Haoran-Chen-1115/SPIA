%%
addpath('../')
addpath('../0_Public')
parameter
filename='../../XDATCAR_PI';

FID=fopen(filename);

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
[~,lines]=system(['wc -l <',filename]);
lines=str2double(lines);
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
MSD = zeros(length(1+NSTART_POS_new:NSKIP_POS_new:NCONF),NTYP,1);
MSD_i = zeros(length(1+NSTART_POS_new:NSKIP_POS_new:NCONF),sum(NITYP),1);
POS_t = zeros(3,NIONS,NBEAD,length(1+NSTART_POS_new:NSKIP_POS_new:NCONF),1);

POSION_INI = POSION_t(:,:,1,1+NSTART_POS_new);

% open position files
NC_loc=0;

for NC=1+NSTART_POS_new:NSKIP_POS_new:NCONF
    NC_loc=NC_loc+1;
    
    DX=0;DY=0;DZ=0;
    for NB=1:NBEAD
        POSION=POSION_t(:,:,NB,NC);
        
        if NC_loc==1
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
        else
            DX_2=POSION(1,:)-POS_t(1,:,NB,NC_loc-1);
            DY_2=POSION(2,:)-POS_t(2,:,NB,NC_loc-1);
            DZ_2=POSION(3,:)-POS_t(3,:,NB,NC_loc-1);
            
            %
            DX_2(DX_2>1/2)=-1+DX_2(DX_2>1/2);
            DY_2(DY_2>1/2)=-1+DY_2(DY_2>1/2);
            DZ_2(DZ_2>1/2)=-1+DZ_2(DZ_2>1/2);
            DX_2(DX_2<-1/2)=1+DX_2(DX_2<-1/2);
            DY_2(DY_2<-1/2)=1+DY_2(DY_2<-1/2);
            DZ_2(DZ_2<-1/2)=1+DZ_2(DZ_2<-1/2);
            
            %
            dx = POS_t(1,:,NB,NC_loc-1) - POSION_INI(1,:);
            dy = POS_t(2,:,NB,NC_loc-1) - POSION_INI(2,:);
            dz = POS_t(3,:,NB,NC_loc-1) - POSION_INI(3,:);
            
            DX_2 = dx + DX_2;
            DY_2 = dy + DY_2;
            DZ_2 = dz + DZ_2;
        end
        
        POS_t(:,:,NB,NC_loc) = POSION_INI + [DX_2;DY_2;DZ_2];
        DX=DX+DX_2;DY=DY+DY_2;DZ=DZ+DZ_2;
    end
    DX=DX/NBEAD;DY=DY/NBEAD;DZ=DZ/NBEAD;
    if NC_loc==1
        POSION_INI=...
            POSION_INI + [DX;DY;DZ];
        DX=zeros(size(DX));
        DY=zeros(size(DY));
        DZ=zeros(size(DZ));
    end
    
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

POS_t2=reshape(mod(POS_t,1),3,NIONS,[]);
POS_t=reshape(A*reshape(POS_t,3,[]),3,NIONS,[]);

save('POS_t_MLFF.mat','POS_t','POS_t2','-v7.3');

%%
load('POS_t_MLFF.mat','POS_t2')
POS_t2=reshape(POS_t2,3,NIONS,NBEAD,[]);

POS_t2_ori=POS_t2;

%%
POSION_new=cell(1,NTYP);
strength=cell(1,NTYP);
res=cell(1,NTYP);
N=cell(1,NTYP);

NIS=0;
for NT=1:NTYP
    POS_t2=POS_t2_ori(:,NIS+1:NIS+NITYP(NT),:,:);
    
    NBins=1000;
    NWidth=5; % Average over +-Nwidth box
    POS_t2(POS_t2>0.955)=POS_t2(POS_t2>0.955)-1;
    POS_t2=reshape(A*reshape(POS_t2,3,[]),3,NITYP(NT),NBEAD,[]);
    
    [N{NT},E,C]=histcn(reshape(POS_t2,3,[])',NBins,NBins,NBins);
    N{NT}=N{NT}(1:NBins,1:NBins,1:NBins);
    N3=N{NT};
    N2=N{NT};
    N{NT}=zeros(size(N{NT}));
    for NX=1:NBins
        box_x = NX-NWidth:NX+NWidth;
        box_x(box_x<1)=box_x(box_x<1)+1000;
        box_x(box_x>1000)=box_x(box_x>1000)-1000;
        
        N{NT}(NX,:,:) = mean(N2(box_x,:,:),1);
    end
    N2=N{NT};
    
    for NY=1:NBins
        box_y = NY-NWidth:NY+NWidth;
        box_y(box_y<1)=box_y(box_y<1)+1000;
        box_y(box_y>1000)=box_y(box_y>1000)-1000;
        
        N{NT}(:,NY,:) = mean(N2(:,box_y,:),2);
    end
    N2=N{NT};
    
    for NZ=1:NBins
        box_z = NZ-NWidth:NZ+NWidth;
        box_z(box_z<1)=box_z(box_z<1)+1000;
        box_z(box_z>1000)=box_z(box_z>1000)-1000;
        
        N{NT}(:,:,NZ) = mean(N2(:,:,box_z),3);
    end
    
    %%
    NBlock_n=ceil(sqrt(NITYP(NT)))*2;
    while mod(1000,NBlock_n)~=0
        NBlock_n=NBlock_n+1;
    end
    NBlock=linspace(0,NBins,NBlock_n+1);
    POSION_new{NT}=zeros(3,NBlock_n^3);
    strength{NT}=zeros(1,NBlock_n^3);
    res{NT}=zeros(1,NBlock_n^3);
    NB=0;
    eps=1E-2;
    
    maxstrength=max(max(max(N{NT})));
    for NZ=1:NBlock_n
        NZ_b=NBlock(NZ)+1:NBlock(NZ+1);
        for NY=1:NBlock_n
            NY_b=NBlock(NY)+1:NBlock(NY+1);
            for NX=1:NBlock_n
                %NB=NB+1;
                NX_b=NBlock(NX)+1:NBlock(NX+1);
                [x,y,z]=ndgrid(C{1}(NX_b),C{2}(NY_b),C{3}(NZ_b));
                x=reshape(x,[],1);
                y=reshape(y,[],1);
                z=reshape(z,[],1);
                if max(reshape(N{NT}(NX_b,NY_b,NZ_b),1,[]))>0.2*maxstrength
                    [gfit,res2]=...
                        gaussfitn([x,y,z],reshape(N{NT}(NX_b,NY_b,NZ_b),[],1),...
                        [],[],[],'Display','off');
                    x_in = gfit{3}(1)<=max(x)+eps && gfit{3}(1)>=min(x)-eps;
                    y_in = gfit{3}(2)<=max(y)+eps && gfit{3}(2)>=min(y)-eps;
                    z_in = gfit{3}(3)<=max(z)+eps && gfit{3}(3)>=min(z)-eps;
                    
                    if x_in && y_in && z_in
                        NB=NB+1;
                        res{NT}(NB)=res2;
                        POSION_new{NT}(:,NB)=gfit{3};
                        strength{NT}(NB)=gfit{2};
                    end
                end
            end
        end
    end
    
    POSION_new{NT}=POSION_new{NT}(:,1:NB);
    strength{NT}=strength{NT}(1:NB);res{NT}=res{NT}(1:NB);
    NIS = NIS + NITYP(NT);
end
save('POS_gauss','POSION_new','NB',...
    'POSION_new','res','strength','N','-v7.3');

%%
% Transform to primitive cell
ia=cell(1,NTYP);
ic=cell(1,NTYP);
times=cell(1,NTYP);
rare_site=cell(1,NTYP);
POS_new_pri=cell(1,NTYP);
POSION_new_pri=cell(1,NTYP);
for NT=1:NTYP
    POSION_new_cry = A\POSION_new{NT};
    POSION_new_cry = mod(POSION_new_cry,1);
    POSION_new_cry(POSION_new_cry<2E-2)=...
        POSION_new_cry(POSION_new_cry<2E-2)+1;
    
    N_sc = N_sc_ex*diag(N_sc_base);
    MAXC = ceil(det(N_sc)^(1/3));
    POSION_new_pri{NT} =  N_sc'*POSION_new_cry;
    POSION_new_pri{NT} = mod(POSION_new_pri{NT},1);
    POSION_new_pri{NT}(POSION_new_pri{NT}<MAXC*2E-2)=...
        POSION_new_pri{NT}(POSION_new_pri{NT}<MAXC*2E-2)+1;
    
    [~,ia{NT},ic{NT}] = ...
        uniquetol(POSION_new_pri{NT}',5E-2,'ByRows',true,'DataScale',1);
    
    times{NT}=accumarray(ic{NT},1).';
    rare_site{NT}=find(times{NT}<1);
    POS_new_pri{NT}=zeros(3,length(ia{NT})-length(rare_site{NT}));
    % Rarely occupied, or just an overfitted position from the orbital
    
    j=0;
    for i=1:length(ia{NT})
        if ~ismember(i,rare_site{NT})
            j=j+1;
            POS_new_pri{NT}(:,j) = mean(POSION_new_pri{NT}(:,ic{NT}==i),2);
        end
    end
    
end

%
FID=fopen([ROOT_DIR,'/BEAD_1_sym/HEAD_1']);
fread(FID,[3,3],'double');
fread(FID,1,'double');
fread(FID,1,'double');
fread(FID,1,'int');fread(FID,1,'int');
fread(FID,1,'int');NITYP_EQ=zeros(1,NTYP);
for NT=1:NTYP
    NITYP_EQ(NT)=fread(FID,1,'int');
end
fclose(FID);
POSION_EQ=cell(1,NTYP);
for NT=1:NTYP
    POSS=fopen(['../../BEAD_1_sym/POS_1_',int2str(NT)]);
    POSION_EQ{NT}=fread(POSS,[3,NITYP_EQ(NT)],'double');
    fclose(POSS);
    
    POSION_EQ{NT} = A*POSION_EQ{NT};
end

SYM_tol=2E-2;
POS_EQ_pri=cell(1,NTYP);
for NT=1:NTYP
    POSION_EQ_cry = A\POSION_EQ{NT};
    POSION_EQ_pri = N_sc'*POSION_EQ_cry;
    POSION_EQ_pri = mod(POSION_EQ_pri,1);
    POSION_EQ_pri(POSION_EQ_pri<SYM_tol)=...
        POSION_EQ_pri(POSION_EQ_pri<SYM_tol)+1;
    POS_EQ_pri{NT}=uniquetol(POSION_EQ_pri',5E-2,'ByRows',true)';
end

save('POS_gauss','POSION_new_pri','POS_new_pri','POS_EQ_pri',...
    'ia','ic','times','rare_site','-append');

%%
if L_SYM
    %% Determine the symmetry operations
    % Here we assume the same symmetry operations
    % of the original primitivie cell is obeyed
    % This is at least true for Li2MgH16 & maybe Cu2Se
    A_prim = A/N_sc'; 

    NITYP_prim = zeros(1,NTYP);
    NITYP_new = zeros(1,NTYP);
    for NT=1:NTYP
        NITYP_prim(NT)=size(POS_EQ_pri{NT},2);
        NITYP_new(NT)=size(POS_new_pri{NT},2);
    end
    
    [IBRAV,N_SYM_TRUE,SYM_OP_TRUE,TRANS_ROT]...
        = Symmetry(A_prim,POS_EQ_pri,NTYP,NITYP_prim);
    
    %% Now symmetrize the positions
    load('POS_gauss.mat','POS_new_pri')
    SYM_tol=2E-2; % The tolerance should be released.
    
    for NT=1:NTYP
        POS_new_pri{NT}(abs(POS_new_pri{NT}-1)<SYM_tol) ...
            = POS_new_pri{NT}(abs(POS_new_pri{NT}-1)<SYM_tol) - 1;
        
        POS_EQ_pri{NT}(abs(POS_EQ_pri{NT}-1)<SYM_tol) ...
            = POS_EQ_pri{NT}(abs(POS_EQ_pri{NT}-1)<SYM_tol) - 1;
    end
    
    % There may be some overall translation
    % Since the mass center of initial (defective) configuration
    % could be different from the perfect one.
    DP=zeros(3,sum(NITYP_new));
    NIS=0;
    for NT=1:NTYP
        [EQ_idx,idx]=ismembertol(...
            POS_EQ_pri{NT}',POS_new_pri{NT}',SYM_tol,'ByRows',true);
        POS_new_pri{NT}=POS_new_pri{NT}(:,idx(idx~=0));
        DP(:,NIS+1:NIS+NITYP_new(NT)) = POS_new_pri{NT} - POS_EQ_pri{NT}(:,EQ_idx);
        NIS=NIS+NITYP_new(NT);
    end
    DP = mean(DP,2);
    
    for NT=1:NTYP
        POS_new_pri{NT}=POS_new_pri{NT}-DP;
    end
    
    POS_SYM=POS_new_pri;
    for NT=1:NTYP
        for NS=2:N_SYM_TRUE
            POS_tmp = mod(SYM_OP_TRUE(:,:,NS)' * POS_new_pri{NT},1);
            POS_tmp=mod(POS_tmp-TRANS_ROT(:,NS),1);
            POS_tmp(abs(POS_tmp-1)<SYM_tol) ...
                = POS_tmp(abs(POS_tmp-1)<SYM_tol) - 1;
            
            [~,locB]=ismembertol(POS_new_pri{NT}',POS_tmp',SYM_tol,...
                'DataScale', 1,'ByRows',true);
            POS_SYM{NT}=POS_SYM{NT}+POS_tmp(:,locB);
        end
        POS_SYM{NT}=POS_SYM{NT}/N_SYM_TRUE;
    end
    POSION_new=POS_SYM;
    
    for NT=1:NTYP
        POSION_new{NT}=POSION_new{NT}+DP;
        POS_new_pri{NT}=POS_new_pri{NT}+DP;
    end
    
else
    load('POS_gauss.mat','POS_new_pri')
    SYM_tol=2E-2; % The tolerance should be released.
    for NT=1:NTYP
        POS_new_pri{NT}(abs(POS_new_pri{NT}-1)<SYM_tol) ...
            = POS_new_pri{NT}(abs(POS_new_pri{NT}-1)<SYM_tol) - 1;
    end
    
    POSION_new=POS_new_pri;
end
save('../POSION_new.mat','POSION_new')
