%% parameter
addpath('../')
addpath('../0_Public')
parameter;
gpuDevice(jobid);
maxNumCompThreads(6);
if NK==NK_BEGIN
    FID_OUT=fopen(['out_C',int2str(NCL),'_G',int2str(jobid),'.file'],'w');
else
    FID_OUT=fopen(['out_C',int2str(NCL),'_G',int2str(jobid),'.file'],'a');
end

%% Read in informations
[A,B,OMEGA,...
    NKPTS,ISPIN,NTYP,NITYP,NCPU,...
    NELE,VKPT,NELE_TYP,POMASS]=...
    RD_HEAD(ROOT_DIR);
%NBANDS = max(ceil(NELE * 0.6),ceil(NELE/2+sum(NITYP)/2));
NBANDS = ceil(NELE/2 * NBANDS_mul_Ef);

%% PAW datasets
[NMAX,RG,SIMPI,LLMAX,WAE_PS,WAE,WPS,...
    LMAX,LMMAXC,CH0,CH1,CH2,CH3,LPS]=RD_PAW(ROOT_DIR,NTYP);

PSPNL = cell(1,NTYP); NPSNL = 100;
PSMAXN = zeros(1,NTYP);
for NT=1:NTYP
    FID=fopen([ROOT_DIR,'/BEAD_1_sym/PSPNL_1_',int2str(NT)]);
    PSMAXN(NT) = fread(FID,1,'double');
    PSPNL{NT}=fread(FID,[NPSNL+1,LMAX(NT)],'double');
    fclose(FID);
end

%% Transformation operator correction (Old version)
if L_OVERL_OLD
    [QPAW3,CQIJ]=OVERL_CORE_1(ROOT_DIR,NTYP,NITYP,LMMAXC,CH0,CH1,CH2,CH3,WAE,WPS,SIMPI);
end

%% Load in equilibrium positions
POSION_EQ=cell(1,NTYP);
for NT=1:NTYP
    % Phase's space cost is huge. Instead read positions.
    POSS=fopen([ROOT_DIR,'/BEAD_1_sym/POS_1_',int2str(NT)]);
    POSION_EQ{NT}=fread(POSS,[3,NITYP(NT)],'double');
    fclose(POSS);
end

% Load in shifted overlap kernal
load('CQIJ3.mat','CQIJ3_div','r','Nr')
% Calculate the lower part of the kernal
for NT1=1:NTYP
    for NT2=1:NT1-1
        r{NT1,NT2}=r{NT2,NT1};
        for i=1:length(r{NT1,NT2})
            R=[0,-r{NT1,NT2}(i),0,r{NT1,NT2}(i)];
            Rotate1=Rotate_SH(R,LMAX,LMMAXC,LPS,NT1);
            Rotate2=Rotate_SH(R,LMAX,LMMAXC,LPS,NT2);
            CQIJ3_div{NT1,NT2}(:,:,i)=...
                Rotate1*CQIJ3_div{NT2,NT1}(:,:,i).'*(Rotate2.');
        end
    end
end

CQIJ3_fun = cell(NTYP,NTYP);
for NT1=1:NTYP
    for NT2=1:NTYP
        CQIJ3_fun{NT1,NT2}=cell(LMMAXC(NT1),LMMAXC(NT2));
        for i=1:LMMAXC(NT1)
            for j=1:LMMAXC(NT2)
                CQIJ3_fun{NT1,NT2}{i,j}=@(x) interp1(r{NT1,NT2},squeeze(CQIJ3_div{NT1,NT2}(i,j,:)),x,'spline');
            end
        end
    end
end

%% Determine the denser k-mesh, and find irreducible k-points
% k-mesh of the supercell
%[VKPT_sym,NKPTS_sym,VKPT_f,~,SYM_OP_TRUE]...
%    =Kpoints(N_k_nscf,A,NTYP,NITYP,POSION_EQ,1);

%SYM_OP_min(:,:,1) = eye(3);
%SYM_OP_min(:,:,2) = -eye(3); % minimal symmetry with inversion symmetry
%[VKPT,NKPTS_nosym,~]...
%    =Kpoints_gsym(VKPT_f,A,SYM_OP_min);

%[~,VKPT_sym_idx] = ismembertol(...
%    VKPT_sym(:,1:3),VKPT(:,1:3),'ByRows',true);

%k_weight = VKPT_sym(:,4)';
%VKPT_sym = VKPT_sym(:,1:3)';
%VKPT = VKPT(:,1:3)';

%% For Non-diagonal supercell, determine the irreducible k-mesh from the diagonal one
%% Information about the primitive cell
FID=fopen([ROOT_DIR,'/BEAD_1_primitive/HEAD_1']);
A_prim=fread(FID,[3,3],'double');
fread(FID,2,'double');
fread(FID,2,'int');
NTYP=fread(FID,1,'int');NITYP_prim=zeros(1,NTYP);
for NT=1:NTYP
    NITYP_prim(NT)=fread(FID,1,'int');
end
fclose(FID);
B_prim=2*pi*inv(A_prim);
OMEGA_prim=det(A);

% Primitive cell ion position
if L_POS_new
    load('../POSION_new.mat')
    POSION_prim=POSION_new;
else
    % Primitive cell ion position
    POSION_prim=cell(1,NTYP);
    for NT=1:NTYP
        POSS=fopen([ROOT_DIR,'/BEAD_1_primitive/POS_1_',int2str(NT)]);
        POSION_prim{NT}=fread(POSS,[3,NITYP_prim(NT)],'double');
        fclose(POSS);
    end
end

%% Diagonal supercell (base)
A_base = A_prim*diag(N_sc_base);
B_base = 2*pi*inv(A_base);

NITYP_base = NITYP_prim * det(diag(N_sc_base));
POSION_base = cell(1,NTYP);
NX = 0:N_sc_base(1)-1;
NY = 0:N_sc_base(2)-1;
NZ = 0:N_sc_base(3)-1;
POS_Nsc = zeros(det(diag(N_sc_base)),3);
POS_Nsc(:,1)= reshape(repmat(NX,[1,N_sc_base(2)*N_sc_base(3)]),1,[]);
POS_Nsc(:,2)= reshape(repmat(NY,[N_sc_base(1),N_sc_base(3)]),1,[]);
POS_Nsc(:,3)= reshape(repmat(NZ,[N_sc_base(1)*N_sc_base(2),1]),1,[]);
for NT=1:NTYP
    POSION_base{NT}(1,:) = reshape(POSION_prim{NT}(1,:)+POS_Nsc(:,1),1,[])/N_sc_base(1);
    POSION_base{NT}(2,:) = reshape(POSION_prim{NT}(2,:)+POS_Nsc(:,2),1,[])/N_sc_base(2);
    POSION_base{NT}(3,:) = reshape(POSION_prim{NT}(3,:)+POS_Nsc(:,3),1,[])/N_sc_base(3);
end

%% Determine the denser full k-mesh of the basic diagonal supercell
eps=1E-5;

if L_path
    KP=cell(size(VK_hsym,1)-1,1);
    for i=1:size(VK_hsym,1)-1
        KP{i}=[...
            linspace(VK_hsym(i,1),VK_hsym(i+1,1),NK_line(i)+1).',...
            linspace(VK_hsym(i,2),VK_hsym(i+1,2),NK_line(i)+1).',...
            linspace(VK_hsym(i,3),VK_hsym(i+1,3),NK_line(i)+1).'];
        if i>1
            KP{i}=KP{i}(2:end,:);
        end
    end
    VKPTS_full = cell2mat(KP).';
    NKPTS_full = size(VKPTS_full,2);
    VKPTS_full(VKPTS_full>0.5+eps) =...
        VKPTS_full(VKPTS_full>0.5+eps) - 1;
else
    NKPTS_full = N_k_nscf(1)*N_k_nscf(2)*N_k_nscf(3);
    
    K1 = (0:N_k_nscf(1)-1)/N_k_nscf(1);
    K2 = (0:N_k_nscf(2)-1)/N_k_nscf(2);
    K3 = (0:N_k_nscf(3)-1)/N_k_nscf(3);
    
    clear VKPTS_full
    VKPTS_full(1,:) = reshape(repmat(K1,1,N_k_nscf(2)*N_k_nscf(3)),1,[]);
    VKPTS_full(2,:) = reshape(repmat(K2,N_k_nscf(1),N_k_nscf(3)),1,[]);
    VKPTS_full(3,:) = reshape(repmat(K3,N_k_nscf(1)*N_k_nscf(2),1),1,[]);
    VKPTS_full(VKPTS_full>0.5+eps) = VKPTS_full(VKPTS_full>0.5+eps) - 1;
end
%% Transform to the non-diagonal supercell
eps=1E-5;
VKPT_sc_ex = N_sc_ex*VKPTS_full;
VKPT_sc_ex = mod(VKPT_sc_ex,1);
VKPT_sc_ex(VKPT_sc_ex>1-eps) = VKPT_sc_ex(VKPT_sc_ex>1-eps) -1;
VKPT_sc_ex = [VKPT_sc_ex(3,:);VKPT_sc_ex(2,:);VKPT_sc_ex(1,:)];
[VKPT_NDsc,~,ND_sym_idx] = uniquetol(VKPT_sc_ex','ByRows',true);
VKPT_NDsc(VKPT_NDsc>0.5+eps) = VKPT_NDsc(VKPT_NDsc>0.5+eps) -1;
VKPT_NDsc = VKPT_NDsc';
VKPT_NDsc = [VKPT_NDsc(3,:);VKPT_NDsc(2,:);VKPT_NDsc(1,:)];

if L_POS_new
    N_sc=N_sc_ex*diag(N_sc_base);
    MAXC=det(N_sc);
    NX=2*MAXC+1;NY=2*MAXC+1;NZ=2*MAXC+1;
    
    POSION_new = cell(1,NTYP);
    NITYP_new = zeros(1,NTYP);
    for NT=1:NTYP
        POSION_ori = POSION_prim{NT};
        
        PX=reshape(repmat(POSION_ori(1,:).',1,NX)+[0:MAXC,-MAXC:-1],1,[]);
        PX=reshape(repmat(...
            reshape(PX...
            ,1,[],NX),1,1,NY*NZ),1,[]);
        
        PY=reshape(repmat(POSION_ori(2,:).',1,NY)+[0:MAXC,-MAXC:-1],1,[]);
        PY=reshape(repmat(...
            reshape(PY...
            ,1,[],NY),1,NX,NZ),1,[]);
        
        PZ=reshape(repmat(POSION_ori(3,:).',1,NZ)+[0:MAXC,-MAXC:-1],1,[]);
        PZ=reshape(repmat(...
            reshape(PZ...
            ,1,[],NZ),1,NX*NY,1),1,[]);
        
        POS_tmp = [PX;PY;PZ];
        POS_tmp = N_sc'\POS_tmp;
        POS_tmp = mod(POS_tmp,1);
        POS_tmp(POS_tmp>1-eps)=POS_tmp(POS_tmp>1-eps)-1;
        POS_tmp = uniquetol(POS_tmp','ByRows',true)';
        POSION_new{NT} = POS_tmp;
        NITYP_new(NT) = length(POS_tmp);
    end
    POSION_EQ = POSION_new;
    NITYP_EQ = NITYP_new;
else
    NITYP_EQ = NITYP;
end

[IBRAV,N_SYM_TRUE,SYM_OP_TRUE,TRANS_ROT]...
    = Symmetry(A,POSION_EQ,NTYP,NITYP_EQ);

[VKPT_sym,NKPTS_sym,~]...
    =Kpoints_gsym(VKPT_NDsc,A,SYM_OP_TRUE);

if ~L_IBZ
    SYM_OP_min(:,:,1) = eye(3);
    SYM_OP_min(:,:,2) = -eye(3);
    [VKPT_sym,NKPTS_sym,~]...
        =Kpoints_gsym(VKPT_NDsc,A,SYM_OP_min);
end

k_weight=VKPT_sym(:,4)';
VKPT_sym=VKPT_sym(:,1:3)';

%% Calculation
%NKPTS=1;
%for NK=1:NKPTS_sym
%for NK=4:4
t3=tic;
%NK_c=int2str(VKPT_sym_idx(NK));
%NK_sym=VKPT_sym_idx(NK);
NK_c=int2str(NK);
NK_sym=NK;
%% Read in necessary informations:
%% Mesh informations
B=gpuArray(B);
NC_c='1';
[G_INDEX,G,G2,NPL,NPLWV,NGX,NGY,NGZ]=...
    RD_INDEX(ROOT_DIR,NC_c,'1',VKPT_sym(:,NK_sym),B,ENCUT);
G_INDEX=gpuArray(G_INDEX);

fprintf(FID_OUT,'Bases number: %d\n',NPL);

% For solving the eigenvalue problem,
% begin with a random trial function.
%NBANDS=NPL;
WAVE_conf = single(rand(NPL,NBANDS) ...
    + 1i*(0.2*rand(NPL,NBANDS)-0.1));
if sum(abs(VKPT_sym(:,NK_sym)).^2)<1E-3
    WAVE_conf(1,:)=rand(1,NBANDS,'single');
end
weight = gather(G(4,:).'*AUTOA^2+1);
WAVE_conf = WAVE_conf./weight;

WAVE_norm = sqrt(sum(abs(WAVE_conf).^2,1));
WAVE_conf = WAVE_conf./WAVE_norm;

%%
for ISP=1:ISPIN
    ISP_c=int2str(ISP);
    %% Overlap matrix and transformation matrix related
    NC_c='1';
    [QPROJ_ORI,TRANS_CORE]=...
        OVERL_CORE_nscf(PSMAXN,PSPNL,NTYP,...
        LMAX,LMMAXC,RG,NMAX,...
        NPL,CH0,CH1,CH2,CH3,LPS,...
        WAE_PS,WAE,WPS,G2,OMEGA);
    
    %% Generate the Transformation matrix for equilibrium configuration
    TRANS_EQ=eye(NPL,'single','gpuArray');
    for NT=1:NTYP
        POSION1=gpuArray(POSION_EQ{NT});
        CREXP=single(exp(CITPI*G_INDEX(1:3,1:NPL).'*POSION1));
        CREXP=conj(CREXP)*CREXP.';
        % For \hat{T}_c=|\tilde{p}_i^a><\tilde{\phi}_i^a|...
        TRANS_EQ=TRANS_EQ+CREXP.*TRANS_CORE{NT};
    end
    TRANS_EQ=gather(TRANS_EQ);
    clear CREXP
    
    %%
    %GREEN = zeros(NBANDS,NBANDS,'single','gpuArray');
    %for NC=jobid*NSKIP+1+NSTART:NSKIP*NGPU:NCONF
    Spec_av = cell(BASIS.NKPTS_red,length(ETA));
    Spec_av2 = cell(BASIS.NKPTS_red,length(ETA));
    
    for NE=1:length(ETA)
        for NK_s=1:BASIS.NKPTS_red
            Spec_av{NK_s,NE}=zeros(size(E_sample),'like',gather(B));
            Spec_av2{NK_s,NE}=zeros(size(E_sample),'like',gather(B));
        end
    end
    
    NCAC = length(N_BEGIN+1:NSKIP:N_END);
    CONF_CAC = N_BEGIN+1:NSKIP:N_END;
    E_t = zeros(NBANDS,NBEAD,NCAC);
    for NCs=1:NCAC
        %        for NC=1:NSKIP*NGPU:NCONF
        NC = CONF_CAC(NCs);
        NC_c=int2str(NC);
        for NB=1:NBEAD
            fprintf(FID_OUT,'Kpoint %d, Configuration %d, Bead %d\n',[NK,NC,NB]);
            NB_c=int2str(NB);
            %% Positions of ions
            POSION=cell(1,NTYP);
            for NT=1:NTYP
                % Phase's space cost is huge. Instead read positions.
                POSS=fopen([ROOT_DIR,'/BEAD_',NC_c,'/POS_',NB_c,'_',int2str(NT)]);
                POSION{NT}=fread(POSS,[3,NITYP(NT)],'double');
                fclose(POSS);
            end
            
            t2=tic;
            % V_proj = |<k_s+\tilde{G_s}|S(R0,R)|psi_NK>|^2
            FID=fopen([ROOT_DIR,'/BEAD_',NC_c,'/HEAD_1']);
            fread(FID,[3,3],'double');
            fread(FID,2,'double');
            fread(FID,3+NTYP,'int');
            NCPU=fread(FID,1,'int');
            fclose(FID);

            % E_conf being the eigenvalues
            [WAVE_conf,E_conf,V_proj,N_notconv] = ...
                Green_c(ROOT_DIR,NC_c,NB_c,ISP_c,...
                NTYP,NITYP,LMAX,LMMAXC,CH0,CH1,CH2,CH3,LPS,...
                NPL,NPLWV,NGX,NGY,NGZ,NCPU,...
                QPROJ_ORI,TRANS_CORE,...
                POSION,G_INDEX,G,HSQDTM,CITPI,...
                POSION_EQ, TRANS_EQ, A, B, ...
                CQIJ3_fun, r, RYTOEV, WAVE_conf, NBANDS, NELE);
            if N_notconv~=0
                fprintf(FID_OUT,'Not converged');
                error('Not converged');
            end
            
            E_t(:,NB,NCs)=gather(E_conf);
            E = E_sample;
            for NE=1:length(ETA)
                G_conf = ETA(NE)./((E_conf-E).^2+ETA(NE)^2)/pi;
                
                % Spectral function should be
                % sum_{NK,G} V_conf(k_p+G_p) delta(E-E_{NK})
                for NK_s=1:BASIS.NKPTS_red
                    Spec = sum(sum(V_proj(BASIS.G_kp{NK_s},:),1).'.*G_conf,1);
                    Spec_av{NK_s,NE}=...
                        gather(Spec_av{NK_s,NE}+Spec);
                end
                
                FID=fopen([ROOT_DIR,'/BEAD_',int2str(NC),'/HEAD_',int2str(NB)]);
                fread(FID,[3,3],'double');
                EFERMI=fread(FID,1,'double');
                fclose(FID);
                
                % With respect to instantaneous Fermi energy
                G_conf = ETA(NE)./((E_conf-EFERMI_Av-E).^2+ETA(NE)^2)/pi;
                for NK_s=1:BASIS.NKPTS_red
                    Spec = sum(sum(V_proj(BASIS.G_kp{NK_s},:),1).'.*G_conf,1);
                    Spec_av2{NK_s,NE}=...
                        gather(Spec_av{NK_s,NE}+Spec);
                end
                
            end
            t2_end=toc(t2);
            fprintf(FID_OUT,'Loop: %6.4f\n',t2_end);
        end
    end
    
    %Spec_av=gather(Spec_av);
    
    fn=['../Spec_',int2str(NK),'_',int2str(jobid),'.mat'];
    save(fn,'Spec_av','Spec_av2','-v7.3')
    save(['../Eigen_',int2str(NK),'_',int2str(jobid),'.mat'],'E_t','-v7.3')
end
t3_end=toc(t3);
fprintf(FID_OUT,'K-loop: %6.4f\n',t3_end);
fprintf(FID_OUT,' \n');
%end

fclose(FID_OUT);
clear
