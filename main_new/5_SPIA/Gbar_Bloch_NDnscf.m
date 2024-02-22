%% parameter
addpath('../')
addpath('../0_Public')
parameter;
gpuDevice(jobid);
maxNumCompThreads(2);
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
NBANDS = ceil(NELE/2 * NBANDS_mul);

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

if L_Liquid
    for NT=1:NTYP
	NITYP_EQ(NT)=0;
        POSION_EQ{NT}=[];
    end
end

%[~,idx]=ismembertol(VKPT_nosym',VKPT_sym','ByRows',true);
%idx=find(~idx);
%VKPT_cal=VKPT_nosym(:,idx);
%NKPTS_cal=length(idx);
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

if LBloch
    %% Load in Bloch basis
    WAVE_gamma=zeros(length(BASIS.G_INDEX_sc));
    Ni=0;
    for NK_s=1:BASIS.NKPTS_red
        for i=Ni+1:Ni+length(BASIS.E_BAND{NK_s})
            WAVE_gamma(BASIS.G_kp{NK_s},BASIS.Band_idx(i))=BASIS.WAVE{NK_s}(:,i-Ni);
        end
        Ni=Ni+length(BASIS.G_kp{NK_s});
    end
    clear BASIS
    BASIS_k=gpuArray(single(WAVE_gamma(:,1:NBANDS)));
else
    % For systems where Bloch bases can diviate
    % heavily from the initial one
    % we must set NBANDS=NPL
    NBANDS=NPL;
    BASIS_k=0;
end
fprintf(FID_OUT,'Bases number: %d\n',NBANDS);
clear WAVE_gamma

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
    if ~L_Liquid
        for NT=1:NTYP
            POSION1=gpuArray(POSION_EQ{NT});
            CREXP=single(exp(CITPI*G_INDEX(1:3,1:NPL).'*POSION1));
            CREXP=conj(CREXP)*CREXP.';
            % For \hat{T}_c=|\tilde{p}_i^a><\tilde{\phi}_i^a|...
            TRANS_EQ=TRANS_EQ+CREXP.*TRANS_CORE{NT};
        end
        TRANS_EQ=gather(TRANS_EQ);
        clear CREXP
    end
    
    %%
    GREEN = zeros(NBANDS,NBANDS,'single','gpuArray');
    %for NC=jobid*NSKIP+1+NSTART:NSKIP*NGPU:NCONF
    NCAC = length(N_BEGIN+1:NSKIP:N_END);
    CONF_CAC = N_BEGIN+1:NSKIP:N_END;
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

            FID=fopen([ROOT_DIR,'/BEAD_',NC_c,'/HEAD_1']);
            fread(FID,[3,3],'double');
            fread(FID,2,'double');
            fread(FID,3+NTYP,'int');
            NCPU=fread(FID,1,'int');
            fclose(FID);

            FHAM = Green_c(ROOT_DIR,NC_c,NB_c,ISP_c,...
                NTYP,NITYP,LMAX,LMMAXC,CH0,CH1,CH2,CH3,LPS,...
                NPL,NPLWV,NGX,NGY,NGZ,NCPU,...
                QPROJ_ORI,TRANS_CORE,...
                POSION,G_INDEX,G,HSQDTM,CITPI,...
                POSION_EQ, TRANS_EQ, BASIS_k, A, B, ...
                CQIJ3_fun, r, kB, TEMP, KOMEGA,EFERMI_Av, LBloch);
            %FHAM=BASIS_k'*TdT0'...
            %    *(BASIS_trial*diag(1./((1i*ETA+EFERMI_Av)...
            %    -E_BAND_now))*BASIS_trial')...
            %    *TdT0*BASIS_k;
            
            %clear E_BAND_now
            GREEN=GREEN+FHAM;
            clear FHAM
            t2_end=toc(t2);
            fprintf(FID_OUT,'Loop: %6.4f\n',t2_end);
        end
    end
    
    GREEN=gather(GREEN);
    
    fn=['../GREEN_',int2str(NK),'_',int2str(jobid),'.mat'];
    save(fn,'GREEN','-v7.3')
    
end
t3_end=toc(t3);
fprintf(FID_OUT,'K-loop: %6.4f\n',t3_end);
fprintf(FID_OUT,' \n');
%end

fclose(FID_OUT);
clear
