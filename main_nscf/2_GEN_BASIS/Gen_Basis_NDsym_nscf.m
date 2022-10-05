% Generate k-points in terms of the primitive cell
addpath('../')
addpath('../0_Public')
parameter;
gpuDevice(1);
maxNumCompThreads(12);
%pool = parpool('local',12,'IdleTimeout',1200);

%% Lattice informations of the primitive cell
FID=fopen([ROOT_DIR,'/BEAD_1_primitive/HEAD_1']);
A=fread(FID,[3,3],'double');
EFERMI=fread(FID,1,'double');
TOTEN=fread(FID,1,'double');
NKPTS=fread(FID,1,'int');ISPIN=fread(FID,1,'int');
NTYP=fread(FID,1,'int');NITYP=zeros(1,NTYP);
for NT=1:NTYP
    NITYP(NT)=fread(FID,1,'int');
end
NCPU=fread(FID,1,'int');
fclose(FID);
if isempty(NCPU)
    NCPU=1;
end
B=2*pi*inv(A);
OMEGA=det(A);

%% Primitive cell ion position and the k-mesh
POSION=cell(1,NTYP);
for NT=1:NTYP
    % Phase's space cost is huge. Instead read positions.
    POSS=fopen([ROOT_DIR,'/BEAD_1_primitive/POS_1_',int2str(NT)]);
    POSION{NT}=fread(POSS,[3,NITYP(NT)],'double');
    fclose(POSS);
end

% k-mesh with only time-reversal symmetry
[VKPTS_f,NKPTS_f]=Full_K(N_k_nscf.*N_sc_base,1);

SYM_OP_min(:,:,1) = eye(3);
SYM_OP_min(:,:,2) = -eye(3); % minimal symmetry with inversion symmetry
[VKPT,NKPTS,VKPTS_nosym_idx]...
    =Kpoints_gsym(VKPTS_f,A,SYM_OP_min);

k_weight = VKPT(:,4)';
VKPT = VKPT(:,1:3)';

%% G_INDEX
FID=fopen([ROOT_DIR,'/BEAD_1_primitive/INDEX_1_1_1']);
NGX=fread(FID,1,'int');NGY=fread(FID,1,'int');NGZ=fread(FID,1,'int');
fclose(FID);

NPLWV=NGX*NGY*NGZ;

G_INDEX_f=zeros(3,NPLWV);
GB3=[0:NGZ/2,-NGZ/2+1:-1];
GB2=[0:NGY/2,-NGY/2+1:-1];
GB1=[0:NGX/2,-NGX/2+1:-1];

G_INDEX_f(1,:)=repmat(GB1,1,NGY*NGZ);
G_INDEX_f(2,:)=reshape(repmat(GB2,NGX,NGZ),1,NPLWV);
G_INDEX_f(3,:)=reshape(repmat(GB3,NGX*NGY,1),1,NPLWV);

NPL_ss=zeros(NKPTS,1); G_INDEX_prim=cell(1,NKPTS);
for NK=1:NKPTS
    G_INDEX_tmp = G_INDEX_f + VKPT(:,NK);
    
    G_f(1:3,:)=gpuArray(B).'*G_INDEX_tmp(1:3,:);
    G_f(4,:)=G_f(1,:).^2+G_f(2,:).^2+G_f(3,:).^2;
    
    idx = find(HSQDTM*G_f(4,:)<ENCUT);
    G_INDEX_prim{NK} = gather([G_INDEX_tmp(:,idx);idx]);
    NPL_ss(NK) = length(G_INDEX_prim{NK});
end

%% Diagonal supercell (base)
A_base = diag(N_sc_base)*A;
B_base = 2*pi*inv(A_base);

NITYP_base = NITYP * det(diag(N_sc_base));
POSION_base = cell(1,NTYP);
NX = 0:N_sc_base(1)-1;
NY = 0:N_sc_base(2)-1;
NZ = 0:N_sc_base(3)-1;
POS_Nsc = zeros(det(diag(N_sc_base)),3);
POS_Nsc(:,1)= reshape(repmat(NX,[1,N_sc_base(2)*N_sc_base(3)]),1,[]);
POS_Nsc(:,2)= reshape(repmat(NY,[N_sc_base(1),N_sc_base(3)]),1,[]);
POS_Nsc(:,3)= reshape(repmat(NZ,[N_sc_base(1)*N_sc_base(2),1]),1,[]);
for NT=1:NTYP
    POSION_base{NT}(1,:) = reshape(POSION{NT}(1,:)+POS_Nsc(:,1),1,[])/N_sc_base(1);
    POSION_base{NT}(2,:) = reshape(POSION{NT}(2,:)+POS_Nsc(:,2),1,[])/N_sc_base(2);
    POSION_base{NT}(3,:) = reshape(POSION{NT}(3,:)+POS_Nsc(:,3),1,[])/N_sc_base(3);
end

%% Determine the denser full k-mesh of the basic diagonal supercell
eps=1E-5;
NKPTS_full = N_k_nscf(1)*N_k_nscf(2)*N_k_nscf(3);

K1 = (0:N_k_nscf(1)-1)/N_k_nscf(1);
K2 = (0:N_k_nscf(2)-1)/N_k_nscf(2);
K3 = (0:N_k_nscf(3)-1)/N_k_nscf(3);

clear VKPTS_full
VKPTS_full(1,:) = reshape(repmat(K1,1,N_k_nscf(2)*N_k_nscf(3)),1,[]);
VKPTS_full(2,:) = reshape(repmat(K2,N_k_nscf(1),N_k_nscf(3)),1,[]);
VKPTS_full(3,:) = reshape(repmat(K3,N_k_nscf(1)*N_k_nscf(2),1),1,[]);
VKPTS_full(VKPTS_full>0.5+eps) = VKPTS_full(VKPTS_full>0.5+eps) - 1;

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

[A_sc,B_sc,OMEGA_sc,...
    ~,~,NTYP_sc,NITYP_sc,...
    NELE_sc,NBANDS_sc,~]=RD_HEAD(ROOT_DIR,NELE_TYP,NBANDS_mul);

POSION_sc=cell(1,NTYP_sc);
for NT=1:NTYP
    POSS=fopen([ROOT_DIR,'/BEAD_1/POS_1_',int2str(NT)]);
    POSION_sc{NT}=fread(POSS,[3,NITYP_sc(NT)],'double');
    fclose(POSS);
end

[IBRAV,N_SYM_TRUE,SYM_OP_TRUE,TRANS_ROT]...
    = Symmetry(A_sc,POSION_sc,NTYP,NITYP);

[VKPT_sym,NKPTS_sym,~]...
    =Kpoints_gsym(VKPT_NDsc,A_sc,SYM_OP_TRUE);

k_weight=VKPT_sym(:,4)';
VKPT_sym=VKPT_sym(:,1:3)';

%% Load in shifted overlap kernal
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

if L_OVERL_OLD
    [QPAW3,CQIJ]=OVERL_CORE_1(ROOT_DIR,NTYP,NITYP,LMMAXC,CH0,CH1,CH2,CH3,WAE,WPS,SIMPI);
end

load('CQIJ3.mat')
CQIJ3_fun = cell(NTYP,NTYP);
for NT1=1:NTYP
    for NT2=NT1:NTYP
        CQIJ3_fun{NT1,NT2}=cell(LMMAXC(NT1),LMMAXC(NT2));
        for i=1:LMMAXC(NT1)
            for j=1:LMMAXC(NT2)
                CQIJ3_fun{NT1,NT2}{i,j}=@(x) interp1(r{NT1,NT2},squeeze(CQIJ3_div{NT1,NT2}(i,j,:)),x,'spline');
            end
        end
    end
end

%% Correspondence to primitive cell k-points
FID=fopen([ROOT_DIR,'/BEAD_1/INDEX_1_1_1']);
NGX_sc=fread(FID,1,'int');NGY_sc=fread(FID,1,'int');NGZ_sc=fread(FID,1,'int');
fclose(FID);

NPLWV_sc=NGX_sc*NGY_sc*NGZ_sc;
    
G_INDEX_sc_f=zeros(3,NPLWV_sc);
GB3_sc=[0:NGZ_sc/2,-NGZ_sc/2+1:-1];
GB2_sc=[0:NGY_sc/2,-NGY_sc/2+1:-1];
GB1_sc=[0:NGX_sc/2,-NGX_sc/2+1:-1];

G_INDEX_sc_f(1,:)=repmat(GB1_sc,1,NGY_sc*NGZ_sc);
G_INDEX_sc_f(2,:)=reshape(repmat(GB2_sc,NGX_sc,NGZ_sc),1,NPLWV_sc);
G_INDEX_sc_f(3,:)=reshape(repmat(GB3_sc,NGX_sc*NGY_sc,1),1,NPLWV_sc);

FID=fopen([ROOT_DIR,'/BEAD_1_sym/HEAD_1']);
A_sc=fread(FID,[3,3],'double');
B_sc=2*pi*inv(A_sc);
fclose(FID);

for NK_sc=1:NKPTS_sym
    tic;
    NK_sc_c=int2str(NK_sc);
    %% G_INDEX_sc
    G_INDEX_tmp = G_INDEX_sc_f + VKPT_sym(:,NK_sc);
    
    G_f = zeros(4,NPLWV_sc,'gpuArray');
    G_f(1:3,:)=gpuArray(B_sc).'*G_INDEX_tmp(1:3,:);
    G_f(4,:)=G_f(1,:).^2+G_f(2,:).^2+G_f(3,:).^2;
    
    idx = find(HSQDTM*G_f(4,:)<ENCUT);
    G_INDEX_sc = gather([G_INDEX_tmp(:,idx);idx]);
    NPL_sc = length(G_INDEX_sc);
    
    %% Correspondence with primitive cell k-mesh
    eps = 1E-5;
    N_sc = diag(N_sc_base)*N_sc_ex;
    G_INDEX_sc_pr=mod(N_sc\G_INDEX_sc(1:3,:),1).';
    G_INDEX_sc_pr(G_INDEX_sc_pr>0.5+eps)=...
        G_INDEX_sc_pr(G_INDEX_sc_pr>0.5+eps)-1;
    % G_INDEX_sc_pr=mod((G_INDEX_sc(1:3,:).'./N_sc),1);
    VKPT_sc_pr = uniquetol(G_INDEX_sc_pr,'ByRows',true);
    VKPT_sc_pr(VKPT_sc_pr>0.5+eps)=...
        VKPT_sc_pr(VKPT_sc_pr>0.5+eps)-1;

    
    [~,VKPT_sp_ia]=ismembertol(...
        VKPT_sc_pr,VKPT.',1E-4,'ByRows',true);
    
    idx=find(VKPT_sp_ia==0);
    
    Factor_sym = ones(length(VKPT_sp_ia),1);
    Factor_sym(idx)=-1;
    
    VKPT_sps=VKPT_sc_pr;
    VKPT_sps(idx,:)=-VKPT_sc_pr(idx,:);
    VKPT_sps(VKPT_sps<-0.5+eps)=VKPT_sps(VKPT_sps<-0.5+eps)+1;
    [~,VKPT_sps_ia]=ismembertol(...
        VKPT_sps,VKPT.',1E-4,'ByRows',true);
    
    VKPTS_idx=VKPT_sps_ia;
    NKPTS_red=length(VKPT_sps_ia);
    G_kp=cell(1,NKPTS_red);
    
    NPL_sum=0;
    for NK=1:NKPTS_red
        NK_c=int2str(VKPT_sps_ia(NK));
        NPL_s=NPL_ss(VKPT_sps_ia(NK));
        G_INDEX_sym=G_INDEX_prim{VKPT_sps_ia(NK)};
        
        G_INDEX_sym(1:3,:)=Factor_sym(NK)*G_INDEX_sym(1:3,:);
        [AA,G_kpidx] = ismembertol(...
            (N_sc\G_INDEX_sc(1:3,:)).',G_INDEX_sym(1:3,:).',1E-4,'ByRows',true);
        
        [~,idx]=sort(G_kpidx(G_kpidx~=0),'ascend');
        G_kpidx=find(AA);
        G_kpidx=G_kpidx(idx);
        
        G_kp{NK}=G_kpidx;
        
        NPL_sum=NPL_sum+NPL_s;
    end

    if NPL_sum~=NPL_sc
        disp('Plane wave number inconsistent')
        return
    end
    
    WAVE=cell(1,NKPTS_red);E_BAND=cell(1,NKPTS_red);
    for NK=1:NKPTS_red
        NK_c=int2str(VKPT_sps_ia(NK));
        NC_c='1_primitive';
        G_INDEX = G_INDEX_prim{VKPT_sps_ia(NK)};
        NPL = length(G_INDEX);
        G = zeros(4,NPL,'gpuArray');
        G2 = zeros(4,NPL,'gpuArray');
        G(1:3,:)=gpuArray(B).'*G_INDEX(1:3,:);
        G(4,:)=G(1,:).^2+G(2,:).^2+G(3,:).^2;

        G2(4,:)=sqrt(G(4,:));
        G2(4,:)=max(G2(4,:),1E-12);
        G2(1:3,:)=bsxfun(@rdivide,G(1:3,:),G2(4,:));

        G=single(G);
        for ISP=1:ISPIN
%             tic;
            ISP_c=int2str(ISP);
            if L_OVERL_OLD
                [QPROJ_ORI,FOVL_CORE,TRANS_CORE,TRANS_c_CORE]=...
                    OVERL_CORE_2_OLD(ROOT_DIR,NTYP,LMAX,LMMAXC,...
                    NC_c,NK_c,ISP_c,RG,NMAX,...
                    NPL,CH0,CH1,CH2,CH3,...
                    WAE_PS,WAE,WPS,G2,CQIJ,QPAW3,OMEGA);
            else
                [QPROJ_ORI,TRANS_CORE]=...
                    OVERL_CORE_nscf(PSMAXN,PSPNL,NTYP,...
                    LMAX,LMMAXC,RG,NMAX,...
                    NPL,CH0,CH1,CH2,CH3,LPS,...
                    WAE_PS,WAE,WPS,G2,OMEGA);
            end
            
            NC_c='1_primitive';
            NB_c='1';
            if L_OVERL_OLD
                [FHAM,TRANS,FOVL]=HAMILT_OLD(ROOT_DIR,NC_c,NB_c,ISP_c,...
                    NTYP,NITYP,LMMAXC,CH0,CH1,CH2,CH3,...
                    QPROJ_ORI,NPL,NPLWV,NGX,NGY,NGZ,...
                    TRANS_CORE,TRANS_c_CORE,FOVL_CORE,...
                    POSION,G_INDEX,G,HSQDTM,CITPI);
            else
                [FHAM,TRANS]=HAMILT(ROOT_DIR,NC_c,NB_c,ISP_c,...
                    NTYP,NITYP,LMAX,LMMAXC,CH0,CH1,CH2,CH3,LPS,...
                    NPL,NPLWV,NGX,NGY,NGZ,NCPU,...
                    QPROJ_ORI,TRANS_CORE,...
                    POSION,G_INDEX,G,HSQDTM,CITPI);
                
                % Rigorous correction for overlap matrix
                FOVL=OVERL(POSION, POSION, TRANS, TRANS, A, B, ...
                    CQIJ3_fun, r, G_INDEX, CITPI, QPROJ_ORI, ...
                    NTYP,NITYP,NPL,LMAX,LMMAXC, ...
                    CH0,CH1,CH2,CH3,LPS);
                FOVL=triu(FOVL,1)+triu(FOVL,1)'+diag(real(diag(FOVL)));
            end
            
            %tic;
            F=chol(FOVL,'lower');
            FHAM = F\FHAM/F';
            FHAM=triu(FHAM,1)+triu(FHAM,1)'+diag(real(diag(FHAM)));
            [WAVE{NK},E_tmp]=eig(FHAM);
            WAVE{NK}=gather(F'\WAVE{NK});
            %toc;
            %[WAVE{NK},E_tmp]=eig(gather(FHAM),gather(FOVL));
            if Factor_sym(NK)==-1
                WAVE{NK}=conj(WAVE{NK});
            end
            E_BAND{NK}=gather(diag(real(E_tmp)));

%             toc;
        end
        
    end
    toc;
    %% Band component analysis
    ss=0;
    for NK=1:NKPTS_red
        ss=ss+length(G_kp{NK});
    end
    
    E_BAND_ALL=[];
    for NK=1:NKPTS_red
        E_BAND_ALL=[E_BAND_ALL;E_BAND{NK}];
    end
    
    [E_BAND_ALL,E_idx]=sort(E_BAND_ALL,'ascend');
    
    Band_kp=zeros(ss,1);
    Band_idx = zeros(1,NPL_sc);
    Ni=0; 
    for NK=1:NKPTS_red
        for i=Ni+1:Ni+length(G_kp{NK})
            Band_kp(E_idx==i)=NK;
            Band_idx(i)=find(E_idx==i);
        end
        Ni=Ni+length(G_kp{NK});
    end
    VKPT_red=VKPT(:,VKPTS_idx);
    
    BASIS(NK_sc).E_BAND=E_BAND;
    BASIS(NK_sc).G_INDEX_sc=G_INDEX_sc;
    BASIS(NK_sc).G_kp=G_kp;
    BASIS(NK_sc).NKPTS_red=NKPTS_red;
    BASIS(NK_sc).VKPTS_idx=VKPTS_idx;
    BASIS(NK_sc).WAVE=WAVE;
    BASIS(NK_sc).Band_kp=Band_kp;
    BASIS(NK_sc).E_BAND_ALL=E_BAND_ALL;
    BASIS(NK_sc).E_idx=E_idx;
    BASIS(NK_sc).Band_idx=Band_idx;
    BASIS(NK_sc).VKPT_red=VKPT_red;
    BASIS(NK_sc).Factor_sym=Factor_sym;
    
    toc;
    clear E_BAND G_INDEX_sc G_kp NKPTS_red VKPTS_idx WAVE Band_kp E_BAND_ALL VKPT_red
end

save('../BASIS_NDsym_nscf.mat','BASIS','VKPT','-v7.3')
