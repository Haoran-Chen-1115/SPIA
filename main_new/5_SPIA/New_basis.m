addpath('../')
addpath('../0_Public')
parameter;
gpuDevice;
if LBloch
    load('BASIS_NDsym_nscf.mat')
end
path='../';

%% Read in informations
[A,B,OMEGA_sc,...
    NKPTS,ISPIN,NTYP,NITYP,NCPU,...
    NELE,VKPT,NELE_TYP,POMASS]=...
    RD_HEAD(ROOT_DIR);
NBANDS = ceil(NELE/2 * NBANDS_mul);

%% Load in equilibrium positions
POSION_EQ=cell(1,NTYP);
for NT=1:NTYP
    % Phase's space cost is huge. Instead read positions.
    POSS=fopen([ROOT_DIR,'/BEAD_1_sym/POS_1_',int2str(NT)]);
    POSION_EQ{NT}=fread(POSS,[3,NITYP(NT)],'double');
    fclose(POSS);
end

%%
FID=fopen([ROOT_DIR,'/BEAD_1_primitive/HEAD_1']);
A_prim=fread(FID,[3,3],'double');
fread(FID,2,'double');
fread(FID,2,'int');
NTYP=fread(FID,1,'int');NITYP_prim=zeros(1,NTYP);
for NT=1:NTYP
    NITYP_prim(NT)=fread(FID,1,'int');
end
fclose(FID);
A_prim=A_prim*diag(N_sc_base);
B_prim=2*pi*inv(A_prim);
OMEGA_prim=det(A_prim);

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

if L_POS_new
    N_sc = N_sc_ex*diag(N_sc_base);
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

%%
%[VKPTS_f,NKPTS_f]=Full_K(N_k_nscf.*N_sc_base,1);
[VKPTS_f,NKPTS_f]=Full_K(N_k_nscf,1);

SYM_OP_min(:,:,1) = eye(3);
SYM_OP_min(:,:,2) = -eye(3); % minimal symmetry with inversion symmetry
[VKPT,NKPTS,VKPTS_nosym_idx]...
    =Kpoints_gsym(VKPTS_f,A,SYM_OP_min);

VKPT = VKPT(:,1:3)';


FID=fopen([ROOT_DIR,'/BEAD_1_primitive/INDEX_1_1_1']);
NGX=fread(FID,1,'int');NGY=fread(FID,1,'int');NGZ=fread(FID,1,'int');
fclose(FID);

NGX=NGX*N_sc_base(1);
NGY=NGY*N_sc_base(2);
NGZ=NGZ*N_sc_base(2);

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
    
    clear G_f
    G_f(1:3,:)=gpuArray(B_prim).'*G_INDEX_tmp(1:3,:);
    G_f(4,:)=G_f(1,:).^2+G_f(2,:).^2+G_f(3,:).^2;
    
    idx = find(HSQDTM*G_f(4,:)<ENCUT);
    G_INDEX_prim{NK} = gather([G_INDEX_tmp(:,idx);idx]);
    NPL_ss(NK) = length(G_INDEX_prim{NK});
end
%%
if ~LBloch
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
    
    A_sc=A; B_sc=B;
end
%%
NK_group = floor(linspace(1,NKPTS_sym+1,NCLUSTER+1));
NK_BEGIN = NK_group(NCL);
NK_END = NK_group(NCL+1)-1;

Gbar_inv_all=cell(1,NKPTS_sym);
NBANDS_ss = NELE * 2;
E_BAND_ALL_NK = zeros(NBANDS_ss,NKPTS_sym);
for NK=NK_BEGIN:NK_END
    NK_sym=NK;
    
    if LBloch
        VKPT_NK=(BASIS(NK_sym).VKPT_red.').*BASIS(NK).Factor_sym;
    end
    load([path,'Gbar_',int2str(NK),'.mat'],'Gbars')
    load('EFERMI.mat')
    
    % Heff = (1i*ETA + EFERMI_Av)*eye(length(Gbars))-(Gi+Gi')/2;
    WAVE_gamma=zeros(length(Gbars));
    Gi_new = zeros(length(Gbars));
    if LBloch
        NKPTS_red = BASIS(NK).NKPTS_red;
    else
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
        
        %% G_INDEX_sc
        NK_sc=1;
        G_INDEX_tmp = G_INDEX_sc_f + VKPT_sym(:,NK);
        
        G_f = zeros(4,NPLWV_sc,'gpuArray');
        G_f(1:3,:)=gpuArray(B_sc).'*G_INDEX_tmp(1:3,:);
        G_f(4,:)=G_f(1,:).^2+G_f(2,:).^2+G_f(3,:).^2;
        
        idx = find(HSQDTM*G_f(4,:)<ENCUT);
        G_INDEX_sc = gather([G_INDEX_tmp(:,idx);idx]);
        NPL_sc = length(G_INDEX_sc);
        
        G=zeros(4,NPL_sc,'like',B_sc);
        %G_INDEX=gpuArray(G_INDEX);
        %G=zeros(4,NPL,'gpuArray');
        %G=zeros(4,NPL);
        G(1:3,:)=B_sc.'*G_INDEX_sc(1:3,:);
        G(4,:)=G(1,:).^2+G(2,:).^2+G(3,:).^2;
        
        clear G2
        G2(4,:)=sqrt(G(4,:));
        G2(4,:)=max(G2(4,:),1E-12);
        G2(1:3,:)=bsxfun(@rdivide,G(1:3,:),G2(4,:));
        
        G=single(G);
        %%
        eps = 1E-5;
        N_sc = diag(N_sc_base)*N_sc_ex;
        G_INDEX_sc_pr=mod(N_sc_ex\G_INDEX_sc(1:3,:),1).';
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
        for NK_s=1:NKPTS_red
            NK_c=int2str(VKPT_sps_ia(NK_s));
            NPL_s=NPL_ss(VKPT_sps_ia(NK_s));
            G_INDEX_sym=G_INDEX_prim{VKPT_sps_ia(NK_s)};
            
            G_INDEX_sym(1:3,:)=Factor_sym(NK_s)*G_INDEX_sym(1:3,:);
            [AA,G_kpidx] = ismembertol(...
                (N_sc_ex\G_INDEX_sc(1:3,:)).',G_INDEX_sym(1:3,:).',1E-4,'ByRows',true);
            
            [~,idx]=sort(G_kpidx(G_kpidx~=0),'ascend');
            G_kpidx=find(AA);
            G_kpidx=G_kpidx(idx);
            
            G_kp{NK_s}=G_kpidx;
            
            NPL_sum=NPL_sum+NPL_s;
        end
        
        if NPL_sum~=NPL_sc
            error('Plane wave number inconsistent')
            %return
        end
        
        %% Generate the Transformation matrix for equilibrium configuration
        [QPROJ_ORI,TRANS_CORE]=...
            OVERL_CORE_nscf(PSMAXN,PSPNL,NTYP,...
            LMAX,LMMAXC,RG,NMAX,...
            NPL_sc,CH0,CH1,CH2,CH3,LPS,...
            WAE_PS,WAE,WPS,G2,OMEGA_sc);
        
        TRANS_EQ=eye(NPL_sc,'single','gpuArray');
        for NT=1:NTYP
            POSION1=gpuArray(POSION_EQ{NT});
            CREXP=single(exp(CITPI*G_INDEX_sc(1:3,1:NPL_sc).'*POSION1));
            CREXP=conj(CREXP)*CREXP.';
            % For \hat{T}_c=|\tilde{p}_i^a><\tilde{\phi}_i^a|...
            TRANS_EQ=TRANS_EQ+CREXP.*TRANS_CORE{NT};
        end
        TRANS_EQ=gather(TRANS_EQ);
        clear CREXP
        
        FOVL_EQ=OVERL(POSION_EQ, POSION_EQ, TRANS_EQ, TRANS_EQ, A_sc, B_sc, ...
            CQIJ3_fun, r, G_INDEX_sc, CITPI, QPROJ_ORI, ...
            NTYP,NITYP_EQ,NPL_sc,LMAX,LMMAXC, ...
            CH0,CH1,CH2,CH3,LPS);
        FOVL_EQ=(FOVL_EQ+FOVL_EQ')/2;
        clear QPROJ_ORI TRANS_CORE
    end
    
    WAVE=cell(1,NKPTS_red);
    E_BAND=cell(1,NKPTS_red);
    if ~LBloch
        Gi_new2=cell(1,NKPTS_red);
    end
    
    for NK_s=1:NKPTS_red
        if LBloch
            idx = find(BASIS(NK).Band_kp==NK_s);
            idx = idx(idx<=length(Gbars));
            NBANDS_s = length(idx);

            Gbars_s = Gbars(idx,idx);
            Gi_s = inv(Gbars_s);

            Heff = EFERMI_Av*eye(length(Gbars_s))-(Gi_s+Gi_s')/2;

            [BASIS_re,E_BAND_re] = eig(Heff);
            WAVE{NK_s} = BASIS(NK).WAVE{NK_s}(:,1:NBANDS_s)*BASIS_re;
            E_BAND{NK_s} =diag(E_BAND_re);

            Gi_new(idx,idx) = BASIS_re'*Gi_s*BASIS_re;
        else
            idx = G_kp{NK_s};
            NBANDS_s=length(idx); % =NPL_s
            
            % Set a set of orthonormal bases
            FOVL = FOVL_EQ(idx,idx);
            [F,~]=eig(FOVL);
            E=real(diag(F'*FOVL*F));
            F=F./sqrt(E.');
            Gbars_s = ...
                F'*Gbars(idx,idx)*F;
            Gi_s = inv(Gbars_s);
            
            % Expand in orthonormal bases
            %Gbars_s = ...
            %    BASIS(NK).WAVE{NK_s}'*Gbars(idx,idx)*BASIS(NK).WAVE{NK_s};
            %Gi_s = inv(Gbars_s);

            Heff = EFERMI_Av*eye(length(Gbars_s))-(Gi_s+Gi_s')/2;
            clear Gbars_s

            [BASIS_re,E_BAND_re] = eig(Heff);
            E_BAND{NK_s} =gather(diag(E_BAND_re));
            WAVE{NK_s} = gather(F*BASIS_re);

            %idx2 = find(BASIS(NK).Band_kp==NK_s);
            Gi_new2{NK_s} = gather(BASIS_re'*Gi_s*BASIS_re);
        end
    end
    clear F E FOVL_EQ BASIS_re E_BAND_re Gi_s Heff
    
    % Band component analysis
    ss=0;
    for NK_s=1:NKPTS_red
        ss=ss+length(G_kp{NK_s});
    end
    
    E_BAND_ALL=[];
    for NK_s=1:NKPTS_red
        E_BAND_ALL=[E_BAND_ALL;E_BAND{NK_s}];
    end
    
    [E_BAND_ALL,E_idx]=sort(E_BAND_ALL,'ascend');
    
    Band_idx = zeros(1,length(Gbars));
    Band_kp=zeros(ss,1);
    Ni=0; 
    for NK_s=1:NKPTS_red
        for i=Ni+1:Ni+length(E_BAND{NK_s})
            Band_idx(i)=find(E_idx==i);
            if ~LBloch
                Band_kp(E_idx==i)=NK;
            end
        end
        Ni=Ni+length(E_BAND{NK_s});
    end
    
    if ~LBloch
        Ni=0;
        for NK_s=1:NKPTS_red
            for i=Ni+1:Ni+length(E_BAND{NK_s})
                Gi_new(Band_idx(Ni+1:Ni+length(E_BAND{NK_s})),...
                    Band_idx(i))=Gi_new2{NK_s}(:,i-Ni);
            end
            Ni=Ni+length(G_kp{NK_s});
        end
    end
    
    if ~LBloch
        VKPT_red=VKPT(:,VKPTS_idx);
        
        BASIS.E_BAND=E_BAND;
        BASIS.G_INDEX_sc=G_INDEX_sc;
        BASIS.G_kp=G_kp;
        BASIS.NKPTS_red=NKPTS_red;
        BASIS.VKPTS_idx=VKPTS_idx;
        BASIS.WAVE=WAVE;
        BASIS.Band_kp=Band_kp;
        BASIS.E_BAND_ALL=E_BAND_ALL;
        BASIS.E_idx=E_idx;
        BASIS.Band_idx=Band_idx;
        BASIS.VKPT_red=VKPT_red;
        BASIS.Factor_sym=Factor_sym;
    else
        BASIS(NK).E_BAND=E_BAND;
        BASIS(NK).WAVE=WAVE;
        BASIS(NK).E_BAND_ALL=E_BAND_ALL;
        BASIS(NK).E_idx=E_idx;
        BASIS(NK).Band_idx=Band_idx;
    end

    Gbar_inv_all{NK} = diag(Gi_new);
    save([path,'/BASIS_NDsym_new_',int2str(NK),'.mat'],'BASIS','VKPT','-v7.3')

    E_BAND_ALL_NK(:,NK) ...
        =BASIS.E_BAND_ALL(1:NBANDS_ss);
end
E_BAND_ALL=E_BAND_ALL_NK;

%%
EFERMI_gauss = EFERMI_Av;
for NK=NK_BEGIN:NK_END
    Gbar_inv = Gbar_inv_all{NK};
    Gbar_inv = Gbar_inv - EFERMI_Av + EFERMI_gauss;
    Gbar_diag = 1./Gbar_inv;

    save([path,'/Gbar_inv_new_',int2str(NK),'.mat'],'Gbar_diag','Gbar_inv')
end

EFERMI_Av = EFERMI_gauss;
if ~exist('../EFERMI_new.mat','file')
    save([path,'/EFERMI_new.mat'],'EFERMI_Av')
end
gpuDevice([]);
