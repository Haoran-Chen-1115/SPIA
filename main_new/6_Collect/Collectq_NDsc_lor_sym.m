%% Initialization
addpath('../')
addpath('../0_Public/')
parameter;
maxNumCompThreads(12);
%gpuDevice;

BETA=kB*TEMP;

tic;
[A,B,OMEGA,...
    NKPTS,ISPIN,NTYP,NITYP,NCPU,...
    NELE,VKPT,NELE_TYP,POMASS]=...
    RD_HEAD(ROOT_DIR);
NBANDS = ceil(NELE/2 * NBANDS_mul);
NIONS=sum(NITYP);

%% Load Basis
tic;
% load('EFERMI.mat')

% q_weight=[];

NBANDS_s=2*NELE;


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
    POSION_prim=cell(1,NTYP);
    for NT=1:NTYP
        % Phase's space cost is huge. Instead read positions.
        POSS=fopen([ROOT_DIR,'/BEAD_1_primitive/POS_1_',int2str(NT)]);
        POSION_prim{NT}=fread(POSS,[3,NITYP_prim(NT)],'double');
        fclose(POSS);
    end
end

%% Load in equilibrium positions
if L_POS_new
    MAXC=det(N_sc_ex);
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
        POS_tmp = N_sc_ex'\POS_tmp;
        POS_tmp = mod(POS_tmp,1);
        POS_tmp(POS_tmp>1-eps)=POS_tmp(POS_tmp>1-eps)-1;
        POS_tmp = uniquetol(POS_tmp','ByRows',true)';
        POSION_new{NT} = POS_tmp;
        NITYP_new(NT) = length(POS_tmp);
    end
    POSION_EQ = POSION_new;
    NITYP_EQ = NITYP_new;
else
    POSION_EQ=cell(1,NTYP);
    for NT=1:NTYP
        % Phase's space cost is huge. Instead read positions.
        POSS=fopen([ROOT_DIR,'/BEAD_1_sym/POS_1_',int2str(NT)]);
        POSION_EQ{NT}=fread(POSS,[3,NITYP(NT)],'double');
        fclose(POSS);
    end
    NITYP_EQ = NITYP;
end

%% Diagonal supercell (base)
A_base = diag(N_sc_base)*A_prim;
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

[IBRAV,N_SYM_TRUE,SYM_OP_TRUE,TRANS_ROT]...
    = Symmetry(A,POSION_EQ,NTYP,NITYP_EQ);

SYM_OP_INV = SYM_INV(SYM_OP_TRUE,A);

[VKPT_sym,NKPTS_sym,~]...
    =Kpoints_gsym(VKPT_NDsc,A,SYM_OP_TRUE);

if ~L_IBZ && ~L_SYM
    SYM_OP_min(:,:,1) = eye(3);
    SYM_OP_min(:,:,2) = -eye(3);
    [VKPT_sym,NKPTS_sym,~]...
        =Kpoints_gsym(VKPT_NDsc,A,SYM_OP_min);
else
    SYM_OP_min(:,:,1) = eye(3);
    SYM_OP_min(:,:,2) = -eye(3);
    [VKPT_nosym,NKPTS_nosym,VKPTS_nosym_idx]...
        =Kpoints_gsym(VKPT_NDsc,A,SYM_OP_min);
    
    k_weight_nosym=VKPT_nosym(:,4)';
    VKPT_nosym=VKPT_nosym(:,1:3)';
end

k_weight=VKPT_sym(:,4)';
VKPT_sym=VKPT_sym(:,1:3)';

if L_IBZ && L_SYM
    [~,VKPTS_sym_idx]=ismembertol(VKPT_sym.',VKPT_nosym.','ByRows',true);
end

%% Indexfor k and k+q
% Now we can calculate \lambda_q,
% which should be equal to \sum_\nu q_weight*\lambda_q\nu
% As we calculated every q point, q_weight = 1/N_q.
NBEAD1 = NBEAD * 16;
tau = (0:(NBEAD-1))/NBEAD;
tau1 = (0:(NBEAD1-1))/NBEAD1;

% Transform W11 into the form of
% W11_w(NBANDS,NBANDS,NK,NQ)
NBANDS_ss = NELE * 2;
E_BAND_ALL = zeros(NBANDS_ss,NKPTS_sym);
ImE_BAND_ALL = zeros(NBANDS_ss,NKPTS_sym);
BETA = kB*TEMP;
ETA = (2*KOMEGA+1)*pi*BETA;
for NK=1:NKPTS_sym
    if ~L_IBZ && L_SYM
        NK_sym=VKPTS_sym_idx(NK);
    else
        NK_sym=NK;
    end
    load(['Gbar_inv_new_',int2str(NK_sym),'.mat'],'Gbar_inv','Gbar_diag')
    
    load('EFERMI_new.mat')
    E_BAND_ALL(:,NK) ...
        = real(EFERMI_Av-Gbar_inv(1:NBANDS_ss));
    ImE_BAND_ALL(:,NK) = -imag(1i*ETA-Gbar_inv(1:NBANDS_ss));
    %E_BAND_ALL(:,NK) ...
    %    =BASIS(NK_sym).E_BAND_ALL(1:NBANDS_ss);
end

%% Determine the q-transfer grid
load('../BASIS_NDsym_new_1.mat')
BASIS=BASIS;

VKPTS=BASIS.VKPT_red.*BASIS.Factor_sym';
DX=reshape(VKPTS(1,:)-VKPTS(1,:)',1,[]);
DY=reshape(VKPTS(2,:)-VKPTS(2,:)',1,[]);
DZ=reshape(VKPTS(3,:)-VKPTS(3,:)',1,[]);
DK=[DX;DY;DZ].*N_sc_base.';
DK=mod(DK,1);
DK(DK-0.5>1E-5)=DK(DK-0.5>1E-5)-1;
[VQPTS,idx,~]=uniquetol(DK.',1E-4,'ByRows',true,'DataScale',1);
NQPTS=size(VQPTS,1);

%% Initialization
Nsigma=80;
Nomega=1;
lambda = cell(Nomega,Nsigma);
lambda_q = cell(Nomega,Nsigma);
lambda2_q = cell(Nomega,Nsigma);
lambda2 = cell(Nomega,Nsigma);

phase_space = zeros(1,Nsigma);
dos_ef = zeros(Nomega,Nsigma);

No=1;
for nsig=1:Nsigma
    lambda_q{No,nsig}=zeros(NQPTS,NBEAD/2+1);
    lambda{No,nsig}=zeros(1,NBEAD/2+1);
    
    lambda2_q{No,nsig}=zeros(NQPTS,NBEAD1/2+1);
    lambda2{No,nsig}=zeros(1,NBEAD1/2+1);
    phase_space(No,nsig)=0;
    
    gamma=0.02*nsig;
    
    % load('../../main_nscf_spec_eig/dos_Av.mat','EFERMI_Av')
    load('../EFERMI.mat')
    EFERMI_lor = EFERMI_Av;
    
    % Calculate band weight and DOS at Fermi energy
    w0lor = @(x,d) d./(x.^2+d.^2)/pi;
    
    wk1 = w0lor( EFERMI_lor - E_BAND_ALL ,gamma);
    
    dos_ef(No,nsig) = 0;
    for NK=1:NKPTS_sym
        dos_ef(No,nsig) = dos_ef(No,nsig) +...
            sum(k_weight(NK)*wk1(:,NK))/sum(k_weight(1:NKPTS_sym));
    end
end

%% Transform to another representation
for NK=1:NKPTS_sym
    W11_Q1 = zeros(NBANDS_ss,NBANDS_ss,NBEAD);
    W11_Q2 = zeros(NBANDS_ss,NBANDS_ss,NBEAD1);
    
    tic;
    
    if ~L_IBZ && L_SYM
        NK_sym=VKPTS_sym_idx(NK);
        NC_c='1';
        [G_INDEX,G,G2,NPL,NPLWV,NGX,NGY,NGZ]=...
            RD_INDEX(ROOT_DIR,'1','1',VKPT_nosym(:,NK_sym),B,ENCUT);
        
        %% Symmetrize Gamma
        G_INDEX_equiv=cell(1,N_SYM_TRUE);
        VKPT_equiv=cell(1,N_SYM_TRUE);
        VKPT_equiv_idx=zeros(1,N_SYM_TRUE);
        VKPT_equiv_factor = ones(1,N_SYM_TRUE);
        G_equiv_idx=cell(1,N_SYM_TRUE);
        
        G_INDEX_o = G_INDEX(1:3,:);% + VKPT_sym(1:3,NK_sym);
        for NS=1:N_SYM_TRUE
            G_INDEX_equiv{NS} = ...
                gather(SYM_OP_INV(:,:,NS)' * G_INDEX_o);
            VKPT_equiv{NS} = SYM_OP_INV(:,:,NS)' * VKPT_nosym(1:3,NK_sym);
            VKPT_equiv{NS} = mod( VKPT_equiv{NS},1);
            VKPT_equiv{NS}( VKPT_equiv{NS}>0.5+eps) =...
                VKPT_equiv{NS}( VKPT_equiv{NS}>0.5+eps) - 1;
            
            [~,VKPT_equiv_idx(NS)] = ismembertol(VKPT_equiv{NS}',...
                VKPT_nosym','ByRows',true);
            
            VKPT_equiv_inv = VKPT_equiv{NS};
            if VKPT_equiv_idx(NS)==0
                VKPT_equiv_factor(NS)=-1;
                
                VKPT_equiv_inv = -VKPT_equiv{NS};
                VKPT_equiv_inv(VKPT_equiv_inv<-0.5+eps)=...
                    VKPT_equiv_inv(VKPT_equiv_inv<-0.5+eps)+1;
                [~,VKPT_equiv_idx(NS)]=ismembertol(...
                    VKPT_equiv_inv',VKPT_nosym',1E-4,'ByRows',true);
            end
            
            [G_INDEX_o_k,~,~,~,~,~,~,~]=...
                RD_INDEX(ROOT_DIR,'1','1',...
                VKPT_equiv_inv,B,ENCUT);
            G_INDEX_o_k = G_INDEX_o_k(1:3,:);% + VKPT_equiv{NS};
            
            G_INDEX_o_k = gather(G_INDEX_o_k * VKPT_equiv_factor(NS));
            [~,G_equiv_idx{NS}] = ismembertol(G_INDEX_equiv{NS}',G_INDEX_o_k'...
                ,'ByRows',true);
            %[~,G_equiv_idx_inv{NS}] = ismembertol(G_INDEX_o_k',G_INDEX_equiv{NS}'...
            %    ,'ByRows',true);
        end
        
        unique_NK = unique(VKPT_equiv_idx);
        NKPTS_unique = length(unique_NK);
        weight_NK = zeros(NKPTS_unique,1);
        for NS=1:N_SYM_TRUE
            [~,idx]=ismember(VKPT_equiv_idx(NS),unique_NK);
            weight_NK(idx)=weight_NK(idx)+1;
        end
        
        if NKPTS_unique~=k_weight(NK)/k_weight_nosym(NK_sym)
            error('Symmetry operations yield in-consistent k meshes')
        end
        
        Gamma_k = zeros([NBANDS,NBANDS,NBEAD],'single');
        DTbars_k = zeros([NBANDS,NBANDS,NBEAD],'single');
        for NS=1:NKPTS_unique
            NK_loc = unique_NK(NS);
            load(['../Gamma_',int2str(NK_loc),'.mat'],'Gamma','DTbars')
            
            Gamma_k = Gamma_k + Gamma*weight_NK(NS)/N_SYM_TRUE;
            DTbars_k = DTbars_k + DTbars*weight_NK(NS)/N_SYM_TRUE;
        end
        Gamma=Gamma_k;
        DTbars=DTbars_k;
        clear DTbars_k
    else
        NK_sym=NK;
        NC_c='1';
        [G_INDEX,G,G2,NPL,NPLWV,NGX,NGY,NGZ]=...
            RD_INDEX(ROOT_DIR,'1','1',VKPT_sym(:,NK_sym),B,ENCUT);
        
        load(['Gamma_',int2str(NK_sym),'.mat'],'Gamma','DTbars')
    end
    
    load(['Gbar_inv_new_',int2str(NK_sym),'.mat'],'Gbar_inv','Gbar_diag')
    
    %% Solve Bethe-Salpeter Equation
    Gamma_k=Gamma-abs(DTbars).^2;
    Gamma_tau=real(fft(Gamma_k(1:NBANDS,1:NBANDS,:),NBEAD,3));
    W11_tau=zeros([NBANDS,NBANDS],'like',Gamma_k);
    for NB=1:NBEAD
        W11_tau(:,:,NB)=(eye(NBANDS,'like',Gamma_k)+Gamma_tau(:,:,NB)*...
            diag(abs(Gbar_diag(1:NBANDS)).^2))\Gamma_tau(:,:,NB);
    end
    W11_w=real(gather(ifft(W11_tau(1:NBANDS_s,1:NBANDS_s,:),NBEAD,3)));
    
    % Interpolate W
    W11_tau_s = gather(W11_tau(1:NBANDS_s,1:NBANDS_s,:));
    
    Wq_tau = W11_tau_s;
    Wq_tau(:,:,end+1) = Wq_tau(:,:,1);
    Wq_tau = permute(Wq_tau,[3,1,2]);
    W11_tau = permute(interp1([tau,1],Wq_tau,tau1,'linear'),[2,3,1]);
    
    W11_w2 = real(gather(ifft(W11_tau,NBEAD1,3)));
    
    % W11_w(NBANDS,NBANDS,NK,NQ)
    for NB = 1:NBEAD
        W11_Q1(:,:,NB) = ...
            gather(W11_w(1:NBANDS_ss,1:NBANDS_ss,NB));
    end
    for NB = 1:NBEAD1
        W11_Q2(:,:,NB) = ...
            W11_w2(1:NBANDS_ss,1:NBANDS_ss,NB);
    end
    
    save(['../W11_',int2str(NK),'.mat'],'W11_Q1','W11_Q2','-v7.3');
    clear W11_tau W11_tau1 W11_w W11_w2 Wq_tau
    toc;
    
    load(['BASIS_NDsym_new_',int2str(NK_sym),'.mat'])
    
    %load(['../W11_',int2str(NK),'.mat'],'W11_Q1','W11_Q2');
    %%
    for nsig=1:Nsigma
        gamma=0.02*nsig;
        
        load('../EFERMI.mat')
        %load('../../main_nscf_spec_eig/dos_Av.mat','EFERMI_Av')
        EFERMI_lor = EFERMI_Av;
        
        % Calculate band weight and DOS at Fermi energy
        w0lor = @(x,d) d./(x.^2+d.^2)/pi;
        
        wk1 = w0lor( EFERMI_lor - E_BAND_ALL ,gamma);
        % Calculate lambda
        VKPTS=BASIS.VKPT_red.*BASIS.Factor_sym';
        DX=reshape(VKPTS(1,:)'-VKPTS(1,:),1,[]);
        DY=reshape(VKPTS(2,:)'-VKPTS(2,:),1,[]);
        DZ=reshape(VKPTS(3,:)'-VKPTS(3,:),1,[]);
        DK=[DX;DY;DZ].*N_sc_base.';
        DK=mod(DK,1);
        DK(DK-0.5>1E-5)=DK(DK-0.5>1E-5)-1;
        [VQPTS_NK,~,idx]=uniquetol(DK.',1E-4,'ByRows',true,'DataScale',1);
        NQPTS_NK=size(VQPTS_NK,1);
        if NQPTS~=NQPTS_NK
            error('q-grid wrong');
        end
        
        Ni=0; Band_idx_all=cell(1,BASIS.NKPTS_red);
        for NK_s=1:BASIS.NKPTS_red
            for i=Ni+1:Ni+length(BASIS.E_BAND{NK_s})
                if BASIS.Band_idx(i)<=NBANDS_s
                    Band_idx_all{NK_s}(i-Ni)=BASIS.Band_idx(i);
                end
            end
            Ni=Ni+length(BASIS.E_BAND{NK_s});
        end
        
        for NQ=1:NQPTS
            idx_q=find(idx==NQ);
            NK_start=mod(idx_q-1,BASIS.NKPTS_red)+1;
            NK_end=ceil(idx_q/BASIS.NKPTS_red);
            
            NK_start_u = unique(NK_start);
            for NK_s=1:length(NK_start_u)
                idx_kq=find(NK_start==NK_start_u(NK_s));
                NK_kq_end=NK_end(idx_kq);
                
                Band_idx_k = Band_idx_all{NK_start_u(NK_s)};
                Band_idx_kq = cell2mat(Band_idx_all(NK_kq_end));
                
                wks1=...
                    w0lor( EFERMI_lor - E_BAND_ALL(Band_idx_k,NK) ...
                    ,gamma);
                wks2=...
                    w0lor( EFERMI_lor - E_BAND_ALL(Band_idx_kq,NK)...
                    ,gamma);
                
                weight = k_weight(NK) * wks1 .* wks2.'...
                    /sum(k_weight(1:NKPTS_sym));
                
                for NB = 1:NBEAD/2+1
                    lambda_q{No,nsig}(NQ,NB) = lambda_q{No,nsig}(NQ,NB) ...
                        + sum(sum(weight .*...
                        W11_Q1(Band_idx_k,Band_idx_kq,NB)))/BETA/dos_ef(No,nsig);
                end
                
                for NB = 1:NBEAD1/2+1
                    lambda2_q{No,nsig}(NQ,NB) = lambda2_q{No,nsig}(NQ,NB) ...
                        + sum(sum(weight .*...
                        W11_Q2(Band_idx_k,Band_idx_kq,NB)))/BETA/dos_ef(No,nsig);
                end
            end
        end
    end
end

rho_2 = zeros(Nomega,Nsigma);
omega2_2 = zeros(Nomega,Nsigma);
for nsig=1:Nsigma
    for NB=1:NBEAD/2+1
        lambda{No,nsig}(NB)=sum(lambda_q{No,nsig}(:,NB));
    end
    
    for NB=1:NBEAD1/2+1
        lambda2{No,nsig}(NB)=sum(lambda2_q{No,nsig}(:,NB));
    end
    
    Ncut=NBEAD/2;
    [rho_2(No,nsig),omega2_2(No,nsig)]=Gbar_Eliashberg(lambda2{No,nsig},Ncut,mustar);
end

lam = zeros(Nomega,Nsigma);
for No=1:Nomega
    for nsig = 1:Nsigma
        lam_no(No,nsig) = lambda{No,nsig}(1);
        lam(No,nsig) = lambda2{No,nsig}(1);
    end
end

save('Result_q_lor','rho_2','omega2_2','mustar','Nsigma',...
    'lambda2','lam','lambda2_q',...
    'lambda','lam_no','lambda_q','EFERMI_lor','dos_ef')
%save('Intermediate','W11_Q1','E_BAND_ALL','k_weight','-v7.3')

%% Determine the transition temperature through interpolating lambda
BETA=kB*TEMP;
Ncut=NBEAD/2;
No=1;
%nsig=12; % For LaH9.6, nsig=12;
%if I_Ef~=1
%    dos_ef_sigma=dos_ef;
%    load('../4_Spec/dos_Av.mat','dos_ef','centers');
%    idx = find(centers==0);
%    dos_ef_all=dos_ef;
%    dos_ef=dos_ef_sigma;
%    [~,nsig] = min(abs(dos_ef-dos_ef_all(idx)));
%end
T_c=zeros(1,Nsigma);
for nsig=1:Nsigma
    lam_f=@(x) interp1(BETA*[0:NBEAD/2],...
        lambda2{No,nsig}(1:NBEAD/2+1),x,'pchip');
    
    rho=100;
    T_up = TEMP+200;
    T_down = TEMP-200;
    while abs(rho)>1E-3
        
        T_tmp = (T_up+T_down)/2;
        BETA_tmp = kB*T_tmp;
        
        nu2 = (BETA_tmp)*[0:NBEAD/2];
        lambda_nu = lam_f(nu2);%/BETA_tmp^2;
        
        %     lambda_nu(1)=lambda2{ns}(1);
        %     lambda_nu(2:end) = lambda_nu(2:end)./([1:NBEAD/2].^2);
        lambda_nu(end+1:end+NBEAD/2-1)=lambda_nu(end-1:-1:2);
        
        
        rho=Gbar_Eliashberg2(lambda_nu,Ncut,mustar,omega2_2(No,nsig)*BETA/BETA_tmp);
        if rho<0
            T_up = T_tmp;
        else
            T_down = T_tmp;
        end
        
        if T_up-T_down<1E-4
            disp('fail to converge')
            break
        end
    end
    T_c(nsig)=T_tmp;
end
save('T_c','T_c','mustar');
