%% parameter
addpath('../')
addpath('../0_Public')
parameter;
gpuDevice(jobid);
maxNumCompThreads(2);
if NK==NK_BEGIN
  FID_OUT=fopen(['out_C',int2str(NCL),'_T',int2str(jobid),'.file'],'w');
else
  FID_OUT=fopen(['out_C',int2str(NCL),'_T',int2str(jobid),'.file'],'a');
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

SYM_OP_INV = SYM_INV(SYM_OP_TRUE,A);

[VKPT_sym,NKPTS_sym,~]...
    =Kpoints_gsym(VKPT_NDsc,A,SYM_OP_TRUE);

k_weight=VKPT_sym(:,4)';
VKPT_sym=VKPT_sym(:,1:3)';

if L_IBZ
    [VKPT_nosym,NKPTS_nosym,~]...
        =Kpoints_gsym(VKPT_NDsc,A,SYM_OP_TRUE);
    k_weight_nosym=VKPT_nosym(:,4)';
    VKPT_nosym=VKPT_nosym(:,1:3)';
else
    SYM_OP_min(:,:,1) = eye(3);
    SYM_OP_min(:,:,2) = -eye(3);
    [VKPT_nosym,NKPTS_nosym,~]...
        =Kpoints_gsym(VKPT_NDsc,A,SYM_OP_min);

    k_weight_nosym=VKPT_nosym(:,4)';
    VKPT_nosym=VKPT_nosym(:,1:3)';
end

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
    RD_INDEX(ROOT_DIR,NC_c,'1',VKPT_nosym(:,NK_sym),B,ENCUT);
G_INDEX=gpuArray(G_INDEX);

if ~L_Liquid
    if ~L_IBZ && L_SYM
	LBloch=false;
    else
        LBloch=true; % Now the new basis set has been determined
    end
    load(['../BASIS_NDsym_new_',int2str(NK_sym),'.mat'],'BASIS')
    %NBANDS=1024;
    %% Load in Bloch basis
    WAVE_gamma=zeros(length(BASIS.G_INDEX_sc),length(BASIS.Band_idx));
    Ni=0;
    for NK_s=1:BASIS.NKPTS_red
        for i=Ni+1:Ni+length(BASIS.E_BAND{NK_s})
            WAVE_gamma(BASIS.G_kp{NK_s},BASIS.Band_idx(i))=BASIS.WAVE{NK_s}(:,i-Ni);
        end
        Ni=Ni+length(BASIS.E_BAND{NK_s});
    end
    clear BASIS
    BASIS_k=gpuArray(single(WAVE_gamma(:,1:NBANDS)));
    clear WAVE_gamma

    load(['Gbar_inv_new_',int2str(NK),'.mat']);
else
    LBloch=false;
    BASIS_k=eye([NPL,NBANDS],'single','gpuArray');
    load(['Gbar_inv_',int2str(NK),'.mat']);
end
Gbar_inv=gpuArray(Gbar_inv(1:NBANDS));

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
    if ~L_IBZ && L_SYM
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
        idx_cal_SYM = find(VKPT_equiv_idx==NK_sym);
    end
    
    fprintf(FID_OUT,'Kpoint: %f %f %f\n',VKPT_nosym(:,NK_sym));
    if ~L_IBZ && L_SYM
        NS_cal=0;
        for NS=1:N_SYM_TRUE
            if VKPT_equiv_idx(NS)==NK
                NS_cal=NS_cal+1;
            end
        end
        fprintf(FID_OUT,'Symmetrized %d times\n',NS_cal);
    end
    %%
    if ~L_IBZ && L_SYM
        FHAM_t=complex(zeros(NPL,NPL,NBEAD,'single'));
    else
        Tbar_F=complex(zeros(NBANDS,NBANDS,NBEAD,'single'));
    end
    
    %Tbar=complex(zeros(NBANDS,NBANDS,'single'));
    %Tbar_F=complex(zeros(NBANDS,NBANDS,'single'));
    Tbar_GFw=zeros(NBANDS,NBANDS,NBEAD,'single');
    Tbar_Fw=complex(zeros(NBANDS,NBANDS,NBEAD,'single'));
    
    %for NC=jobid*NSKIP+1+NSTART:NSKIP*NGPU:NCONF
    NCAC = length(N_BEGIN+1:NSKIP:N_END);
    CONF_CAC = N_BEGIN+1:NSKIP:N_END;
    for NCs=1:NCAC
        %        for NC=1:NSKIP*NGPU:NCONF
        NC = CONF_CAC(NCs);
        NC_c=int2str(NC);
        for NB=1:NBEAD
            fprintf(FID_OUT,'Kpoint %d, Configuration %d, Bead %d\n',[NK_sym,NC,NB]);
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
            
            if ~L_IBZ && L_SYM
                FHAM_t(:,:,NB)=gather(FHAM);
            else
                FHAM=diag(Gbar_inv)*FHAM*diag(Gbar_inv)...
                    -diag(Gbar_inv);
                Tbar_F(:,:,NB)=gather(FHAM(1:NBANDS,1:NBANDS));
            end
            %clear E_BAND_now 
            clear FHAM
            t2_end=toc(t2);
            fprintf(FID_OUT,'Loop: %6.4f\n',t2_end);
        end
        
        %%
        t4=tic;
        
        %BASIS_k=gather(BASIS_k);
        %Gbar_inv=gather(Gbar_inv);
        %%
        if ~L_IBZ && L_SYM
            FHAM_t = ifft(FHAM_t, NBEAD, 3);
            Gi = BASIS_k*diag(Gbar_inv);
        
            Tbar_Fw=gpuArray(Tbar_Fw);
            Tbar_GFw=gpuArray(Tbar_GFw);
            for NB=1:NBEAD
                FHAM = gpuArray(FHAM_t(:,:,NB));
                for NS=1:N_SYM_TRUE
                    if VKPT_equiv_idx(NS)==NK
                        tic;
                        Tbar_F = ...
                            Gi'*FHAM(G_equiv_idx{NS},G_equiv_idx{NS})*Gi;
                        if NB==1
                            Tbar_F = Tbar_F - diag(Gbar_inv);
                        end
                        %toc;
                        Tbar_Fw(:,:,NB)=(Tbar_Fw(:,:,NB)+Tbar_F/NS_cal);
                        %toc;
                        Tbar_GFw(:,:,NB)=(Tbar_GFw(:,:,NB)+abs(Tbar_F).^2/NS_cal);
                        toc;
                    end
                end
            end
            Tbar_Fw=gather(Tbar_Fw);
            Tbar_GFw=gather(Tbar_GFw);
            clear FHAM Gi
        else
            % Tbar used for checking that Tbar -> 0
            %Tbar=Tbar+sum(Tbar_F,3)./NBEAD;
            % Turn to frequency domain. Quasi-static approximation applied.
            Tbar_F=ifft(Tbar_F, NBEAD, 3);
            Tbar_Fw=gather(Tbar_Fw+Tbar_F);
            % Gamma (Scattering Amplitude; or the flucatuation of T matrix)
            Tbar_GFw=gather(Tbar_GFw+abs(Tbar_F).^2);
        end
        t4_end=toc(t4);
        fprintf(FID_OUT,'NC-gather: %6.4f\n',t4_end);
    end
    
    if jobid==1
        save(['../FERMI_Surface_',int2str(NK_sym),'.mat'],'NBANDS')
        save(['../INDEX_',int2str(NK_sym),'.mat'],'G','G_INDEX')
    end

    Tbar_GFw=gather(Tbar_GFw);Tbar_Fw=gather(Tbar_Fw);
    fn=['../Tbar_GFw_',int2str(NK_sym),'_',int2str(jobid),'.mat'];
    save(fn,'Tbar_GFw','Tbar_Fw','-v7.3')
    
end
t3_end=toc(t3);
fprintf(FID_OUT,'K-loop: %6.4f\n',t3_end);
fprintf(FID_OUT,' \n');
%end

fclose(FID_OUT);
clear
