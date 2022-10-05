%% parameter
addpath('../')
addpath('../0_Public')
parameter;
gpuDevice;
maxNumCompThreads(9);
if NK==NK_BEGIN
  FID_OUT=fopen(['out_C',int2str(NCL),'_T',int2str(jobid),'.file'],'w');
else
  FID_OUT=fopen(['out_C',int2str(NCL),'_T',int2str(jobid),'.file'],'a');
end

%% Read in informations
[A,B,OMEGA,...
    NKPTS,ISPIN,NTYP,NITYP,NCPU,...
    NELE,NBANDS,VKPT]=RD_HEAD(ROOT_DIR,NELE_TYP,NBANDS_mul);

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
POSION_prim=cell(1,NTYP);
for NT=1:NTYP
    % Phase's space cost is huge. Instead read positions.
    POSS=fopen([ROOT_DIR,'/BEAD_1_primitive/POS_1_',int2str(NT)]);
    POSION_prim{NT}=fread(POSS,[3,NITYP_prim(NT)],'double');
    fclose(POSS);
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
    = Symmetry(A,POSION_EQ,NTYP,NITYP);

[VKPT_sym,NKPTS_sym,~]...
    =Kpoints_gsym(VKPT_NDsc,A,SYM_OP_TRUE);

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

    %NBANDS=1024;
    %% Load in Bloch basis
    WAVE_gamma=zeros(length(BASIS.G_INDEX_sc));
    Ni=0;
    for NK_s=1:BASIS.NKPTS_red
        for i=Ni+1:Ni+length(BASIS.G_kp{NK_s})
            WAVE_gamma(BASIS.G_kp{NK_s},BASIS.Band_idx(i))=BASIS.WAVE{NK_s}(:,i-Ni);
        end
        Ni=Ni+length(BASIS.G_kp{NK_s});
    end
    clear BASIS
    BASIS_k=gpuArray(single(WAVE_gamma(:,1:NBANDS)));
    clear WAVE_gamma

    load(['Gbar_inv_',int2str(NK),'.mat']);
    Gbar_inv=gpuArray(Gbar_inv(1:NBANDS));
    
    %%
    for ISP=1:ISPIN
        ISP_c=int2str(ISP);
        %% Overlap matrix and transformation matrix related
        NC_c='1';
        [QPROJ_ORI,TRANS_CORE]=...
            OVERL_CORE_nscf(PSMAXN,PSPNL,NTYP,...
            LMAX,LMMAXC,RG,NMAX,...
            NPL,CH0,CH1,CH2,CH3,...
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
        Tbar=complex(zeros(NBANDS,NBANDS,'single','gpuArray'));
        Tbar_F=complex(zeros(NBANDS,NBANDS,'single','gpuArray'));
        Tbar_GFw=zeros(NBANDS,NBANDS,'single','gpuArray');
        Tbar_Fw=complex(zeros(NBANDS,NBANDS,'single','gpuArray'));
        
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
                FHAM = Green_c(ROOT_DIR,NC_c,NB_c,ISP_c,...
                    NTYP,NITYP,LMAX,LMMAXC,CH0,CH1,CH2,CH3,LPS,...
                    NPL,NPLWV,NGX,NGY,NGZ,NCPU,...
                    QPROJ_ORI,TRANS_CORE,...
                    POSION,G_INDEX,G,HSQDTM,CITPI,...
                    POSION_EQ, TRANS_EQ, BASIS_k, A, B, ...
                    CQIJ3_fun, r, kB, TEMP, KOMEGA,EFERMI_Av);
                %FHAM=BASIS_k'*TdT0'...
                %    *(BASIS_trial*diag(1./((1i*ETA+EFERMI_Av)...
                %    -E_BAND_now))*BASIS_trial')...
                %    *TdT0*BASIS_k;
                %FHAM=Gbar_inv*FHAM*Gbar_inv-Gbar_inv;
                FHAM=diag(Gbar_inv)*FHAM*diag(Gbar_inv)...
                    -diag(Gbar_inv);
                Tbar_F(:,:,NB)=FHAM(1:NBANDS,1:NBANDS);
                %clear E_BAND_now 
                clear FHAM
                t2_end=toc(t2);
                fprintf(FID_OUT,'Loop: %6.4f\n',t2_end);
            end
            % Tbar used for checking that Tbar -> 0
            Tbar=Tbar+sum(Tbar_F,3)./NBEAD;
            % Turn to frequency domain. Quasi-static approximation applied.
            Tbar_F=ifft(Tbar_F, NBEAD, 3);
            Tbar_Fw=Tbar_Fw+Tbar_F;
            % Gamma (Scattering Amplitude; or the flucatuation of T matrix)
            Tbar_GFw=Tbar_GFw+abs(Tbar_F).^2;
        end
        
        if jobid==1
            save(['../FERMI_Surface_',int2str(NK),'.mat'],'NBANDS')
            save(['../INDEX_',int2str(NK),'.mat'],'G','G_INDEX')
        end

        Tbar=gather(Tbar);Tbar_GFw=gather(Tbar_GFw);Tbar_Fw=gather(Tbar_Fw);
        fn=['../Tbar_GFw_',int2str(NK),'_',int2str(jobid),'.mat'];
        save(fn,'Tbar_GFw','Tbar','Tbar_Fw','-v7.3')
        
    end
    t3_end=toc(t3);
    fprintf(FID_OUT,'K-loop: %6.4f\n',t3_end);
    fprintf(FID_OUT,' \n');
%end

fclose(FID_OUT);
