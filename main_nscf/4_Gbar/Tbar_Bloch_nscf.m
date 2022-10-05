%% parameter
addpath('../')
addpath('../0_Public')
parameter;
gpuDevice(jobid);
if NK==1
  FID_OUT=fopen(['out_T',int2str(jobid),'.file'],'w');
else
  FID_OUT=fopen(['out_T',int2str(jobid),'.file'],'a');
end

%% Read in informations
[A,B,OMEGA,...
    NKPTS,ISPIN,NTYP,NITYP,...
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
[VKPT_sym,NKPTS_sym,VKPT_f,~,SYM_OP_TRUE]...
    =Kpoints(N_k_nscf,A,NTYP,NITYP,POSION_EQ,1);

SYM_OP_min(:,:,1) = eye(3);
SYM_OP_min(:,:,2) = -eye(3); % minimal symmetry with inversion symmetry
[VKPT,NKPTS_nosym,~]...
    =Kpoints_gsym(VKPT_f,A,SYM_OP_min);

[~,VKPT_sym_idx] = ismembertol(...
    VKPT_sym(:,1:3),VKPT(:,1:3),'ByRows',true);

k_weight = VKPT_sym(:,4)';
VKPT_sym = VKPT_sym(:,1:3)';
VKPT = VKPT(:,1:3)';

%% Calculation
%NKPTS=1;
%for NK=1:NKPTS_sym
%for NK=4:4
    t3=tic;
    NK_c=int2str(VKPT_sym_idx(NK));
    NK_sym=VKPT_sym_idx(NK);
    %NK_c=int2str(NK);
    %NK_sym=NK;
    %% Read in necessary informations:
    %% Mesh informations
    B=gpuArray(B);
    NC_c='1';
    [G_INDEX,G,G2,NPL,NPLWV,NGX,NGY,NGZ]=...
        RD_INDEX(ROOT_DIR,NC_c,'1',VKPT(:,NK_sym),B,ENCUT);

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
    BASIS_k=gpuArray(WAVE_gamma(:,1:NBANDS));
    load(['Gbar_inv_',int2str(NK),'.mat']);
    Gbar_inv=gpuArray(diag(Gbar_inv(1:NBANDS)));
    
    %% Initialization of wavefunction prediction
    ICALL=zeros(NBANDS,1);ICALL_out=-1;
    %BASIS_trial = BASIS_k;
    NIS=0;
    for NT=1:NTYP
        POS_L1(:,NIS+1:NIS+NITYP(NT)) = POSION_EQ{NT};
        NIS = NIS + NITYP(NT);
    end
    
    %%
    for ISP=1:ISPIN
        ISP_c=int2str(ISP);
        %% Overlap matrix and transformation matrix related
        if L_OVERL_OLD
            [QPROJ_ORI,FOVL_CORE,TRANS_CORE,TRANS_c_CORE]=...
                OVERL_CORE_2_OLD(ROOT_DIR,NTYP,LMAX,LMMAXC,...
                NC_c,NK_c,ISP_c,RG,NMAX,...
                NPL,CH0,CH1,CH2,CH3,...
                WAE_PS,WAE,WPS,G2,CQIJ,QPAW3,OMEGA);
        else
            [QPROJ_ORI,RES]=...
                OVERL_CORE_nscf(PSMAXN,PSPNL,NTYP,...
                LMAX,LMMAXC,RG,NMAX,...
                NPL,CH0,CH1,CH2,CH3,...
                WAE_PS,WAE,WPS,G2,OMEGA);
        end

        %% Generate the Transformation matrix for equilibrium configuration
        TRANS_EQ=single(eye(NPL,'gpuArray'));
        for NT=1:NTYP
            POSION1=gpuArray(POSION_EQ{NT});
            CREXP=single(exp(CITPI*G_INDEX(1:3,1:NPL).'*POSION1));
            CREXP=conj(CREXP)*CREXP.';
            % For \hat{T}_c=|\tilde{p}_i^a><\tilde{\phi}_i^a|...
            TRANS_CORE=gather((RES{NT}*QPROJ_ORI{NT}.')/sqrt(OMEGA));
            TRANS_EQ=TRANS_EQ+CREXP.*TRANS_CORE;
        end
        TRANS_EQ=gather(TRANS_EQ);
        clear CREXP TRANS_CORE
        
        %%
        Tbar=complex(zeros(NBANDS,NBANDS,'single','gpuArray'));
        Tbar_F=complex(zeros(NBANDS,NBANDS,'single','gpuArray'));
        Tbar_GFw=zeros(NBANDS,NBANDS,'single','gpuArray');
        Tbar_Fw=complex(zeros(NBANDS,NBANDS,'single','gpuArray'));
        
        %for NC=jobid*NSKIP+1+NSTART:NSKIP*NGPU:NCONF
        for NC=N_BEGIN+1:NSKIP:N_END
        %for NC=4001:NSKIP*NGPU:NCONF
            NC_c=int2str(NC);
            for NB=1:NBEAD
                %[NK,NC,NB]
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
                
                %% Calculate Hamiltonian, Overlap Matrix and Transformation Matrix
                t2=tic;
                if L_OVERL_OLD
                    [FHAM,TRANS,FOVL]=HAMILT_OLD(ROOT_DIR,NC_c,NB_c,ISP_c,...
                        NTYP,NITYP,LMMAXC,CH0,CH1,CH2,CH3,...
                        QPROJ_ORI,NPL,NPLWV,NGX,NGY,NGZ,...
                        TRANS_CORE,TRANS_c_CORE,FOVL_CORE,...
                        POSION,G_INDEX,G,HSQDTM,CITPI);
                else
                    [FHAM,TRANS]=HAMILT(ROOT_DIR,NC_c,NB_c,ISP_c,...
                        NTYP,NITYP,LMMAXC,CH0,CH1,CH2,CH3,...
                        NPL,NPLWV,NGX,NGY,NGZ,...
                        QPROJ_ORI,RES,OMEGA,...
                        POSION,G_INDEX,G,HSQDTM,CITPI);
                    
                    % Rigorous correction for overlap matrix
                    FOVL=OVERL(POSION, POSION, TRANS, TRANS, A, B, ...
                        CQIJ3_fun, r, G_INDEX, CITPI, QPROJ_ORI, ...
                        NTYP,NITYP,NPL,LMAX,LMMAXC, ...
                        CH0,CH1,CH2,CH3,LPS);
                    FOVL=triu(FOVL,1)+triu(FOVL,1)'+diag(real(diag(FOVL)));
                end           

                %% Inverse Green's function under pseudo plane waves
                %FHAM=gpuArray(FHAM);FOVL=gpuArray(FOVL);
                % Green's function
                BETA=kB*TEMP;
                ETA=(2*KOMEGA+1)*pi*BETA;
                FHAM=(1i*ETA+EFERMI_Av)*FOVL-FHAM;
                clear FOVL
                %% Expand with Bloch basis
                TdT0=OVERL(POSION, POSION_EQ, TRANS, TRANS_EQ, A, B, ...
                    CQIJ3_fun, r, G_INDEX, CITPI, QPROJ_ORI, ...
                    NTYP,NITYP,NPL,LMAX,LMMAXC, ...
                    CH0,CH1,CH2,CH3,LPS);
                
                %% Predict the wavefunction
                NIS=0;
                for NT=1:NTYP
                    POS_now(:,NIS+1:NIS+NITYP(NT)) = POSION{NT};
                    NIS = NIS + NITYP(NT);
                end
%                 Wave_pred;
%                 ICALL_out=ICALL_out+1;
                
                %% Calculate the lowest eigen-pairs
                %t1=tic;
                %[BASIS_trial,E_BAND_now,DESUM]=...
                %    Eigen_DAV_fun(G,FHAM,FOVL,...
                %    BASIS_trial,NBANDS,NSIM,HSQDTM,NELE,EDIFF);
                %BASIS_trial=single(BASIS_trial);
                %E_BAND_now=single(E_BAND_now);
                %if abs(DESUM)>EDIFF
                %    fprintf(FID_OUT,'Configuration %d, Bead %d\n',[NC,NB]);
                %    fprintf(FID_OUT,'Eigen_DAV: %6.4f\n',E_BAND_now(1));
                %    fprintf(FID_OUT,'DESUM: %6.4f\n',DESUM);
                %end
                %toc(t1)
                %fprintf(FID_OUT,'Eigen_DAV: %6.4f\n',t1_end);
                
                %% Green's function
                FHAM=BASIS_k'*TdT0'*(FHAM\TdT0)*BASIS_k;
                %FHAM=BASIS_k'*TdT0'...
                %    *(BASIS_trial*diag(1./((1i*ETA+EFERMI_Av)...
                %    -E_BAND_now))*BASIS_trial')...
                %    *TdT0*BASIS_k;
                FHAM=Gbar_inv*FHAM*Gbar_inv-Gbar_inv;
                Tbar_F(:,:,NB)=FHAM(1:NBANDS,1:NBANDS);
                %clear E_BAND_now 
                clear FHAM FOVL
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
        save(fn,'Tbar_GFw','Tbar','Tbar_Fw')
        
    end
    t3_end=toc(t3);
    fprintf(FID_OUT,'K-loop: %6.4f\n',t3_end);
    fprintf(FID_OUT,' \n');
%end

fclose(FID_OUT);
