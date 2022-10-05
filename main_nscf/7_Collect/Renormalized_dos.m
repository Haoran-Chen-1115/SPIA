addpath('../')
addpath('../0_Public')
parameter;
gpuDevice(2);
load('BASIS_NDsym_nscf.mat')

%% Read in informations
[A,B,OMEGA,...
    NKPTS,ISPIN,NTYP,NITYP,NCPU,...
    NELE,NBANDS,VKPT]=RD_HEAD(ROOT_DIR,NELE_TYP,NBANDS_mul);

%% Load in equilibrium positions
POSION_EQ=cell(1,NTYP);
for NT=1:NTYP
    % Phase's space cost is huge. Instead read positions.
    POSS=fopen([ROOT_DIR,'/BEAD_1_sym/POS_1_',int2str(NT)]);
    POSION_EQ{NT}=fread(POSS,[3,NITYP(NT)],'double');
    fclose(POSS);
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

% %% Determine the denser k-mesh, and find irreducible k-points
% % k-mesh of the supercell
% [VKPT_sym,NKPTS_sym,VKPT_f,~,SYM_OP_TRUE]...
%     =Kpoints(N_k_nscf,A,NTYP,NITYP,POSION_EQ,1);
% N_SYM_TRUE = length(SYM_OP_TRUE);
% 
% SYM_OP_INV = SYM_INV(SYM_OP_TRUE,A);
% 
% SYM_OP_min(:,:,1) = eye(3);
% SYM_OP_min(:,:,2) = -eye(3); % minimal symmetry with inversion symmetry
% [VKPT,NKPTS_nosym,~]...
%     =Kpoints_gsym(VKPT_f,A,SYM_OP_min);
% 
% [~,VKPT_sym_idx] = ismembertol(...
%     VKPT_sym(:,1:3),VKPT(:,1:3),'ByRows',true);
% 
% k_weight = VKPT_sym(:,4)';
% VKPT_sym = VKPT_sym(:,1:3)';
% VKPT = VKPT(:,1:3)';

%%
NBAND_s=2*NELE;
E_BAND_ALL = zeros(NBANDS,NKPTS_sym);
for NK=1:NKPTS_sym
    NK_sym=NK;
    %NK_sym=VKPT_sym_idx(NK);
    
    VKPT_NK=(BASIS(NK_sym).VKPT_red.').*BASIS(NK).Factor_sym;

%     load(['Gbar_inv_',int2str(NK)])
%     load('EFERMI.mat')
%     E_BAND_re = real(EFERMI_Av-Gbar_inv);
    load(['Gbar_',int2str(NK),'.mat'])
    Gbars=gpuArray(Gbars);
    Gbar_inv = inv(Gbars);
    Gbar_inv = -(Gbar_inv + Gbar_inv')/2;
    Gbar_inv = eig(Gbar_inv);
    load('EFERMI.mat')
    E_BAND_re = real(EFERMI_Av+Gbar_inv);
    E_BAND_ALL(:,NK) = gather(E_BAND_re(1:NBANDS));
end

%% Determine the Fermi energy using Gaussian smearing
NBANDS_s = 2*NELE;

degauss = 0.02*HSQDTM/(AUTOA^2); % Gaussian smearing in unit of eV
% Mainly to make a comparison with traditional calculation
sqrtpm1 = 1/sqrt(pi);
w1gauss = @(x,d) erfc(-x/d)/2 + x.*exp(-min(x.^2/d^2,200)) * sqrtpm1/2;

E_up=-1E8; E_low=1E8;
eps = 1E-4;
for NK = 1:NKPTS_sym
    NK_sym=NK;
    
    E_up = max([E_up;BASIS(NK_sym).E_BAND_ALL(1:NBANDS_s)]);
    E_low = min([E_low;BASIS(NK_sym).E_BAND_ALL(1:NBANDS_s)]);
end

for NI = 1:1000
    Ef = (E_up + E_low)/2;
    sumk = 0;
    for NK = 1:NKPTS_sym
        sumk = sumk +...
            sum(w1gauss(Ef-E_BAND_ALL(:,NK),degauss))...
            *k_weight(NK)/sum(k_weight(1:NKPTS_sym));
    end
    
    if abs(sumk-NELE/2) < eps
        EFERMI_re = Ef;
        break
    elseif sumk-NELE/2 < -eps
        E_low = Ef;
    else
        E_up = Ef;
    end
end
% Calculate band weight and DOS at Fermi energy
w0gauss = @(x,d) exp(-min(x.^2/d^2,200)) * sqrtpm1;
wk1 = w0gauss( EFERMI_re - E_BAND_ALL ,degauss)/degauss;

dos_ef = 0;
for NK=1:NKPTS_sym
    dos_ef = dos_ef + sum(k_weight(NK)*wk1(:,NK))/sum(k_weight);
end

%%
E_BAND_ori=zeros(NBANDS,NKPTS_sym);
for NK=1:NKPTS_sym
    NK_sym=NK;
    
    E_BAND_ori(:,NK) = BASIS(NK_sym).E_BAND_ALL(1:NBANDS);
end

E_up=-1E8; E_low=1E8;
eps = 1E-4;
for NK = 1:NKPTS_sym
    NK_sym=NK;
    
    E_up = max([E_up;BASIS(NK_sym).E_BAND_ALL(1:NBANDS_s)]);
    E_low = min([E_low;BASIS(NK_sym).E_BAND_ALL(1:NBANDS_s)]);
end

for NI = 1:1000
    Ef = (E_up + E_low)/2;
    sumk = 0;
    for NK = 1:NKPTS_sym
        sumk = sumk +...
            sum(w1gauss(Ef-E_BAND_ori(:,NK),degauss))...
            *k_weight(NK)/sum(k_weight(1:NKPTS_sym));
    end
    
    if abs(sumk-NELE/2) < eps
        EFERMI_ori = Ef;
        break
    elseif sumk-NELE/2 < -eps
        E_low = Ef;
    else
        E_up = Ef;
    end
    
end
% Calculate band weight and DOS at Fermi energy
w0gauss = @(x,d) exp(-min(x.^2/d^2,200)) * sqrtpm1;
wk1 = w0gauss( EFERMI_ori - E_BAND_ori ,degauss)/degauss;

dos_ef_ori = 0;
for NK=1:NKPTS_sym
    dos_ef_ori = dos_ef_ori + sum(k_weight(NK)*wk1(:,NK))/sum(k_weight);
end

save('dos2_r.mat','k_weight','NKPTS_sym','VKPT_sym','degauss','OMEGA',...
    'E_BAND_ori','EFERMI_ori','dos_ef_ori',...
    'E_BAND_ALL','EFERMI_re','dos_ef');