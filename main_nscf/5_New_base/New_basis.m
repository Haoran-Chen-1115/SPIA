addpath('../')
addpath('../0_Public')
parameter;
%gpuDevice(1);
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

SYM_OP_INV = SYM_INV(SYM_OP_TRUE,A);

[VKPT_sym,NKPTS_sym,~]...
    =Kpoints_gsym(VKPT_NDsc,A,SYM_OP_TRUE);

k_weight=VKPT_sym(:,4)';
VKPT_sym=VKPT_sym(:,1:3)';

%%
Gbar_inv_all=cell(1,NKPTS_sym);
for NK=1:NKPTS_sym
    NK_sym=NK;
    
    VKPT_NK=(BASIS(NK_sym).VKPT_red.').*BASIS(NK).Factor_sym;
    load(['Gbar_',int2str(NK),'.mat'])
    load('EFERMI.mat')
    
    % Heff = (1i*ETA + EFERMI_Av)*eye(length(Gbars))-(Gi+Gi')/2;
    WAVE=cell(1,BASIS(NK).NKPTS_red);E_BAND=cell(1,BASIS(NK).NKPTS_red);
    WAVE_gamma=zeros(length(Gbars));
    Gi_new = zeros(length(Gbars));
    for NK_s=1:BASIS(NK).NKPTS_red
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
            idx = BASIS(NK).G_kp{NK_s};
            NBANDS_s=length(idx); % =NPL_s

            % Expand in orthonormal bases
            Gbars_s = ...
                BASIS(NK).WAVE{NK_s}'*Gbars(idx,idx)*BASIS(NK).WAVE{NK_s};
            Gi_s = inv(Gbars_s);

            Heff = EFERMI_Av*eye(length(Gbars_s))-(Gi_s+Gi_s')/2;

            [BASIS_re,E_BAND_re] = eig(Heff);
            WAVE{NK_s} = BASIS(NK).WAVE{NK_s}*BASIS_re;
            E_BAND{NK_s} =diag(E_BAND_re);

            idx2 = find(BASIS(NK).Band_kp==NK_s);
            Gi_new(idx2,idx2) = BASIS_re'*Gi_s*BASIS_re;
        end
    end
    
    % Band component analysis
    ss=0;
    for NK_s=1:BASIS(NK).NKPTS_red
        ss=ss+length(BASIS(NK).G_kp{NK_s});
    end
    
    E_BAND_ALL=[];
    for NK_s=1:BASIS(NK).NKPTS_red
        E_BAND_ALL=[E_BAND_ALL;E_BAND{NK_s}];
    end
    
    [E_BAND_ALL,E_idx]=sort(E_BAND_ALL,'ascend');
    
    Band_idx = zeros(1,length(Gbars));
    Ni=0; 
    for NK_s=1:BASIS(NK).NKPTS_red
        for i=Ni+1:Ni+length(E_BAND{NK_s})
            Band_idx(i)=find(E_idx==i);
        end
        Ni=Ni+length(E_BAND{NK_s});
    end
    
    BASIS(NK).E_BAND=E_BAND;
    BASIS(NK).WAVE=WAVE;
    BASIS(NK).E_BAND_ALL=E_BAND_ALL;
    BASIS(NK).E_idx=E_idx;
    BASIS(NK).Band_idx=Band_idx;

    Gbar_inv_all{NK} = diag(Gi_new);
end

save('../BASIS_NDsym_new.mat','BASIS','VKPT','-v7.3')

%%
NBANDS_ss = NELE * 2;
E_BAND_ALL = zeros(NBANDS_ss,NKPTS_sym);
for NK=1:NKPTS_sym
    NK_c=int2str(NK);
    NK_sym=NK;
    
    E_BAND_ALL(:,NK) ...
        =BASIS(NK_sym).E_BAND_ALL(1:NBANDS_ss);
end

nisg=2;
degauss = 0.01*nisg*HSQDTM/(AUTOA^2); % Gaussian smearing in unit of eV
% Mainly to make a comparison with traditional calculation
sqrtpm1 = 1/sqrt(pi);
w1gauss = @(x,d) erfc(-x/d)/2 + x.*exp(-min(x.^2/d^2,200)) * sqrtpm1/2;

E_up=-1E8; E_low=1E8;
eps = 1E-4;
for NK = 1:NKPTS_sym
    NK_sym=NK;
    
    E_up = max([E_up;BASIS(NK_sym).E_BAND_ALL(1:NBANDS_ss)]);
    E_low = min([E_low;BASIS(NK_sym).E_BAND_ALL(1:NBANDS_ss)]);
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
        EFERMI_gauss = Ef;
        break
    elseif sumk-NELE/2 < -eps
        E_low = Ef;
    else
        E_up = Ef;
    end
    
end

%%
for NK=1:NKPTS_sym
    Gbar_inv = Gbar_inv_all{NK};
    Gbar_inv = Gbar_inv - EFERMI_Av + EFERMI_gauss;
    Gbar_diag = 1./Gbar_inv;

    save(['../Gbar_inv_new_',int2str(NK),'.mat'],'Gbar_diag','Gbar_inv')
end

EFERMI_Av = EFERMI_gauss;
save('../EFERMI_new.mat','EFERMI_Av')

