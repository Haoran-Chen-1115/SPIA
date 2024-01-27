addpath('../')
addpath('../0_Public')
parameter;
%gpuDevice(1);
load('BASIS_NDsym_new.mat')
path='../';

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
KP=cell(5,1);
% W-L
KP{1}=[linspace(0.5,0.5,6).',linspace(0.25,0.5,6).',linspace(0.75,0.5,6).'];

% L-G
KP{2}=[linspace(0.5,0,11).',linspace(0.5,0,11).',linspace(0.5,0,11).'];
KP{2}=KP{2}(2:end,:);

% G-X
KP{3}=[linspace(0,0.5,11).',linspace(0,0,11).',linspace(0,0.5,11).'];
KP{3}=KP{3}(2:end,:);
%KP{3}=[linspace(0,0.5,16).',linspace(0.5,0,16).',linspace(0,0.5,16).'];
%KP{3}=KP{3}(2:8,:);
%KP{3}(end+1,:)=[0.25,0.25,0.25];

% X-W
KP{4}=[linspace(0.5,0.5,6).',linspace(0,0.25,6).',linspace(0.5,0.75,6).'];
KP{4}=KP{4}(2:end,:);
%KP{4}=[linspace(0.5,0,16).',linspace(0.5,0,16).',linspace(0.5,0,16).'];
%KP{4}=KP{4}(9:end,:);

% W-K
KP{5}=[linspace(0.5,0.375,6).',linspace(0.25,0.375,6).',linspace(0.75,0.75,6).'];
KP{5}=KP{5}(3:end,:);

VKPT_Band=cell2mat(KP);
ll=zeros(5,1);xl=cell(1,5);
ll(1)=sqrt(sum(((KP{1}(end,:)-KP{1}(1,:))*B).^2));
xl{1}=linspace(0,ll(1),size(KP{1},1));
LL=ll(1);
for i=2:5
    ll(i)=sqrt(sum(((KP{i-1}(end,:)-KP{i}(end,:))*B).^2));
    xl{i}=linspace(0,ll(i),size(KP{i},1)+1)+LL;
    xl{i}=xl{i}(2:end);
    LL=LL+ll(i);
end
xl=cell2mat(xl);

xt=zeros(1,6);
for i=2:6
    xt(i)=ll(i-1)+xt(i-1);
end
xt=xt/LL;

% VKPT_Band(13,:)=[0.5,0.5,0.5];

NBAND_s=NELE/det(diag(N_sc_base));
E_BAND_ori=zeros(NBAND_s,length(VKPT_Band));
E_BAND_ori(:,:)=NaN;
eps=1E-3;

VKPT_Band = mod(VKPT_Band,1);
VKPT_Band(VKPT_Band>0.5+eps) = VKPT_Band(VKPT_Band>0.5+eps) - 1;

%%
E_BAND=zeros(NBAND_s,length(VKPT_Band));
E_BAND_ALL = zeros(NBANDS,NKPTS_sym);

E_BAND(:,:)=NaN;
for NK=1:NKPTS_sym
    NK_sym=NK;
    
    VKPT_NK=(BASIS(NK_sym).VKPT_red.').*BASIS(NK).Factor_sym;

    load([path,'Gbar_inv_new_',int2str(NK)])
    load([path,'EFERMI_new.mat'])
    E_BAND_re = real(EFERMI_Av-Gbar_inv);
    E_BAND_ALL(:,NK) = E_BAND_re(1:NBANDS);
    E_BAND_s = cell(1,BASIS(NK_sym).NKPTS_red);
    Ni=0;
    for NK_s=1:BASIS(NK_sym).NKPTS_red
        for i=Ni+1:Ni+length(BASIS(NK_sym).G_kp{NK_s})
            if BASIS(NK_sym).Band_idx(i)<=NBANDS
                E_BAND_s{NK_s}(i-Ni) = E_BAND_re(BASIS(NK_sym).Band_idx(i));
            end
        end
        Ni=Ni+length(BASIS(NK_sym).G_kp{NK_s});
    end


    for NS=1:N_SYM_TRUE
        VKPT_SYM = (N_sc_ex\ (SYM_OP_INV(:,:,NS)' * N_sc_ex*VKPT_NK.')).';
        VKPT_SYM = mod(VKPT_SYM,1);
        VKPT_SYM(VKPT_SYM>0.5+eps) = VKPT_SYM(VKPT_SYM>0.5+eps) - 1;

        [~,VKPT_ia]=ismembertol(...
            VKPT_Band,VKPT_SYM,1E-3,...
            'ByRows',true);
        idx=find(VKPT_ia~=0);
        VKPT_ia=VKPT_ia(VKPT_ia~=0);
        for i=1:length(VKPT_ia)
            E_BAND(:,idx(i))=E_BAND_s{VKPT_ia(i)}(1:NBAND_s);
            E_BAND_ori(:,idx(i))=BASIS(NK).E_BAND{VKPT_ia(i)}(1:NBAND_s);
        end
        
        [~,VKPT_ia]=ismembertol(...
            VKPT_Band,-VKPT_SYM,1E-3,...
            'ByRows',true);
        idx=find(VKPT_ia~=0);
        VKPT_ia=VKPT_ia(VKPT_ia~=0);
        for i=1:length(VKPT_ia)
            E_BAND(:,idx(i))=E_BAND_s{VKPT_ia(i)}(1:NBAND_s);
            E_BAND_ori(:,idx(i))=BASIS(NK).E_BAND{VKPT_ia(i)}(1:NBAND_s);
        end
    end
%     VKPT_NK=(BASIS(NK).VKPT_red.').*BASIS(NK).Factor_sym;
%     VKPT_NK(VKPT_NK+0.5<1E-5)=VKPT_NK(VKPT_NK+0.5<1E-5)+1;
%     [~,VKPT_ia]=ismembertol(...
%         VKPT_sym,VKPT_NK,...
%         'ByRows',true);
%     idx=find(VKPT_ia~=0);
%     VKPT_ia=VKPT_ia(VKPT_ia~=0);
%     for i=1:length(VKPT_ia)
%         E_BAND(:,idx(i))=BASIS(NK).E_BAND{VKPT_ia(i)}(1:6);
%     end  
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
    
%%
E_BAND_ALL=zeros(NBANDS,NKPTS_sym);
for NK=1:NKPTS_sym
    NK_sym=NK;
    
    E_BAND_ALL(:,NK) = BASIS(NK_sym).E_BAND_ALL(1:NBANDS);
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
            sum(w1gauss(Ef-E_BAND_ALL(:,NK),degauss))...
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
save([path,'band'],'xl','LL','xt','E_BAND_ori','EFERMI_ori','E_BAND','EFERMI_re')
%%
load([path,'band.mat'])
idx2=find(~isnan(E_BAND_ori(1,:)));
figure;plot(xl(idx2)/LL,E_BAND_ori(:,idx2)-EFERMI_ori,'Color','k')
hold on
%plot(xl(idx2)/LL,E_BAND(:,idx2)-EFERMI_re,'--','Color','r')

% ylim([-26,10])
%yline(0,':')
plot([0 1],[0 0],'--k')
xticks(xt);
xticklabels({'W','L','\Gamma','X','W','K'})
yy=ylim;
for i=2:length(xt)
    plot([xt(i) xt(i)],ylim,'-k','LineWidth',1);
end
