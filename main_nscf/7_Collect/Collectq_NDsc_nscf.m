%% Initialization
clear
addpath('../')
addpath('../0_Public/')
parameter;
maxNumCompThreads(30);
%gpuDevice;

BETA=kB*TEMP;

tic;
[A,B,OMEGA,...
    NKPTS,ISPIN,NTYP,NITYP,NCPU,...
    NELE,NBANDS,VKPT]=RD_HEAD(ROOT_DIR,NELE_TYP,NBANDS_mul);
NIONS=sum(NITYP);

%% Load in equilibrium positions
POSION_EQ=cell(1,NTYP);
for NT=1:NTYP
    % Phase's space cost is huge. Instead read positions.
    POSS=fopen([ROOT_DIR,'/BEAD_1_sym/POS_1_',int2str(NT)]);
    POSION_EQ{NT}=fread(POSS,[3,NITYP(NT)],'double');
    fclose(POSS);
end

%% Load Basis
tic;
load('BASIS_NDsym_new.mat')
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

%% Index for k and k+q
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
for NK=1:NKPTS_sym
    NK_c=int2str(NK);
    NK_sym=NK;
    
    load(['Gbar_inv_new_',int2str(NK),'.mat'],'Gbar_inv','Gbar_diag')
    load('EFERMI_new.mat')
    E_BAND_ALL(:,NK) ...
        = real(EFERMI_Av-Gbar_inv(1:NBANDS_ss));
    %E_BAND_ALL(:,NK) ...
    %    =BASIS(NK_sym).E_BAND_ALL(1:NBANDS_ss);
end

%% Transform to another representation
W11_Q1 = zeros(NBANDS_ss,NBANDS_ss,NKPTS_sym,NBEAD);
W11_Q2 = zeros(NBANDS_ss,NBANDS_ss,NKPTS_sym,NBEAD1);
for NK=1:NKPTS_sym
    tic;
    NK_c=int2str(NK);
    NK_sym=NK;
    NC_c='1';
    [G_INDEX,G,G2,NPL,NPLWV,NGX,NGY,NGZ]=...
        RD_INDEX(ROOT_DIR,'1','1',VKPT_sym(:,NK_sym),B,ENCUT);
    
    %% Solve Bethe-Salpeter Equation
    load(['Gamma_',int2str(NK),'.mat'],'Gamma','DTbars')
    load(['Gbar_inv_new_',int2str(NK),'.mat'],'Gbar_inv','Gbar_diag')
    
    % Gamma_k=gpuArray(Gamma)-abs(gpuArray(DTbars)).^2;
    Gamma_k=Gamma-abs(DTbars).^2;
    Gamma_tau=real(fft(Gamma_k(1:NBANDS,1:NBANDS,:),NBEAD,3));
    W11_tau=zeros(NBANDS,NBANDS,'single');%,'gpuArray');
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
        W11_Q1(:,:,NK,NB) = ...
            gather(W11_w(1:NBANDS_ss,1:NBANDS_ss,NB));
    end
    for NB = 1:NBEAD1
        W11_Q2(:,:,NK,NB) = ...
            W11_w2(1:NBANDS_ss,1:NBANDS_ss,NB);
    end
    toc;
end

%% Determine the q-transfer grid
VKPTS=BASIS(1).VKPT_red.*BASIS(1).Factor_sym';
DX=reshape(VKPTS(1,:)-VKPTS(1,:)',1,[]);
DY=reshape(VKPTS(2,:)-VKPTS(2,:)',1,[]);
DZ=reshape(VKPTS(3,:)-VKPTS(3,:)',1,[]);
DK=[DX;DY;DZ].*N_sc_base.';
DK=mod(DK,1);
DK(DK-0.5>1E-5)=DK(DK-0.5>1E-5)-1;
[VQPTS,idx,~]=uniquetol(DK.',1E-4,'ByRows',true,'DataScale',1);
NQPTS=size(VQPTS,1);

%% Determine the Fermi energy using Gaussian smearing
Nsigma = 10;
lambda = cell(1,Nsigma); lambda_q = cell(1,Nsigma);
lambda2_q = cell(1,Nsigma); lambda2 = cell(1,Nsigma);

phase_space = zeros(1,Nsigma); dos_ef = zeros(1,Nsigma);
rho_2 = zeros(1,Nsigma); omega2_2 = zeros(1,Nsigma);
EFERMI_gauss=zeros(1,Nsigma);
for nisg=1:Nsigma
    degauss = 0.01*nisg*HSQDTM/(AUTOA^2); % Gaussian smearing in unit of eV
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
            EFERMI_gauss(nisg) = Ef;
            break
        elseif sumk-NELE/2 < -eps
            E_low = Ef;
        else
            E_up = Ef;
        end
        
    end
    
    % Calculate band weight and DOS at Fermi energy
    w0gauss = @(x,d) exp(-min(x.^2/d^2,200)) * sqrtpm1;
    wk1 = w0gauss( EFERMI_gauss(nisg) - E_BAND_ALL ,degauss)/degauss;
    
    dos_ef(nisg) = 0;
    for NK=1:NKPTS_sym
        dos_ef(nisg) = dos_ef(nisg) +...
            sum(k_weight(NK)*wk1(:,NK))/sum(k_weight(1:NKPTS_sym));
    end
    
    % Calculate lambda
    lambda_q{nisg}=zeros(NQPTS,NBEAD/2+1);
    lambda{nisg}=zeros(1,NBEAD/2+1);
    
    lambda2_q{nisg}=zeros(NQPTS,NBEAD1/2+1);
    lambda2{nisg}=zeros(1,NBEAD1/2+1);
    phase_space(nisg)=0;
    for NK=1:NKPTS_sym
        VKPTS=BASIS(NK).VKPT_red.*BASIS(NK).Factor_sym';
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
        
        Ni=0; Band_idx_all=cell(1,BASIS(1).NKPTS_red);
        for NK_s=1:BASIS(1).NKPTS_red
            for i=Ni+1:Ni+length(BASIS(NK).E_BAND{NK_s})
                if BASIS(NK).Band_idx(i)<=NBANDS_s
                    Band_idx_all{NK_s}(i-Ni)=BASIS(NK).Band_idx(i);
                end
            end
            Ni=Ni+length(BASIS(NK).E_BAND{NK_s});
        end
        
        for NQ=1:NQPTS
            idx_q=find(idx==NQ);
            NK_start=mod(idx_q-1,BASIS(1).NKPTS_red)+1;
            NK_end=ceil(idx_q/BASIS(1).NKPTS_red);
            
            NK_start_u = unique(NK_start);
            for NK_s=1:length(NK_start_u)
                idx_kq=find(NK_start==NK_start_u(NK_s));
                NK_kq_end=NK_end(idx_kq);
                
                Band_idx_k = Band_idx_all{NK_start_u(NK_s)};
                Band_idx_kq = cell2mat(Band_idx_all(NK_kq_end));
                
                wks1=...
                    w0gauss( EFERMI_gauss(nisg) - E_BAND_ALL(Band_idx_k,NK) ...
                    ,degauss)/degauss;
                wks2=...
                    w0gauss( EFERMI_gauss(nisg) - E_BAND_ALL(Band_idx_kq,NK)...
                    ,degauss)/degauss;
                weight = k_weight(NK) * wks1 .* wks2.'...
                    /sum(k_weight(1:NKPTS_sym));
                
                for NB = 1:NBEAD/2+1
                    lambda_q{nisg}(NQ,NB) = lambda_q{nisg}(NQ,NB) ...
                        + sum(sum(weight .*...
                        W11_Q1(Band_idx_k,Band_idx_kq,NK,NB)))/BETA/dos_ef(nisg);
                end
                
                for NB = 1:NBEAD1/2+1
                    lambda2_q{nisg}(NQ,NB) = lambda2_q{nisg}(NQ,NB) ...
                        + sum(sum(weight .*...
                        W11_Q2(Band_idx_k,Band_idx_kq,NK,NB)))/BETA/dos_ef(nisg);
                end
            end
        end
    end
    
    for NB=1:NBEAD/2+1
        lambda{nisg}(NB)=sum(lambda_q{nisg}(:,NB));
    end
    
    for NB=1:NBEAD1/2+1
        lambda2{nisg}(NB)=sum(lambda2_q{nisg}(:,NB));
    end
    Ncut=NBEAD/2;
    [rho_2(nisg),omega2_2(nisg)]=Gbar_Eliashberg(lambda2{nisg},Ncut,mustar);
end

lam = zeros(Nsigma,1);
for nisg = 1:Nsigma
    lam_no(nisg) = lambda{nisg}(1);
    lam(nisg) = lambda2{nisg}(1);
end

save('Result_q_re','rho_2','omega2_2','mustar','Nsigma',...
    'lambda2','lam','lambda2_q',...
    'lambda','lam_no','lambda_q')
%save('Intermediate','W11_Q1','E_BAND_ALL','k_weight','-v7.3')

%% Determine the transition temperature through interpolating lambda
BETA=kB*TEMP;
Ncut=NBEAD/2;
ns=2;
lam_f=@(x) interp1(BETA^2*[0:NBEAD/2].^2,BETA^2*[0:NBEAD/2].^2.*lambda2{ns}(1:NBEAD/2+1),x,'pchip');

rho=100;
T_up = TEMP+50;
T_down = TEMP-50;
while abs(rho)>1E-3
    
    T_tmp = (T_up+T_down)/2;
    BETA_tmp = kB*T_tmp;
    
    nu2 = (BETA_tmp)^2*[0:NBEAD/2].^2;
    lambda_nu = lam_f(nu2)/BETA_tmp^2;
    
    lambda_nu(1)=lambda2{ns}(1);
    lambda_nu(2:end) = lambda_nu(2:end)./([1:NBEAD/2].^2);
    lambda_nu(end+1:end+NBEAD/2-1)=lambda_nu(end-1:-1:2);
    
    
    rho=Gbar_Eliashberg2(lambda_nu,Ncut,mustar,omega2_2(ns)*BETA/BETA_tmp);
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
T_c=T_tmp
