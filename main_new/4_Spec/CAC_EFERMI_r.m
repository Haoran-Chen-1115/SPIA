%% Parameter
addpath('../')
addpath('../0_Public/')
%load('BASIS_NDsym_struct.mat')
parameter;
path='../';

%%
[A,B,OMEGA,...
    NKPTS,ISPIN,NTYP,NITYP,NCPU,...
    NELE,VKPT,NELE_TYP,POMASS]=...
    RD_HEAD(ROOT_DIR);
%NBANDS = max(ceil(NELE * 0.6),ceil(NELE/2+sum(NITYP)/2));
NBANDS = ceil(NELE/2 * NBANDS_mul_Ef);

%% Load in equilibrium positions
POSION_EQ=cell(1,NTYP);
for NT=1:NTYP
    % Phase's space cost is huge. Instead read positions.
    POSS=fopen([ROOT_DIR,'/BEAD_1_sym/POS_1_',int2str(NT)]);
    POSION_EQ{NT}=fread(POSS,[3,NITYP(NT)],'double');
    fclose(POSS);
end

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

if L_POS_new
    load('../POSION_new.mat')
    POSION_prim=POSION_new;
else
    % Primitive cell ion position
    POSION_prim=cell(1,NTYP);
    for NT=1:NTYP
        % Phase's space cost is huge. Instead read positions.
        POSS=fopen([ROOT_DIR,'/BEAD_1_primitive/POS_1_',int2str(NT)]);
        POSION_prim{NT}=fread(POSS,[3,NITYP_prim(NT)],'double');
        fclose(POSS);
    end
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

%% Determine the k-path in the IBZ of the base diagonal supercell
eps=1E-5;
 
if L_path
    error('Uniform grid must be used to calculate DOS');
end
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
    [IBRAV,N_SYM_TRUE,SYM_OP_TRUE,TRANS_ROT]...
        = Symmetry(A,POSION_new,NTYP,NITYP_new);
else
    [IBRAV,N_SYM_TRUE,SYM_OP_TRUE,TRANS_ROT]...
        = Symmetry(A,POSION_EQ,NTYP,NITYP);
end
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

SYM_OP_INV = SYM_INV(SYM_OP_TRUE,A);
%% Correspondence between k-points in different cells
CONF_CAC = NSTART+1:NSKIP:NCONF;
NCONF_CAC = length(CONF_CAC);
CONF_idx = floor(linspace(0,NCONF_CAC,NGPU+1));

E_all=zeros(NBANDS,NBEAD,NCONF_CAC,NKPTS_sym);

for NK=1:NKPTS_sym
    
    %NCAC_s=0;
    load(['../Eigen_',int2str(NK),'.mat'],'E_t','NCAC_s');
    E_all(:,:,:,NK)=E_t;
    %for jobid=1:NGPU
    %    N_BEGIN = CONF_CAC(CONF_idx(jobid)+1)-1;
    %    N_END = CONF_CAC(CONF_idx(jobid+1));
    %    NCAC=length(N_BEGIN+1:NSKIP:N_END);
    %    
    %    load(['../Eigen_',int2str(NK),'_',int2str(jobid)],'E_t');
    %    E_all(:,:,NCAC_s+1:NCAC_s+NCAC,NK)=E_t;
    %    NCAC_s=NCAC_s+NCAC;
    %end
end
E_all_sq=reshape(E_all,[],NKPTS_sym);

%%
%nisg=1;
%degauss = 0.01*nisg*HSQDTM/(AUTOA^2); % Gaussian smearing in unit of eV
degauss = 0.2;
% Mainly to make a comparison with traditional calculation
sqrtpm1 = 1/sqrt(pi);
w1gauss = @(x,d) erfc(-x/d)/2 + x.*exp(-min(x.^2/d^2,200)) * sqrtpm1/2;

eps = 1E-4;
E_up=max(max(E_all_sq));
E_low=min(min(E_all_sq));
    
for NI = 1:1000
    Ef = (E_up + E_low)/2;
    sumk = 0;
    for NK = 1:NKPTS_sym
        sumk = sumk +...
            sum(w1gauss(Ef-E_all_sq(:,NK),degauss))...
            *k_weight(NK)/sum(k_weight(1:NKPTS_sym));
    end
    sumk=sumk/NCAC_s/NBEAD;
    
    if abs(sumk-NELE/2) < eps
        EFERMI_gauss = Ef;
        break
    elseif sumk-NELE/2 < -eps
        E_low = Ef;
    else
        E_up = Ef;
    end
    
end
EFERMI_Av=EFERMI_gauss;

%%
NBANDS_s = 2*NELE;

degauss=0.2;
degauss2=degauss/2; % interval

%Ndown = ceil(abs(min(min(E_all_sq))-EFERMI_Av+degauss2/2)/degauss2);
Ndown = ceil(abs(-10-EFERMI_Av+degauss2/2)/degauss2);
centers = -Ndown*degauss2:degauss2:20;

% Calculate band weight and DOS at Fermi energy
w0gauss = @(x,d) exp(-min(x.^2/d^2,200)) * sqrtpm1;
dos_ef=zeros(1,length(centers));

for Ne = 1:length(centers)
    wk1 = w0gauss( EFERMI_Av +centers(Ne) - E_all_sq ,degauss)/degauss;

    for NK=1:NKPTS_sym
        dos_ef(Ne) = dos_ef(Ne) + sum(k_weight(NK)*wk1(:,NK))/sum(k_weight);
    end
    dos_ef(Ne)=dos_ef(Ne)/NCAC_s/NBEAD;
end

save('dos_Av','EFERMI_Av','NBEAD','NCAC_s','OMEGA','centers','degauss','dos_ef')
save('../EFERMI.mat','EFERMI_Av')
