%% Parameter
addpath('../')
addpath('../0_Public/')
parameter;
scriptG='Gbar_Bloch_NDnscf';
scriptT='Tbar_Bloch_NDnscf';
jbatch=cell(1,NGPU);

[A,B,OMEGA,...
    NKPTS,ISPIN,NTYP,NITYP,NCPU,...
    NELE,VKPT,NELE_TYP,POMASS]=...
    RD_HEAD(ROOT_DIR);

%% Load in equilibrium positions
POSION_EQ=cell(1,NTYP);
for NT=1:NTYP
    % Phase's space cost is huge. Instead read positions.
    POSS=fopen([ROOT_DIR,'/BEAD_1_sym/POS_1_',int2str(NT)]);
    POSION_EQ{NT}=fread(POSS,[3,NITYP(NT)],'double');
    fclose(POSS);
end

% Calculate Average Fermi Energy
load('EFERMI.mat');

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

if L_path
    KP=cell(size(VK_hsym,1)-1,1);
    for i=1:size(VK_hsym,1)-1
        KP{i}=[...
            linspace(VK_hsym(i,1),VK_hsym(i+1,1),NK_line(i)+1).',...
            linspace(VK_hsym(i,2),VK_hsym(i+1,2),NK_line(i)+1).',...
            linspace(VK_hsym(i,3),VK_hsym(i+1,3),NK_line(i)+1).'];
        if i>1
            KP{i}=KP{i}(2:end,:);
        end
    end
    VKPTS_full = cell2mat(KP).';
    NKPTS_full = size(VKPTS_full,2);
    VKPTS_full(VKPTS_full>0.5+eps) =...
        VKPTS_full(VKPTS_full>0.5+eps) - 1;
else
    NKPTS_full = N_k_nscf(1)*N_k_nscf(2)*N_k_nscf(3);
    
    K1 = (0:N_k_nscf(1)-1)/N_k_nscf(1);
    K2 = (0:N_k_nscf(2)-1)/N_k_nscf(2);
    K3 = (0:N_k_nscf(3)-1)/N_k_nscf(3);
    
    clear VKPTS_full
    VKPTS_full(1,:) = reshape(repmat(K1,1,N_k_nscf(2)*N_k_nscf(3)),1,[]);
    VKPTS_full(2,:) = reshape(repmat(K2,N_k_nscf(1),N_k_nscf(3)),1,[]);
    VKPTS_full(3,:) = reshape(repmat(K3,N_k_nscf(1)*N_k_nscf(2),1),1,[]);
    VKPTS_full(VKPTS_full>0.5+eps) = VKPTS_full(VKPTS_full>0.5+eps) - 1;
end
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

%[~,idx]=ismembertol(VKPT_nosym',VKPT_sym','ByRows',true);
%idx=find(~idx);
%VKPT_cal=VKPT_nosym(:,idx);
%NKPTS_cal=length(idx);

disp([int2str(NKPTS_sym),' k points in all:']);
for i=1:NKPTS_sym
    disp([int2str(i),':(',num2str(VKPT_sym(1,i)),',',num2str(VKPT_sym(2,i)),',',num2str(VKPT_sym(3,i)),')']);
end
%% Collect Gbar
%Gbars=cell(1,NKPTS_sym);
%Gbar_diag=cell(1,NKPTS_sym);
%Gbar_inv=cell(1,NKPTS_sym);

if LBloch
    load('BASIS_NDsym_nscf.mat','BASIS')
else
    BASIS=0;
end


NK_group = floor(linspace(1,NKPTS_sym+1,NCLUSTER+1));
NK_BEGIN = NK_group(NCL);
NK_END = NK_group(NCL+1)-1;
patchJobStorageLocation;
for NK=NK_BEGIN:NK_END
%for NK=4:4
  tic;
  NC_c='1';
  %NK_c=int2str(VKPT_sym_idx(NK));
  %NK_sym=VKPT_sym_idx(NK);
  NK_c=int2str(NK);
  NK_sym=NK;
  [G_INDEX,G,G2,NPL,NPLWV,NGX,NGY,NGZ]=...
      RD_INDEX(ROOT_DIR,NC_c,'1',VKPT_sym(:,NK_sym),B,ENCUT);

  NBANDS=NPL;

  %% Dispatch
  %CONF = floor(linspace(NSTART,NCONF,NGPU+1));
  CONF_CAC = NSTART+1:NSKIP:NCONF;
  NCONF_CAC = length(CONF_CAC);
  CONF_idx = floor(linspace(0,NCONF_CAC,NGPU+1));
  %load('BASIS_NDsym_nscf.mat','BASIS')
  %BASIS=BASIS(NK_sym);
  for jobid=1:NGPU
      %N_BEGIN = CONF(jobid);
      %N_END = CONF(jobid+1);
      N_BEGIN = CONF_CAC(CONF_idx(jobid)+1)-1;
      N_END = CONF_CAC(CONF_idx(jobid+1));
      if LBloch
          workspace = struct('jobid', jobid,...
              'EFERMI_Av',EFERMI_Av,'NCL',NCL,...
              'NK_BEGIN',NK_BEGIN,'NK',NK,...
              'N_BEGIN',N_BEGIN,'N_END',N_END,'BASIS',BASIS(NK_sym));
      else
          workspace = struct('jobid', jobid,...
              'EFERMI_Av',EFERMI_Av,'NCL',NCL,...
              'NK_BEGIN',NK_BEGIN,'NK',NK,...
              'N_BEGIN',N_BEGIN,'N_END',N_END,'BASIS',BASIS);
      end 
      jbatch{jobid}=batch(scriptG,'Workspace',workspace);
  end

  for jobid=1:NGPU
      wait(jbatch{jobid})
  end

  Gbars=zeros(NBANDS);
  for jobid=1:NGPU
      fn=['GREEN_',int2str(NK),'_',int2str(jobid),'.mat'];
      load(fn)
      Gbars=Gbars+GREEN;
  end

  N_CAC=0;
  %CONF = floor(linspace(NSTART,NCONF,NGPU+1));
  CONF_CAC = NSTART+1:NSKIP:NCONF;
  NCONF_CAC = length(CONF_CAC);
  CONF_idx = floor(linspace(0,NCONF_CAC,NGPU+1));
  for jobid = 1:NGPU
      %N_BEGIN = CONF(jobid);
      %N_END = CONF(jobid+1);
      N_BEGIN = CONF_CAC(CONF_idx(jobid)+1)-1;
      N_END = CONF_CAC(CONF_idx(jobid+1));

      N_CAC = N_CAC + length(N_BEGIN+1:NSKIP:N_END);
  end

  Gbars=Gbars./(NBEAD*N_CAC);
  Gbars=gather(Gbars);
  Gbar_diag=diag(Gbars);
  Gbar_inv=1./Gbar_diag;
  toc;

  save(['../Gbar_',int2str(NK),'.mat'],'Gbars','EFERMI_Av','VKPT','VKPT_sym','k_weight','-v7.3')
  save(['../Gbar_inv_',int2str(NK),'.mat'],'Gbar_inv','Gbar_diag','VKPT','VKPT_sym','k_weight','-v7.3')

  clear Gbars
  for jobid=1:NGPU
    delete(['../GREEN_',int2str(NK),'_',int2str(jobid),'.mat'])
  end
  
  save(['../Gbar_OK_',int2str(NK),'.mat'],'jobid')

end
gpuDevice([])
