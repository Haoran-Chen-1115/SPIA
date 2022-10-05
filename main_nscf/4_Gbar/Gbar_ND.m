%% Parameter
addpath('../')
addpath('../0_Public/')
parameter;
scriptG='Gbar_Bloch_NDnscf';
scriptT='Tbar_Bloch_NDnscf';
jbatch=cell(1,NGPU);

%if exist('Gbar.mat','file')==2
%    for i=1:NGPU
%        delete(['GREEN_',int2str(i),'.mat']);
%    end
%    delete('Gbar.mat');delete('Gbar_inv.mat');
%end

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

%% Collect Gbar
%Gbars=cell(1,NKPTS_sym);
%Gbar_diag=cell(1,NKPTS_sym);
%Gbar_inv=cell(1,NKPTS_sym);

load('BASIS_NDsym_nscf.mat','BASIS')
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
  CONF = floor(linspace(NSTART,NCONF,NGPU+1));
  %load('BASIS_NDsym_nscf.mat','BASIS')
  %BASIS=BASIS(NK_sym);
  for jobid=1:NGPU
      N_BEGIN = CONF(jobid);
      N_END = CONF(jobid+1);
      workspace = struct('jobid', jobid, 'EFERMI_Av',EFERMI_Av,'NCL',NCL,'NK_BEGIN',NK_BEGIN,'NK',NK,'N_BEGIN',N_BEGIN,'N_END',N_END,'BASIS',BASIS(NK_sym));
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
  CONF = floor(linspace(NSTART,NCONF,NGPU+1));
  for jobid = 1:NGPU
      N_BEGIN = CONF(jobid);
      N_END = CONF(jobid+1);
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

end

if NK==0
   %% Dispatch
  CONF = floor(linspace(NSTART,NCONF,NGPU+1));
  for jobid=1:NGPU
      N_BEGIN = CONF(jobid);
      N_END = CONF(jobid+1);
      workspace = struct('jobid', jobid, 'EFERMI_Av',EFERMI_Av,'NCL',NCL,'NK_BEGIN',NK_BEGIN,'NK',NK,'N_BEGIN',N_BEGIN,'N_END',N_END,'BASIS',BASIS(NK_sym));
      jbatch{jobid}=batch(scriptT,'Workspace',workspace);
  end

  for jobid=1:NGPU
      wait(jbatch{jobid})
  end
 
  load(['FERMI_Surface_',int2str(NK),'.mat'])
  Gamma=zeros(NBANDS,NBANDS,NBEAD,'single');
  Tbars=zeros(NBANDS,NBANDS,'single');
  DTbars=zeros(NBANDS,NBANDS,NBEAD,'single');

  N_CAC=0;
  CONF = floor(linspace(NSTART,NCONF,NGPU+1));
  for jobid = 1:NGPU
      N_BEGIN = CONF(jobid);
      N_END = CONF(jobid+1);
      N_CAC = N_CAC + length(N_BEGIN+1:NSKIP:N_END);
  end

  for jobid=1:NGPU
      fn=['Tbar_GFw_',int2str(NK),'_',int2str(jobid),'.mat'];
      load(fn)
      Gamma=Gamma+Tbar_GFw;
      Tbars=Tbars+Tbar;
      DTbars=DTbars+Tbar_Fw;
  end
  Gamma=Gamma./N_CAC;
  Tbars=Tbars./N_CAC;
  DTbars=DTbars./N_CAC;
  toc;

  save(['../Tbar_',int2str(NK),'.mat'],'Tbars','EFERMI_Av','VKPT','VKPT_sym','k_weight','-v7.3')
  save(['../Gamma_',int2str(NK),'.mat'],'Gamma','DTbars','EFERMI_Av','VKPT','VKPT_sym','k_weight','-v7.3')

  clear Gamma DTbars Tbars

  for jobid=1:NGPU
    delete(['../Tbar_GFw_',int2str(NK),'_',int2str(jobid),'.mat'])
  end
end
