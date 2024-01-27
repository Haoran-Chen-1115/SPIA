%% Parameter
addpath('../')
addpath('../0_Public/')
parameter;
if exist('LCAC_Ef','var')
    L_path=~LCAC_Ef;
end

scriptG='Spec_cac';
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
if L_path || I_Ef==1
    load('EFERMI.mat');
else
    load('EFERMI_old.mat');
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

%% Determine the k-path in the IBZ of the base diagonal supercell
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
    POSION_EQ = POSION_new;
    NITYP_EQ = NITYP_new;
else
    NITYP_EQ = NITYP;
end

[IBRAV,N_SYM_TRUE,SYM_OP_TRUE,TRANS_ROT]...
    = Symmetry(A,POSION_EQ,NTYP,NITYP_EQ);
    
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
end

%% For divide k points in the IBZ
% if L_path
%     min_space1=1;min_space2=1;min_space3=1;
%     for i=1:size(VK_hsym,1)-1
%         if abs(VK_hsym(i,1)-VK_hsym(i+1,1))>eps
%             min_space1=min(min_space1,...
%                 abs((VK_hsym(i,1)-VK_hsym(i+1,1))/NK_line));
%         end
%         if abs(VK_hsym(i,2)-VK_hsym(i+1,2))>eps
%             min_space2=min(min_space2,...
%                 abs((VK_hsym(i,2)-VK_hsym(i+1,2))/NK_line));
%         end
%         if abs(VK_hsym(i,3)-VK_hsym(i+1,3))>eps
%             min_space3=min(min_space3,...
%                 abs((VK_hsym(i,3)-VK_hsym(i+1,3))/NK_line));
%         end
%     end
%     NKPTS_1 = 1/min_space1;
%     NKPTS_2 = 1/min_space2;
%     NKPTS_3 = 1/min_space3;
%     [VKPTS_f,NKPTS_f]=Full_K([NKPTS_1,NKPTS_2,NKPTS_3],1);
if ~L_path
    [VKPTS_f,NKPTS_f]=Full_K(N_k_nscf.*N_sc_base,1);
    
    SYM_OP_min(:,:,1) = eye(3);
    SYM_OP_min(:,:,2) = -eye(3); % minimal symmetry with inversion symmetry
    [VKPT,NKPTS,VKPTS_nosym_idx]...
        =Kpoints_gsym(VKPTS_f,A,SYM_OP_min);
    
    VKPT = VKPT(:,1:3)';
end
FID=fopen([ROOT_DIR,'/BEAD_1_primitive/INDEX_1_1_1']);
NGX=fread(FID,1,'int');NGY=fread(FID,1,'int');NGZ=fread(FID,1,'int');
fclose(FID);
NPLWV=NGX*NGY*NGZ;

G_INDEX_f=zeros(3,NPLWV);
if mod(NGZ,2)==0
    GB3=[0:NGZ/2,-NGZ/2+1:-1];
else
    GB3=[0:(NGZ-1)/2,-(NGZ+1)/2+1:-1];
end
if mod(NGY,2)==0
    GB2=[0:NGY/2,-NGY/2+1:-1];
else
    GB2=[0:(NGY-1)/2,-(NGY+1)/2+1:-1];
end
if mod(NGX,2)==0
    GB1=[0:NGX/2,-NGX/2+1:-1];
else
    GB1=[0:(NGX-1)/2,-(NGX+1)/2+1:-1];
end

G_INDEX_f(1,:)=repmat(GB1,1,NGY*NGZ);
G_INDEX_f(2,:)=reshape(repmat(GB2,NGX,NGZ),1,NPLWV);
G_INDEX_f(3,:)=reshape(repmat(GB3,NGX*NGY,1),1,NPLWV);

if ~L_path
    NPL_ss=zeros(NKPTS,1); G_INDEX_prim=cell(1,NKPTS);
    for NK=1:NKPTS
        G_INDEX_tmp = G_INDEX_f + VKPT(:,NK);
        
        G_f = zeros(4,NPLWV,'like',B);
        G_f(1:3,:)= B_prim.'*G_INDEX_tmp(1:3,:);
        G_f(4,:)=G_f(1,:).^2+G_f(2,:).^2+G_f(3,:).^2;
        
        idx = find(HSQDTM*G_f(4,:)<ENCUT);
        G_INDEX_prim{NK} = gather([G_INDEX_tmp(:,idx);idx]);
        NPL_ss(NK) = length(G_INDEX_prim{NK});
    end
end
%%
if ~LBloch
    FID=fopen([ROOT_DIR,'/BEAD_1/INDEX_1_1_1']);
    NGX_sc=fread(FID,1,'int');NGY_sc=fread(FID,1,'int');NGZ_sc=fread(FID,1,'int');
    fclose(FID);
    
    NPLWV_sc=NGX_sc*NGY_sc*NGZ_sc;
    
    G_INDEX_sc_f=zeros(3,NPLWV_sc);
    if mod(NGZ_sc,2)==0
        GB3_sc=[0:NGZ_sc/2,-NGZ_sc/2+1:-1];
    else
        GB3_sc=[0:(NGZ_sc-1)/2,-(NGZ_sc+1)/2+1:-1];
    end
    if mod(NGY_sc,2)==0
        GB2_sc=[0:NGY_sc/2,-NGY_sc/2+1:-1];
    else
        GB2_sc=[0:(NGY_sc-1)/2,-(NGY_sc+1)/2+1:-1];
    end
    if mod(NGX_sc,2)==0
        GB1_sc=[0:NGX_sc/2,-NGX_sc/2+1:-1];
    else
        GB1_sc=[0:(NGX_sc-1)/2,-(NGX_sc+1)/2+1:-1];
    end
    
    G_INDEX_sc_f(1,:)=repmat(GB1_sc,1,NGY_sc*NGZ_sc);
    G_INDEX_sc_f(2,:)=reshape(repmat(GB2_sc,NGX_sc,NGZ_sc),1,NPLWV_sc);
    G_INDEX_sc_f(3,:)=reshape(repmat(GB3_sc,NGX_sc*NGY_sc,1),1,NPLWV_sc);
    
    A_sc=A; B_sc=B;
end

for NK=1:NKPTS_sym
    NK_sym=NK;
    
    if LBloch
        VKPT_NK=(BASIS(NK_sym).VKPT_red.').*BASIS(NK).Factor_sym;
    end

    if LBloch
        NKPTS_red = BASIS(NK).NKPTS_red;
    else
        %% G_INDEX_sc
        NK_sc=1;
        G_INDEX_tmp = G_INDEX_sc_f + VKPT_sym(:,NK);
        
        G_f = zeros(4,NPLWV_sc,'like',B_sc);
        G_f(1:3,:)=B_sc.'*G_INDEX_tmp(1:3,:);
        G_f(4,:)=G_f(1,:).^2+G_f(2,:).^2+G_f(3,:).^2;
        
        idx = find(HSQDTM*G_f(4,:)<ENCUT);
        G_INDEX_sc = gather([G_INDEX_tmp(:,idx);idx]);
        NPL_sc = length(G_INDEX_sc);
        
        G=zeros(4,NPL_sc,'like',B_sc);
        %G_INDEX=gpuArray(G_INDEX);
        %G=zeros(4,NPL,'gpuArray');
        %G=zeros(4,NPL);
        G(1:3,:)=B_sc.'*G_INDEX_sc(1:3,:);
        G(4,:)=G(1,:).^2+G(2,:).^2+G(3,:).^2;
        
        clear G2
        G2(4,:)=sqrt(G(4,:));
        G2(4,:)=max(G2(4,:),1E-12);
        G2(1:3,:)=bsxfun(@rdivide,G(1:3,:),G2(4,:));
        
        G=single(G);
        %%
        eps = 1E-5;
        N_sc = diag(N_sc_base)*N_sc_ex;
        G_INDEX_sc_pr=mod(N_sc\G_INDEX_sc(1:3,:),1).';
        G_INDEX_sc_pr(G_INDEX_sc_pr>0.5+eps)=...
            G_INDEX_sc_pr(G_INDEX_sc_pr>0.5+eps)-1;
        % G_INDEX_sc_pr=mod((G_INDEX_sc(1:3,:).'./N_sc),1);
        VKPT_sc_pr = uniquetol(G_INDEX_sc_pr,'ByRows',true);
        VKPT_sc_pr(VKPT_sc_pr>0.5+eps)=...
            VKPT_sc_pr(VKPT_sc_pr>0.5+eps)-1;
        
        if L_path
            VKPT=VKPT_sc_pr.';NKPTS=size(VKPT,2);
            NPL_ss=zeros(NKPTS,1); G_INDEX_prim=cell(1,NKPTS);
            for NK_s=1:NKPTS
                G_INDEX_tmp = G_INDEX_f + VKPT(:,NK_s);
                
                G_f = zeros(4,NPLWV,'like',B);
                G_f(1:3,:)= B_prim.'*G_INDEX_tmp(1:3,:);
                G_f(4,:)=G_f(1,:).^2+G_f(2,:).^2+G_f(3,:).^2;
                
                idx = find(HSQDTM*G_f(4,:)<ENCUT);
                G_INDEX_prim{NK_s} = gather([G_INDEX_tmp(:,idx);idx]);
                NPL_ss(NK_s) = length(G_INDEX_prim{NK_s});
            end
            % The variables are meaningless when calculating band structure
            VKPT_sps_ia=1:NKPTS;
            VKPTS_idx=1:NKPTS;
            Factor_sym=ones(NKPTS,1);
            NKPTS_red=length(VKPT_sps_ia);
        else
            [~,VKPT_sp_ia]=ismembertol(...
                VKPT_sc_pr,VKPT.',1E-4,'ByRows',true);
            
            idx=find(VKPT_sp_ia==0);
            
            Factor_sym = ones(length(VKPT_sp_ia),1);
            Factor_sym(idx)=-1;
            
            VKPT_sps=VKPT_sc_pr;
            VKPT_sps(idx,:)=-VKPT_sc_pr(idx,:);
            VKPT_sps(VKPT_sps<-0.5+eps)=VKPT_sps(VKPT_sps<-0.5+eps)+1;
            [~,VKPT_sps_ia]=ismembertol(...
                VKPT_sps,VKPT.',1E-4,'ByRows',true);
            
            VKPTS_idx=VKPT_sps_ia;
            NKPTS_red=length(VKPT_sps_ia);
        end
        
        G_kp=cell(1,NKPTS_red);
        NPL_sum=0;
        for NK_s=1:NKPTS_red
            NK_c=int2str(VKPT_sps_ia(NK_s));
            NPL_s=NPL_ss(VKPT_sps_ia(NK_s));
            G_INDEX_sym=G_INDEX_prim{VKPT_sps_ia(NK_s)};
            
            G_INDEX_sym(1:3,:)=Factor_sym(NK_s)*G_INDEX_sym(1:3,:);
            [AA,G_kpidx] = ismembertol(...
                (N_sc\G_INDEX_sc(1:3,:)).',G_INDEX_sym(1:3,:).',1E-4,'ByRows',true);
            
            [~,idx]=sort(G_kpidx(G_kpidx~=0),'ascend');
            G_kpidx=find(AA);
            G_kpidx=G_kpidx(idx);
            
            G_kp{NK_s}=G_kpidx;
            
            NPL_sum=NPL_sum+NPL_s;
        end
        
        if NPL_sum~=NPL_sc
            disp('Plane wave number inconsistent')
            return
        end
    end

    if ~LBloch
        VKPT_red=VKPT(:,VKPTS_idx);
        
        BASIS(NK).G_INDEX_sc=G_INDEX_sc;
        BASIS(NK).G_kp=G_kp;
        BASIS(NK).NKPTS_red=NKPTS_red;
        BASIS(NK).VKPTS_idx=VKPTS_idx;
        BASIS(NK).VKPT_red=VKPT_red;
        BASIS(NK).Factor_sym=Factor_sym;
    end
end

if ~LBloch && NCL==1
    save('../BASIS_NDsym_struct.mat','BASIS')
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

  [~,~,~,...
      ~,~,~,~,~,...
      NELE,~,~,~]=...
      RD_HEAD(ROOT_DIR);
  NBANDS = ceil(NELE/2 * NBANDS_mul_Ef);

  %% Dispatch
  %CONF = floor(linspace(NSTART,NCONF,NGPU+1));
  CONF_CAC = NSTART+1:NSKIP:NCONF;
  NCONF_CAC = length(CONF_CAC);
  CONF_idx = floor(linspace(0,NCONF_CAC,NGPU+1));

  %load('BASIS_NDsym_nscf.mat','BASIS')
  %BASIS=BASIS(NK_sym);
  E_sample = EFERMI_Av-100:0.001:EFERMI_Av+70;
  BETA = kB*TEMP;
  ETA = [(0.0002:0.0002:0.001)*RYTOEV,(0.01:0.01:0.1)*RYTOEV,(2*(0:2:16)*pi+1)*BETA];
  %ETA = [(0.0002:0.0002:0.001)*RYTOEV,(0.01:0.01:0.1)*RYTOEV,(2*(4:2:16)*pi+1)*BETA];
  for jobid=1:NGPU
      %N_BEGIN = CONF(jobid);
      %N_END = CONF(jobid+1);
      N_BEGIN = CONF_CAC(CONF_idx(jobid)+1)-1;
      N_END = CONF_CAC(CONF_idx(jobid+1));
      
      workspace = struct('jobid', jobid,...
          'EFERMI_Av',EFERMI_Av,'NCL',NCL,...
          'NK_BEGIN',NK_BEGIN,'NK',NK,...
          'N_BEGIN',N_BEGIN,'N_END',N_END,'BASIS',BASIS(NK_sym),...
          'E_sample',E_sample,'ETA',ETA);

      jbatch{jobid}=batch(scriptG,'Workspace',workspace);
  end

  for jobid=1:NGPU
      wait(jbatch{jobid})
  end

  Spec_all = cell(BASIS(NK_sym).NKPTS_red,length(ETA));
  Spec_all2 = cell(BASIS(NK_sym).NKPTS_red,length(ETA));
  for NE=1:length(ETA)
      for NK_s=1:BASIS(NK_sym).NKPTS_red
          Spec_all{NK_s,NE}=zeros(size(E_sample),'like',gather(B));
          Spec_all2{NK_s,NE}=zeros(size(E_sample),'like',gather(B));
      end
  end
  for jobid=1:NGPU
      fn=['Spec_',int2str(NK),'_',int2str(jobid),'.mat'];
      load(fn)
      for NE=1:length(ETA)
          for NK_s=1:BASIS(NK_sym).NKPTS_red
              Spec_all{NK_s,NE}=Spec_all{NK_s,NE}+Spec_av{NK_s,NE};
              Spec_all2{NK_s,NE}=Spec_all2{NK_s,NE}+Spec_av2{NK_s,NE};
          end
      end
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

  for NE=1:length(ETA)
      for NK_s=1:BASIS(NK_sym).NKPTS_red
        Spec_all{NK_s,NE}=Spec_all{NK_s,NE}/NBEAD/N_CAC;
        Spec_all2{NK_s,NE}=Spec_all2{NK_s,NE}/NBEAD/N_CAC;
      end
  end
  toc;

  save(['../Spectral_',int2str(NK),'.mat'],...
      'Spec_all','Spec_all2','EFERMI_Av','ETA','E_sample','-v7.3')

  for jobid=1:NGPU
    delete(['../Spec_',int2str(NK),'_',int2str(jobid),'.mat'])
  end

  CONF_CAC = NSTART+1:NSKIP:NCONF;
  NCONF_CAC = length(CONF_CAC);
  CONF_idx = floor(linspace(0,NCONF_CAC,NGPU+1));
  
  E_all=zeros(NBANDS,NBEAD,NCONF_CAC);
  NCAC_s=0;
  for jobid=1:NGPU
      N_BEGIN = CONF_CAC(CONF_idx(jobid)+1)-1;
      N_END = CONF_CAC(CONF_idx(jobid+1));
      NCAC=length(N_BEGIN+1:NSKIP:N_END);
      
      load(['../Eigen_',int2str(NK),'_',int2str(jobid)],'E_t');
      E_all(:,:,NCAC_s+1:NCAC_s+NCAC)=E_t;
      NCAC_s=NCAC_s+NCAC;
  end
  E_t=E_all;

  save(['../Eigen_',int2str(NK),'.mat'],'E_t','NCAC_s','-v7.3')

  for jobid=1:NGPU
    delete(['../Eigen_',int2str(NK),'_',int2str(jobid),'.mat'])
  end

end
