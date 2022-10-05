function [G_INDEX,G,G2,NPL,NPLWV,NGX,NGY,NGZ]=RD_INDEX(ROOT_DIR,NC_c,NK_c,VKPT,B,ENCUT)
%% READ IN BASIC INFORMATIONS
FID=fopen([ROOT_DIR,'/BEAD_',NC_c,'/INDEX_1_',NK_c,'_1']);
NGX=fread(FID,1,'int');NGY=fread(FID,1,'int');NGZ=fread(FID,1,'int');
% NPL_o=fread(FID,1,'int');
% G_INDEX_o=fread(FID,[4,NPL_o],'int')+[VKPT;0];
fclose(FID);

NPLWV=NGX*NGY*NGZ;
%% FULL G MESH, AND APPLY ENERGY CUTOFF 
RYTOEV=13.605826;AUTOA=0.529177249;HSQDTM=RYTOEV*AUTOA*AUTOA;

G_INDEX_f=zeros(3,NPLWV);
GB3=[0:NGZ/2,-NGZ/2+1:-1];
GB2=[0:NGY/2,-NGY/2+1:-1];
GB1=[0:NGX/2,-NGX/2+1:-1];

G_INDEX_f(1,:)=repmat(GB1,1,NGY*NGZ);
G_INDEX_f(2,:)=reshape(repmat(GB2,NGX,NGZ),1,NPLWV);
G_INDEX_f(3,:)=reshape(repmat(GB3,NGX*NGY,1),1,NPLWV);
G_INDEX_f = G_INDEX_f + VKPT;

G_f(1:3,:)=B.'*G_INDEX_f(1:3,:);
G_f(4,:)=G_f(1,:).^2+G_f(2,:).^2+G_f(3,:).^2;

idx = find(HSQDTM*G_f(4,:)<ENCUT);
G_INDEX = [G_INDEX_f(:,idx);idx];
NPL = length(G_INDEX);

%%
G=zeros(4,NPL,'like',B);
%G_INDEX=gpuArray(G_INDEX);
%G=zeros(4,NPL,'gpuArray');
%G=zeros(4,NPL);
G(1:3,:)=B.'*G_INDEX(1:3,:);
G(4,:)=G(1,:).^2+G(2,:).^2+G(3,:).^2;

G2(4,:)=sqrt(G(4,:));
G2(4,:)=max(G2(4,:),1E-12);
G2(1:3,:)=bsxfun(@rdivide,G(1:3,:),G2(4,:));

G=single(G);

end
