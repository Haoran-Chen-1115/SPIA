%% Calulate T_0^\dagger T
%% Initialization
addpath('../')
addpath('../0_Public')
%parameter;
parpool(10,'IdleTimeout', 2880)

%% Pseudo basis we care about
% Read in information about grid and the index for calculated G vectors
% It should be the same for all beads
FID=fopen('../../BEAD_1/INDEX_1_1_1');
NGX=fread(FID,1,'int');NGY=fread(FID,1,'int');NGZ=fread(FID,1,'int');
NPL=fread(FID,1,'int');
G_INDEX=fread(FID,[4,NPL],'int');
fclose(FID);

NPLWV=NGX*NGY*NGZ;
% GREEN=complex(zeros(NPL,NPL,'single','gpuArray'));
% GREEN=zeros(NPL,NPL,'gpuArray')+1i*zeros(NPL,NPL,'gpuArray');

% Everything here except EFERMI and TOTEN should be the same for all beads
FID=fopen('../../BEAD_1/HEAD_1');
A=fread(FID,[3,3],'double');
EFERMI=fread(FID,1,'double');
TOTEN=fread(FID,1,'double');
NKPTS=fread(FID,1,'int');ISPIN=fread(FID,1,'int');
NTYP=fread(FID,1,'int');NITYP=zeros(1,NTYP);
for NT=1:NTYP
    NITYP(NT)=fread(FID,1,'int');
end
fclose(FID);
B=2*pi*inv(A);
OMEGA=det(A);

% Force NKPTS=1
NKPTS=1;

% B=gpuArray(B);G_INDEX=gpuArray(G_INDEX);
% G=zeros(4,NPL,'gpuArray');
B=gpuArray(B);G_INDEX=gpuArray(G_INDEX);
G=zeros(4,NPL,'gpuArray');
G(1:3,:)=B.'*G_INDEX(1:3,:);
G(4,:)=G(1,:).^2+G(2,:).^2+G(3,:).^2;

clear G2
G2(4,:)=sqrt(G(4,:));
G2(4,:)=max(G2(4,:),1E-12);
G2(1:3,:)=bsxfun(@rdivide,G(1:3,:),G2(4,:));

G=single(G);

%% Calculate T_core
% Information obtained from POTCAR, which is the same for all configuartions
tic;
[NMAX,RG,SIMPI,LLMAX,WAE_PS,WAE,WPS,...
    LMAX,LMMAXC,CH0,CH1,CH2,CH3,LPS]=RD_PAW('../../',NTYP);

toc;

%% Calculate T^\daggerT
%% Method 1 in k-space
% GEN_TRANS_CORE_recp;

%% Method 2 in real space
GEN_CORE_real_s;

