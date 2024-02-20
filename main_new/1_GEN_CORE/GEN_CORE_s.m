%% Calulate T_0^\dagger T
%% Initialization
addpath('../')
addpath('../0_Public')
%parameter;
parpool(10,'IdleTimeout', 2880)

%% Pseudo basis we care about
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

% Read in information about grid and the index for calculated G vectors
% It should be the same for all beads
[G_INDEX,G,G2,NPL,NPLWV,NGX,NGY,NGZ]=RD_INDEX(ROOT_DIR,'1','1',[0;0;0],B,ENCUT)

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

