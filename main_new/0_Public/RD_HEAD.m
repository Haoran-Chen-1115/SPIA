function [A,B,OMEGA,...
    NKPTS,ISPIN,NTYP,NITYP,NCPU,...
    NELE,VKPT,NELE_TYP,POMASS]=...
    RD_HEAD(ROOT_DIR)

% Lattice
FID=fopen([ROOT_DIR,'/BEAD_1/HEAD_1']);
A=fread(FID,[3,3],'double');
fread(FID,1,'double');
fread(FID,1,'double');
NKPTS=fread(FID,1,'int');ISPIN=fread(FID,1,'int');
NTYP=fread(FID,1,'int');NITYP=zeros(1,NTYP);
for NT=1:NTYP
    NITYP(NT)=fread(FID,1,'int');
end
NCPU=fread(FID,1,'int');
fclose(FID);
if isempty(NCPU)
    NCPU=1;
end
B=2*pi*inv(A);
OMEGA=det(A);

% k-points
FID=fopen([ROOT_DIR,'/BEAD_1/VKPT']);
VKPT=fread(FID,[3,NKPTS],'double');
VKPT(abs(VKPT)<1E-15)=0;
fclose(FID);

FID=fopen([ROOT_DIR,'/BEAD_1/ELECT_1']);
NELE_TYP=zeros(1,NTYP);
POMASS=zeros(1,NTYP);
for NT=1:NTYP
    NELE_TYP(NT)=fread(FID,1,'double');
end
for NT=1:NTYP
    POMASS(NT)=fread(FID,1,'double');
end

% electrons and bands
NELE = sum(NITYP.*NELE_TYP);
end
