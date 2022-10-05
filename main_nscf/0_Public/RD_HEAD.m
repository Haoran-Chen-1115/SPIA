function [A,B,OMEGA,...
    NKPTS,ISPIN,NTYP,NITYP,NCPU,...
    NELE,NBANDS,VKPT]=RD_HEAD(ROOT_DIR,NELE_TYP,NBANDS_mul)

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

% electrons and bands
NELE = sum(NITYP.*NELE_TYP);
NBANDS = ceil(NELE/2 * NBANDS_mul);

% k-points
FID=fopen([ROOT_DIR,'/BEAD_1/VKPT']);
VKPT=fread(FID,[3,NKPTS],'double');
VKPT(abs(VKPT)<1E-15)=0;
fclose(FID);

end
