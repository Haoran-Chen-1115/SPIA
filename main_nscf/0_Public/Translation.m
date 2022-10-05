function [Ts_true,Ns_true]=Translation(ROOT_DIR,N_sc)

FID=fopen([ROOT_DIR,'/BEAD_1_primitive/HEAD_1']);
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

N_trans = det(N_sc);
% In case that there appears supercell types like 1*1*N_trans
T0 = 0:N_trans-1;
clear Ts
Ts(1,:) = reshape(repmat(T0,1,N_trans*N_trans),1,[]);
Ts(2,:) = reshape(repmat(T0,N_trans,N_trans),1,[]);
Ts(3,:) = reshape(repmat(T0,N_trans*N_trans,1),1,[]);
Ts = Ts.';

% There must be N_trans different translations
% they together give a set of configurations 
% that guarantee the periodic boundary condition of the primitive cell is preserved
Ts_true = zeros(1,3);
Ns_true = 1;
Nsc_inv = inv(N_sc);
eps = 1E-4;
for Ns = 2:N_trans^3
    Ts_tmp =  mod(Ts(Ns,:)*Nsc_inv,1);
    Ts_tmp(abs(Ts_tmp-1)<1E-4)=Ts_tmp(abs(Ts_tmp-1)<1E-4)-1;
    
    if ~ismembertol(Ts_tmp,Ts_true,'ByRows',true)
        Ns_true = Ns_true+1;
        Ts_true=[Ts_true;Ts_tmp];
    end
end

if Ns_true~=N_trans
    disp('Translation number incorrect')
end