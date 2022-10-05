function [BASIS_trial,E_tmp,DESUM]=Eigen_DAV_fun(G,FHAM,FOVL,BASIS_k,...
    NBANDS,NSIM,HSQDTM,NELE,EDIFF)
%% Kinetic energy
TKIN = diag(HSQDTM*G(4,:));

%% Generate the subspace
% Use the equilibrium wavefunction to be the trial function
BASIS_trial_in=BASIS_k;
NPL=size(BASIS_k,1);
NBLOCK=ceil(NBANDS/NSIM);
if mod(NBANDS,NSIM)~=0
  NSIM_inner=[repmat(NSIM,1,NBLOCK-1),mod(NBANDS,NSIM)];
else
  NSIM_inner=repmat(NSIM,1,NBLOCK);
end


% NBANDS=ceil(NBANDS/NSIM)*NSIM;
EBREAK=EDIFF/NBANDS/NSIM;

% Orthonormalization
FHAM_basis = BASIS_trial_in' * FHAM * BASIS_trial_in;
FHAM_basis = diag(real(diag(FHAM_basis)))+...
    triu(FHAM_basis,1)+triu(FHAM_basis,1)';

FOVL_basis = BASIS_trial_in' * FOVL * BASIS_trial_in;
FOVL_basis = diag(real(diag(FOVL_basis)))+...
    triu(FOVL_basis,1)+triu(FOVL_basis,1)';

[basis_ini,E_tmp]=eig(gather(FHAM_basis),gather(FOVL_basis));
% basis_ini=basis_ini...
%                 ./(sqrt(diag(real(basis_ini'*FOVL_basis*basis_ini))).');
BASIS_trial_in = BASIS_trial_in * basis_ini;
clear FHAM_basis FOVL_basis basis_ini

% E_pre=diag(BASIS_trial_in'*FHAM*BASIS_trial_in);
% [~,idx]=sort(real(E_pre),'ascend');
% BASIS_trial_in=BASIS_trial_in(:,idx);

err=1E9;
NITER = 2; % maximum ineer iteration loop
NELM = 20; % maximum outer iteration loop
E=zeros(NBANDS,1);
% while err > EBREAK

%%
ESUM=1E9;
for i=1:NELM
    %% Handle with blocks
    %tic;
    FOVL_pre_l = BASIS_trial_in' * FOVL;
    
    NBL=1;
    for j=1:NBLOCK
        basis_sub=zeros(NPL,NSIM_inner(j)*NITER,'gpuArray');
        basis_sub(:,1:NSIM_inner(j))=BASIS_trial_in(:,NBL:NBL+NSIM_inner(j)-1);
        E_BAND_sub = 1E9;
        for ITER=1:NITER
%             tic;
            basis_sub_iter=basis_sub(:,1:NSIM_inner(j)*ITER);
            E_BAND_sub_last = E_BAND_sub;
            FOVL_sub_r=FOVL*basis_sub_iter;
            FHAM_sub_r=FHAM*basis_sub_iter;
            
            FOVL_sub=basis_sub_iter'*FOVL_sub_r;
            FHAM_sub=basis_sub_iter'*FHAM_sub_r;
%             toc;
            % Preconditioning
            RES_VEC=basis_sub_iter(:,(ITER-1)*NSIM_inner(j)+1:ITER*NSIM_inner(j));
%             EKIN = 1.5*sum(TKIN*abs(RES_VEC).^2,1);
%             EKIN(EKIN<3)=3;
%             X = diag(TKIN)./EKIN;
%             X2 = 27+X.*(18+X.*(12+8.*X));
%             PRECON = X2./(X2+16*X.^4);
%             PRECON = PRECON./EKIN*2;
%             toc;
            % Subspace
            if ~isempty(find(~isfinite(FOVL_sub)))
                DESUM=1E9;
                return
            end
            FHAM_sub=diag(real(diag(FHAM_sub)))+...
                triu(FHAM_sub,1)+triu(FHAM_sub,1)';
            FOVL_sub=diag(real(diag(FOVL_sub)))+...
                triu(FOVL_sub,1)+triu(FOVL_sub,1)';
            
            [WAVE_sub,E_BAND_sub]=eig(gather(FHAM_sub),gather(FOVL_sub));
            E_BAND_sub=diag(E_BAND_sub);
            E_BAND_sub=E_BAND_sub(1:NSIM_inner(j));
            WAVE_sub=WAVE_sub(:,1:NSIM_inner(j));
%             WAVE_sub=WAVE_sub...
%                 ./(sqrt(diag(real(WAVE_sub'*gather(FOVL_sub)*WAVE_sub))).');
%             toc;
            DEMAX = max(abs(E_BAND_sub-E_BAND_sub_last));
            clear basis_sub_iter
            
            if (abs(DEMAX)<EBREAK || ITER==NITER)
                break
            end
%             toc;
            % Residual vector
            RES_VEC=FHAM_sub_r*WAVE_sub-FOVL_sub_r*WAVE_sub*diag(E_BAND_sub);
%             toc;
            % Preconditioning
            EKIN = 1.5*sum(TKIN*abs(RES_VEC).^2,1);
            EKIN(EKIN<3)=3;
            X = diag(TKIN)./EKIN;
            X2 = 27+X.*(18+X.*(12+8.*X));
            PRECON = X2./(X2+16*X.^4);
            PRECON = PRECON./EKIN*2;
             
            
            % Apply Preconditioning
            basis_inter=PRECON.*RES_VEC;
            basis_inter=basis_inter./(sqrt(diag(real(basis_inter'*FOVL*basis_inter))).');
%             toc;
            % Orthonormal to all other bands
            FOVL_pre = FOVL_pre_l * basis_inter;
            basis_inter = basis_inter - BASIS_trial_in*FOVL_pre;
            
            basis_sub(:,ITER*NSIM_inner(j)+1:(ITER+1)*NSIM_inner(j)) = basis_inter;
%             toc;
            % Slower than directly calculate eigens with CPU
            % T_chol=chol(FOVL_sub,'lower');
            % FHAM_sub2=T_chol\FHAM_sub/(T_chol');
            % [WAVE_sub2,E_BAND_sub2]=eig(FHAM_sub2);
            % WAVE_sub2=(T_chol')\WAVE_sub2;
        end
        BASIS_trial_in(:,NBL:NBL+NSIM_inner(j)-1)=...
            basis_sub(:,1:NSIM_inner(j)*ITER)*WAVE_sub(:,1:NSIM_inner(j));
        FOVL_pre_l(NBL:NBL+NSIM_inner(j)-1,:)=BASIS_trial_in(:,NBL:NBL+NSIM_inner(j)-1)'*FOVL;
        E(NBL:NBL+NSIM_inner(j)-1)=E_BAND_sub;
        clear basis_sub WAVE_sub
        
%         CHAM(:,NBL:NBL+NSIM_inner(j)-1)=BASIS_trial_in'...
%             *(FHAM*BASIS_trial_in(:,NBL:NBL+NSIM_inner(j)-1)...
%             -FOVL*BASIS_trial_in(:,NBL:NBL+NSIM_inner(j)-1)*diag(E_BAND_sub));
        
        NBL=NBL+NSIM_inner(j);
    end
%     [E,idx]=sort(E,'ascend');
%     basis_k=basis_k(:,idx);
    %toc;
    ESUM_last = ESUM;
    ESUM = 2*sum(E(1:NELE/2));
    DESUM = ESUM-ESUM_last;
    
    if abs(DESUM)<EDIFF
        BASIS_trial = BASIS_trial_in;
        E_tmp=single(E);
        E_tmp(1)
        DESUM
        break
    end
    
    %% Subspace rotation
    CHAM = BASIS_trial_in'*(FHAM*BASIS_trial_in-FOVL*BASIS_trial_in*diag(E));
    NBL=1;
    for j=1:NBLOCK
        CHAM(NBL:NBL+NSIM_inner(j)-1,NBL:NBL+NSIM_inner(j)-1)=diag(E(NBL:NBL+NSIM_inner(j)-1));
        NBL=NBL+NSIM_inner(j);
    end
    CHAM = triu(CHAM,1)+triu(CHAM,1)'+diag(real(diag(CHAM)));
    [CHAM,E_tmp] = eig(CHAM);
    E_tmp=diag(E_tmp);
    BASIS_trial_in=BASIS_trial_in*CHAM;
    clear CHAM
    %toc;
    
    if i==NELM
        E_tmp(1)
        DESUM
    end
    
    BASIS_trial = BASIS_trial_in;
end


%end
