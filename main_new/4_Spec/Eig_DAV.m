function [WAVE_conf,E_conf,N_notconv]=...
    Eig_DAV(HAM,OVL,WAVE_conf,B,NBANDS,NELE,NPL,RYTOEV)
%% Davidson iteration for determing eigenfunctions
DAV_DIM = 2; % maximum number of bases
NELM = 20; % maximum outer iteration loop

NBANDS_max = NBANDS*DAV_DIM;

% In Ry.
conv_thr=1E-8;
ethr=max(1E-4,0.1*min(1E-2,conv_thr/NELE));
%empty_ethr= max(ethr*5,1E-5);
empty_ethr=ethr; % same accuracy for all bands

%% Diagonal part of the Hamiltonian, used for preconditioning
h_diag=real(diag(HAM).');
s_diag=real(diag(OVL).');

%% Iteration
% Use the original Davidson algorithm instead of block one
ethr_all = [ethr*ones(ceil(NELE/2)+10,1);...
    empty_ethr*ones(NBANDS-ceil(NELE/2)-10,1)]*RYTOEV;

avg_iter=0;
for j=1:5
    %%
    WAVE_in = gpuArray(zeros(NPL,NBANDS_max,'like',WAVE_conf));
    Hpsi = gpuArray(zeros(NPL,NBANDS_max,'like',WAVE_conf));
    Spsi = gpuArray(zeros(NPL,NBANDS_max,'like',WAVE_conf));
    WAVE_in(:,1:NBANDS) = WAVE_conf;
    NBANDS_cac=NBANDS;
    
    %% Hamiltonian in the subspace
    %Hpsi=H_on_psi();
    Hpsi(:,1:NBANDS_cac)=...
        HAM*WAVE_in(:,1:NBANDS_cac);
    Spsi(:,1:NBANDS_cac)=...
        OVL*WAVE_in(:,1:NBANDS_cac);
    %H_on_psi(WAVE_in,SV,CDIJ,GPROJ,...
    %    G,G_INDEX,NPLWV,NGX,NGY,NGZ,NPL,...
    %    NTYP,NITYP,LMMAXC,NBANDS_cac);
    
    FHAM=WAVE_in(:,1:NBANDS_cac)'*Hpsi(:,1:NBANDS_cac);
    FOVL=WAVE_in(:,1:NBANDS_cac)'*Spsi(:,1:NBANDS_cac);
    if j==1
        FHAM = (FHAM + FHAM')/2;
        FOVL = (FOVL + FOVL')/2;
        F=chol(FOVL,'lower');
        FHAM = F\FHAM/F';
        FHAM=(FHAM+FHAM')/2;
        [basis_in,E_tmp]=eig(FHAM);
        basis_in=F'\basis_in;
        E_old=diag(real(E_tmp));
    else
        basis_in=eye(NBANDS_cac,'like',WAVE_in);
        E_old=real(diag(FHAM));
    end
    
    E_in=E_old;
    
    % Optimize the subspaces
    Conv=false(1,NBANDS);
    max_daviter=40;
    for dav_iter=1:max_daviter
        t_in=tic;
        basis_in=[basis_in(:,~Conv),...
            basis_in(:,NBANDS+1:end),...
            basis_in(:,Conv)];
        E_in = [E_in(~Conv);...
            E_in(NBANDS+1:end);...
            E_in(Conv)];
        N_notconv = length(find(~Conv));
        E_notconv = E_in(1:N_notconv);
        
        % Residual vector
        Hpsi(:,1:NBANDS_cac)=...
            Hpsi(:,1:NBANDS_cac)*basis_in;
        WAVE_in(:,1:NBANDS_cac) =...
            WAVE_in(:,1:NBANDS_cac)*basis_in;
        Spsi(:,1:NBANDS_cac) =...
            Spsi(:,1:NBANDS_cac)*basis_in;
        
        % RES,residual vector
        WAVE_in(:,NBANDS_cac+1:NBANDS_cac+ N_notconv)...
            = Hpsi(:,1:N_notconv)...
            -Spsi(:,1:N_notconv)*diag(E_notconv);
        
        % Preconditioning, version in QE
        x = h_diag - E_notconv.*s_diag;
        x = 0.5*(x+1+sqrt(1+(x-1).^2));
        WAVE_in(:,NBANDS_cac+1:NBANDS_cac+ N_notconv)...
            = WAVE_in(:,NBANDS_cac+1:NBANDS_cac+ N_notconv)...
            ./x.';
        
        RNORM = sqrt(sum(...
            abs(WAVE_in(:,NBANDS_cac+1:NBANDS_cac+ N_notconv)).^2 ...
            ,1));
        WAVE_in(:,NBANDS_cac+1:NBANDS_cac+ N_notconv)...
            = WAVE_in(:,NBANDS_cac+1:NBANDS_cac+ N_notconv)./RNORM;
        %if ~isempty(G0_idx)
        %    RES(G0_idx)=real(RES(G0_idx));
        %end
        
        % Hamiltonian acts on the new bases
        % Hres
        Hpsi(:,NBANDS_cac+1:NBANDS_cac+ N_notconv)...
            = HAM * WAVE_in(:,NBANDS_cac+1:NBANDS_cac+ N_notconv);
        Spsi(:,NBANDS_cac+1:NBANDS_cac+ N_notconv)...
            = OVL * WAVE_in(:,NBANDS_cac+1:NBANDS_cac+ N_notconv);
        
        % Hamiltonian in the new subspace
        %WAVE_in(:,NBANDS_cac+1:NBANDS_cac+ N_notconv) = RES;
        %Hpsi(:,NBANDS_cac+1:NBANDS_cac+ N_notconv) = Hres;
        NBANDS_cac = NBANDS_cac + N_notconv;
        
        FHAM=WAVE_in(:,1:NBANDS_cac)'*Hpsi(:,1:NBANDS_cac);
        FOVL=WAVE_in(:,1:NBANDS_cac)'*Spsi(:,1:NBANDS_cac);
        
        % Diagonalization
        FHAM = (FHAM + FHAM')/2;
        FOVL = (FOVL + FOVL')/2;
        
        %[F,E]=eig(FOVL);
        %for i=1:length(F)
        %    F(:,i)=F(:,i)/sqrt(E(i,i));
        %end
        %FHAM = F'*FHAM*F;
        F=chol(FOVL,'lower');
        FHAM = F\FHAM/F';
        FHAM=(FHAM+FHAM')/2;
        [basis_in,E_tmp]=eig(FHAM);
        basis_in=F'\basis_in;
        E_tmp=diag(real(E_tmp));
        
        %WAVE_trial = WAVE_trial*basis_in;
        E_in = E_tmp;
        DE=abs(real(E_in(1:NBANDS)-E_old));
        Conv=DE<ethr_all;
        
        N_notconv = length(find(~Conv));
        
        E_old=E_in(1:NBANDS);
        
        if N_notconv==0 || dav_iter==max_daviter
            WAVE_conf=gather(...
                WAVE_in(:,1:NBANDS_cac)...
                *basis_in(:,1:NBANDS));
            
            E_conf=E_in(1:NBANDS);
            break
        end
        
        if NBANDS_cac+N_notconv>NBANDS_max
            WAVE_in(:,1:NBANDS)=...
                WAVE_in(:,1:NBANDS_cac)...
                *basis_in(:,1:NBANDS);
            WAVE_in(:,NBANDS+1:end)=0;
            
            Hpsi(:,1:NBANDS)=...
                Hpsi(:,1:NBANDS_cac)*basis_in(:,1:NBANDS);
            Hpsi(:,NBANDS+1:end)=0;
            
            Spsi(:,1:NBANDS)=...
                Spsi(:,1:NBANDS_cac)*basis_in(:,1:NBANDS);
            Spsi(:,NBANDS+1:end)=0;
            
            NBANDS_cac=NBANDS;
            basis_in=eye(NBANDS_cac,'like',WAVE_in);
            
            E_in=E_in(1:NBANDS);
        end
        t_in_end=toc(t_in);
        %fprintf(FID_OUT,'step %d: %6.4f s   ',[avg_iter+dav_iter,t_in_end]);
        %if mod(avg_iter+dav_iter,5)==0
        %   fprintf(FID_OUT,'\n');
        %end
    end
    avg_iter=avg_iter+dav_iter;
    if N_notconv==0
        break
    end
end
%fprintf(FID_OUT,'\n');
clear Hpsi basis_in FHAM FOVL

end
