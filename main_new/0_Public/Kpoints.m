%% Now determine the irreducible k-mesh/q-mesh
function ...
    [VKPTS_SYM,NKPTS,VKPTS_full,VKPTS_SYM_idx,SYM_OP_TRUE]...
    =Kpoints(N_KP,A,NTYP,NITYP,POSION,IFLAG)
%% Symmerty operations
eps = 1E-5;
% [IBRAV,N_SYM_ALL,SYM_OP_ALL,TRANS_ROT]...
%     = Symmetry(A,POSION,NTYP,NITYP);
% 
% idx_op = ismembertol(TRANS_ROT.',[0,0,0],'ByRows',true);
% SYM_OP_TRUE = SYM_OP_ALL(:,:,idx_op);
% N_SYM_TRUE = length(find(idx_op));
[IBRAV,N_SYM_TRUE,SYM_OP_TRUE,TRANS_ROT]...
    = Symmetry(A,POSION,NTYP,NITYP);
% Take inverse ones
SYM_OP_INV = SYM_INV(SYM_OP_TRUE,A);

%POS_fake=cell(1);
%POS_fake{1} = [0;0;0]; NTYP_fake = 1; NITYP_fake=1;
%B = inv(A).';
%A_fake = B.*N_KP;
%[IBRAV_k,N_SYM_TRUE,SYM_OP_INV,TRANS_ROT_INV]...
%    = Symmetry(A_fake,POS_fake,NTYP_fake,NITYP_fake);

%% Generate the full k-mesh
NKPTS_full = N_KP(1)*N_KP(2)*N_KP(3);


K1 = (0:N_KP(1)-1)/N_KP(1);
K2 = (0:N_KP(2)-1)/N_KP(2);
K3 = (0:N_KP(3)-1)/N_KP(3);

clear VKPTS_full
if IFLAG == 1 % k-mesh (consistent with VASP)
    VKPTS_full(1,:) = reshape(repmat(K1,1,N_KP(2)*N_KP(3)),1,[]);
    VKPTS_full(2,:) = reshape(repmat(K2,N_KP(1),N_KP(3)),1,[]);
    VKPTS_full(3,:) = reshape(repmat(K3,N_KP(1)*N_KP(2),1),1,[]);
    VKPTS_full(VKPTS_full>0.5+eps) = VKPTS_full(VKPTS_full>0.5+eps) - 1;
elseif IFLAG == 2 % q-mesh (consistent with QE)
    VKPTS_full(1,:) = reshape(repmat(K1,N_KP(2)*N_KP(3),1),1,[]);
    VKPTS_full(2,:) = reshape(repmat(K2,N_KP(3),N_KP(1)),1,[]);
    VKPTS_full(3,:) = reshape(repmat(K3,1,N_KP(1)*N_KP(2)),1,[]);
    VKPTS_full(VKPTS_full>0.5+eps) = VKPTS_full(VKPTS_full>0.5+eps) - 1;
end
%% Find irreducible k-points
NKPTS = 1; K_weight = 1;
VKPTS_SYM = zeros(3,1); % Gamma point must be a irreducible k point
VKPTS_SYM_idx = zeros(1,NKPTS_full); 
VKPTS_SYM_idx(1) = 1; % Index for irreducible k points.
for NK = 2:NKPTS_full
    VKPT_SYM = zeros(3,N_SYM_TRUE);
    for NS = 1:N_SYM_TRUE
        VKPT_SYM(:,NS) = SYM_OP_INV(:,:,NS)' * VKPTS_full(:,NK);
    end
    VKPT_SYM = mod(VKPT_SYM,1);
    VKPT_SYM(VKPT_SYM>0.5+eps) = VKPT_SYM(VKPT_SYM>0.5+eps) - 1;
    if any(ismembertol(VKPT_SYM',VKPTS_SYM','ByRows',true))
        idx=find(ismembertol(VKPTS_SYM',VKPT_SYM','ByRows',true));
        K_weight(idx) = K_weight(idx) + 1;
        VKPTS_SYM_idx(NK) = idx;
    else
        NKPTS = NKPTS + 1;
        VKPTS_SYM = [VKPTS_SYM, VKPTS_full(:,NK)];
        % VKPTS_SYM(:,NKPTS) = VKPTS_full(:,NK);
        K_weight(NKPTS) = 1;
        VKPTS_SYM_idx(NK) = NKPTS;
    end
end

VKPTS_SYM = [VKPTS_SYM',K_weight'];
end