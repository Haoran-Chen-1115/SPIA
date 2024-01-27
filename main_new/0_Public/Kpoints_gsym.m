%% Now determine the irreducible k-mesh/q-mesh under given symmetry
function ...
    [VKPTS_SYM,NKPTS,VKPTS_SYM_idx]...
    =Kpoints_gsym(VKPTS_full,A,SYM_OP_TRUE)
%% Symmerty operations
eps = 1E-5;
N_SYM_TRUE = size(SYM_OP_TRUE,3);
% Take inverse ones
SYM_OP_INV = SYM_INV(SYM_OP_TRUE,A);

%% Generate the full k-mesh
NKPTS_full = size(VKPTS_full,2);

%% Find irreducible k-points
% NKPTS = 1; K_weight = 1;
% VKPTS_SYM = zeros(3,1); % Gamma point must be a irreducible k point
% VKPTS_SYM_idx = zeros(1,NKPTS_full); 
% VKPTS_SYM_idx(1) = 1; % Index for irreducible k points.
% for NK = 2:NKPTS_full
NKPTS = 0; K_weight = 0;
VKPTS_SYM = [];zeros(3,1); % Gamma point must be a irreducible k point
VKPTS_SYM_idx = zeros(1,NKPTS_full); 
%VKPTS_SYM_idx(1) = 1; % Index for irreducible k points.
for NK = 1:NKPTS_full
    VKPT_SYM = zeros(3,N_SYM_TRUE);
    for NS = 1:N_SYM_TRUE
        VKPT_SYM(:,NS) = SYM_OP_INV(:,:,NS)' * VKPTS_full(:,NK);
    end
    VKPT_SYM = mod(VKPT_SYM,1);
    VKPT_SYM(VKPT_SYM>0.5+eps) = VKPT_SYM(VKPT_SYM>0.5+eps) - 1;
    if NK~=1
        if any(ismembertol(VKPT_SYM',VKPTS_SYM','ByRows',true))
            idx=find(ismembertol(VKPTS_SYM',VKPT_SYM','ByRows',true));
            K_weight(idx) = K_weight(idx) + 1;
            VKPTS_SYM_idx(NK) = idx;
        else
            NKPTS = NKPTS + 1;
            VKPTS_SYM = [VKPTS_SYM, VKPTS_full(:,NK)];
            K_weight(NKPTS) = 1;
            VKPTS_SYM_idx(NK) = NKPTS;
        end
    else
        NKPTS = NKPTS + 1;
        VKPTS_SYM = [VKPTS_SYM, VKPTS_full(:,NK)];
        K_weight(NKPTS) = 1;
        VKPTS_SYM_idx(NK) = NKPTS;
    end
end

VKPTS_SYM = [VKPTS_SYM',K_weight'];
end