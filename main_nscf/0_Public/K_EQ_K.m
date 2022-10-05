function idx = K_EQ_K(VKPT,VKPTS,SYM_OP_INV)
% The function is to find the equivalent k-point of VKPT in the set VKPTS
eps = 1E-5;
NQPTS = size(VKPT,2);
N_SYM_TRUE = size(SYM_OP_INV,3);
% for NQ = 1:NQPTS
NQ=1;
    VKPT_SYM = zeros(3,N_SYM_TRUE);
    for NS = 1:N_SYM_TRUE
        VKPT_SYM(:,NS) = SYM_OP_INV(:,:,NS)' * VKPT(:,NQ);
    end
    VKPT_SYM = mod(VKPT_SYM,1);
    VKPT_SYM(VKPT_SYM>0.5+eps) = VKPT_SYM(VKPT_SYM>0.5+eps) - 1;
    VKPT_SYM(VKPT_SYM<-0.5+eps) = VKPT_SYM(VKPT_SYM<-0.5+eps) + 1;
    [idx_tmp,idx_local]=ismembertol(VKPT_SYM',VKPTS',eps,'ByRows',true);
    idx_tmp = find(idx_tmp);
    SYM_OP = SYM_OP_INV(:,:,idx_tmp(1));
%     if length(idx_tmp)==1
%         SYM_OP = SYM_OP_INV(:,:,idx_local(idx_tmp));
%         break
%     end
% end

idx = zeros(1,NQPTS);
for NQ=1:NQPTS
    VKPT_SYM = SYM_OP' * VKPT(:,NQ);
    VKPT_SYM = mod(VKPT_SYM,1);
    VKPT_SYM(VKPT_SYM>0.5+eps) = VKPT_SYM(VKPT_SYM>0.5+eps) - 1;
    VKPT_SYM(VKPT_SYM<-0.5+eps) = VKPT_SYM(VKPT_SYM<-0.5+eps) + 1;
    
    idx(NQ)=find(ismembertol(VKPTS',VKPT_SYM','ByRows',true));
end

end