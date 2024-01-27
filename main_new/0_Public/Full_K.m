function [VKPTS_f,NKPTS_f]=Full_K(N_KP,IFLAG)
NKPTS_f = N_KP(1)*N_KP(1)*N_KP(1);

K1 = (0:N_KP(1)-1)/N_KP(1);
K2 = (0:N_KP(2)-1)/N_KP(2);
K3 = (0:N_KP(3)-1)/N_KP(3);

eps = 1E-5;
if IFLAG == 1 % k-mesh (consistent with VASP)
    VKPTS_f(1,:) = reshape(repmat(K1,1,N_KP(2)*N_KP(3)),1,[]);
    VKPTS_f(2,:) = reshape(repmat(K2,N_KP(1),N_KP(3)),1,[]);
    VKPTS_f(3,:) = reshape(repmat(K3,N_KP(1)*N_KP(2),1),1,[]);
    VKPTS_f(VKPTS_f>0.5+eps) = VKPTS_f(VKPTS_f>0.5+eps) - 1;
elseif IFLAG == 2 % q-mesh (consistent with QE)
    VKPTS_f(1,:) = reshape(repmat(K1,N_KP(2)*N_KP(3),1),1,[]);
    VKPTS_f(2,:) = reshape(repmat(K2,N_KP(3),N_KP(1)),1,[]);
    VKPTS_f(3,:) = reshape(repmat(K3,1,N_KP(1)*N_KP(2)),1,[]);
    VKPTS_f(VKPTS_f>0.5+eps) = VKPTS_f(VKPTS_f>0.5+eps) - 1;
end
end