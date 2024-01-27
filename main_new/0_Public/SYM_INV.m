function SYM_OP_INV = SYM_INV(SYM_OP_TRUE,A)
B = inv(A);
N_SYM_TRUE = size(SYM_OP_TRUE,3);
SYM_OP_tmp = zeros(size(SYM_OP_TRUE));
for NS = 1:N_SYM_TRUE
    SYM_OP_tmp(:,:,NS) = (B*B') * SYM_OP_TRUE(:,:,NS) * (A'*A);
end
SYM_OP_INV=round(SYM_OP_tmp);
if any(reshape(abs(SYM_OP_INV-SYM_OP_tmp)>1E-5,1,[]))
    error('Inverse operation wrong!')
end
end