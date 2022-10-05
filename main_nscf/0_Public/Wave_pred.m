%% Determine the alignment parameter \alpha and \beta
ICALL_out=ICALL_out+1;
if ICALL_out<2
    
elseif ICALL_out<3
    % Linear extrapolation
    DP_1=mod(POS_now-POS_L1+1.5,1)-0.5;
    DP_2=mod(POS_L1-POS_L2+1.5,1)-0.5;
    
    A0=reshape(DP_1,[],1)'*reshape(DP_1,[],1);
    A1=-2*reshape(DP_1,[],1)'*reshape(DP_2,[],1);
    A2=reshape(DP_2,[],1)'*reshape(DP_2,[],1);
    alpha=-A1/2/A2;
    beta=0;
else
    DP_1=POS_now-POS_L1;
    DP_2=POS_L1-POS_L2;
    DP_3=POS_L2-POS_L3;
    
    A0=reshape(DP_1,[],1)'*reshape(DP_1,[],1); % A0= |Δr_(n+1)|^2
    A1=-2*reshape(DP_1,[],1)'*reshape(DP_2,[],1); % A1= -2Δr_n·Δr_(n+1)
    A2=reshape(DP_2,[],1)'*reshape(DP_2,[],1); % A2= |Δr_n|^2
    
    B1=-2*reshape(DP_1,[],1)'*reshape(DP_3,[],1); % B1= -2Δr_(n-1)·Δr_(n+1)
    B2=reshape(DP_3,[],1)'*reshape(DP_3,[],1); % B2= |Δr_(n-1)|^2
    
    AB=2*reshape(DP_2,[],1)'*reshape(DP_3,[],1); % AB= 2Δr_(n-1)·Δr_n
    
    if abs((A0+A2+A1)/A0)<1E-4 && abs((A0+B2+B1)/A0)<1E-4
        % parrallel position change
        alpha=2;
        beta=-1;
    else
        % minimize ||αΔr_n + βΔr_(n-1) - Δr_(n+1)||^2
        % = α^2|Δr_n|^2 + β^2|Δr_(n-1)|^2 + |Δr_(n+1)|^2
        %    + 2αβΔr_n·Δr_(n-1）
        %    - 2αΔr_n·Δr_(n+1）- 2βΔr_（n-1）·Δr_(n+1）
        alpha = (2*A1*B2-B1*AB)/(AB^2-4*A2*B2);
        beta  = (2*B1*A2-A1*AB)/(AB^2-4*A2*B2);
        
    end
    
end

%% Next we will minimize Σ_m|ψ_m(t_n) - ψ_m(t_(n-1))|
% for a continous change in wavefunction.

% Calculate matrix U_nm=<m|T^† T'|n'>
% In VASP, where an approximate form is calculated
% U_nm≈<m|n'>+Σ_ij (Σ_a<m|p_i^a>) Q_ij (Σ_a'<p_j'^a'|n'>).
% This is to assume that the displacement is very small.
% I'll calculate the accurate one here.
D_BASIS_L1=zeros(NPL,NBANDS);
ICALL = min(ICALL,ICALL_out-1);

if ICALL_out>=2
    TRANS_L1=eye(NPL,'single','gpuArray');
    TRANS_L2=eye(NPL,'single','gpuArray');
%     TdTp=eye(NPL,'single','gpuArray');
    NIS=0;
    POSION_1=cell(1,NTYP);POSION_2=cell(1,NTYP);
    for NT=1:NTYP
        POSION_1{NT}=POS_L1(:,NIS+1:NIS+NITYP(NT));
        POSION_2{NT}=POS_L2(:,NIS+1:NIS+NITYP(NT));
        
        CREXP_L1=single(exp(-CITPI*G_INDEX(1:3,1:NPL).'...
            *gpuArray(POSION_1{NT})));
        CREXP_L2=single(exp(-CITPI*G_INDEX(1:3,1:NPL).'...
            *gpuArray(POSION_2{NT})));
        
        %             TRANS_L2 = TRANS_L1;
        TRANS_L2=TRANS_L2+(CREXP_L2*CREXP_L2').*TRANS_CORE{NT};
        TRANS_L1=TRANS_L1+(CREXP_L1*CREXP_L1').*TRANS_CORE{NT};
        
%         TdTp=TdTp+(CREXP_L2*CREXP_L1').*FOVL_CORE{NT};
        
        NIS=NIS+NITYP(NT);
    end
    %     TdTp = TRANS_L2'+ TRANS_L1 + TdTp;
    %     clear TRANS_c_L1 TRANS_c_L2
    
    TdTp=OVERL(POSION_1, POSION_2, TRANS_L1, TRANS_L2, A, ...
        CQIJ3_fun, r, G_INDEX, CITPI, QPROJ_ORI, ...
        NTYP,NITYP,NPL,LMMAXC, ...
        CH0,CH1,CH2,CH3);
    clear TRANS_L2 TRANS_L1
    
    U = BASIS_L1' * TdTp * BASIS_L2;
    
    UU = U'*U; % The matrix we are to diagonize.

    UU=triu(UU,1)'+triu(UU,1)+diag(real(diag(UU)));
    
    [AA,CosT2]=eig(UU);
    AA=AA';
    CosT2=diag(CosT2);
    if ~isempty(find(CosT2<0))
        ICALL(:)=0;
    end
    
    CosT=abs(sqrt(complex(CosT2)));
%     CosT=sqrt(CosT2);
    AA_prime = diag(1./CosT) * AA * U';
    
    AA_prime = (AA'*AA_prime).'; % inv(A')=A'^H;
    
    %% Rotate wavefunctions and perform alignment
    BASIS_L2 = BASIS_L2 * AA_prime;
    D_BASIS_L2 = D_BASIS_L2 * AA_prime;
    
    U = diag(BASIS_L2' * TdTp * BASIS_L1);
    clear TdTp
    idx_cross = find(abs(U)<0.9);
    
    ICALL = ICALL + 1;
    ICALL(idx_cross) = 1;
    
    D_BASIS_L1 = BASIS_L1-BASIS_L2;
    D_BASIS_L1(:,idx_cross) = 0;
    
    % extrapolate
    idx_interp2 = find(ICALL>2);
    idx_interp1 = find(ICALL==2);
    
    alpha0 = zeros(1,NBANDS); alpha0(idx_interp1) = 1; alpha0(idx_interp2) = alpha;
    beta0  = zeros(1,NBANDS); beta0(idx_interp2) = beta;
    BASIS_trial = BASIS_L1 + alpha0.*D_BASIS_L1 + beta0.*D_BASIS_L2;
    
    POS_L3 = POS_L2;
else
    ICALL = 1;
    BASIS_trial = BASIS_L1;
end

BASIS_L2 = BASIS_L1;
D_BASIS_L2 = D_BASIS_L1;
POS_L2 = POS_L1;
POS_L1 = POS_now;