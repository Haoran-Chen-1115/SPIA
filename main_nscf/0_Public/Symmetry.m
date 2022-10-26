%% Determine the symmetric operations of the system
function ...
    [IBRAV,N_SYM_TRUE,SYM_OP_TRUE,TRANS_ROT] ...
    = Symmetry(A,POSION,NTYP,NITYP)
% Determine the correct Bravais lattice type, because the specific axis is
% already fixed in Bravais
% Test all possible inversion and invertion

% Sometimes the original axis does not provide the highest symmetry
IBRAV_ORI = Bravais(A,inv(A));
IBRAV = IBRAV_ORI;

% P = perms([1,2,3]);
% for NP = 1:length(P)
%     for N1 = [-1,1]
%         for N2 = [-1,1]
%             for N3 = [-1,1]
%                 A_TRANS_tmp = zeros(3,3);
%                 A_TRANS_tmp(P(NP,1),1) = N1;
%                 A_TRANS_tmp(P(NP,2),2) = N2;
%                 A_TRANS_tmp(P(NP,3),3) = N3;
%                 
%                 if abs(det(A_TRANS_tmp)-1)<1E-10
%                     A_prime = A * A_TRANS_tmp';
%                     B_prime = inv(A_prime);
%                     IBRAV_tmp = Bravais(A_prime,B_prime);
%                     if IBRAV_tmp >= IBRAV
%                         % A higher index means higher symmetry in this case,
%                         % because body-centered cases never mix with
%                         % face-centered or others for the same lattice system
%                         IBRAV = IBRAV_tmp;
%                         A_new = A_prime;
%                         A_TRANS = A_TRANS_tmp;
%                     end
%                 end
%             end
%         end
%     end
% end
% For the algorithm above,
% if IBRAV = IBRAV_ORI, the A_TRANS = eye(3) must have been returned.
%A=A_new;
for N9=-2:2
    for N8=-2:2
        for N7=-2:2
            for N6=-2:2
                for N5=-2:2
                    ID1=N5*N9-N6*N8;
                    for N4=-2:2
                        ID2=N6*N7-N4*N9;
                        ID3=N4*N8-N5*N7;
                        for N3=-2:2
                            ID4=N3*ID3;
                            for N2=-2:2
                                ID5=N2*ID2+ID4;
                                for N1=-2:2
                                    if ((N1*ID1+ID5)==1)
                                        N=[N1,N2,N3;N4,N5,N6;N7,N8,N9];
                                        A_prime = A*N';
                                        B_prime = inv(A_prime);
                                        IBRAV_tmp = Bravais(A_prime,B_prime);
                                        if IBRAV_tmp >= IBRAV
                                            % A higher index means higher symmetry in this case,
                                            % because body-centered cases never mix with
                                            % face-centered or others for the same lattice system
                                            IBRAV = IBRAV_tmp;
                                            A_new = A_prime;
                                            A_TRANS = N;
                                        end
                                        
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%% First the pure point group symmetry of the Bravais lattice
% I follow the description of 
% https://www.cryst.ehu.es/cryst/get_point_genpos.html

% I found that in hexagonal and cubic cases, generators are different from 
% that in VASP, which should be checked some day

% In principle, the maximum symmetry should be set here

% For each Bravais lattice:
% using Seitz symbol
if IBRAV == 1
    % IBRAV = 1, Triclinic, C_i
    BRAV_GEN(:,:,1) = [-1,0,0;0,-1,0;0,0,-1];  % Inverse operation
    
elseif IBRAV == 2
    % IBRAV = 2, Primitive Monoclinic, C_2h
    % Rotation axis a2
    BRAV_GEN(:,:,1) = [-1,0,0;0,1,0;0,0,-1]; % 2_{010}
    BRAV_GEN(:,:,2) = [-1,0,0;0,-1,0;0,0,-1];  % -1
    
elseif IBRAV == 3
    % IBRAV = 3, Base-centered Monoclinic, C_2h
    % Rotation axis a2 - a1 = b*y
    BRAV_GEN(:,:,1) = [0,-1,0;-1,0,0;0,0,-1]; % 2_{1-10}
    BRAV_GEN(:,:,2) = [-1,0,0;0,-1,0;0,0,-1];  % -1
    
elseif IBRAV == 4
    % IBRAV = 4, Primitive Orthorhombic, D_2h
    BRAV_GEN(:,:,1) = [-1,0,0;0,-1,0;0,0,1]; % 2_{001}
    BRAV_GEN(:,:,2) = [-1,0,0;0,1,0;0,0,-1]; % 2_{010}
    BRAV_GEN(:,:,3) = [-1,0,0;0,-1,0;0,0,-1]; % -1

elseif IBRAV == 5
    % IBRAV = 5, Base-centered Orthorhombic, D_2h
    BRAV_GEN(:,:,1) = [-1,0,0;0,-1,0;0,0,1]; % 2_{001}
    BRAV_GEN(:,:,2) = [0,-1,0;-1,0,0;0,0,-1]; % 2_{010}:a1->-a2,a2->-a1,a3->-a3
    BRAV_GEN(:,:,3) = [-1,0,0;0,-1,0;0,0,-1]; % -1
    
elseif IBRAV == 6
    % IBRAV = 6, Body-centered Orthorhombic, D_2h
    BRAV_GEN(:,:,1) = [0,1,0;1,0,0;-1,-1,-1]; % 2_{001}:a1->a2,a2->a1,a3->-a1-a2-a3
    BRAV_GEN(:,:,2) = [0,0,1;-1,-1,-1;1,0,0]; % 2_{010}:a1->a3,a2->-a1-a2-a3,a3->a1
    BRAV_GEN(:,:,3) = [-1,0,0;0,-1,0;0,0,-1]; % -1
    
elseif IBRAV == 7
    % IBRAV = 7, Face-centered Orthorhombic, D_2h
    BRAV_GEN(:,:,1) = [0,1,-1;1,0,-1;0,0,-1]; % 2_{001}
    BRAV_GEN(:,:,2) = [0,-1,1;0,-1,0;1,-1,0]; % 2_{010}
    BRAV_GEN(:,:,3) = [-1,0,0;0,-1,0;0,0,-1]; % -1
    
elseif IBRAV == 8
    % IBRAV = 8, Primitive Tetragonal, D_4h
    BRAV_GEN(:,:,1) = [-1,0,0;0,-1,0;0,0,1]; % 2_{001}
    BRAV_GEN(:,:,2) = [0,-1,0;1,0,0;0,0,1]; % 4^{+}_{001}
    BRAV_GEN(:,:,3) = [-1,0,0;0,1,0;0,0,-1]; % 2_{010}
    BRAV_GEN(:,:,4) = [-1,0,0;0,-1,0;0,0,-1]; % -1
    
elseif IBRAV == 9
    % IBRAV = 9, Body-centered Tetragonal, D_4h
    BRAV_GEN(:,:,1) = [0,1,0;1,0,0;-1,-1,-1]; % 2_{001}
    BRAV_GEN(:,:,2) = [1,1,1;0,0,-1;-1,0,0]; % 4^{+}_{001}
    BRAV_GEN(:,:,3) = [0,0,1;-1,-1,-1;1,0,0]; % 2_{010}
    BRAV_GEN(:,:,4) = [-1,0,0;0,-1,0;0,0,-1]; % -1

elseif IBRAV == 10
    % IBRAV = 10, Rhombohedral
    % In rhombohedral axes
    BRAV_GEN(:,:,1) = [0,0,1;1,0,0;0,1,0]; % 3^{+}_{111}
    BRAV_GEN(:,:,2) = [0,-1,0;-1,0,0;0,0,-1]; % 2_{1-10}
    BRAV_GEN(:,:,3) = [-1,0,0;0,-1,0;0,0,-1]; % -1

elseif IBRAV == 11
    % IBRAV = 11, Hexagonal, D_6h
    BRAV_GEN(:,:,1) = [0,-1,0;1,-1,0;0,0,1]; % 3^{+}_{001}
    BRAV_GEN(:,:,2) = [-1,0,0;0,-1,0;0,0,1]; % 2_{001}
    BRAV_GEN(:,:,3) = [0,1,0;1,0,0;0,0,-1]; % 2_{110}
    BRAV_GEN(:,:,4) = [-1,0,0;0,-1,0;0,0,-1]; % -1
    % Different from that in VASP, where only three generators are used
    % Check some day

elseif IBRAV == 12
    % IBRAV = 12, Primitive Cubic, O_h
    BRAV_GEN(:,:,1) = [-1,0,0;0,-1,0;0,0,1]; % 2_{001}
    BRAV_GEN(:,:,2) = [-1,0,0;0,1,0;0,0,-1]; % 2_{010}
    BRAV_GEN(:,:,3) = [0,0,1;1,0,0;0,1,0]; % 3^{+}_{111}
    BRAV_GEN(:,:,4) = [0,1,0;1,0,0;0,0,-1]; % 2_{110}
    BRAV_GEN(:,:,5) = [-1,0,0;0,-1,0;0,0,-1]; % -1

elseif IBRAV == 13
    % IBRAV = 13, Body-centered Cubic, O_h
    BRAV_GEN(:,:,1) = [0,1,0;1,0,0;-1,-1,-1]; % 2_{001}
    BRAV_GEN(:,:,2) = [0,0,1;-1,-1,-1;1,0,0]; % 2_{010}
    BRAV_GEN(:,:,3) = [0,0,1;1,0,0;0,1,0]; % 3^{+}_{111}
    BRAV_GEN(:,:,4) = [-1,0,0;0,-1,0;1,1,1]; % 2_{110}
    BRAV_GEN(:,:,5) = [-1,0,0;0,-1,0;0,0,-1]; % -1
    
elseif IBRAV == 14
    % IBRAV = 14, Face-centered Cubic, O_h
    BRAV_GEN(:,:,1) = [0,1,-1;1,0,-1;0,0,-1]; % 2_{001}
    BRAV_GEN(:,:,2) = [0,-1,1;0,-1,0;1,-1,0]; % 2_{010}
    BRAV_GEN(:,:,3) = [0,0,1;1,0,0;0,1,0]; % 3^{+}_{111}
    BRAV_GEN(:,:,4) = [-1,0,1;0,-1,1;0,0,1]; % 2_{110}
    BRAV_GEN(:,:,5) = [-1,0,0;0,-1,0;0,0,-1]; % -1
end

N_GEN = size(BRAV_GEN,3);

% %% Exchange element according to the transformation matrix
% for NG = 1:N_GEN
%     BRAV_GEN(:,:,NG) = A_TRANS\BRAV_GEN(:,:,NG)/A_TRANS';
% end

%% Now, generate all possible symmetry operations using the generators.
% Cover every multiplication of the generators
SYM_OP(:,:,1) = eye(3);
N_SYM = 1;

for NG = 1:N_GEN
    BRAV_SYM = BRAV_GEN(:,:,NG);
    
    % In case that the operation already exists
    if ismember(reshape(BRAV_SYM,1,9),reshape(SYM_OP,9,[])','rows')
        continue
    end
    
    GEN_order = 0;
    while 1
        GEN_order = GEN_order + 1;
        if ismember(reshape(BRAV_SYM,1,9),reshape(eye(3),[],9),'rows')
            break
        end
        
        BRAV_SYM = BRAV_GEN(:,:,NG) * BRAV_SYM;
    end
    
    N_SYM_tmp = N_SYM; 
    for NS1 = 1:N_SYM
        BRAV_SYM = SYM_OP(:,:,NS1);
        
        for NO = 1:GEN_order - 1
            BRAV_SYM = BRAV_SYM * BRAV_GEN(:,:,NG);
            
            for NS2 = 1:N_SYM
                BRAV_SYM = BRAV_SYM * SYM_OP(:,:,NS2);
                if ~ismember(reshape(BRAV_SYM,1,9),reshape(SYM_OP,9,[])','rows')
                    N_SYM_tmp = N_SYM_tmp + 1;
                    SYM_OP(:,:,N_SYM_tmp) =  BRAV_SYM;
                end
            end
            
        end
        
    end
    N_SYM = N_SYM_tmp;
end


%% Corresponding operation in the original lattice
for NS = 1:N_SYM
    SYM_OP(:,:,NS) = A_TRANS\SYM_OP(:,:,NS)*A_TRANS;
end

%% Check whether the symmetry operation is allowed
% for the Bravais lattice with basis

% Act on all the atoms in the primitive cell.
% I'm not sure whether it works for a supercell or not.
SYM_tol = 1E-4;
POS_tmp = POSION;
for NT=1:NTYP
    POS_tmp{NT}(abs(POS_tmp{NT}-1)<SYM_tol)...
        = POS_tmp{NT}(abs(POS_tmp{NT}-1)<SYM_tol) - 1;
end

% Sort with a tolerance
POS_tmp = sort_POS(POS_tmp,NTYP,SYM_tol);
POS_ORI = cell2mat(POS_tmp);

% NIS = 0;
% for NT=1:NTYP
%     POS_ORI(:,NIS+1:NIS+NITYP(NT)) = uniquetol(POS_tmp{NT}',1E-8,'ByRows',true)';
%     NIS = NIS + NITYP(NT);
% end

N_SYM_TRUE = 1;
TRANS_ROT(:,N_SYM_TRUE) = [0;0;0];
SYM_OP_TRUE(:,:,N_SYM_TRUE) = SYM_OP(:,:,1);
for NS = 2:N_SYM
    LABORT = false;
    POS_SYM = cell(1,NTYP);
    for NT = 1:NTYP
        POS_SYM{NT} = mod(A_TRANS'\SYM_OP(:,:,NS)' * POSION{NT},1);
        POS_SYM{NT}(abs(POS_SYM{NT}-1)<SYM_tol) ...
            = POS_SYM{NT}(abs(POS_SYM{NT}-1)<SYM_tol) - 1;
        POS_SYM{NT} = A_TRANS'*POS_SYM{NT};
    end
    
    % Check whether translation will take back the ions
    [NITYP_B,NT_B] = max(NITYP);
    
    % First try if the operation is a pure point group operation
    T1 = [0;0;0];
    POS_TRY = cell(1,NTYP);
    for NT = 1:NTYP
        POS_TRY{NT} = mod(POS_SYM{NT} - T1,1);
        POS_TRY{NT}(abs(POS_TRY{NT}-1)<SYM_tol) ...
            = POS_TRY{NT}(abs(POS_TRY{NT}-1)<SYM_tol) - 1;
    end
    POS_TRY = sort_POS(POS_TRY,NTYP,SYM_tol);
    POS_TRY = cell2mat(POS_TRY);
    
    if MAT_EQ(POS_TRY, POS_ORI, SYM_tol)
        N_SYM_TRUE = N_SYM_TRUE + 1;
        TRANS_ROT(:,N_SYM_TRUE) = T1;
        SYM_OP_TRUE(:,:,N_SYM_TRUE) = SYM_OP(:,:,NS);
        LABORT = true;
    end
    
    % Then space group operation
    if ~LABORT
        for NS_s = 1:N_SYM_TRUE
            T1 = TRANS_ROT(:,NS_s);
            POS_TRY = cell(1,NTYP);
            for NT = 1:NTYP
                POS_TRY{NT} = mod(POS_SYM{NT} - T1,1);
                POS_TRY{NT}(abs(POS_TRY{NT}-1)<SYM_tol) ...
                    = POS_TRY{NT}(abs(POS_TRY{NT}-1)<SYM_tol) - 1;
            end
            POS_TRY = sort_POS(POS_TRY,NTYP,SYM_tol);
            POS_TRY = cell2mat(POS_TRY);
            
            if MAT_EQ(POS_TRY, POS_ORI, SYM_tol)
                N_SYM_TRUE = N_SYM_TRUE + 1;
                TRANS_ROT(:,N_SYM_TRUE) = T1;
                SYM_OP_TRUE(:,:,N_SYM_TRUE) = SYM_OP(:,:,NS);
                LABORT = true;
                break
            end
        end
    end
    
    if ~LABORT
        for NI = 1:NITYP_B
            T1 = POS_SYM{NT_B}(:,NI)-POS_tmp{NT_B}(:,1);
            POS_TRY = cell(1,NTYP);
            for NT = 1:NTYP
                POS_TRY{NT} = mod(POS_SYM{NT} - T1,1);
                POS_TRY{NT}(abs(POS_TRY{NT}-1)<SYM_tol) ...
                    = POS_TRY{NT}(abs(POS_TRY{NT}-1)<SYM_tol) - 1;
            end
            POS_TRY = sort_POS(POS_TRY,NTYP,SYM_tol);
            POS_TRY = cell2mat(POS_TRY);
            
            if MAT_EQ(POS_TRY, POS_ORI, SYM_tol)
                N_SYM_TRUE = N_SYM_TRUE + 1;
                TRANS_ROT(:,N_SYM_TRUE) = T1;
                SYM_OP_TRUE(:,:,N_SYM_TRUE) = SYM_OP(:,:,NS);
                break
            end
            
        end
    end
    
end 
end

function LEQ = MAT_EQ(A,B,tol)
    LEQ = all(reshape(abs(A-B)<tol ,1,[]));
end

function POS_tmp = sort_POS(POS_tmp,NTYP,tol)
for NT=1:NTYP
    % a1 is the fastest index
    [~,idx1] = sort(POS_tmp{NT}(1,:),'ascend');
    POS_tmp{NT} = POS_tmp{NT}(:,idx1);
    [P1,~,P1_IC]=uniquetol(POS_tmp{NT}(1,:),tol);
    
    % sort over a2
    for N2 = 1:length(P1)
        POS_tmp2 = POS_tmp{NT}(:,P1_IC==N2);
        [~,idx2] = sort(POS_tmp2(2,:),'ascend');
        POS_tmp{NT}(:,P1_IC==N2) = POS_tmp2(:,idx2);
    end
    
    % sort over a3
    [P2,~,P2_IC]=uniquetol(POS_tmp{NT}(1:2,:)',tol,'ByRows',true);
    for N3 = 1:size(P2,1)
        POS_tmp2 = POS_tmp{NT}(:,P2_IC==N3);
        [~,idx2] = sort(POS_tmp2(3,:),'ascend');
        POS_tmp{NT}(:,P2_IC==N3) = POS_tmp2(:,idx2);
    end
end

end
