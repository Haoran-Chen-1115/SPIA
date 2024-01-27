function [V,E,V_proj,N_notconv]=Green_c(ROOT_DIR,NC_c,NB_c,ISP_c,...
    NTYP,NITYP,LMAX,LMMAXC,CH0,CH1,CH2,CH3,LPS,...
    NPL,NPLWV,NGX,NGY,NGZ,NCPU,...
    QPROJ_ORI,TRANS_CORE,...
    POSION,G_INDEX,G,HSQDTM,CITPI,...
    POSION_EQ, TRANS_EQ, A, B, ...
    CQIJ3_fun, r, RYTOEV, V, NBANDS, NELE)
%% Calculate Hamiltonian, Overlap Matrix and Transformation Matrix
[FHAM,TRANS]=HAMILT(ROOT_DIR,NC_c,NB_c,ISP_c,...
    NTYP,NITYP,LMAX,LMMAXC,CH0,CH1,CH2,CH3,LPS,...
    NPL,NPLWV,NGX,NGY,NGZ,NCPU,...
    QPROJ_ORI,TRANS_CORE,...
    POSION,G_INDEX,G,HSQDTM,CITPI);

% Rigorous correction for overlap matrix
FOVL=OVERL(POSION, POSION, TRANS, TRANS, A, B, ...
    CQIJ3_fun, r, G_INDEX, CITPI, QPROJ_ORI, ...
    NTYP,NITYP,NPL,LMAX,LMMAXC, ...
    CH0,CH1,CH2,CH3,LPS);
%FOVL=triu(FOVL,1)+triu(FOVL,1)'+diag(real(diag(FOVL)));
FOVL = (FOVL'+FOVL)/2;
FHAM = (FHAM'+FHAM)/2;

%% Eigenvalues
[V,E,N_notconv]=Eig_DAV(FHAM,FOVL, V, B, NBANDS, NELE, NPL, RYTOEV);

%% Expand with equilibrium ones
TdT0=OVERL(POSION, POSION_EQ, TRANS, TRANS_EQ, A, B, ...
    CQIJ3_fun, r, G_INDEX, CITPI, QPROJ_ORI, ...
    NTYP,NITYP,NPL,LMAX,LMMAXC, ...
    CH0,CH1,CH2,CH3,LPS);

V_proj = abs(TdT0'*V).^2;
end

function [FHAM,TRANS]=HAMILT(ROOT_DIR,NC_c,NB_c,ISP_c,...
    NTYP,NITYP,LMAX,LMMAXC,CH0,CH1,CH2,CH3,LPS,...
    NPL,NPLWV,NGX,NGY,NGZ,NCPU,...
    QPROJ_ORI,TRANS_CORE,...
    POSION,G_INDEX,G,HSQDTM,CITPI)

if NCPU==1
    SVV=fopen([ROOT_DIR,'/BEAD_',NC_c,'/SV_',NB_c,'_',ISP_c]);
    SV=reshape(fread(SVV,NPLWV,'double'),NGX,NGY,NGZ);
    fclose(SVV);
else
    SV = cell(1,NCPU);SV_idx=cell(1,NCPU);
    for NU=1:NCPU
        %
        NROW=NGZ;
        NPX=0;
        for N2=1:NGX
            NODE_TARGET=mod(N2-1,NCPU)+1;
            if NODE_TARGET==NU
                NPX=NPX+1;
                SV_idx{NU}(NPX)=N2;
            end
        end
        NCOL = NPX*NGY;
        NP = NCOL*NROW;
        
        %
        NU_c=int2str(NU);
        SVV=fopen([ROOT_DIR,'/BEAD_',NC_c,'/SV_',NB_c,'_',ISP_c,'_',NU_c]);
        SV{NU}=reshape(fread(SVV,NP,'double'),NGZ,NPX,NGY);
        fclose(SVV);
        SV{NU}=reshape(permute(SV{NU},[3,1,2]),[],1); 
    end
    SV_idx=cell2mat(SV_idx);
    SV=reshape(cell2mat(SV.'),NGY,NGZ,NGX);
    SV(:,:,SV_idx)=SV;
    SV=permute(SV,[3,1,2]);
end

%% Local Part and Pseudo Kinetic Energy
% Calculate local potential
%                 SV=gpuArray(SV);
SV=gpuArray(single(SV));
SVF=reshape(fftn(SV)/NPLWV,1,NPLWV);

FIND=mod(bsxfun(@minus,G_INDEX(1,:).',G_INDEX(1,:))+NGX,NGX)+1+...
    mod(bsxfun(@minus,G_INDEX(2,:).',G_INDEX(2,:))+NGY,NGY)*NGX+...
    mod(bsxfun(@minus,G_INDEX(3,:).',G_INDEX(3,:))+NGZ,NGZ)*NGX*NGY;
FHAM=SVF(FIND);
FHAM=FHAM+diag(G(4,:)*HSQDTM);
clear SVF FIND KIN

%% Nonlocal Part, Overlap and Transformation Matrix
%                 FOVL=single(eye(NPL,'gpuArray'));
TRANS=eye(NPL,'single','gpuArray');
%                 TRANS_c=single(eye(NPL,'gpuArray'));
%                 FOVL=eye(NPL,'gpuArray');
%                 TRANS=eye(NPL,'gpuArray');
%                 TRANS_c=eye(NPL,'gpuArray');

for NT=1:NTYP
    NT_c=int2str(NT);
    
    if NCPU==1
        CDIJJ=fopen([ROOT_DIR,'/BEAD_',NC_c,'/CDIJ_',NB_c,'_',ISP_c,'_',NT_c]);
        % Non local pseudopotential strength and Overlap Strength
        %                     for NT1=1:NT-1 % Because of my miss when implementing vasp output
        %                         fread(CDIJJ,...
        %                             [LMMAXC(NT)*LMMAXC(NT),NITYP(NT1)],'double');
        %                     end
        CDIJ=reshape(fread(CDIJJ,...
            [LMMAXC(NT)*LMMAXC(NT),NITYP(NT)],'double'),LMMAXC(NT),LMMAXC(NT),NITYP(NT));
        
        fclose(CDIJJ);
    else
        CDIJ=cell(1,NCPU);
        for NU=1:NCPU
            NU_c=int2str(NU);
            
            CDIJJ=fopen([ROOT_DIR,'/BEAD_',NC_c,'/CDIJ_',NB_c,'_',ISP_c,'_',NT_c,'_',NU_c]);
            CDIJ{NU}=fread(CDIJJ,...
                [LMMAXC(NT)*LMMAXC(NT)*NITYP(NT),1],'double');
            NWRITTEN = length(CDIJ{NU});
            CDIJ{NU}(end+1:end+LMMAXC(NT)*LMMAXC(NT)*NITYP(NT)-NWRITTEN,1)=0;
            
            fclose(CDIJJ);
        end
        CDIJ=reshape(cell2mat(CDIJ),LMMAXC(NT),LMMAXC(NT),NITYP(NT),NCPU);
        CDIJ=sum(CDIJ,4);
    end
    % Calculate Non local part of Hamiltonian
    % QPROJ_ORI{NT}=gpuArray(QPROJ_ORI{NT});
    %                     CDIJ=gpuArray(CDIJ);
    CDIJ=gpuArray(single(CDIJ));
    
    % Calculate phase
    %                     POSION1=single(gpuArray(POSION{NT}));
    POSION1=gpuArray(POSION{NT});
    CREXP=single(exp(CITPI*G_INDEX(1:3,1:NPL).'*POSION1));
    
    % 'Phasefactor' ignored in spherical Bessel function
    CQFAK=zeros(1,LMMAXC(NT),'single','gpuArray');
    %                     CQFAK=zeros(1,LMMAXC(NT),'gpuArray');
    NCHM=1;
    LCH=LPS{NT};
    for NL=1:LMAX(NT)
        if LCH(NL)==0
            CQFAK(NCHM)=1;
            NCHM=NCHM+1;
        elseif LCH(NL)==1
            CQFAK(NCHM:NCHM+2)=1i;
            NCHM=NCHM+3;
        elseif LCH(NL)==2
            CQFAK(NCHM:NCHM+4)=-1;
            NCHM=NCHM+5;
        elseif LCH(NL)==3
            CQFAK(NCHM:NCHM+6)=-1i;
            NCHM=NCHM+7;
        end
    end
    % CQFAK(1:CH0{NT}(2))=1;
    % CQFAK(CH0{NT}(2)+1:CH0{NT}(2)+3*CH1{NT}(2))=1i;
    % CQFAK(CH0{NT}(2)+3*CH1{NT}(2)+1:CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2))=-1;
    % CQFAK(CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2)+1:...
    %     CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2)+7*CH3{NT}(2))=-1i;
    
    QPROJ=bsxfun(@times,QPROJ_ORI{NT},reshape(CQFAK,1,LMMAXC(NT)));
    
    QPROJ2=bsxfun(@times,...
        reshape(QPROJ(1:NPL,:),NPL,LMMAXC(NT),1),...
        reshape(CREXP(1:NPL,:),NPL,1,NITYP(NT)));
    QPROJ21=permute(QPROJ2,[2,1,3]);
    
    D_Q1=reshape(sum(bsxfun(@times,...
        reshape(CDIJ,LMMAXC(NT),1,LMMAXC(NT),NITYP(NT)),...
        reshape(QPROJ21,LMMAXC(NT),NPL,1,NITYP(NT))),1),NPL,LMMAXC(NT),NITYP(NT));
    
    QPROJ21=conj(reshape(QPROJ2,NPL,LMMAXC(NT)*NITYP(NT)));
    D_Q1=reshape(D_Q1,NPL,LMMAXC(NT)*NITYP(NT)).';
    
    FHAM=FHAM+QPROJ21*D_Q1;
    
    clear QPROJ21 QPROJ2 QPROJ D_Q1
    
    % Calculate T
    % Transformation matrix
    CREXP=conj(CREXP)*CREXP.';
    %TRANS_CORE=(RES{NT}*QPROJ_ORI{NT}.')/sqrt(OMEGA);
    TRANS=TRANS+CREXP.*TRANS_CORE{NT};
    % For \hat{T}_c=|\tilde{p}_i^a><\tilde{\phi}_i^a|...
    %                     TRANS_c=TRANS_c+CREXP.*TRANS_c_CORE{NT};
    % Overlap matrix
    %                     FOVL=FOVL+CREXP.*FOVL_CORE{NT};
    
    %TRANS=gather(TRANS);
    clear CREXP
end

% Replace S by S'=S+(\hat{T}-\hat{T}_c+h.c.)
%                 TRANS_c=TRANS-TRANS_c;
%                 FOVL=FOVL+TRANS_c+TRANS_c';
%                 clear TRANS_c

%FHAM=triu(FHAM,1)+triu(FHAM,1)'+diag(real(diag(FHAM)));
%                 FOVL=triu(FOVL)+triu(FOVL,1)';

end

function TdT0=OVERL(POS_1, POS_2, TRANS_1, TRANS_2, A, B, ...
    CQIJ3_fun, r, G_INDEX, CITPI, QPROJ_ORI, ...
    NTYP,NITYP_1, NPL,LMAX,LMMAXC, ...
    CH0,CH1,CH2,CH3,LPS)

CQIJ_DIV=cell(NTYP,NTYP);
idx_phase=cell(NTYP,NTYP);
TdT0=zeros(NPL,NPL,'single','gpuArray');
%TdT0_2=zeros(NPL,NPL,'single','gpuArray');
CQFAK_=cell(1,NTYP);
for NT=1:NTYP
    CQFAK_{NT}=zeros(1,LMMAXC(NT),'single');
    %                     CQFAK=zeros(1,LMMAXC(NT),'gpuArray');
    NCHM=1;
    LCH=LPS{NT};
    for NL=1:LMAX(NT)
        if LCH(NL)==0
            CQFAK_{NT}(NCHM)=1;
            NCHM=NCHM+1;
        elseif LCH(NL)==1
            CQFAK_{NT}(NCHM:NCHM+2)=1i;
            NCHM=NCHM+3;
        elseif LCH(NL)==2
            CQFAK_{NT}(NCHM:NCHM+4)=-1;
            NCHM=NCHM+5;
        elseif LCH(NL)==3
            CQFAK_{NT}(NCHM:NCHM+6)=-1i;
            NCHM=NCHM+7;
        end
    end

    % CQFAK_{NT}(1:CH0{NT}(2))=1;
    % CQFAK_{NT}(CH0{NT}(2)+1:CH0{NT}(2)+3*CH1{NT}(2))=1i;
    % CQFAK_{NT}(CH0{NT}(2)+3*CH1{NT}(2)+1:CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2))=-1;
    % CQFAK_{NT}(CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2)+1:...
    %     CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2)+7*CH3{NT}(2))=-1i;
%     CQFAK_{NT}(CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2)+1:...
%         CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2)+7*CH3{NT}(2))=-1i;
end

for NT1=1:NTYP
    for NT2=1:NTYP
        % For now, the positions are stored in a 'direct'
        % manner
        DX=-POS_1{NT1}(1,:).'+POS_2{NT2}(1,:);
        DY=-POS_1{NT1}(2,:).'+POS_2{NT2}(2,:);
        DZ=-POS_1{NT1}(3,:).'+POS_2{NT2}(3,:);
        
        % This is only right for a cuboid lattice
        DX(DX>1/2)=-1+DX(DX>1/2);
        DY(DY>1/2)=-1+DY(DY>1/2);
        DZ(DZ>1/2)=-1+DZ(DZ>1/2);
        DX(DX<-1/2)=1+DX(DX<-1/2);
        DY(DY<-1/2)=1+DY(DY<-1/2);
        DZ(DZ<-1/2)=1+DZ(DZ<-1/2);
        
        DX_2=A(1,1)*DX+A(1,2)*DY+A(1,3)*DZ;
        DY_2=A(2,1)*DX+A(2,2)*DY+A(2,3)*DZ;
        DZ_2=A(3,1)*DX+A(3,2)*DY+A(3,3)*DZ;
        
        DX=DX_2;DY=DY_2;DZ=DZ_2;
        
        DR=sqrt(DX.^2+DY.^2+DZ.^2);
        %if NT1==NT2
        %    DR=triu(DR)...
        %        +triu((r{NT1,NT2}(end)+1)*ones(NITYP_1(NT1)),1).';
        %end
        %                         DR(DR<1E-8)=1E-10;
        idx=find(DR<r{NT1,NT2}(end));
        idx_phase{NT1,NT2}(1,:)=mod(idx-1,NITYP_1(NT1))+1;
        idx_phase{NT1,NT2}(2,:)=ceil(idx/NITYP_1(NT1));
        % index for calculate overall phase factor
        %         [NI2,NI2_ia,NI2_ic]=unique(idx_phase{NT1,NT2}(2,:));
        
        DP=[DX(idx),DY(idx),DZ(idx),DR(idx)];
        NP=length(idx);
        CQIJ_DIV{NT1,NT2}=zeros(LMMAXC(NT1),LMMAXC(NT2),NP);
        
        CQIJ_interp=zeros(LMMAXC(NT1),LMMAXC(NT2),NP);
        for i=1:LMMAXC(NT1)
            for j=1:LMMAXC(NT2)
                CQIJ_interp(i,j,:) = CQIJ3_fun{NT1,NT2}{i,j}(DP(:,4));
            end
        end
        
        %if NT1==NT2
        %    for i=1:NP
        %        R=DP(i,:);
        %        Rotate=Rotate_SH(R,LMAX,LMMAXC,LPS,NT1);
        %        CQIJ_DIV{NT1,NT2}(:,:,i)=...
        %            Rotate*CQIJ_interp(:,:,i)*(Rotate.');
        %    end
        %else
            for i=1:NP
                R=DP(i,:);
                Rotate1=Rotate_SH(R,LMAX,LMMAXC,LPS,NT1);
                Rotate2=Rotate_SH(R,LMAX,LMMAXC,LPS,NT2);
                CQIJ_DIV{NT1,NT2}(:,:,i)=...
                    Rotate1*CQIJ_interp(:,:,i)*(Rotate2.');
            end
        %end
        
        % Different from that in the previous part,
        % the phase factor does not sum over all ions,
        % because the kernals are different.
        POSION_NT1=POS_1{NT1}(:,idx_phase{NT1,NT2}(1,:));
        CREXP_first=single(exp(-CITPI*G_INDEX(1:3,1:NPL).'*POSION_NT1));
        CREXP_second=single(exp(-1i*G_INDEX(1:3,1:NPL).'*B*(DP(:,1:3).')));
        CREXP_second=conj(CREXP_first.*CREXP_second);
        
        %toc;
        
        QPROJ_NT2=QPROJ_ORI{NT2}.*reshape(CQFAK_{NT2},1,LMMAXC(NT2)); % NPL,LMMAXC(2)
        QPROJ2_NT2=reshape(QPROJ_NT2,NPL,LMMAXC(NT2),1)...
            .*reshape(CREXP_second,NPL,1,NP); % NPL,LMMAXC(2),NP
        QPROJ2_NT2=permute(QPROJ2_NT2,[2,1,3]); % LMMAXC(2),NPL,NP
        
        %Q_Q1=multiprod(CQIJ_DIV{NT1,NT2},QPROJ2_NT2); % LMMAXC(1),NPL,NP
        Q_Q1=pagefun(@mtimes,CQIJ_DIV{NT1,NT2},QPROJ2_NT2); % LMMAXC(1),NPL,NP
        Q_Q1=reshape(permute(Q_Q1,[1,3,2]),LMMAXC(NT1)*NP,NPL); % LMMAXC(1)*NP,NPL
        
        QPROJ_NT1=QPROJ_ORI{NT1}.*reshape(conj(CQFAK_{NT1}),1,LMMAXC(NT1)); % NPL,LMMAXC(1)
        QPROJ2_NT1=reshape(QPROJ_NT1,NPL,LMMAXC(NT1),1)...
            .*reshape(CREXP_first,NPL,1,NP); % NPL,LMMAXC(1),NP
        QPROJ2_NT1=reshape(QPROJ2_NT1,NPL,LMMAXC(NT1)*NP); % NPL,LMMAXC(1)*NP
        
        TdT0 = TdT0 + QPROJ2_NT1*Q_Q1; % NPL,NPL
        
        %for i=1:NP
        %   TdT0_CORE=...
        %       bsxfun(@times,QPROJ_ORI{NT1},reshape(conj(CQFAK_{NT1}),1,LMMAXC(NT1)))...
        %       *CQIJ_DIV{NT1,NT2}(:,:,i)...
        %       *bsxfun(@times,QPROJ_ORI{NT2},reshape(CQFAK_{NT2},1,LMMAXC(NT2))).';
        %   %                             CREXP=conj(CREXP_first(:,i).*CREXP_second(:,i))*CREXP_first(:,i).';
        %   TdT0_2=TdT0_2+...
        %       CREXP_first(:,i)*CREXP_second(:,i).'...
        %       .*TdT0_CORE;
        %end
        %toc;
    end
end

% This is the T^\dagger T_0 under two sets of
% configuration-dependent pseudo plane wave basis
TdT0=TRANS_1'+TRANS_2-eye(NPL,'single','gpuArray')+TdT0;
end