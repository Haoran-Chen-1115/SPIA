function TdT0=OVERL(POS_1, POS_2, TRANS_1, TRANS_2, A, B, ...
    CQIJ3_fun, r, G_INDEX, CITPI, QPROJ_ORI, ...
    NTYP,NITYP,NPL,LMAX,LMMAXC, ...
    CH0,CH1,CH2,CH3,LPS)

CQIJ_DIV=cell(NTYP,NTYP);
idx_phase=cell(NTYP,NTYP);
TdT0=zeros(NPL,NPL,'single','gpuArray');
%TdT0_2=zeros(NPL,NPL,'single','gpuArray');
CQFAK_=cell(1,NTYP);
for NT=1:NTYP
    % 'Phasefactor' ignored in spherical Bessel function
    CQFAK_{NT}=zeros(1,LMMAXC(NT),'single');
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
    
%     CQFAK_{NT}=single(zeros(1,LMMAXC(NT)));
%     %                     CQFAK=zeros(1,LMMAXC(NT),'gpuArray');
%     CQFAK_{NT}(1:CH0{NT}(2))=1;
%     CQFAK_{NT}(CH0{NT}(2)+1:CH0{NT}(2)+3*CH1{NT}(2))=1i;
%     CQFAK_{NT}(CH0{NT}(2)+3*CH1{NT}(2)+1:CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2))=-1;
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
        %        +triu((r{NT1,NT2}(end)+1)*ones(NITYP(NT1)),1).';
        %end
        %                         DR(DR<1E-8)=1E-10;
        idx=find(DR<r{NT1,NT2}(end));
        idx_phase{NT1,NT2}(1,:)=mod(idx-1,NITYP(NT1))+1;
        idx_phase{NT1,NT2}(2,:)=ceil(idx/NITYP(NT1));
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
