function [FHAM,TRANS,FOVL]=HAMILT_OLD(ROOT_DIR,NC_c,NB_c,ISP_c,...
    NTYP,NITYP,LMMAXC,CH0,CH1,CH2,CH3,...
    QPROJ_ORI,NPL,NPLWV,NGX,NGY,NGZ,...
    TRANS_CORE,TRANS_c_CORE,FOVL_CORE,...
    POSION,G_INDEX,G,HSQDTM,CITPI)
SVV=fopen([ROOT_DIR,'/BEAD_',NC_c,'/SV_',NB_c,'_',ISP_c]);
SV=reshape(fread(SVV,NPLWV,'double'),NGX,NGY,NGZ);
fclose(SVV);

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
FOVL=single(eye(NPL,'gpuArray'));
TRANS=single(eye(NPL,'gpuArray'));
TRANS_c=single(eye(NPL,'gpuArray'));
%                 FOVL=eye(NPL,'gpuArray');
%                 TRANS=eye(NPL,'gpuArray');
%                 TRANS_c=eye(NPL,'gpuArray');

for NT=1:NTYP
    NT_c=int2str(NT);
    
    CDIJJ=fopen([ROOT_DIR,'/BEAD_',NC_c,'/CDIJ_',NB_c,'_',ISP_c,'_',NT_c]);
    % Non local pseudopotential strength and Overlap Strength
    %                     for NT1=1:NT-1 % Because of my miss when implementing vasp output
    %                         fread(CDIJJ,...
    %                             [LMMAXC(NT)*LMMAXC(NT),NITYP(NT1)],'double');
    %                     end
    CDIJ=reshape(fread(CDIJJ,...
        [LMMAXC(NT)*LMMAXC(NT),NITYP(NT)],'double'),LMMAXC(NT),LMMAXC(NT),NITYP(NT));
    
    fclose(CDIJJ);
    % Calculate Non local part of Hamiltonian
    
    % QPROJ_ORI{NT}=gpuArray(QPROJ_ORI{NT});
    %                     CDIJ=gpuArray(CDIJ);
    CDIJ=gpuArray(single(CDIJ));
    
    % Calculate phase
    %                     POSION1=single(gpuArray(POSION{NT}));
    POSION1=gpuArray(POSION{NT});
    CREXP=single(exp(CITPI*G_INDEX(1:3,1:NPL).'*POSION1));
    
    % 'Phasefactor' ignored in spherical Bessel function
    CQFAK=single(zeros(1,LMMAXC(NT),'gpuArray'));
    %                     CQFAK=zeros(1,LMMAXC(NT),'gpuArray');
    CQFAK(1:CH0{NT}(2))=1;
    CQFAK(CH0{NT}(2)+1:CH0{NT}(2)+3*CH1{NT}(2))=1i;
    CQFAK(CH0{NT}(2)+3*CH1{NT}(2)+1:CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2))=-1;
    CQFAK(CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2)+1:...
        CH0{NT}(2)+3*CH1{NT}(2)+5*CH2{NT}(2)+7*CH3{NT}(2))=-1i;
    
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
    TRANS=TRANS+CREXP.*TRANS_CORE{NT};
    % For \hat{T}_c=|\tilde{p}_i^a><\tilde{\phi}_i^a|...
    TRANS_c=TRANS_c+CREXP.*TRANS_c_CORE{NT};
    % Overlap matrix
    FOVL=FOVL+CREXP.*FOVL_CORE{NT};
    
    %TRANS=gather(TRANS);
    clear CREXP
end

% Replace S by S'=S+(\hat{T}-\hat{T}_c+h.c.)
TRANS_c=TRANS-TRANS_c;
FOVL=FOVL+TRANS_c+TRANS_c';
%                 clear TRANS_c

FHAM=triu(FHAM,1)+triu(FHAM,1)'+diag(real(diag(FHAM)));
FOVL=triu(FOVL,1)+triu(FOVL,1)'+diag(real(diag(FOVL)));

end
