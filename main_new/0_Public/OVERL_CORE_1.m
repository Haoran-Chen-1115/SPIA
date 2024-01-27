function [QPAW3,CQIJ]=OVERL_CORE_1(ROOT_DIR,NTYP,NITYP,LMMAXC,CH0,CH1,CH2,CH3,WAE,WPS,SIMPI)
% i.e. <\tilde{\phi_i^a|\phi_j^a>-<\tilde{\phi}_i^a|\tilde{\phi}_j^a>
QPAW3=cell(1,NTYP);CQIJ=cell(1,NTYP);

for NT=1:NTYP
    NCH=0;NCHM=0;
    QPAW3{NT}=zeros(LMMAXC(NT),'gpuArray');
    for NL=1:CH0{NT}(2)
        NCH=NCH+1;NCHM=NCHM+1;
        QPAW3{NT}(NCHM,1:CH0{NT}(2))=sum((WPS{NT}(:,NCH).*WAE{NT}(:,1:CH0{NT}(2))-...
            WPS{NT}(:,NCH).*WPS{NT}(:,1:CH0{NT}(2))).*SIMPI{NT},1);
    end
    for NL=1:CH1{NT}(2)
        LL=CH0{NT}(2)+1:CH0{NT}(2)+CH1{NT}(2);
        NCH=NCH+1;
        for NM=1:3
            NCHM=NCHM+1;
            QPAW3{NT}(NCHM,NCHM-3*(NL-1):3:NCHM+3*(CH1{NT}(2)-NL))=...
                sum((WPS{NT}(:,NCH).*WAE{NT}(:,LL)-...
                WPS{NT}(:,NCH).*WPS{NT}(:,LL)).*SIMPI{NT},1);
        end
    end
    for NL=1:CH2{NT}(2)
        LL=CH0{NT}(2)+CH1{NT}(2)+1:CH0{NT}(2)+CH1{NT}(2)+CH2{NT}(2);
        NCH=NCH+1;
        for NM=1:5
            NCHM=NCHM+1;
            QPAW3{NT}(NCHM,NCHM-5*(NL-1):5:NCHM+5*(CH2{NT}(2)-NL))=...
                sum((WPS{NT}(:,NCH).*WAE{NT}(:,LL)-...
                WPS{NT}(:,NCH).*WPS{NT}(:,LL)).*SIMPI{NT},1);
        end
    end
    for NL=1:CH3{NT}(2)
        LL=CH0{NT}(2)+CH1{NT}(2)+CH2{NT}(2)+1:CH0{NT}(2)+CH1{NT}(2)+CH2{NT}(2)+CH3{NT}(2);
        NCH=NCH+1;
        for NM=1:7
            NCHM=NCHM+1;
            QPAW3{NT}(NCHM,NCHM-7*(NL-1):7:NCHM+7*(CH3{NT}(2)-NL))=...
                sum((WPS{NT}(:,NCH).*WAE{NT}(:,LL)-...
                WPS{NT}(:,NCH).*WPS{NT}(:,LL)).*SIMPI{NT},1);
        end
    end
    QPAW3{NT}=gpuArray(single(QPAW3{NT}));
    
    CQIJJ=fopen([ROOT_DIR,'/BEAD_1/CQIJ_1_1_',int2str(NT)]);
    CQIJ_ORI=reshape(fread(CQIJJ,...
        [LMMAXC(NT)*LMMAXC(NT),NITYP(NT)],'double'),LMMAXC(NT),LMMAXC(NT),NITYP(NT));
    CQIJ{NT}=gpuArray(CQIJ_ORI(:,:,1));
    fclose(CQIJJ);
end

end
