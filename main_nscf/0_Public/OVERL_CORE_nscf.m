function [QPROJ_ORI,RES]=...
    OVERL_CORE_nscf(PSMAXN,PSPNL,NTYP,...
    LMAX,LMMAXC,RG,NMAX,...
    NPL,CH0,CH1,CH2,CH3,LPS,...
    WAE_PS,WAE,WPS,G2,OMEGA)

NPSNL = 100;
RES=cell(1,NTYP);
QPROJ_L = cell(1,NTYP); QPROJ_ORI = cell(1,NTYP);
for NT=1:NTYP
    %% Calculate <\psi_\bm{k}|\phi_i^a> around each type of ion
    % PROJJ=fopen([ROOT_DIR,'/BEAD_',NC_c,'/PROJ_1_',NK_c,'_',ISP_c,'_',int2str(NT)]);
    % % LMAX(NT)=fread(PROJJ,1,'int');LMMAXC(NT)=fread(PROJJ,1,'int');
    % % CH0{NT}=fread(PROJJ,2,'int');CH1{NT}=fread(PROJJ,2,'int');
    % % CH2{NT}=fread(PROJJ,2,'int');CH3{NT}=fread(PROJJ,2,'int');
    % fread(PROJJ,10,'int');
    % QPROJ_ORI{NT}=fread(PROJJ,[NPL,LMMAXC(NT)],'double');
    % fclose(PROJJ);
    %% Directly construct QPROJ
    LCH = LPS{NT};
    %     LCH = [zeros(1,CH0{NT}(2)),...
    %         ones(1,CH1{NT}(2)),...
    %         2*ones(1,CH2{NT}(2)),...
    %         3*ones(1,CH3{NT}(2))];
    G_grid = linspace(0,PSMAXN(NT),NPSNL+1);
    G_grid = [-G_grid(2),G_grid(1:end-1)];
    QPROJ_L{NT} = zeros(NPL,LMAX(NT),'gpuArray');
    QPROJ_ORI{NT} = zeros(NPL,LMMAXC(NT),'gpuArray');
    NCH = 0;
    for NL = 1:LMAX(NT)
        %                 PSPNL_fun = @(x) interp1(G_grid, PSPNL{NT}(:,NL), x, 'cubic');
        %                 QPROJ_L{NT}(:,NL) =...
        %                     PSPNL_fun(gather(G2(4,:)'))/sqrt(OMEGA)/2/sqrt(pi);
        idx = ceil(G2(4,:)'*NPSNL/PSMAXN(NT));
        REM = mod(G2(4,:)'*NPSNL/PSMAXN(NT),1);
        V1 = PSPNL{NT}(idx,NL);
        V2 = PSPNL{NT}(idx+1,NL);
        V3 = PSPNL{NT}(idx+2,NL);
        V4 = PSPNL{NT}(idx+3,NL);
        T0 = V2;
        T1=((6*V3)-(2*V1)-(3*V2)-V4)/6;
        T2=(V1+V3-(2*V2))/2;
        T3=(V4-V1+(3*(V2-V3)))/6;
        QPROJ_L{NT}(:,NL) =...
            (T0+REM.*(T1+REM.*(T2+REM.*T3)))/sqrt(OMEGA);
        
        if LCH(NL)==0
            NCH = NCH + 1;
            QPROJ_ORI{NT}(:,NCH)=...
                QPROJ_L{NT}(:,NL)/2/sqrt(pi);
        elseif LCH(NL)==1
            NCH = NCH + 3;
            FAK2 = sqrt(3/pi)/2;
            Y1=FAK2*[G2(2,:);G2(3,:);G2(1,:)].';
            QPROJ_ORI{NT}(:,NCH-2:NCH)=...
                QPROJ_L{NT}(:,NL).*Y1;
        elseif LCH(NL)==2
            NCH = NCH + 5;
            FAK3 = sqrt(15/pi)/4;
            Y2=FAK3*[2*G2(1,:).*G2(2,:);...
                2*G2(3,:).*G2(2,:);...
                (3*G2(3,:).^2-1)/sqrt(3);...
                2*G2(3,:).*G2(1,:);...
                (G2(1,:).^2-G2(2,:).^2)].';
            QPROJ_ORI{NT}(:,NCH-4:NCH)=...
                QPROJ_L{NT}(:,NL).*Y2;
        elseif LCH(NL)==3
            NCH = NCH + 7;
            FAK3 = [sqrt(35/2/pi)/4;sqrt(105/pi)/2;sqrt(21/2/pi)/4;...
                sqrt(7/pi)/4;sqrt(21/2/pi)/4;sqrt(105/pi)/4;sqrt(35/2/pi)/4];
            Y3=(FAK3.*[(3*G2(1,:).^2-G2(2,:).^2).*G2(2,:);...
                G2(1,:).*G2(2,:).*G2(3,:);...
                G2(2,:).*(5*G2(3,:).^2-1);...
                G2(3,:).*(5*G2(3,:).^2-3);...
                G2(1,:).*(5*G2(3,:).^2-1);...
                (G2(1,:).^2-G2(2,:).^2).*G2(3,:);...
                (G2(1,:).^2-3*G2(2,:).^2).*G2(1,:)]).';
            QPROJ_ORI{NT}(:,NCH-6:NCH)=...
                QPROJ_L{NT}(:,NL).*Y3;
        end
        
    end
    
    NR=find(WAE_PS{NT}(:,1)==0,1,'first');
    GR=bsxfun(@times,reshape(G2(4,:),NPL,1),reshape(RG{NT},1,NMAX{NT}));
    Wi=cell(1,LMAX(NT));
    for NL=1:LMAX(NT)
        Wi{NL}=@(x) interp1(gather(RG{NT}),gather(WAE_PS{NT}(:,NL).'./x),x,'spline').*x;
    end
    
    WAEi=cell(1,LMAX(NT));WPSi=cell(1,LMAX(NT));
    for NL=1:LMAX(NT)
        WAEi{NL}=@(x) interp1(RG{NT},WAE{NT}(:,NL).',x,'spline').*x;
        WPSi{NL}=@(x) interp1(RG{NT},WPS{NT}(:,NL).',x,'spline').*x;
    end
    
    % Grids and spherical Bessel functions
    xR=@(x) bsxfun(@times,reshape(G2(4,:),[],1),reshape(x,1,[]));
    j0i=@(x) (G2(4,:)<1E-10).' + sin(xR(x))...
        ./xR(x).*(G2(4,:)>1E-10).';
    j1i=@(x) (sin(xR(x))./(xR(x).^2)...
        -cos(xR(x))./xR(x))...
        .*(G2(4,:)>1E-10).';
    j2i=@(x) (-j0i(x)+3*j1i(x)./xR(x))...
        .*(G2(4,:)>1E-10).';
    j3i=@(x) (-j1i(x)+5*j2i(x)./xR(x))...
        .*(G2(4,:)>1E-10).';
    
    RES{NT}=zeros(NPL,LMMAXC(NT),'gpuArray');
    NCH=0;NCHM=0;
    for NL=1:LMAX(NT)
        % L=0
        %for NL=1:CH0{NT}(2)
        if LCH(NL)==0
            NCH=NCH+1;NCHM=NCHM+1;
            
            FAK1 = 4*pi/(2*sqrt(pi));
            
            j=@(x) gather(j0i(x).*Wi{NCH}(x));
            RES{NT}(:,NCHM)=FAK1*integral(j,...
                gather(min(RG{NT})),gather(RG{NT}(NR)),'ArrayValued',true,'RelTol',1E-12);
            %end
            
            % L=1
            %if CH1{NT}(2)~=0
        elseif LCH(NL)==1
            FAK2 = 2*pi*sqrt(3/pi);
            Y1=FAK2*[G2(2,:);G2(3,:);G2(1,:)].';
            
            %tic;
            %for NL=1:CH1{NT}(2)
            NCH=NCH+1;
            j=@(x) gather(j1i(x).*Wi{NCH}(x));
            for NM=1:3
                NCHM=NCHM+1;
                RES{NT}(:,NCHM)=...
                    gather(Y1(:,NM)).*integral(j,...
                    gather(min(RG{NT})),gather(RG{NT}(NR)),'ArrayValued',true,'RelTol',1E-12);
            end
            %end
            %toc;
            %end
            
            % L=2
            %if CH2{NT}(2)~=0
        elseif LCH(NL)==2
            FAK3 = pi*sqrt(15/pi);
            Y2=FAK3*[2*G2(1,:).*G2(2,:);...
                2*G2(3,:).*G2(2,:);...
                (3*G2(3,:).^2-1)/sqrt(3);...
                2*G2(3,:).*G2(1,:);...
                (G2(1,:).^2-G2(2,:).^2)].';
            
            %for NL=1:CH2{NT}(2)
            NCH=NCH+1;
            j=@(x) gather(j2i(x).*Wi{NCH}(x));
            for NM=1:5
                NCHM=NCHM+1;
                RES{NT}(:,NCHM)=...
                    gather(Y2(:,NM)).*integral(j,...
                    gather(min(RG{NT})),gather(RG{NT}(NR)),'ArrayValued',true,'RelTol',1E-12);
            end
            %end
            %end
            
            % L=3
            %if CH3{NT}(2)~=0
        elseif LCH(NL)==3
            FAK3 = [pi*sqrt(35/2/pi);2*pi*sqrt(105/pi);pi*sqrt(21/2/pi);...
                pi*sqrt(7/pi);pi*sqrt(21/2/pi);pi*sqrt(105/pi);pi*sqrt(35/2/pi)];
            Y3=(FAK3.*[(3*G2(1,:).^2-G2(2,:).^2).*G2(2,:);...
                G2(1,:).*G2(2,:).*G2(3,:);...
                G2(2,:).*(5*G2(3,:).^2-1);...
                G2(3,:).*(5*G2(3,:).^2-3);...
                G2(1,:).*(5*G2(3,:).^2-1);...
                (G2(1,:).^2-G2(2,:).^2).*G2(3,:);...
                (G2(1,:).^2-3*G2(2,:).^2).*G2(1,:)]).';
            
            %for NL=1:CH3{NT}(2)
            NCH=NCH+1;
            j=@(x) gather(j3i(x).*Wi{NCH}(x));
            for NM=1:7
                NCHM=NCHM+1;
                RES{NT}(:,NCHM)=...
                    gather(Y3(:,NM)).*integral(j,...
                    gather(min(RG{NT})),gather(RG{NT}(NR)),'ArrayValued',true,'RelTol',1E-12);
            end
            %end
        end
    end
end

for NT=1:NTYP
    RES{NT}=gpuArray(single(RES{NT}));
    QPROJ_ORI{NT}=gpuArray(single(QPROJ_ORI{NT}));
    % RES{NT}=gpuArray(RES{NT});
    % QPROJ_ORI{NT}=gpuArray(QPROJ_ORI{NT});
end

%TRANS_CORE=cell(1,NTYP);
%for NT=1:NTYP
%    TRANS_CORE{NT}=gather((RES{NT}*QPROJ_ORI{NT}.')/sqrt(OMEGA));
%end

end
