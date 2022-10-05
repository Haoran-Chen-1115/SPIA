CQIJ3_div=cell(NTYP,NTYP);

NR=zeros(1,NTYP);Wi=cell(1,NTYP);
for NT=1:NTYP
    NR(NT)=gather(find(WAE_PS{NT}(:,1)==0,1,'first'));
    %     GR=bsxfun(@times,reshape(G2_AE(4,:),NPL_AE,1),reshape(RG{NT},1,NMAX{NT}));
    Wi{NT}=cell(1,LMAX(NT));
    for NL=1:LMAX(NT)
        Wi{NT}{NL}=@(x)...
            interp1(gather(RG{NT}),gather(WAE_PS{NT}(:,NL)./(RG{NT}.^2)),x,'spline')...
            .*(x<RG{NT}(NR(NT)));
    end
end

D_log=zeros(NTYP,NTYP);
r=cell(NTYP,NTYP);
for NT1=1:NTYP
    for NT2=NT1:NTYP
        %     r=gather(linspace(0,2*RG{NT}(NR),100));
%         D_log(NT1,NT2)=1+2*(min(gather(RG{NT1}(2)/RG{NT1}(1)),...
%             gather(RG{NT2}(2)/RG{NT2}(1)))...
%             -1);
%         r{NT1,NT2}=[0,(gather(RG{NT1}(NR(NT1)))...
%             +gather(RG{NT2}(NR(NT2))))*exp((-120:0)*log(D_log(NT1,NT2)))];
        r1=linspace(0,gather(RG{NT1}(NR(NT1))),80);
        r2=linspace(gather(RG{NT1}(NR(NT1))),...
            gather(RG{NT1}(NR(NT1)))+gather(RG{NT2}(NR(NT2))),41);
        r{NT1,NT2}=[r1,r2(2:end)];
        %if NT1==NT2
        %    r{NT1,NT2}=[r1,r2(2:end)];
        %else
        %    r{NT1,NT2} = r2(22:end);
        %end
    end
end
%%
time=cell(NTYP,NTYP);
for NT1=1:NTYP
    for NT2=NT1:NTYP
        t2=tic;
        time{NT1,NT2}=zeros(LMMAXC(NT1),length(r{NT1,NT2}));
        CQIJ3_div{NT1,NT2}=zeros(LMMAXC(NT1),LMMAXC(NT2),length(r{NT1,NT2}));
        CQIJ3_tmp1=zeros(LMMAXC(NT1),LMMAXC(NT2),length(r{NT1,NT2}));
        CQIJ3_tmp2=zeros(LMMAXC(NT1),LMMAXC(NT2),length(r{NT1,NT2}));
        Nr_inter=120;
        %if NT1 == NT2
        %    Nr_inter=120;
        %else
        %    Nr_inter=20;
        %end
        
        parfor Nr=1:Nr_inter
            t3=tic;
            R=[0,0,r{NT1,NT2}(Nr)];
            R(4)=sqrt(sum(R.^2,2));
            % Integral interval
%             rmin=max(0,gather(R(4)-RG{NT2}(NR(NT2))));
%             r_min_fun =@(phi,the) max(0,gather(2*R(4)*sin(the)...
%                 -sqrt(4*R(4)^2*sin(the).^2 ...
%                 +RG{NT2}(NR(NT2))^2-R(4)^2)));
%             rmax=gather(RG{NT1}(NR(NT1)));
%             sin_the1=max(gather((RG{NT1}(NR(NT1))^2+R(4)^2-RG{NT2}(NR(NT2))^2)...
%                 /(2*RG{NT1}(NR(NT1))*R(4))),-1);
%             the1=max(0,real(asin(sin_the1)));
%             the_min=the1; the_max=pi-the1;
%             phi_min=0;
%             phi_max=2*pi;
            
            if R(4)<=RG{NT2}(NR(NT2))
                r_min_fun = 0;
                rmax = gather(RG{NT1}(NR(NT1)));
                r2 = gather(RG{NT2}(NR(NT2)));
                %r_max_fun = @(phi,the) min(rmax,...
                %    R(4)*cos(the) + sqrt(r2^2-R(4)^2*sin(the).^2));
                r_max_fun = rmax;
                phi_min = 0;
                phi_max = 2*pi;
                the_min = 0;
                the_max = pi;
            else
                r2 = gather(RG{NT2}(NR(NT2)));
                rmax = gather(RG{NT1}(NR(NT1)));
                %r_min_fun = @(phi,the) min(rmax,...
                %    max(R(4)*cos(the) - sqrt(r2^2-R(4)^2*sin(the).^2),0));
                r_min_fun = 0;
                r_max_fun = rmax;
                phi_min = 0;
                phi_max = 2*pi;
                the_min = 0;
                the_max = pi;
                %the_max = asin(min(r2/R(4),1));
            end
                
            %% Define the integrand at first
            RR=@(r,the,phi)...
                sqrt((r.*sin(the).*cos(phi)-R(1)).^2 ...
                +(r.*sin(the).*sin(phi)-R(2)).^2 ...
                +(r.*cos(the)-R(3)).^2);
            XX=@(r,the,phi) r.*sin(the).*cos(phi)-R(1);
            YY=@(r,the,phi) r.*sin(the).*sin(phi)-R(2);
            ZZ=@(r,the,phi) r.*cos(the)-R(3);
            
            Wi_fun=cell(1,LMAX(NT1));
            NCHM1=0;
            for NCH1 = 1:LMAX(NT1)
                
                if LPS{NT1}(NCH1) == 0
                    % L=0
                    FAK0 = 1/(2*sqrt(pi));
                    for NM=1:1
                        NCHM1=NCHM1+1;
                        Wi_fun{NCHM1}=@(r,the,phi)...
                            FAK0*Wi{NT1}{NCH1}(r);
                    end
                    
                elseif LPS{NT1}(NCH1) == 1
                    % L=1
                    FAK1 = 1/2*sqrt(3/pi);
                    Y1=cell(1,3);
                    
                    Y1{1} = @(the,phi) FAK1*sin(the).*sin(phi);
                    Y1{2} = @(the,phi) FAK1*cos(the);
                    Y1{3} = @(the,phi) FAK1*sin(the).*cos(phi);
                    for NM=1:3
                        NCHM1=NCHM1+1;
                        Wi_fun{NCHM1}=@(r,the,phi)...
                            Y1{NM}(the,phi).*Wi{NT1}{NCH1}(r);
                    end
                elseif LPS{NT1}(NCH1) == 2
                    % L=2
                    FAK2 = 1/4*sqrt(15/pi);
                    Y2=cell(1,5);
                    
                    Y2{1}=@(the,phi) FAK2*sin(the).^2.*sin(2*phi);
                    Y2{2}=@(the,phi) FAK2*2*sin(the).*cos(the).*sin(phi);
                    Y2{3}=@(the,phi) FAK2/sqrt(3)*(3*cos(the).^2-1);
                    Y2{4}=@(the,phi) FAK2*2*sin(the).*cos(the).*cos(phi);
                    Y2{5}=@(the,phi) FAK2*sin(the).^2.*cos(2*phi);
                    
                    for NM=1:5
                        NCHM1=NCHM1+1;
                        Wi_fun{NCHM1}=@(r,the,phi)...
                            Y2{NM}(the,phi).*Wi{NT1}{NCH1}(r);
                    end
                elseif LPS{NT1}(NCH1)==3
                    % L=3
                    FAK3 = [1/4*sqrt(35/2/pi); 1/4*sqrt(105/pi);...
                        1/4*sqrt(21/2/pi); 1/4*sqrt(7/pi);...
                        1/4*sqrt(21/2/pi); 1/4*sqrt(105/pi);...
                        1/4*sqrt(35/2/pi)];
                    
                    Y3=cell(1,7);
                    
                    Y3{1}=@(the,phi) FAK3(1)*sin(the).^3.*sin(3*phi);
                    Y3{2}=@(the,phi) FAK3(2)*sin(the).^2.*cos(the).*sin(2*phi);
                    Y3{3}=@(the,phi) FAK3(3)*sin(the).*(5*cos(the).^2-1).*sin(phi);
                    Y3{4}=@(the,phi) FAK3(4)*(5*cos(the).^3-3*cos(the));
                    Y3{5}=@(the,phi) FAK3(5)*sin(the).*(5*cos(the).^2-1).*cos(phi);
                    Y3{6}=@(the,phi) FAK3(6)*sin(the).^2.*cos(the).*cos(2*phi);
                    Y3{7}=@(the,phi) FAK3(7)*sin(the).^3.*cos(3*phi);
                    
                    for NM=1:7
                        NCHM1=NCHM1+1;
                        Wi_fun{NCHM1}=@(r,the,phi)...
                            Y3{NM}(the,phi).*Wi{NT1}{NCH1}(r);
                    end
                end
            end
            
            Wi_fun_div=cell(1,LMAX(NT2));
            NCHM2=0;
            for NCH2 = 1:LMAX(NT2)
                
                if LPS{NT2}(NCH2) == 0
                    % L=0
                    FAK0 = 1/(2*sqrt(pi));
                    for NM=1:1
                        NCHM2=NCHM2+1;
                        Wi_fun_div{NCHM2}=@(r,the,phi)...
                            FAK0*Wi{NT2}{NCH2}(RR(r,the,phi));
                    end
                elseif LPS{NT2}(NCH2)==1
                    % L=1
                    FAK1 = 1/2*sqrt(3/pi);
                    Y1_div=cell(1,3);
                    Y1_div{1} = @(r,the,phi) FAK1*...
                        YY(r,the,phi)...
                        ./RR(r,the,phi);
                    Y1_div{2} = @(r,the,phi) FAK1*...
                        ZZ(r,the,phi)...
                        ./RR(r,the,phi);
                    Y1_div{3} = @(r,the,phi) FAK1*...
                        XX(r,the,phi)...
                        ./RR(r,the,phi);
                    
                    for NM=1:3
                        NCHM2=NCHM2+1;
                        Wi_fun_div{NCHM2}=@(r,the,phi)...
                            Y1_div{NM}(r,the,phi).*Wi{NT2}{NCH2}(RR(r,the,phi));
                    end
                elseif LPS{NT2}(NCH2)==2
                    % L=2
                    FAK2 = 1/4*sqrt(15/pi);
                    Y2_div=cell(1,5);
                    Y2_div{1}=@(r,the,phi) FAK2...
                        *2*XX(r,the,phi).*YY(r,the,phi)...
                        ./(RR(r,the,phi).^2);
                    Y2_div{2}=@(r,the,phi) FAK2...
                        *2*ZZ(r,the,phi).*YY(r,the,phi)...
                        ./(RR(r,the,phi).^2);
                    Y2_div{3}=@(r,the,phi) FAK2/sqrt(3)...
                        *(3*ZZ(r,the,phi).^2 ...
                        ./(RR(r,the,phi).^2)-1);
                    Y2_div{4}=@(r,the,phi) FAK2...
                        *2*XX(r,the,phi).*ZZ(r,the,phi)...
                        ./(RR(r,the,phi).^2);
                    Y2_div{5}=@(r,the,phi) FAK2...
                        *(XX(r,the,phi).^2 - YY(r,the,phi).^2) ...
                        ./(RR(r,the,phi).^2);
                    
                    for NM=1:5
                        NCHM2=NCHM2+1;
                        Wi_fun_div{NCHM2}=@(r,the,phi)...
                            Y2_div{NM}(r,the,phi).*Wi{NT2}{NCH2}(RR(r,the,phi));
                    end
                elseif LPS{NT2}(NCH2)==3
                    % L=3
                    FAK3 = [1/4*sqrt(35/2/pi); 1/4*sqrt(105/pi);...
                        1/4*sqrt(21/2/pi); 1/4*sqrt(7/pi);...
                        1/4*sqrt(21/2/pi); 1/4*sqrt(105/pi);...
                        1/4*sqrt(35/2/pi)];
                    
                    Y3_div=cell(1,7);
                    
                    Y3_div{1}=@(r,the,phi) FAK3(1)...
                        *(3*XX(r,the,phi).^2 -YY(r,the,phi).^2).*YY(r,the,phi) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{2}=@(r,the,phi) FAK3(2)...
                        *2*XX(r,the,phi).*YY(r,the,phi).*ZZ(r,the,phi) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{3}=@(r,the,phi) FAK3(3)...
                        *YY(r,the,phi).*(5*ZZ(r,the,phi).^2-RR(r,the,phi).^2) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{4}=@(r,the,phi) FAK3(4)...
                        *ZZ(r,the,phi).*(5*ZZ(r,the,phi).^2-3*RR(r,the,phi).^2) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{5}=@(r,the,phi) FAK3(5)...
                        *XX(r,the,phi).*(5*ZZ(r,the,phi).^2-RR(r,the,phi).^2) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{6}=@(r,the,phi) FAK3(6)...
                        *(XX(r,the,phi).^2-YY(r,the,phi).^2).*ZZ(r,the,phi) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{7}=@(r,the,phi) FAK3(7)...
                        *(XX(r,the,phi).^2 -3*YY(r,the,phi).^2).*XX(r,the,phi) ...
                        ./(RR(r,the,phi).^3);
                    
                    for NM=1:7
                        NCHM2=NCHM2+1;
                        Wi_fun_div{NCHM2}=@(r,the,phi)...
                            Y3_div{NM}(r,the,phi).*Wi{NT2}{NCH2}(RR(r,the,phi));
                    end
                end
            end
            
            CQIJ3_tmp=zeros(LMMAXC(NT1),LMMAXC(NT2));
            if NT1==NT2
                for NL1=1:LMMAXC(NT1)
                    t1=tic;
                    for NL2=1:LMMAXC(NT2)
                        Integrand2=@(phi,the,r)...
                            gather(Wi_fun_div{NL2}(r,the,phi)...
                            .*Wi_fun{NL1}(r,the,phi).*(r.^2).*sin(the));                        
                        % Convergency is greatly improved by determing the
                        % integral interval first
                        CQIJ3_tmp(NL1,NL2)=...
                            integral3(Integrand2,phi_min,phi_max...
                            ,the_min,the_max,r_min_fun,r_max_fun,...
                            'AbsTol',1E-10,'RelTol',1E-6); 
                    end
                    %time{NT1,NT2}(NL1,Nr)=toc(t1);
                    %disp([NT1,NT2,Nr,r{NT1,NT2}(Nr),NL1,toc(t1)])
                end
            else
                for NL1=1:LMMAXC(NT1)
                    t1=tic;
                    for NL2=1:LMMAXC(NT2)
                        % t2=tic;
                        Integrand2=@(phi,the,r)...
                            gather(Wi_fun_div{NL2}(r,the,phi)...
                            .*Wi_fun{NL1}(r,the,phi).*(r.^2).*sin(the));
                        CQIJ3_tmp(NL1,NL2)=...
                            integral3(Integrand2,phi_min,phi_max...
                            ,the_min,the_max,r_min_fun,r_max_fun...
                            ,'AbsTol',1E-10,'RelTol',1E-6); 
                    end
                    %time{NT1,NT2}(NL1,Nr)=toc(t1);
                    %disp([NT1,NT2,Nr,r{NT1,NT2}(Nr),NL1,toc(t1)])
                end
            end
            CQIJ3_tmp1(:,:,Nr)=CQIJ3_tmp;
            disp(['Loop time for',...
                '(',int2str(NT1),',',int2str(NT2),')',...
                'at (',int2str(Nr),',',num2str(r{NT1,NT2}(Nr)),'A ) is',...
                num2str(toc(t3)),'s'])
            
        end

        % For distant overlaps, use 'iterated' methods.
        parfor Nr=Nr_inter+1:length(r{NT1,NT2})
            t3=tic;
            R=[0,r{NT1,NT2}(Nr),0];
            R(4)=sqrt(sum(R.^2,2));
            % Integral interval
            
            if R(4)<=RG{NT2}(NR(NT2))
                r_min_fun = 0;
                rmax = gather(RG{NT1}(NR(NT1)));
                r_max_fun = @(phi,the) min(RG{NT1}(NR(NT1)),...
                    R(4)*cos(the) + sqrt(RG{NT2}(NR(NT2))^2-R(4)^2*sin(the).^2));
                phi_min = 0;
                phi_max = 2*pi;
                the_min = 0;
                the_max = pi;
            else
                r_min_fun = @(phi,the) @(phi,the) min(RG{NT1}(NR(NT1)),...
                    R(4)*cos(the) - sqrt(RG{NT2}(NR(NT2))^2-R(4)^2*sin(the).^2));
                rmax = gather(RG{NT1}(NR(NT1)));
                r_max_fun = RG{NT1}(NR(NT1));
                phi_min = 0;
                phi_max = 2*pi;
                the_min = 0;
                the_max = asin(min(RG{NT2}(NR(NT2))/R(4),1));
            end
%             rmin=max(0,gather(R(4)-RG{NT2}(NR(NT2))));
%             r_min_fun =@(phi,the) max(0,gather(2*R(4)*sin(the)...
%                 -sqrt(4*R(4)^2*sin(the).^2 ...
%                 +RG{NT2}(NR(NT2))^2-R(4)^2)));
%             rmax=gather(RG{NT1}(NR(NT1)));
%             sin_the1=max(gather((RG{NT1}(NR(NT1))^2+R(4)^2-RG{NT2}(NR(NT2))^2)...
%                 /(2*RG{NT1}(NR(NT1))*R(4))),-1);
%             the1=max(0,real(asin(sin_the1)));
%             the_min=the1; the_max=pi-the1;
%             phi_min=real(asin(sin_the1));
%             phi_max=pi-real(asin(sin_the1));
            %% Define the integrand at first
            RR=@(r,the,phi)...
                sqrt((r.*sin(the).*cos(phi)-R(1)).^2 ...
                +(r.*sin(the).*sin(phi)-R(2)).^2 ...
                +(r.*cos(the)-R(3)).^2);
            XX=@(r,the,phi) r.*sin(the).*cos(phi)-R(1);
            YY=@(r,the,phi) r.*sin(the).*sin(phi)-R(2);
            ZZ=@(r,the,phi) r.*cos(the)-R(3);
            
            Wi_fun=cell(1,LMAX(NT1));
            NCHM1=0;
            for NCH1 = 1:LMAX(NT1)
                
                if LPS{NT1}(NCH1) == 0
                    % L=0
                    FAK0 = 1/(2*sqrt(pi));
                    for NM=1:1
                        NCHM1=NCHM1+1;
                        Wi_fun{NCHM1}=@(r,the,phi)...
                            FAK0*Wi{NT1}{NCH1}(r);
                    end
                    
                elseif LPS{NT1}(NCH1) == 1
                    % L=1
                    FAK1 = 1/2*sqrt(3/pi);
                    Y1=cell(1,3);
                    
                    Y1{1} = @(the,phi) FAK1*sin(the).*sin(phi);
                    Y1{2} = @(the,phi) FAK1*cos(the);
                    Y1{3} = @(the,phi) FAK1*sin(the).*cos(phi);
                    for NM=1:3
                        NCHM1=NCHM1+1;
                        Wi_fun{NCHM1}=@(r,the,phi)...
                            Y1{NM}(the,phi).*Wi{NT1}{NCH1}(r);
                    end
                elseif LPS{NT1}(NCH1) == 2
                    % L=2
                    FAK2 = 1/4*sqrt(15/pi);
                    Y2=cell(1,5);
                    
                    Y2{1}=@(the,phi) FAK2*sin(the).^2.*sin(2*phi);
                    Y2{2}=@(the,phi) FAK2*2*sin(the).*cos(the).*sin(phi);
                    Y2{3}=@(the,phi) FAK2/sqrt(3)*(3*cos(the).^2-1);
                    Y2{4}=@(the,phi) FAK2*2*sin(the).*cos(the).*cos(phi);
                    Y2{5}=@(the,phi) FAK2*sin(the).^2.*cos(2*phi);
                    
                    for NM=1:5
                        NCHM1=NCHM1+1;
                        Wi_fun{NCHM1}=@(r,the,phi)...
                            Y2{NM}(the,phi).*Wi{NT1}{NCH1}(r);
                    end
                elseif LPS{NT1}(NCH1)==3
                    % L=3
                    FAK3 = [1/4*sqrt(35/2/pi); 1/4*sqrt(105/pi);...
                        1/4*sqrt(21/2/pi); 1/4*sqrt(7/pi);...
                        1/4*sqrt(21/2/pi); 1/4*sqrt(105/pi);...
                        1/4*sqrt(35/2/pi)];
                    
                    Y3=cell(1,7);
                    
                    Y3{1}=@(the,phi) FAK3(1)*sin(the).^3.*sin(3*phi);
                    Y3{2}=@(the,phi) FAK3(2)*sin(the).^2.*cos(the).*sin(2*phi);
                    Y3{3}=@(the,phi) FAK3(3)*sin(the).*(5*cos(the).^2-1).*sin(phi);
                    Y3{4}=@(the,phi) FAK3(4)*(5*cos(the).^3-3*cos(the));
                    Y3{5}=@(the,phi) FAK3(5)*sin(the).*(5*cos(the).^2-1).*cos(phi);
                    Y3{6}=@(the,phi) FAK3(6)*sin(the).^2.*cos(the).*cos(2*phi);
                    Y3{7}=@(the,phi) FAK3(7)*sin(the).^3.*cos(3*phi);
                    
                    for NM=1:7
                        NCHM1=NCHM1+1;
                        Wi_fun{NCHM1}=@(r,the,phi)...
                            Y3{NM}(the,phi).*Wi{NT1}{NCH1}(r);
                    end
                end
            end
            
            Wi_fun_div=cell(1,LMAX(NT2));
            NCHM2=0;
            for NCH2 = 1:LMAX(NT2)
                
                if LPS{NT2}(NCH2) == 0
                    % L=0
                    FAK0 = 1/(2*sqrt(pi));
                    for NM=1:1
                        NCHM2=NCHM2+1;
                        Wi_fun_div{NCHM2}=@(r,the,phi)...
                            FAK0*Wi{NT2}{NCH2}(RR(r,the,phi));
                    end
                elseif LPS{NT2}(NCH2)==1
                    % L=1
                    FAK1 = 1/2*sqrt(3/pi);
                    Y1_div=cell(1,3);
                    Y1_div{1} = @(r,the,phi) FAK1*...
                        YY(r,the,phi)...
                        ./RR(r,the,phi);
                    Y1_div{2} = @(r,the,phi) FAK1*...
                        ZZ(r,the,phi)...
                        ./RR(r,the,phi);
                    Y1_div{3} = @(r,the,phi) FAK1*...
                        XX(r,the,phi)...
                        ./RR(r,the,phi);
                    
                    for NM=1:3
                        NCHM2=NCHM2+1;
                        Wi_fun_div{NCHM2}=@(r,the,phi)...
                            Y1_div{NM}(r,the,phi).*Wi{NT2}{NCH2}(RR(r,the,phi));
                    end
                elseif LPS{NT2}(NCH2)==2
                    % L=2
                    FAK2 = 1/4*sqrt(15/pi);
                    Y2_div=cell(1,5);
                    Y2_div{1}=@(r,the,phi) FAK2...
                        *2*XX(r,the,phi).*YY(r,the,phi)...
                        ./(RR(r,the,phi).^2);
                    Y2_div{2}=@(r,the,phi) FAK2...
                        *2*ZZ(r,the,phi).*YY(r,the,phi)...
                        ./(RR(r,the,phi).^2);
                    Y2_div{3}=@(r,the,phi) FAK2/sqrt(3)...
                        *(3*ZZ(r,the,phi).^2 ...
                        ./(RR(r,the,phi).^2)-1);
                    Y2_div{4}=@(r,the,phi) FAK2...
                        *2*XX(r,the,phi).*ZZ(r,the,phi)...
                        ./(RR(r,the,phi).^2);
                    Y2_div{5}=@(r,the,phi) FAK2...
                        *(XX(r,the,phi).^2 - YY(r,the,phi).^2) ...
                        ./(RR(r,the,phi).^2);
                    
                    for NM=1:5
                        NCHM2=NCHM2+1;
                        Wi_fun_div{NCHM2}=@(r,the,phi)...
                            Y2_div{NM}(r,the,phi).*Wi{NT2}{NCH2}(RR(r,the,phi));
                    end
                elseif LPS{NT2}(NCH2)==3
                    % L=3
                    FAK3 = [1/4*sqrt(35/2/pi); 1/4*sqrt(105/pi);...
                        1/4*sqrt(21/2/pi); 1/4*sqrt(7/pi);...
                        1/4*sqrt(21/2/pi); 1/4*sqrt(105/pi);...
                        1/4*sqrt(35/2/pi)];
                    
                    Y3_div=cell(1,7);
                    
                    Y3_div{1}=@(r,the,phi) FAK3(1)...
                        *(3*XX(r,the,phi).^2 -YY(r,the,phi).^2).*YY(r,the,phi) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{2}=@(r,the,phi) FAK3(2)...
                        *2*XX(r,the,phi).*YY(r,the,phi).*ZZ(r,the,phi) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{3}=@(r,the,phi) FAK3(3)...
                        *YY(r,the,phi).*(5*ZZ(r,the,phi).^2-RR(r,the,phi).^2) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{4}=@(r,the,phi) FAK3(4)...
                        *ZZ(r,the,phi).*(5*ZZ(r,the,phi).^2-3*RR(r,the,phi).^2) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{5}=@(r,the,phi) FAK3(5)...
                        *XX(r,the,phi).*(5*ZZ(r,the,phi).^2-RR(r,the,phi).^2) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{6}=@(r,the,phi) FAK3(6)...
                        *(XX(r,the,phi).^2-YY(r,the,phi).^2).*ZZ(r,the,phi) ...
                        ./(RR(r,the,phi).^3);
                    Y3_div{7}=@(r,the,phi) FAK3(7)...
                        *(XX(r,the,phi).^2 -3*YY(r,the,phi).^2).*XX(r,the,phi) ...
                        ./(RR(r,the,phi).^3);
                    
                    for NM=1:7
                        NCHM2=NCHM2+1;
                        Wi_fun_div{NCHM2}=@(r,the,phi)...
                            Y3_div{NM}(r,the,phi).*Wi{NT2}{NCH2}(RR(r,the,phi));
                    end
                end
            end
            CQIJ3_tmp=zeros(LMMAXC(NT1),LMMAXC(NT2));
            if NT1==NT2
                for NL1=1:LMMAXC(NT1)
                    t1=tic;
                    for NL2=1:LMMAXC(NT2)
                        Integrand2=@(phi,the,r)...
                            gather(Wi_fun_div{NL2}(r,the,phi)...
                            .*Wi_fun{NL1}(r,the,phi).*(r.^2).*sin(the));                        
                        % Convergency is greatly improved by determing the
                        % integral interval first
                        CQIJ3_tmp(NL1,NL2)=...
                            integral3(Integrand2,phi_min,phi_max...
                            ,the_min,the_max,r_min_fun,r_max_fun...
                            ,'AbsTol',1E-9,'RelTol',1E-5,'Method','iterated'); 
                    end
                    %time{NT1,NT2}(NL1,Nr)=toc(t1);
                    %disp([NT1,NT2,Nr,r{NT1,NT2}(Nr),NL1,toc(t1)])
                end
            else
                for NL1=1:LMMAXC(NT1)
                    t1=tic;
                    for NL2=1:LMMAXC(NT2)
                        t2=tic;
                        Integrand2=@(phi,the,r)...
                            gather(Wi_fun_div{NL2}(r,the,phi)...
                            .*Wi_fun{NL1}(r,the,phi).*(r.^2).*sin(the));
                        CQIJ3_tmp(NL1,NL2)=...
                            integral3(Integrand2,phi_min,phi_max...
                            ,the_min,the_max,r_min_fun,r_max_fun...
                            ,'AbsTol',1E-9,'RelTol',1E-5,'Method','iterated'); 
                    end
                    %time{NT1,NT2}(NL1,Nr)=toc(t1);
                    %disp([NT1,NT2,Nr,r{NT1,NT2}(Nr),NL1,toc(t1)])
                end
            end
            CQIJ3_tmp2(:,:,Nr)=CQIJ3_tmp;
            disp(['Loop time for',...
                '(',int2str(NT1),',',int2str(NT2),')',...
                'at (',int2str(Nr),',',num2str(r{NT1,NT2}(Nr)),'A ) is',...
                num2str(toc(t3)),'s'])
            
        end
        
        CQIJ3_div{NT1,NT2}=CQIJ3_tmp1+CQIJ3_tmp2;
        disp(['Loop time for',...
            '(',int2str(NT1),',',int2str(NT2),')',...
            'is',...
            num2str(toc(t2)),'s'])
    end
end

for NT=1:NTYP
    for Nr=1:length(r{NT,NT})
        CQIJ3_div{NT,NT}(:,:,Nr)=sign(CQIJ3_div{NT,NT}(:,:,Nr)).*...
            abs(triu(CQIJ3_div{NT,NT}(:,:,Nr),1).'+triu(CQIJ3_div{NT,NT}(:,:,Nr),0));
    end 
end

for NT1=1:NTYP
    for NT2 = NT1:NTYP
        for Nr=1:length(r{NT1,NT2})
            R = [0,0,r{NT1,NT2}(Nr),r{NT1,NT2}(Nr)];
            Rotate1=Rotate_SH(R,LMAX,LMMAXC,LPS,NT1);
            Rotate2=Rotate_SH(R,LMAX,LMMAXC,LPS,NT2);
            CQIJ3_div{NT1,NT2}(:,:,Nr)=...
                Rotate1.'*CQIJ3_div{NT1,NT2}(:,:,Nr)*Rotate2;
        end
    end
end


save('../CQIJ3','CQIJ3_div','r','Nr')
