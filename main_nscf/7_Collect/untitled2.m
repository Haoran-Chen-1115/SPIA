addpath ../0_Public
parameter

Nsigma=10;
lambda=cell(1,Nsigma);
lambda2=cell(1,Nsigma);lam=zeros(1,Nsigma);
rho_2=zeros(1,Nsigma);omega2_2=zeros(1,Nsigma);
for ns=1:Nsigma
    load('Result_q.mat','lambda2_q')
    lambda2_q1=lambda2_q{ns}(1,:);
    lambda2_q2=lambda2_q{ns}(2,:);
    load('../../../size_3_right/main_nscf/7_Collect/Result_q.mat','lambda2_q')
    lambda2_q1p=lambda2_q{ns}(1,:);
    lambda2_q3=lambda2_q{ns}(2,:);
    
    lambda2{ns}=(lambda2_q1+lambda2_q1p)/2/4+lambda2_q2/4*6+lambda2_q3/4;
    
    Ncut=NBEAD/2;
    [rho_2(ns),omega2_2(ns)]=Gbar_Eliashberg(lambda2{ns},Ncut,mustar);
    lam(ns)=lambda2{ns}(1);
end

for ns=1:Nsigma
    load('Result_q.mat','lambda_q')
    lambda_q1=lambda_q{ns}(1,:);
    lambda_q2=lambda_q{ns}(2,:);
    load('../../../size_3_right/main_nscf/7_Collect/Result_q.mat','lambda2_q')
    lambda_q1p=lambda_q{ns}(1,:);
    lambda_q3=lambda_q{ns}(2,:);
    
    lambda{ns}=(lambda_q1+lambda_q1p)/2/4+lambda_q2/4*6+lambda_q3/4;
end
save('Result.mat','Nsigma','lam','lambda2','omega2_2','rho_2','lambda');

%%
ns=2;
BETA=kB*TEMP;
Ncut=NBEAD/2;
lam_f=@(x) interp1(BETA^2*[0:NBEAD/2].^2,BETA^2*[0:NBEAD/2].^2.*lambda2{ns}(1:NBEAD/2+1),x,'pchip');

rho=100;
T_up = TEMP+50;
T_down = TEMP-50;
while abs(rho)>1E-3
    
    T_tmp = (T_up+T_down)/2;
    BETA_tmp = kB*T_tmp;
    
    nu2 = (BETA_tmp)^2*[0:NBEAD/2].^2;
    lambda_nu = lam_f(nu2)/BETA_tmp^2;
    
    lambda_nu(1)=lambda2{ns}(1);
    lambda_nu(2:end) = lambda_nu(2:end)./([1:NBEAD/2].^2);
    lambda_nu(end+1:end+NBEAD/2-1)=lambda_nu(end-1:-1:2);
    
    
    rho=Gbar_Eliashberg2(lambda_nu,Ncut,mustar,omega2_2(ns)*BETA/BETA_tmp);
    if rho<0
        T_up = T_tmp;
    else
        T_down = T_tmp;
    end

    if T_up-T_down<1E-4
        disp('fail to converge')
        break
    end
end
T_c=T_tmp