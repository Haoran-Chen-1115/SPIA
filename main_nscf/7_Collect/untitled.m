addpath ../0_Public
parameter

load('Result_q.mat')
ns=2;
lambda2_q1=lambda2_q{ns}(1,:);
lambda2_q2=lambda2_q{ns}(2,:);
load('../../../size_3_right/main_nscf/7_Collect/Result_q.mat')
lambda2_q1p=lambda2_q{ns}(1,:);
lambda2_q3=lambda2_q{ns}(2,:);

lambda2=(lambda2_q1+lambda2_q1p)/2/4+lambda2_q2/4*6+lambda2_q3/4;
mustar=0.12;

Ncut=NBEAD/2;
[rho_2,omega2_2]=Gbar_Eliashberg(lambda2,Ncut,mustar);

BETA=kB*TEMP;
Ncut=NBEAD/2;
lam_f=@(x) interp1(BETA^2*[0:NBEAD/2].^2,BETA^2*[0:NBEAD/2].^2.*lambda2(1:NBEAD/2+1),x,'pchip');

rho=100;
T_up = TEMP+50;
T_down = TEMP-50;
while abs(rho)>1E-3
    
    T_tmp = (T_up+T_down)/2;
    BETA_tmp = kB*T_tmp;
    
    nu2 = (BETA_tmp)^2*[0:NBEAD/2].^2;
    lambda_nu = lam_f(nu2)/BETA_tmp^2;
    
    lambda_nu(1)=lambda2(1);
    lambda_nu(2:end) = lambda_nu(2:end)./([1:NBEAD/2].^2);
    lambda_nu(end+1:end+NBEAD/2-1)=lambda_nu(end-1:-1:2);
    
    
    rho=Gbar_Eliashberg2(lambda_nu,Ncut,mustar,omega2_2*BETA/BETA_tmp);
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