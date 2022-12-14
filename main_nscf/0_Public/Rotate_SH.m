function Rotate=Rotate_SH(R,LMAX,LMMAXC,LPS,NT)
%% Rotate Spherical Harmonics for calculating T0^\daggerT
%% Calculate the angles
if R(4)==0
    Rotate=eye(LMMAXC(NT));
    return
end

cos_b=R(3)/R(4);
sin_b=sqrt(1-cos_b^2);
if sin_b<1E-8
    cos_a=1;
    sin_a=0;
else
    cos_a=R(1)/R(4)/sin_b;
    sin_a=R(2)/R(4)/sin_b;
end

cos_2a=cos_a^2-sin_a^2;
sin_2a=2*sin_a*cos_a;
cos_2b=cos_a^2-sin_a^2;
sin_2b=2*sin_a*cos_a;

cos_3a=-3*cos_a+4*cos_a^3;
sin_3a=3*sin_a-4*sin_a^3;
cos_3b=-3*cos_b+4*cos_b^3;
sin_3b=3*sin_b-4*sin_b^3;

%% The rotation matrices
Rotate=zeros(LMMAXC(NT),LMMAXC(NT));
NCHM=0;
for NCH = 1:LMAX(NT)
    if LPS{NT}(NCH)==0
        % L=0
        NCHM=NCHM+1;
        Rotate(NCHM,NCHM)=1;
    elseif LPS{NT}(NCH)==1
        % L=1
        NCHM=NCHM+3;
        Rotate(NCHM-2:NCHM,NCHM-2:NCHM)=...
            [[sin_a*sin_b,-sin_a*cos_b,-cos_a];...
            [cos_b,           sin_b,          0         ];...
            [cos_a*sin_b,-cos_a*cos_b, sin_a]];
    elseif LPS{NT}(NCH)==2
        NCHM=NCHM+5;
        Rotate(NCHM-4:NCHM,NCHM-4:NCHM)=...
            [[-cos_2a*sin_b, -sin_2a*sin_2b/2,      sqrt(3)/2*sin_2a*cos_b^2,  cos_2a*cos_b,-sin_2a*(1+sin_b^2)/2];...
            [-cos_a*cos_b,  -sin_a*cos_2b,        -sqrt(3)*sin_a*sin_2b/2,   -cos_a*sin_b, -sin_a*sin_2b/2      ];...
            [ 0,             sqrt(3)*sin_2b/2,    (1-3*cos_2b)/4,              0,           -sqrt(3)*cos_b^2/2  ];...
            [ sin_a*cos_b,  -cos_a*cos_2b,        -sqrt(3)*cos_a*sin_2b/2,    sin_a*sin_b, -cos_a*sin_2b/2      ];...
            [ sin_2a*sin_b, -cos_2a*sin_2b/2,      sqrt(3)*cos_2a*cos_b^2/2, -sin_2a*cos_b, cos_2a*(-3+cos_2b)/4]];
    elseif LPS{NT}(NCH)==3
        NCHM=NCHM+7;
        Rotate(NCHM-6:NCHM,NCHM-6:NCHM)=...
            [[ 1/8*sin_3a*sin_b*(-7+cos_2b),             sqrt(3/2)/2*cos_3a*sin_2b,          sqrt(15)/4*sin_3a*sin_b*cos_b^2,...
            -sqrt(5/2)/2*sin_a*(1+2*cos_2a)*cos_b^3,  -sqrt(15)/4*cos_3a*cos_b^2,...
            -sqrt(3/2)/8*sin_3a*(-5*cos_b+cos_3b),    -1/8*cos_a*(-1+2*cos_2a)*(-5+3*cos_2b)];...
            [ sqrt(3/2)/8*sin_2a*(-5*cos_b+cos_3b),     cos_2a*cos_2b,                      sqrt(5/2)/4*sin_2a*cos_b*(-1+3*cos_2b),...
            sqrt(15)*sin_a*cos_a*sin_b*cos_b^2,       sqrt(5/2)/2*cos_2a*sin_2b,...
            1/4*sin_2a*sin_b*(-1+3*cos_2b),           sqrt(3/2)/2*cos_2a*sin_2b];...
            [-sqrt(15)/4*sin_a*sin_b*cos_b^2,          -sqrt(5/2)*cos_a*sin_b*cos_b,        1/16*sin_a*(sin_b-15*sin_3b),...
            -sqrt(3/2)/8*sin_a*(cos_b-5*cos_3b),       1/8*cos_a*(-3+5*cos_2b),...
            sqrt(5/2)/4*sin_a*cos_b*(-1+3*cos_2b),    sqrt(15)/4*cos_a*cos_b^2];...
            [-sqrt(5/2)/2*cos_b^3,                      0,                                  sqrt(3/2)/8*(cos_b-5*cos_3b),...
            1/8*(3*sin_b-5*sin_3b),                   0,...
            -sqrt(15)/2*sin_b*cos_b^2,                 0];...
            [-sqrt(15)/4*cos_a*sin_b*cos_b^2,           sqrt(5/2)*sin_a*sin_b*cos_b,        1/16*cos_a*(sin_b-15*sin_3b),...
            -sqrt(3/2)/8*cos_a*(cos_b-5*cos_3b),      -1/8*sin_a*(-3+5*cos_2b),...
            sqrt(5/2)/4*cos_a*cos_b*(-1+3*cos_2b),    -sqrt(15)/4*sin_a*cos_b^2];...
            [ sqrt(3/2)/4*cos_2a*cos_b*(-3+cos_2b),    -sin_2a*cos_2b,                      sqrt(5/2)/8*cos_2a*(cos_b+3*cos_3b),...
            sqrt(15)/2*cos_2a*sin_b*cos_b^2,         -sqrt(10)*sin_a*cos_a*sin_b*cos_b,...
            1/4*cos_2a*sin_b*(-1+3*cos_2b),          -sqrt(6)*sin_a*cos_a*sin_b*cos_b];...
            [ 1/8*cos_3a*sin_b*(-7+cos_2b),            -sqrt(3/2)/2*sin_3a*sin_2b,          sqrt(15)/4*cos_3a*sin_b*cos_b^2,...
            -sqrt(5/2)/2*cos_a*(-1+2*cos_2a)*cos_b^3,  sqrt(15)/4*sin_3a*cos_b^2,...
            -sqrt(3/2)/8*cos_3a*(-5*cos_b+cos_3b),     1/8*sin_a*(1+2*cos_2a)*(-5+3*cos_2b)]];
    end
end

% % L=0
% for NL=1:CH0{NT}(2)
%     NCH=NCH+1;NCHM=NCHM+1;
%     Rotate(NCHM,NCHM)=1;
% end
% 
% % L=1
% if CH1{NT}(2)~=0
%     for NL=1:CH1{NT}(2)
%         NCH=NCH+1;NCHM=NCHM+3;
%         Rotate(NCHM-2:NCHM,NCHM-2:NCHM)=...
%             [[sin_a*sin_b,-sin_a*cos_b,-cos_a];...
%             [cos_b,           sin_b,          0         ];...
%             [cos_a*sin_b,-cos_a*cos_b, sin_a]];
%     end
% end
% 
% %% Hard-coded
% % L=2
% if CH2{NT}(2)~=0
%     for NL=1:CH2{NT}(2)
%         NCH=NCH+1;NCHM=NCHM+5;
%         Rotate(NCHM-4:NCHM,NCHM-4:NCHM)=...
%             [[-cos_2a*sin_b, -sin_2a*sin_2b/2,      sqrt(3)/2*sin_2a*cos_b^2,  cos_2a*cos_b,-sin_2a*(1+sin_b^2)/2];...
%             [-cos_a*cos_b,  -sin_a*cos_2b,        -sqrt(3)*sin_a*sin_2b/2,   -cos_a*sin_b, -sin_a*sin_2b/2      ];...
%             [ 0,             sqrt(3)*sin_2b/2,    (1-3*cos_2b)/4,              0,           -sqrt(3)*cos_b^2/2  ];...
%             [ sin_a*cos_b,  -cos_a*cos_2b,        -sqrt(3)*cos_a*sin_2b/2,    sin_a*sin_b, -cos_a*sin_2b/2      ];...
%             [ sin_2a*sin_b, -cos_2a*sin_2b/2,      sqrt(3)*cos_2a*cos_b^2/2, -sin_2a*cos_b, cos_2a*(-3+cos_2b)/4]];
%         %             [[(sin_a^2-cos_a^2)*sin_b,-2*sin_a*cos_a*sin_b*cos_b,...
%         %               sqrt(3)*sin_a*cos_a*cos_b^2,...
%         %               (cos_a^2-sin_a^2)*cos_b,-sin_a*cos_a*(2-cos_b^2)],...
%         %              [-cos_a*cos_b,sin_a*(sin_b^2-cos_b^2),...
%         %               -sqrt(3)*sin_a*sin_b*cos_b,...
%         %               -cos_a*sin_b,-sin_a*sin_b*cos_b],...
%         %              [0,sqrt(3)*sin_a*cos_a,...
%         %               sin_b^2-cos_b^2/2,...
%         %               0,-sqrt(3)*cos_b^2/2],...
%         %              [sin_a*cos_b,cos_a*(sin_b^2-cos_b^2),...
%         %               -sqrt(3)*cos_a*sin_b*cos_b,...
%         %               sin_a*cos_b,-cos_a*sin_b*cos_b],...
%         %              [2*sin_a*cos_a*sin_b,(sin_a^2-cos_a^2)*sin_b*cos_b,...
%         %               sqrt(3)*(cos_a^2-sin_a^2)*cos_b^2/2,...
%         %               -2*sin_a*cos_a*cos_b,(sin_a^2-cos_a^2)*(1-cos_b^2/2)]];
%     end
% end
% 
% % L=3
% if CH3{NT}(2)~=0
%     for NL=1:CH3{NT}(2)
%         NCH=NCH+1;NCHM=NCHM+7;
%         Rotate(NCHM-6:NCHM,NCHM-6:NCHM)=...
%             [[ 1/8*sin_3a*sin_b*(-7+cos_2b),             sqrt(3/2)/2*cos_3a*sin_2b,          sqrt(15)/4*sin_3a*sin_b*cos_b^2,...
%             -sqrt(5/2)/2*sin_a*(1+2*cos_2a)*cos_b^3,  -sqrt(15)/4*cos_3a*cos_b^2,...
%             -sqrt(3/2)/8*sin_3a*(-5*cos_b+cos_3b),    -1/8*cos_a*(-1+2*cos_2a)*(-5+3*cos_2b)];...
%             [ sqrt(3/2)/8*sin_2a*(-5*cos_b+cos_3b),     cos_2a*cos_2b,                      sqrt(5/2)/4*sin_2a*cos_b*(-1+3*cos_2b),...
%             sqrt(15)*sin_a*cos_a*sin_b*cos_b^2,       sqrt(5/2)/2*cos_2a*sin_2b,...
%             1/4*sin_2a*sin_b*(-1+3*cos_2b),           sqrt(3/2)/2*cos_2a*sin_2b];...
%             [-sqrt(15)/4*sin_a*sin_b*cos_b^2,          -sqrt(5/2)*cos_a*sin_b*cos_b,        1/16*sin_a*(sin_b-15*sin_3b),...
%             -sqrt(3/2)/8*sin_a*(cos_b-5*cos_3b),       1/8*cos_a*(-3+5*cos_2b),...
%             sqrt(5/2)/4*sin_a*cos_b*(-1+3*cos_2b),    sqrt(15)/4*cos_a*cos_b^2];...
%             [-sqrt(5/2)/2*cos_b^3,                      0,                                  sqrt(3/2)/8*(cos_b-5*cos_3b),...
%             1/8*(3*sin_b-5*sin_3b),                   0,...
%             -sqrt(15)/2*sin_b*cos_b^2,                 0];...
%             [-sqrt(15)/4*cos_a*sin_b*cos_b^2,           sqrt(5/2)*sin_a*sin_b*cos_b,        1/16*cos_a*(sin_b-15*sin_3b),...
%             -sqrt(3/2)/8*cos_a*(cos_b-5*cos_3b),      -1/8*sin_a*(-3+5*cos_2b),...
%             sqrt(5/2)/4*cos_a*cos_b*(-1+3*cos_2b),    -sqrt(15)/4*sin_a*cos_b^2];...
%             [ sqrt(3/2)/4*cos_2a*cos_b*(-3+cos_2b),    -sin_2a*cos_2b,                      sqrt(5/2)/8*cos_2a*(cos_b+3*cos_3b),...
%             sqrt(15)/2*cos_2a*sin_b*cos_b^2,         -sqrt(10)*sin_a*cos_a*sin_b*cos_b,...
%             1/4*cos_2a*sin_b*(-1+3*cos_2b),          -sqrt(6)*sin_a*cos_a*sin_b*cos_b];...
%             [ 1/8*cos_3a*sin_b*(-7+cos_2b),            -sqrt(3/2)/2*sin_3a*sin_2b,          sqrt(15)/4*cos_3a*sin_b*cos_b^2,...
%             -sqrt(5/2)/2*cos_a*(-1+2*cos_2a)*cos_b^3,  sqrt(15)/4*sin_3a*cos_b^2,...
%             -sqrt(3/2)/8*cos_3a*(-5*cos_b+cos_3b),     1/8*sin_a*(1+2*cos_2a)*(-5+3*cos_2b)]];
%     end
% end

%% Recursion relation
% Rotate_recursion;

% %% Calculate CQIJ3
% CQIJ3_interp=zeros(LMMAXC);
% for i=1:LMMAXC
%     for j=1:LMMAXC
%         CQIJ3_interp(i,j)=CQIJ3_fun{i,j}(R(4));
%     end
% end
% CQIJ3_Rotate=Rotate*CQIJ3_interp*(Rotate.');
end