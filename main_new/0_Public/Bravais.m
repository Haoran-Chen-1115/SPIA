function IBRAV = Bravais(A,B)
%% For determing the summing weight and the equivilant k-points,
%  there are three steps
%  1. Determine the lattice type
%  2. Determine the symmetric operations
%  3. Act the symmetric operations on each k-point

%% A subroutine to determine the type of the lattice
% Follow the description of https://doi.org/10.1016/j.commatsci.2017.01.017
% Here we are in fact only deal with the supercell, 
% instead of the primitive cell.

% IBRAV =
% 1, Tricilinic
% 2, Primitive Monoclinic
% 3, Base-centered Monoclinic
% 4, Primitive Orthorhombic
% 5, Base-centered Orthorhombic
% 6, Body-centered Orthorhombic
% 7, Face-centered Orthorhombic
% 8, Primitive Tetragonal
% 9, Body-centered Tetragonal
% 10, Rhombohedral
% 11, Hexagonal
% 12, Primitive Cubic
% 13, Body-centered Cubic
% 14, Face-centered Cubic

a = zeros(1,3);
b = zeros(1,3);
for i = 1:3
    a(i) = norm(A(:,i));
    b(i) = norm(B(i,:));
end

cos12 = A(:,1)'*A(:,2)/a(1)/a(2);
cos23 = A(:,2)'*A(:,3)/a(2)/a(3);
cos13 = A(:,3)'*A(:,1)/a(3)/a(1);
cos12_23 = (A(:,1)+A(:,2))'*(A(:,2)+A(:,3));
cos13_23 = (A(:,1)+A(:,3))'*(A(:,2)+A(:,3));
cos12_13 = (A(:,1)+A(:,2))'*(A(:,1)+A(:,3));

a1p2 = norm(A(:,1)+A(:,2));
a1p3 = norm(A(:,1)+A(:,3));
a2p3 = norm(A(:,2)+A(:,3));
a1m2 = norm(A(:,1)-A(:,2));
a1m3 = norm(A(:,1)-A(:,3));
a2m3 = norm(A(:,2)-A(:,3));

eps = 1E-5;
IBRAV=0;
if abs(cos12-cos23)<eps && abs(cos12-cos13)<eps
    if abs(a(1)-a(2))<eps && abs(a(1)-a(3))<eps
        if abs(cos12)<eps
            IBRAV=12; % Primitive cubic
        elseif abs(cos12+1/3)<eps
            IBRAV = 13; % Body-centered cubic
        elseif abs(cos12-1/2)<eps
            IBRAV = 14; % Face-centered cubic
        else
            IBRAV = 10; % Rhobohedral
        end
    elseif abs(cos12)<eps
        if abs(a(1)-a(2))<eps
            IBRAV=8; % Primitive tetragonal
        elseif a(3)-a(2)>eps && a(2)-a(1)>eps % require a(3)>a(2)>a(1)
            IBRAV = 4; % Primitive orthorhombic
        end
    end
elseif abs(cos13-cos23)<eps 
    if abs(cos13)<eps
        if abs(a(1)-a(2))<eps
            % The default specific axis is a3
            if abs(cos12+1/2)<eps
                IBRAV = 11; % Hexagoanl
            elseif cos12<-eps % Require cos12<0
                IBRAV = 5; % Base-centered orthorhombic
            end
        elseif cos12<-eps && a(1)-a(2)>eps % Require cos12<0 and a(1)>a(2)
            % The default specific axis is a2
            IBRAV = 2; % Primitive monoclinic
        end
    else
        if abs(a(1)-a(2))<eps && abs(a(1)-a(3))<eps...
                && abs(cos12_23)<eps && abs(cos12_13)<eps && abs(cos13_23)<eps...
                && abs(a1p3-a2p3)<eps
            % The default specific axis is a3
            IBRAV = 9; % Body-centered tetragonal
        elseif abs(a(1)-a(2))<eps && cos13<-eps && cos23<-eps
            IBRAV = 3; % Base-centered monoclinic
        end
    end
else
    if abs(a(1)-a(2))<eps && abs(a(1)-a(3))<eps...
            && abs(cos12_23)<eps && abs(cos12_13)<eps && abs(cos13_23)<eps...
            && a1p2-a1p3>eps && a1p3-a2p3>eps
        IBRAV = 6; % Body-centered orthorhombic
    elseif abs(a(1)-a2m3)<eps && abs(a(2)-a1m3)<eps && abs(a(3)-a1m2)<eps ...
            && norm(A(:,1)+A(:,2)-A(:,3))-norm(A(:,1)+A(:,3)-A(:,2)) > eps ...
            && norm(A(:,1)+A(:,3)-A(:,2))-norm(A(:,2)+A(:,3)-A(:,1)) > eps
        IBRAV = 7; % Face-centered orthorhombic
    elseif cos12>cos13 && cos13>cos23 && cos23>eps
        IBRAV = 1; % Triclinic
    end
end

end
