function inside = D(x) 
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: D.m
%--------------------------------------------------------------------------
% Description: Jump set
% Return 0 if outside of D, and 1 if inside D
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

global gamma_0 gamma_1 V0 V1 c_0 c_10  

% state
z1 = x(1);
z2 = x(2);
q = x(3);

CalcL = CalculateL(z1);

if(q == 0)
    V0 = gamma_0*(CalcL - CalculateLStar()) + (1/2)*z2^2;
elseif(q == 1)
    V1 = gamma_1*(CalcL - CalculateLStar()) + (1/2)*z2^2;
end

if (q == 0 && V0 >= c_0)||(q == 1 && V1 <= c_10)  
    inside = 1;
else
    inside = 0;
end

end