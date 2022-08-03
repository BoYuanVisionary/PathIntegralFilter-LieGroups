function [outputArg1,outputArg2] = ILQR(algebra,obe,TI,sigma_B)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
inertia_tensor = diag([1,1.11,1.3]);
cx1 = (algebra+2*obe*algebra)*TI/sigma_B^2;
cx2 = (1+2*obe)*TI/sigma_B^2;
cu1 = 0;
cu2 = TI;
cux = 0;
fx1 = 2*algebra*TI;
fu1 = inertia_tensor^(-1)*TI;

N = zeros(6);

A = fx1;
B = fu1;
Q = cx2 / 2;
R = cu2 / 2;
N()
end

