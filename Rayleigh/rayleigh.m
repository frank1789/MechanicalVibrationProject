function [U,  freqs] = rayleigh(pOptimValue, pInputdata)
% rayleigh compute vibration modes return as modes and compute
% naturals frequencies return as freqs.
% Passed:
%   as argument pOptimValue -> mass retrived from optimztion
%   as argument pInputdata -> stiffness of spring from readme
% instaziate local variable from Optimalvalue and data
% m1 = pOptimValue(1);
% m2 = pOptimValue(2);
% m3 = pOptimValue(3);
%
% k1 = pInputdata.stiffness.k1;
% k2 = pInputdata.stiffness.k2;
% k3 = pInputdata.stiffness.k3;
%
% assemble mass matrix [M]
% M = [m1 0 0; 0 m2 0; 0 0 m3];
%
% assemble stifness matrix [K]
% K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3];
%
% instaziate return variable
% freqs = zeros(3,1);
% define simbolic component of vector modes u
% syms u1 u2
% 
% define the Rayleigh's quotient R(X) = w^2 = (x.' * [K] * x) / (x.' * [K] * x)
% 
% rayquot = ([1; u1; u2].' * K * [1; u1; u2])  /  ([1; u1; u2].' * M * [1; u1; u2]);
% 
% remove the imaginary part
% assume(u1, 'real');
% assume(u2, 'real');
% 
% extract the vector mode shape U
% compute the mode shape
% simbolicModeShape = vpasolve([diff(rayquot,u1) == 0, diff(rayquot,u2) ==0],[u1 u2]);
% U = double([[1 simbolicModeShape.u1(1) simbolicModeShape.u2(1) ].' ...
%     [1 simbolicModeShape.u1(2) simbolicModeShape.u2(2) ].' ...
%     [1 simbolicModeShape.u1(3) simbolicModeShape.u2(3) ].']);
% 
% compute the fundamental frequency
% computefreq = @(u1,u2)([1; u1; u2].' * K * [1; u1; u2])  /  ([1; u1; u2].' * M * [1; u1; u2]);
% 
% for i = 1:3
%     freqs(i) = sqrt(computefreq(U(2,i),U(3,i)));
% end

% instaziate local variable from Optimalvalue and data
m1 = pOptimValue(1);
m2 = pOptimValue(2);
m3 = pOptimValue(3);

k1 = pInputdata.stiffness.k1;
k2 = pInputdata.stiffness.k2;
k3 = pInputdata.stiffness.k3;

% assemble mass matrix [M]
M = [m1 0 0; 0 m2 0; 0 0 m3];

% assemble stifness matrix [K]
K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3];

% instaziate return value freqs
freqs = zeros(length(M),1);

% define simbolic component of vector modes u
syms u1 u2

% define the Rayleigh's quotient R(X) = w^2 = (x.' * [K] * x) / (x.' * [K] * x)

rayquot = ([1; u1; u2].' * K * [1; u1; u2])  /  ([1; u1; u2].' * M * [1; u1; u2]);

% remove the imaginary part
assume(u1, 'real');
assume(u2, 'real');

% extract the vector mode shape U
% compute the mode shape
simbolicModeShape = vpasolve([diff(rayquot,u1) == 0, diff(rayquot,u2) ==0],[u1 u2]);
U = double([[1 simbolicModeShape.u1(1) simbolicModeShape.u2(1) ].' ...
    [1 simbolicModeShape.u1(2) simbolicModeShape.u2(2) ].' ...
    [1 simbolicModeShape.u1(3) simbolicModeShape.u2(3) ].']);

% compute the fundamental frequency
computefreq = @(u1,u2)([1; u1; u2].' * K * [1; u1; u2])  /  ([1; u1; u2].' * M * [1; u1; u2]);

for i = 1:3
    freqs(i) = sqrt(computefreq(U(2,i),U(3,i)));
end

end