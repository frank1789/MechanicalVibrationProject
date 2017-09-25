function [ freqs, modes ] = getNaturalFrequencies(pOptimValue, pInputdata)
% getNaturalFrequencies compute vibration modes return as modes and compute
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
% % assemble mass matrix [M]
% M = [m1 0 0; 0 m2 0; 0 0 m3];
%
% assemble stifness matrix [K]
% K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3];
%
% instaziate return variable
% freqs = zeros(3,1);
%
% produces a diagonal matrix D of generalized eigenvalues and a full matrix V
% [ V, D ] = eig(K, M);
%
% compute the frequencies
% for i = 1:length(M)
%     freqs(i) = sqrt(D(i,i));
%     % compute the modes
%     V(1:3,i) = V(1:3,i)./V(1,i);
% end
% 
% return modes
% modes = V;

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

% instaziate return variable
freqs = zeros(3,1);

% produces a diagonal matrix D of generalized eigenvalues and a full matrix V
[ V, D ] = eig(K, M);

% compute the frequencies
for i = 1:length(M)
    freqs(i) = sqrt(D(i,i));
    % compute the modes eignvectors
    V(1:3,i) = V(1:3,i)./V(1,i);
end

% return modes
modes = V;

end