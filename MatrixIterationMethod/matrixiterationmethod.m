function [freq,eig] = matrixiterationmethod(pOptimValue, pInputdata)

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
A = K\M;

% this function, Miter.m and Mult.m should be in the same folder,and set
% the Matalb path to this folder following 6 lines contain
% problem-dependent data
n = length(A);
nvec = length(A);
xm=[1 0 0;0 1 0; 0 0 1];
eps=0.00001;
xs=[1 1 1];
%end of problem-dependent data

[freq,eig,b,c,xx] = miter(A, xs, n, nvec, xm, eps);

% fprintf('Solution of eigenvalue peoblem by \n');
% fprintf('matrix iteration method \n\n');
% fprintf('Natural frequencies: \n\n');
% fprintf('    %0.5f      %0.5f       %0.5f \n\n', freq(1),freq(2),freq(3));
% fprintf('Mode shapes (Columnwise): \n\n');
% for i=1:nvec
%     fprintf('    %0.5f      %0.5f       %0.5f \n', eig(i,1),eig(i,2),eig(i,3));
% end

end