function [ C ] = getdampingmatrix(pOptimalvalue, pInputdata)

% initialize stiffness
k1 = pInputdata.stiffness.k1;
k2 = pInputdata.stiffness.k2;
k3 = pInputdata.stiffness.k3;

% mass
m1= pOptimalvalue(1);
m2= pOptimalvalue(2);
m3= pOptimalvalue(3);

% proportinl value alpha and beta
alpha = pOptimalvalue(5);
beta  = pOptimalvalue(6);

% mass matrix [M]
M = [m1 0 0; 0 m2 0; 0 0 m3];   

% assemble stifness matrix [K]
K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3];

% proportional damping [C] = alpha * [M] + beta * [K[
C = alpha * M + beta * K;
end
