function [ A ] = symbolicequation()
%SYMBOLICEQUATION Summary of this function goes here
%   Detailed explanation goes here

% define symbolic variable
syms m1 m2 m3 c1 c2 c3 k1 k2 k3 x1__ x2__ x3__ x1_ x2_ x3_ sx1 sx2 sx3 

% assemble mass matrix [M]
M = [m1 0 0; 0 m2 0; 0 0 m3];
display(M);

% assemble damping matric [C]
C = [c1 0 0; 0 c2 0; 0 0 c3];
display(C);

% assemble stifness matrix [K]
K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3];
display(K);

% compute matrix [G(s)] -> A
syms s
A = M * s^2 + C*s + K;
display(A)

% compute det([G(s)])
denom = det(A);
display(denom);

% compute inverse inv([G(s)])
num = inv(A);
display(num);

end

