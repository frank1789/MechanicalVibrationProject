function [M,K,freqs,modes] = getNaturalFrequencies( pOptimValue, pInputdata)

m1 = pOptimValue(1);
m2 = pOptimValue(2);
m3 = pOptimValue(3);

c1 = pOptimValue(4);
c2 = pOptimValue(5);
c3 = pOptimValue(6);

k1 = pInputdata.stiffness.k1;
k2 = pInputdata.stiffness.k2;
k3 = pInputdata.stiffness.k3;

% assemble mass matrix [M]
M = [m1 0 0; 0 m2 0; 0 0 m3];
display(M);

% assemble damping matric [C]
C = [c1 0 0; 0 c2 0; 0 0 c3];
display(C);

% assemble stifness matrix [K]
K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3];
display(K);


freqs = zeros(3,1);


%M = [m1,0;0,m2];
%K = [k1+k2,-k2;-k2,k2+k3];
[V,D] = eig(K,M);
for i = 1:length(M)
    freqs(i) = sqrt(D(i,i));
end
modes = V;
 
end