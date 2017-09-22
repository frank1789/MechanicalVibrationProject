function [ err ]=errormio(x0, pForce, pInputdata, pDataset)

% initialize local temporary error variable
e1 = 0;
e2 = 0;
e3 = 0;

% initialize stiffness
k1 = pInputdata.stiffness.k1;
k2 = pInputdata.stiffness.k2;
k3 = pInputdata.stiffness.k3;

% I.C. Initial Conditions pass as argument
%{
gain----------------------------+
damper-----------------------+  |
damper--------------------+  |  |
damper-----------------+  |  |  |
mass---------+--+--+   |  |  |  |
             |  |  |   |  |  |  |
estimated = [m1 m2 m3 c1 c2 c3 gain];
%}
% divide Inititial Condition vector in local variable:
% mass
m1= x0(1);
m2= x0(2);
m3= x0(3);

% dampers
c1 = x0(4);
c2 = x0(5);
c3 = x0(6);

% gain value
%gain = x0(7);

% assemble transfer function
s = tf('s');
M = [m1 0 0; 0 m2 0; 0 0 m3];   % mass matrix [M]
C = [c1 0 0; 0 c2 0; 0 0 c3];   % assemble damping matric [C]
K = [ k1   -k1  0;              % assemble stifness matrix [K]
     -k1 k2+k1 -k2;
       0   -k2  k2+k3];
A = M * s^2 + C*s + K;
G = tf(inv(A));                 % Transfer fuction

% perform simulation
[YS] = lsim(G, pForce, pDataset.time.t);

% compute the difference beetwen referende displacement(pInputdata) and calculated(YS) 
for j = 2:length(pDataset.Displacement.x1)
    e1 = e1 + (YS(j,1) - pDataset.Displacement.x1(j)).^2;
    e2 = e2 + (YS(j,2) - pDataset.Displacement.x2(j)).^2;
    e3 = e3 + (YS(j,3) - pDataset.Displacement.x3(j)).^2;
end

% compute root mean square on error
err = rms(e1 + e2 + e3);

% free memory
clear k1 k2 k3 e1 e2 e3
end