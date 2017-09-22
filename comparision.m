function [ YS ] = comparision(pInputdata, pForce, pOptimvalue, pDataset, pGenplot)


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
m1= pOptimvalue(1);
m2= pOptimvalue(2);
m3= pOptimvalue(3);

% dampers
c1 = pOptimvalue(4);
c2 = pOptimvalue(5);
c3 = pOptimvalue(6);

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

if pGenplot == true
    %comparision plot
    figure();
    hold on
    plot(pDataset.time.t, YS(:,1));
    plot(pDataset.time.t, pDataset.Displacement.x1);
    hold off
    legend('Dispalcent x1','Optimvalue x1');
    grid on
    
    figure();
    hold on
    plot(pDataset.time.t, YS(:,2));
    plot(pDataset.time.t, pDataset.Displacement.x2);
    grid on;
    hold off
    legend('Dispalcent x2','Optimvalue x2');
    
    figure();
    hold on
    plot(pDataset.time.t,YS(:,3));
    plot(pDataset.time.t, pDataset.Displacement.x3);
    grid on;
    hold off
    legend('Dispalcent x3','Optimvalue x3');
end
end
