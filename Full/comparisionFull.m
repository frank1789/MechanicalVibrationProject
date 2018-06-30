function [ YS, residual ] = comparisionFull(pInputdata, pForce, pOptimvalue, pDataset)

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
gain = pOptimvalue(7);

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
YS = lsim(G, gain * pForce, pDataset.time.t);

% compute residual
residual = [pDataset.Displacement.x1 - YS(:,1)  ...
            pDataset.Displacement.x2 - YS(:,2) ...
            pDataset.Displacement.x3 - YS(:,3)];

% plot the result with estimated parameters
figure();
plot(pDataset.time.t, YS(:,1), ...
    pDataset.time.t, pDataset.Displacement.x1, ...
    pDataset.time.t, residual(:,1), 'k-.');
legend('Dispalcent x1','Optimvalue x1','residual');
xlabel({'Time','(s)'})
ylabel({'Displacement','(m)'})
grid on
saveas(gcf,'residualfull1','epsc')

figure();
plot(pDataset.time.t, YS(:,2), ...
    pDataset.time.t, pDataset.Displacement.x2, ...
    pDataset.time.t, residual(:,2), 'k-.');
legend('Dispalcent x2','Optimvalue x2','residual');
grid on
xlabel({'Time','(s)'})
ylabel({'Displacement','(m)'})
saveas(gcf,'residualfull2','epsc')

figure();
plot(pDataset.time.t, YS(:,3), ...
    pDataset.time.t, pDataset.Displacement.x3, ...
    pDataset.time.t, residual(:,3), 'k-.');
legend('Dispalcent x3','Optimvalue x3','residual');
grid on
xlabel({'Time','(s)'})
ylabel({'Displacement','(m)'})
saveas(gcf,'residualfull3','epsc')
end
