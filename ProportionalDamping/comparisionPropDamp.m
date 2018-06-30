function [ YS, residual ] = comparisionPropDamp(pInputdata, pForce, pOptimvalue, pDataset)

% initialize stiffness
k1 = pInputdata.stiffness.k1;
k2 = pInputdata.stiffness.k2;
k3 = pInputdata.stiffness.k3;

% I.C. Initial Conditions pass as argument
%{
beta ----------------------------+
alpha ---------------------+     |
gain ------------------+   |     |
mass --------+--+--+   |   |     |
             |  |  |   |   |     |
           [m1 m2 m3 gain alpha beta];
%}
% divide Inititial Condition vector in local variable:
% mass
m1= pOptimvalue(1);
m2= pOptimvalue(2);
m3= pOptimvalue(3);

% gain value
gain = pOptimvalue(4);

% proportinl value alpha and beta
alpha = pOptimvalue(5);
beta  = pOptimvalue(6);

% assemble transfer function
s = tf('s');
M = [m1 0 0; 0 m2 0; 0 0 m3];   % mass matrix [M]
% assemble stifness matrix [K]
K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3];

% proportional damping [C] = alpha * [M] + beta * [K[
C = alpha * M + beta * K;

A = M * s^2 + C * s + K;
G = tf(inv(A));                 % Transfer fuction

% perform simulation
YS = lsim(G, gain * pForce, pDataset.time.t);

% compute residual
residual = [pDataset.Displacement.x1 - YS(:,1)  ...
            pDataset.Displacement.x2 - YS(:,2) ...
            pDataset.Displacement.x3 - YS(:,3)];

% plot the result with estimated parameters
figure();
hold on
plot(pDataset.time.t, YS(:,1), ...
    pDataset.time.t, pDataset.Displacement.x1, ...
    pDataset.time.t, residual(:,1), 'k');
hold off
leg = legend('{System $x_{1}(t)$}','{Optimvalue $x_{1}(t)$}','residual');
set(leg,'Interpreter','latex');
grid on
xlabel({'Time','(s)'})
ylabel({'Displacement','(m)'})
saveas(gcf,'residualpropdamp1','epsc')
clear leg;

figure();
hold on
plot(pDataset.time.t, YS(:,2), ...
    pDataset.time.t, pDataset.Displacement.x2, ...
    pDataset.time.t, residual(:,2), 'k');
hold off
leg = legend('{System $x_{2}(t)$}','{Optimvalue $x_{2}(t)$}','residual');
set(leg,'Interpreter','latex');
grid on
xlabel({'Time','(s)'})
ylabel({'Displacement','(m)'})
saveas(gcf,'residualpropdamp2','epsc')
clear leg;

figure();
hold on
plot(pDataset.time.t, YS(:,3), ...
    pDataset.time.t, pDataset.Displacement.x3, ...
    pDataset.time.t, residual(:,3), 'k');
hold off
leg = legend('{System $x_{3}(t)$}','{Optimvalue $x_{3}(t)$}','residual');
set(leg,'Interpreter','latex');
grid on
xlabel({'Time','(s)'})
ylabel({'Displacement','(m)'})
saveas(gcf,'residualpropdamp3','epsc')
end
