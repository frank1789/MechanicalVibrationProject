function [freq,eig] = matrixiterationmethod(pOptimValue, pInputdata)
%errormio function to evaluate the simultate response of the linear system
%by lsim. Pass argument
%x0 = initical condition,
%pForce = input force,
%pInputdata = costant from readme file,
%pDataset = value from recorded data
%
%Return value err = error of trhee displacement
%
%How to works:
% initialize local temporary error variable
% e1 = 0;
% e2 = 0;
% e3 = 0;
%
% initialize stiffness
% k1 = pInputdata.stiffness.k1;
% k2 = pInputdata.stiffness.k2;
% k3 = pInputdata.stiffness.k3;
%
% I.C. Initial Conditions pass as argument
%
% gain ---------------------------+
% damper c3--------------------+  |
% damper c2-----------------+  |  |
% damper c1--------------+  |  |  |
% mass---------+--+--+   |  |  |  |
%              |  |  |   |  |  |  |
% estimated = [m1 m2 m3 c1 c2 c3 gain];
%
% divide Inititial Condition vector in local variable:
% % mass
% m1= x0(1);
% m2= x0(2);
% m3= x0(3);
%
% % dampers
% c1 = x0(4);
% c2 = x0(5);
% c3 = x0(6);
%
% gain value
% gain = x0(7);
%
% assemble transfer function
% s = tf('s');
% M = [m1 0 0; 0 m2 0; 0 0 m3];   % mass matrix [M]
% C = [c1 0 0; 0 c2 0; 0 0 c3];   % assemble damping matric [C]
% K = [ k1   -k1  0;              % assemble stifness matrix [K]
%      -k1 k2+k1 -k2;
%        0   -k2  k2+k3];
% A = M * s^2 + C*s + K;
% G = tf(inv(A));                 % Transfer fuction
%
% perform simulation
% [YS] = lsim(G, pForce, pDataset.time.t);
%
% compute the difference beetwen referende displacement(pDataset) and calculated(YS)
% for j = 2:length(pDataset.Displacement.x1)
%     e1 = e1 + (YS(j,1) - pDataset.Displacement.x1(j)).^2;
%     e2 = e2 + (YS(j,2) - pDataset.Displacement.x2(j)).^2;
%     e3 = e3 + (YS(j,3) - pDataset.Displacement.x3(j)).^2;
% end
%
% % compute root mean square on error
% err = rms(e1 + e2 + e3);

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
% display(A);
%===============================================================================
%
%Program 10
%Main program which calls the function MITER
%
%================================================================================
%Run "Program10" in MATLAB command window. Progrm10.m, Miter.m and Mult.m
%should be in the same folder,and set the Matalb path to this folder
%following 6 lines contain problem-dependent data
n = length(A);
nvec = length(A);
% d=[1.0 1.0 1.0; 1.0 2.0 2.0; 1.0 2.0 3.0];
xm=[1 0 0;0 1 0; 0 0 1];
eps=0.00001;
xs=[1 1 1];
%end of problem-dependent data
[freq,eig,b,c,xx]=miter(A,xs,n,nvec,xm,eps);
% fprintf('Solution of eigenvalue peoblem by \n');
% fprintf('matrix iteration method \n\n');
% fprintf('Natural frequencies: \n\n');
% fprintf('    %0.5f      %0.5f       %0.5f \n\n', freq(1),freq(2),freq(3));
% fprintf('Mode shapes (Columnwise): \n\n');
% for i=1:nvec
%     fprintf('    %0.5f      %0.5f       %0.5f \n', eig(i,1),eig(i,2),eig(i,3));
% end

end




