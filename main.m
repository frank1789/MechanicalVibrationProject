close all;
clear;
clc;

% Load folder data
addpath('data_3DOFsystem');
file = dir(fullfile('data_3DOFsystem', '*.mat'));
% load data set: "data_steps"
load data_steps

% number count per encoder revolution is costant
Inputdata.number_count_encoder = 0.0706/16000;

% initialize stiffness
Inputdata.stiffness.k1 = 800;  %[N/m]
Inputdata.stiffness.k2 = 800;  %[N/m]
Inputdata.stiffness.k3 = 400;  %[N/m]

% evaluate displacement and store in struct Displacement
data_steps.Displacement.x1 = x1 * Inputdata.number_count_encoder;
data_steps.Displacement.x2 = x2 * Inputdata.number_count_encoder;
data_steps.Displacement.x3 = x3 * Inputdata.number_count_encoder;

% store original data in the previous struct
data_steps.voltage.v = v;        %[V]
data_steps.time.t = t;           %[s]

% remove generic data
clear x1 x2 x3 t v

% plot figure
figure();
hold on
plot(data_steps.time.t, data_steps.Displacement.x1)
plot(data_steps.time.t, data_steps.Displacement.x2)
plot(data_steps.time.t, data_steps.Displacement.x3)
grid()
legend('Displacement x1', 'Displacement x2','Displacement x3')
hold off

% compute symbolic equation [G(s)]
[ A ] = symbolicequation();

%% use step response to verify the ratio beetween the stiffnesses of the springs.
% compute a new estimation for the voltage-to-force coefficient

% define gain_ x from encoder
gain_x = Inputdata.number_count_encoder;

%unkonw ?
gain_v = 5.250; 

%unkonw ?
gain_tot = (gain_x/gain_v); % ????

% %% Compute the force F
%                                                  % costant from readme
% Inputdata.k_a = 2;                               % servo amp gain
% Inputdata.k_t = 0.1;                             % Servo motor Torque
% Inputdata.k_mp = 26.25;                          % Motor Pinioon pitch radius inverse 
% F = zeros(length(data_steps.time.t), 3);         % initialize F's vector
% f = (Inputdata.k_a * Inputdata.k_t * Inputdata.k_mp) * data_steps.voltage.v;   % tmp vector 
% F(:,1) = f;                                      % matri of forces

%% compute a new estimation for the voltage-to-force coefficient
% syms gv 
% gainvolt = K\[gv; 0; 0];
% 
% F = F * gainvolt;


%% compute transfer function
% s = tf('s');
% 
% denom = 1/det(A);
% num = inv([ m1*s^2 + c1*s + k1,                     -k1,                       0
%                            -k1, m2*s^2 + c2*s + k1 + k2,                     -k2
%                              0,                     -k2, m3*s^2 + c3*s + k2 + k3]);
% % transfer fuction
% G = tf(num, denom);
% display(G);
% [YS,TS,XS] = lsim(G,F,Inputdata.time.t);
% 
% % comparision plot
% figure();
% plot(Inputdata.time.t,YS);
% grid on;
% legend('tf x1','tf x2','tf x3');



%% Optimization estimated parameters
%{ 
 file data impulses.mat: use the impulse response to identify the parameters.
 Choose between response to impulsive force or response to initial conditions:
 in the first case, due to the approximation of the force estimation, you
 can consider the voltage-to-force coefficient as one of the parameters to
 be estimated.

 - repeat the estimation of your parameters by assuming proportional damping.
 How many parameters have to be identified?

 -in order to improve the estimation of the parameters, the springs can be
 detached and the masses can be blocked at their equilibrium position:
 describe a possible strategy for the estimation of the parameters by
 studying the behaviour of three different 1 DOF systems.
%}
% load dataset
load data_impulses.mat

% evaluate displacement and store in struct Displacement
data_impulses.Displacement.x1 = x1 * Inputdata.number_count_encoder;
data_impulses.Displacement.x2 = x2 * Inputdata.number_count_encoder;
data_impulses.Displacement.x3 = x3 * Inputdata.number_count_encoder;

% store original data in the previous struct
data_impulses.voltage.v = v;        %[V]
data_impulses.time.t = t;           %[s]

% compute the force
[ F ] = computeforce( data_impulses );

% remove generic data
clear x1 x2 x3 t v

%{
 I.C. Initial Conditions
gain ---------------------------+
damper c3--------------------+  |
damper c2-----------------+  |  |
damper c1--------------+  |  |  |
mass --------+--+--+   |  |  |  |     
             |  |  |   |  |  |  |
       x0 = [m1 m2 m3 c1 c2 c3  gain];
%}
x0 = [1 1 1 1 1 1];

%{ 
define lower bound and upper bound
gain  ------------+
damper 3--------+ |
damper 2------+ | |
damper 1----+ | | |
mass -+-+-+ | | | |
      | | | | | | |
%}
LB = [1 1 1 0 0 0];
UB = [2 2 2 5 5 5];

objfun = @(x0)errormio(x0, F, Inputdata, data_impulses);
[xfmincon] = fmincon(objfun,x0,[],[],[],[],LB,UB);
display(xfmincon);

% %% fminsearch
%options = optimset('PlotFcns',@optimplotfval);
%f = @(x0)errormio(x0, F, Inputdata, data_impulses);
%[xfminsearch] = fminsearch(f,x0,options)
% %% plot comparison
[ data_impulses.YS ] = comparision(Inputdata, F, xfmincon, data_impulses);
%comparision(Inputdata, F, xfminsearch);


%% compare the result with "goodnessOfFit" function
% reference data
ref_displacement = [data_impulses.Displacement.x1 data_impulses.Displacement.x2 data_impulses.Displacement.x3];

% cost function Normalized root mean square error, where, ? indicates the
% 2-norm of a vector. fit is a row vector of length N and i = 1,...,N, where
% N is the number of channels.
cost_func = 'NRMSE';

fit = goodnessOfFit(data_impulses.YS, ref_displacement, cost_func);
display(fit);
