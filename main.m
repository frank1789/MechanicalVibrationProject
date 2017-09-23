close all;
clear;
clc;

% Load folder data
addpath('data_3DOFsystem', 'Full', 'ProportionalDamping');
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
[ A, K ] = symbolicequation();

% compute a new estimation for the voltage-to-force coefficient
%{
file data steps.mat: use the step response to verify the ratio between the
stiffnesses of the springs. Compute a new estimation for the 
voltage-to-force coefficient.
%}

[ gain_v, Ratio_k3_k2, Ratio_k3_k1 ] = computeSteadyStateRatioStiff(K, data_steps, Inputdata);

geterrorspring(gain_v, Ratio_k3_k2, Ratio_k3_k1, Inputdata, data_steps);

% Optimization estimated parameters
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
problem.x0 = [1.1 1.1 1.01 0.01 0.01 0.01 gain_v];

%{ 
define lower bound and upper bound
gain  --------------------+
damper 3----------------+ |
damper 2--------------+ | |
damper 1------------+ | | |
mass ---------+-+-+ | | | |
              | | | | | | |
%}
problem.lb = [1 1 1 0 0 0 0];
problem.ub = [2 2 2 5 5 5 7];

problem.options = optimoptions('fmincon','Display','iter', ...
                               'Algorithm','sqp','PlotFcn',@optimplotx);
problem.solver = 'fmincon';
problem.objective = @(x0)estimateFull(x0, F, Inputdata, data_impulses);
full.x = fmincon(problem);
display(full.x);

% plot comparison and residual
[ full.YS, full.residuals ] = comparisionFull(Inputdata, F, full.x, data_impulses);

% compare the result with "goodnessOfFit" function
% reference data
ref_displacement = [data_impulses.Displacement.x1 ...
                    data_impulses.Displacement.x2 ...
                    data_impulses.Displacement.x3];

% cost function Normalized root mean square error, where, ? indicates the
% 2-norm of a vector. fit is a row vector of length N and i = 1,...,N, where
% N is the number of channels.
cost_func = 'NRMSE';
full.fit = goodnessOfFit(full.YS, ref_displacement, cost_func);
fprintf('Response full method comparision: %.2f%%\t%.2f%%\t%.2f%%\n', full.fit(1) * 100, full.fit(2) * 100, full.fit(3) * 100);
close all

% Proportinal damping
clear problem
%{
 I.C. Initial Conditions
beta ----------------------------+
alpha ---------------------+     |
gain ------------------+   |     |
mass --------+--+--+   |   |     |
             |  |  |   |   |     |
       x0 = [m1 m2 m3 gain alpha beta];
%}
problem.x0 = [1.1 1.1 1.1 gain_v 1 1];

%{
define lower bound and upper bound
beta ----------------------------+
alpha ---------------------+     |
gain ------------------+   |     |
mass --------+--+--+   |   |     |
             |  |  |   |   |     |
           [m1 m2 m3 gain alpha beta];
%}
problem.lb = [1 1 1 0 0 0];
problem.ub = [2 2 2 7 10 10];
problem.options = optimoptions('fmincon','Display','iter', ...
                               'Algorithm','sqp','PlotFcn',@optimplotx);
problem.solver = 'fmincon';
problem.objective = @(x0)estimatePropDamp(x0, F, Inputdata, data_impulses);
prodamping.x = fmincon(problem);
display(prodamping.x);

% plot comparison and residual
[ prodamping.YS, prodamping.residual ] = comparisionPropDamp(Inputdata, F, prodamping.x, data_impulses);

% compare the result with "goodnessOfFit" function
% reference data
ref_displacement = [data_impulses.Displacement.x1 ...
                    data_impulses.Displacement.x2 ...
                    data_impulses.Displacement.x3];

% cost function Normalized root mean square error, where, ? indicates the
% 2-norm of a vector. fit is a row vector of length N and i = 1,...,N, where
% N is the number of channels.
cost_func = 'NRMSE';
fit = goodnessOfFit(prodamping.YS, ref_displacement, cost_func);
fprintf('Response proportional damping comparision: %.2f%%\t%.2f%%\t%.2f%%\n', fit(1) * 100, fit(2) * 100, fit(3) * 100);
close all


