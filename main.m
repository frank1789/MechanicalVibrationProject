close all;
clear;
clc;

% Load folder data
addpath('data_3DOFsystem', 'Full', 'ProportionalDamping', ... 
    'ModalAnalisys', 'Rayleigh', 'MatrixIterationMethod', ...
    'Sinesweep');
file = dir(fullfile('data_3DOFsystem', '*.mat'));
% load data set: "data_steps"
load data_steps

% number count per encoder revolution is costant
Inputdata.number_count_encoder = 0.0706 / 16000;

% initialize stiffness
Inputdata.stiffness.k1 = 800;  %[N/m]
Inputdata.stiffness.k2 = 800;  %[N/m]
Inputdata.stiffness.k3 = 400;  %[N/m]

% nomianl gain volatge
Inputdata.voltage.v = 5.25;  %[v]

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

[ gain_v, Ratio_k3_k2, Ratio_k3_k1 ] = computeSteadyStateRatioStiff(K, ...
    data_steps, Inputdata);

[ g_v_per, R31_per, R32_per ] = geterrorspring(gain_v, Ratio_k3_k2, ...
    Ratio_k3_k1, Inputdata);
fprintf('Stiffnesses ratios and voltage-to-force coefficients results ') 
fprintf('from analysis:\n\tgain ratio:\t%.5f %%\n\tratio k3/k1:\t %.5f %%\n\tratio k3/k2:\t %.5f %%\n\n', g_v_per, R31_per, R32_per);

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
[ F ] = computeforce(data_impulses);

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
fprintf('Response full method comparision: %.2f %%\t%.2f %%\t%.2f %%\n', full.fit(1:3) * 100);
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
fprintf('The matrix [C] is:\n');
fprintf('\t|% .5f % .5f % .5f |\n', getdampingmatrix(prodamping.x, Inputdata).');

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
prodamping.fit = goodnessOfFit(prodamping.YS, ref_displacement, cost_func);
fprintf('Response proportional damping comparision: %.2f %%\t%.2f %%\t%.2f %%\n', prodamping.fit(1:3) * 100);
close all

% Modal Analisys
% full method
[full.freqs,full.mode] = getNaturalFrequencies(full.x, Inputdata);
fprintf('\nFull method:\n frequencies %.5f Hz, %.5f Hz, %.5f Hz\n modes:\n', full.freqs(1:3));
fprintf('\t|% .5f % .5f % .5f |\n', full.mode.')

% proportinal method
[prodamping.freqs,prodamping.mode] = getNaturalFrequencies(prodamping.x, Inputdata);
fprintf('\nProportional damping:\n frequencies %.5f Hz, %.5f Hz, %.5f Hz\n modes:\n', prodamping.freqs(1:3));
fprintf('\t|% .5f % .5f % .5f |\n', prodamping.mode.')

%% MODAL ANALYSIS
%{
 use Rayleigh quotient and Matrix Iteration Method to estimate the modes of
 the un-damped system. Compare the results with the ones of the eigenvalue
 problem.
%}
% RAYLEIGHT
fprintf('\nMODAL ANALYSIS - Rayleigh:');

% full method
[full.rayleigh.mode,  full.rayleigh.freqs] = rayleigh(full.x, Inputdata);
fprintf('\nFull method:\n frequencies %.5f Hz, %.5f Hz, %.5f Hz\n modes:\n', full.rayleigh.freqs(1:3));
fprintf('\t|% .5f % .5f % .5f |\n', full.rayleigh.mode.');

% proportianl damping
[prodamping.rayleigh.mode,  prodamping.rayleigh.freqs] = rayleigh(prodamping.x, Inputdata);
fprintf('\nProportional damping:\n frequencies %.5f Hz, %.5f Hz, %.5f Hz\n modes:\n', prodamping.rayleigh.freqs(1:3));
fprintf('\t|% .5f % .5f % .5f |\n', prodamping.rayleigh.mode.');

% MATRIX ITERATION METHOD
% full method
fprintf('\nMODAL ANALYSIS - Matrix Iteration Method:');
[full.mim.freqs, full.mim.mode] = matrixiterationmethod(full.x, Inputdata);
fprintf('\nFull method:\n frequencies %.5f Hz, %.5f Hz, %.5f Hz\n modes:\n', full.mim.freqs(1:3));
fprintf('\t|% .5f % .5f % .5f |\n', full.mim.mode.');

% proportianl damping
[prodamping.mim.mode,  prodamping.mim.freqs] = rayleigh(prodamping.x, Inputdata);
fprintf('\nProportional damping:\n frequencies %.5f Hz, %.5f Hz, %.5f Hz\n modes:\n', prodamping.mim.freqs(1:3));
fprintf('\t|% .5f % .5f % .5f |\n', prodamping.mim.mode.');

%% LAPLACE - transfer function
%{
 use Laplace transform to plot the transfer functions between the applied
 force and the positions of the degrees of freedom.
%}

[full.G] = gettranferfunc(full.x, Inputdata);
[prodamping.G] = gettranferfunc(prodamping.x, Inputdata);
display(full.G);
display(prodamping.G);

for i = 1:3
    figure();
    opts = bodeoptions;
opts.Grid = 'on';
opts.MagUnits = 'abs';
    bodeplot(full.G(:,i), prodamping.G(:,i), logspace(0,2,1000), opts);
    grid on;
    namefile = ['bodediagram' int2str(i)];
    legend('free damping', 'proportinal damping','Location','west');
    saveas(gcf,namefile,'epsc')
end

%{
 repeat the previous operations for the proportional damping case (and
 compare the results with the ones of the generic damping case). Use the
 modes of the proportional damping case to write an analytical expression
 for the configurations thanks to modal decomposition.
%}



close all
%% SINE SWEEP
% SINE SWEEP - slow
load data_sine_sweep_slow
sinesweep.slow.displacement.x{1} = x1;
sinesweep.slow.displacement.x{2} = x2;
sinesweep.slow.displacement.x{3} = x3;
sinesweep.slow.v = v;
sinesweep.slow.t = t;
clear x1 x2 x3 t v;
sinesweep.slow.estimtf = cell(1,3);
for i = 1:length(sinesweep.slow.displacement.x)
    [sinesweep.slow.estimtf{i}] = getEstimatedtf(i, sinesweep.slow.t, sinesweep.slow.displacement.x{i}, sinesweep.slow.v, 'slow');
    graphicssinesweep(i, 'sinesweepslow', sinesweep.slow.estimtf{i});
end

% SINE SWEEP - fast
load data_sine_sweep_fast.mat

sinesweep.fast.displacement.x{1} = x1;
sinesweep.fast.displacement.x{2} = x2;
sinesweep.fast.displacement.x{3} = x3;
sinesweep.fast.v = v;
sinesweep.fast.t = t;

clear x1 x2 x3 t v;
sinesweep.fast.estimtf= cell(1,3);
for j = 1:length(sinesweep.fast.displacement.x)
    [sinesweep.fast.estimtf{j}] = getEstimatedtf(j, sinesweep.fast.t, sinesweep.fast.displacement.x{j}, sinesweep.fast.v, 'fast');
    graphicssinesweep(j, 'sinesweepfast', sinesweep.fast.estimtf{j});
end

% generate plot comparision
for k = 1:length(sinesweep.fast.estimtf)
    graphicssinesweep(k, 'sinecompare', sinesweep.slow.estimtf{k}, sinesweep.fast.estimtf{k});
end

close all; clear i k;

%% compute fourier trasform

getFouriertrasform(sinesweep.slow.v, sinesweep.slow.t, 'slow')
getFouriertrasform(sinesweep.fast.v, sinesweep.fast.t, 'fast')
close all
