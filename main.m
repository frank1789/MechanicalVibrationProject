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

% evaluate displacement and store in struct Displacement
Inputdata.Displacement.x1 = x1 * Inputdata.number_count_encoder;
Inputdata.Displacement.x2 = x2 * Inputdata.number_count_encoder;
Inputdata.Displacement.x3 = x3 * Inputdata.number_count_encoder;

% store original data in the previous struct
Inputdata.voltage.v = v;
Inputdata.time.t = t;

% save the data and reinitilaize the workspace
save('mydata.mat', 'Inputdata')
clear 
load mydata

% plot figure
figure();
hold on
plot(Inputdata.time.t, Inputdata.Displacement.x1)
plot(Inputdata.time.t, Inputdata.Displacement.x2)
plot(Inputdata.time.t, Inputdata.Displacement.x3)
grid()
legend('Displacement x1', 'Displacement x2','Displacement x3')
hold off

% compute symbolic equation [G(s)]
[ A ] = symbolicequation();

%% compute a new estimation for the voltage-to-force coefficient

% define gain_ x from encoder
gain_x = Inputdata.number_count_encoder;

%unkonw ?
gain_v = 5.250; 

%unkonw ?
gain_tot = (gain_x/gain_v); % ????

%% search the value for theta 

% assign numerical value
m1 = 1.2;
m2 = 1.3;
m3 = 1.4;

k1 = 800;
k2 = k1;
k3 = 400;

c1 = 3;
c2 = 2;
c3 = 2;

gain = 1;

% costant force parameters

%servo amp gain
k_a = 2;

% Servo motor Torque
k_t = 0.1;

% Motor Pinioon pitch radius inverse
k_mp = 26.25;

% force 
F = zeros(length(Inputdata.time.t), 3);
f = (k_a * k_t * k_mp) * Inputdata.voltage.v;

F(:,1) = f;

%% compute a new estimation for the voltage-to-force coefficient
% syms gv 
% gainvolt = K\[gv; 0; 0];
% 
% F = F * gainvolt;


%% compute transfer function
s = tf('s');

denom = 1/det(A);
num = inv([ m1*s^2 + c1*s + k1,                     -k1,                       0
                           -k1, m2*s^2 + c2*s + k1 + k2,                     -k2
                             0,                     -k2, m3*s^2 + c3*s + k2 + k3]);
% transfer fuction
G = tf(num, denom);
display(G);
[YS,TS,XS] = lsim(G,F,Inputdata.time.t);

% comparision plot
figure();
plot(Inputdata.time.t,YS);
grid on;
legend('tf x1','tf x2','tf x3');

% compare the result with "goodnessOfFit" function
% reference data
ref_displacement = [Inputdata.Displacement.x1 Inputdata.Displacement.x2 Inputdata.Displacement.x3]; 

% cost function Normalized root mean square error, where, ? indicates the
% 2-norm of a vector. fit is a row vector of length N and i = 1,...,N, where
% N is the number of channels.
cost_func = 'NRMSE';

fit = goodnessOfFit(YS,ref_displacement,cost_func);
display(fit);

%% Optimization estimated parameters
%{
% I.C. Initial Conditions
gain-------------------------------------+
stiff k3------------------------------+  |
stiff k2---------------------------+  |  |
stiff k1------------------------+  |  |  |
damper-----------------------+  |  |  |  |
damper--------------------+  |  |  |  |  |
damper-----------------+  |  |  |  |  |  |
mass---------+--+--+   |  |  |  |  |  |  |     
             |  |  |   |  |  |  |  |  |  |
estimated = [m1 m2 m3 c1 c2 c3 k1 k2 k3 gain];
%}
estimated = [0.3 1.4 .13 2 2 2 800 800 400];

%{ 
define lower bound
gain--------------+
damper3---------+ |
damper2-------+ | |
damper1-----+ | | |
mass--+-+-+ | | | |     
%}
LB = [0 0 0 0 0 0 0 0 0];

% upper bound
UB = [2 2 2 5 5 5 800 800 400];


% set initial condition
x0 = estimated;

%% fmincon with boundaries
f = @(x)errormio(x0, A, F, Inputdata);
[xfmincon ] = fmincon(f,x0,[],[],[],[],LB,UB)

%% fminsearch
options = optimset('PlotFcns',@optimplotfval);
f = @(x)errormio(x0, A, F, Inputdata);
[xfminsearch] = fminsearch(f,x0,options)
