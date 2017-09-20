close all;
clear;
clc;

% Load folder data
addpath('data_3DOFsystem');
file = dir(fullfile('data_3DOFsystem', '*.mat'));
load data_steps

% number count per encoder revolution is costant
number_count_encoder = 0.0706/16000;

% evaluate displacement and store in struct Displacement
Displacement.x1 = x1 * number_count_encoder;
Displacement.x2 = x2 * number_count_encoder;
Displacement.x3 = x3 * number_count_encoder;

% plot figure
figure();
hold on
plot(t, Displacement.x1)
plot(t, Displacement.x2)
plot(t, Displacement.x3)
grid()
legend('x1', 'x2','x3')
hold off

% compute symbolic equation
[ A ] = symbolicequation();

%% compute a new estimation for the voltage-to-force coefficient

% define gain_ x from encoder
gain_x = number_count_encoder;

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

% costant force parameters

%servo amp gain
k_a = 2;

% Servo motor Torque
k_t = 0.1;

% Motor Pinioon pitch radius inverse
k_mp = 26.25;

% force 
F = zeros(length(t), 3);
f = (k_a * k_t * k_mp) * v;

F(:,1) = f;

%% compute a new estimation for the voltage-to-force coefficient
% syms gv 
% gainvolt = K\[gv; 0; 0];
% 
% F = F * gainvolt;


% compute transfer function
s = tf('s');

denom = 1/det(A);
num = inv([ m1*s^2 + c1*s + k1,                     -k1,                       0
                           -k1, m2*s^2 + c2*s + k1 + k2,                     -k2
                             0,                     -k2, m3*s^2 + c3*s + k2 + k3]);
% transfer fuction
G = tf(num, denom);
display(G);
[YS,TS,XS] = lsim(G,F,t);

% comparision plot
figure();
plot(t,YS);
grid on;
legend('tf x1','tf x2','tf x3');

% compare the result with "goodnessOfFit" function
% reference data
ref_displacement = [Displacement.x1 Displacement.x2 Displacement.x3]; 

% cost function Normalized root mean square error, where, ? indicates the
%2-norm of a vector. fit is a row vector of length N and i = 1,...,N, where
%N is the number of channels.
cost_func = 'NRMSE';

fit = goodnessOfFit(YS,ref_displacement,cost_func);

%%
% sys = tfest(x1,3);
% data1 = [ Displacement.x1 Displacement.x2 Displacement.x3];
% data = iddata(data1,v,0.005)
%compare(data,YS)


%% Obtain the measured output.

% load iddata1 z1
% yref = z1.y;
%%
%z1 is an iddata object containing measured input/output data. z1.y is the measured output.

%Obtain the estimated output.

% sys = tfest(Displacment.x1,2);
% y_sim = sim(sys,z1(:,[],:));
%sys is a second-order transfer function estimated using the measured input/output data. y is the output estimated using sys and the measured input.

%Calculate the goodness of the fit between the measured and estimated outputs.


y = y_sim.y;
fit = goodnessOfFit(y,yref,cost_func);
%The goodness of fit is calculated using the normalized root mean square error as the cost function.

%Alternatively, you can use compare to calculate the goodness of fit:

opt = compareOptions('InitialCondition','z');
compare(z1,sys,opt);

%% Dummy values for mech quantities (AW trial 1)

P = [1.5 1.4 1.2 0.01 0.001 3 2 2 6.3];

%% Define dummy system with previous dummy mech quantities

%[A, B, C, D] = ss_mech_linear(P);

sys_dum = ss(A,B,C,D);


%%
options = optimset('PlotFcns',@optimplotfval);
%fun = @(x)100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
fun = @(v, t, Displacement, P)estimedparameters(v, t, Displacement, P);
%x0 = [-1.2,1];
%[x, fval] = fminsearch(fun,options)

%%


%rng default % for reproducibility
% tdata = 0:0.1:10;
% ydata = 40*exp(-0.5*tdata) + randn(size(tdata));
% fun = @(x)sseval(x,tdata,ydata);
% 
% x0 = rand(2,1);
% bestx = fminsearch(fun,x0)
% 
% 
% A = bestx(1);
% lambda = bestx(2);
% yfit = A*exp(-lambda*tdata);
% plot(tdata,ydata,'*');
% hold on
% plot(tdata,yfit,'r');
% xlabel('tdata')
% ylabel('Response Data and Curve')
% title('Data and Best Fitting Exponential Curve')
% legend('Data','Fitted Curve')
% hold off
