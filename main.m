close all;
clear;
clc;

% Load folder data
addpath('data_3DOFsystem');
file = dir(fullfile('data_3DOFsystem', '*.mat'));

load data_steps
x1;

% instanziate dipslacement 
displacement = zeros(length(x1), 1);

figure();
plot(t, x1*(0.0706/16000))
grid()

%% simbolic matrix
syms m1 m2 m3 c1 c2 c3 k1 k2 k3 x1__ x2__ x3__ x1_ x2_ x3_ sx1 sx2 sx3 c12 c23

% mass matrix [M]
M = [m1 0 0; 0 m2 0; 0 0 m3];

% damping matric [C]
C = [c1 0 0; 0 c2 0; 0 0 c3];

% stifness matrix [K]
K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3];

% compute matrix [G(s)]
syms s
A = M * s^2 + C*s + K;

% compute det([G(s)])
det = det(A);

% compute inverse inv([G(s)])
num = inv(A);

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

% compute transfer function
s = tf('s');

den = 1/(k1*k2*k3 + c1*k1*k2*s + c1*k1*k3*s + c2*k1*k2*s + c1*k2*k3*s + c2*k1*k3*s + c3*k1*k2*s + c1*c2*c3*s^3 + c1*c2*k2*s^2 + c1*c3*k1*s^2 + c1*c2*k3*s^2 + c1*c3*k2*s^2 + c2*c3*k1*s^2 + c1*c2*m3*s^4 + c1*c3*m2*s^4 + c2*c3*m1*s^4 + c1*k1*m3*s^3 + c1*k2*m2*s^3 + c2*k2*m1*s^3 + c3*k1*m1*s^3 + c1*k2*m3*s^3 + c1*k3*m2*s^3 + c2*k1*m3*s^3 + c2*k3*m1*s^3 + c3*k1*m2*s^3 + c3*k2*m1*s^3 + c1*m2*m3*s^5 + c2*m1*m3*s^5 + c3*m1*m2*s^5 + k1*k2*m1*s^2 + k1*k2*m2*s^2 + k1*k3*m1*s^2 + k1*k2*m3*s^2 + k1*k3*m2*s^2 + k2*k3*m1*s^2 + k1*m1*m3*s^4 + k2*m1*m2*s^4 + k1*m2*m3*s^4 + k2*m1*m3*s^4 + k3*m1*m2*s^4 + m1*m2*m3*s^6);
num = inv([[ m1*s^2 + c1*s + k1,                     -k1,                       0]
[                -k1, m2*s^2 + c2*s + k1 + k2,                     -k2]
[                  0,                     -k2, m3*s^2 + c3*s + k2 + k3]]);

% transfer fuction
G = tf(num,den)
figure();
[YS,TS,XS] = lsim(G,F,t);


% compare plot
figure(1)
hold on
plot(t,YS)
hold off
figure();
plot(TS, YS)
% force are wrong?
%%

%% Dummy values for mech quantities (AW trial 1)

P = [1.2 1.3 1.4 0 0 3 2 2 800 800 800 600 0 0 ];

%% Define dummy system with previous dummy mech quantities

[A, B, C, D] = ss_mech_linear(P);

sys_dum = ss(A,B,C,D);


%%
options = optimset('PlotFcns',@optimplotfval);
fun = @(x)100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
x0 = [-1.2,1];
x = fminsearch(fun,x0,options)

%%


rng default % for reproducibility
tdata = 0:0.1:10;
ydata = 40*exp(-0.5*tdata) + randn(size(tdata));
fun = @(x)sseval(x,tdata,ydata);

x0 = rand(2,1);
bestx = fminsearch(fun,x0)


A = bestx(1);
lambda = bestx(2);
yfit = A*exp(-lambda*tdata);
plot(tdata,ydata,'*');
hold on
plot(tdata,yfit,'r');
xlabel('tdata')
ylabel('Response Data and Curve')
title('Data and Best Fitting Exponential Curve')
legend('Data','Fitted Curve')
hold off

function sse = sseval(x,tdata,ydata)
A = x(1);
lambda = x(2);
sse = sum((ydata - A*exp(-lambda*tdata)).^2);
end