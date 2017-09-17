close all;
clear;
clc;

% Load folder data
addpath('data_3DOFsystem');
file = dir(fullfile('data_3DOFsystem', '*.mat'));

load data_steps
x1;
displacement = [length(x1)];
displacement(1)=0;
%% convertion encoder counts to displacement
    for i = 2:length(x1)
       displacement(i,1) = 0.0706 * ((x1(i) - x1(i-1))/16000);
    end


figure();
plot(t, displacement)
grid()

%% simbolic matrix
syms m1 m2 m3 c1 c2 c3 k1 k2 k3 x1__ x2__ x3__ x1_ x2_ x3_ x1 x2 x3

M = [m1 0 0; 0 m2 0; 0 0 m3]
C = [c1 0 0; 0 c2 0; 0 0 c3]
K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3]

% compute matrix [G(s)]
syms s
A = M * s^2 + C*s + K

% compute det([G(s)])
det = det(A)

% compute inverse inv([G(s)])
num = inv(A)


%% search the value for theta 

% assign numerical value
m1 = 1.5;
m2 = 1.5;
m3 = 1.5;

k1 = 800;
k2 = k1;
k3 = 400;

c1 = 4;
c2 = 4;
c3 = 4;

% costant force parameters
k_a = 2;
k_t = 0.1;
k_mp =1/26.25;

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
lsim(G,F,t);


% compare plot
o = lsim(G,F,t);
figure(1)
hold on
plot(t,o(:,1))

% force are wrong?
%%


