clear;
clc;
close all;
global A w;

PRINT = true;
% PRINT = false;

%% 2nd order system

%Simulation time
tfinal = 500;

%Reference
w = [2 3];
A = zeros(1,length(w));
A(1) = 10;
A(2) = 25;

%Observer
A0 = tf([1 1],[1]);

%Monopoli multiplier
Li = [1 1.5];

%--------------- First set of parameters ----------------

%P
Np_1 = [1];
Dp_1 = [1 -2 1];
kp_1 = 1;
P_1 = tf(kp_1*Np_1,Dp_1);
gP_1 = length(Dp_1) - 1; % system order

%Pm
Nm_1 = [1];
Dm_1 = [1 3 2];
km_1 = 1;
Pm_1 = tf(km_1*Nm_1,Dm_1);
gPm_1 = length(Dm_1)-1;

%Initial conditions
y0_1  = 0;

%Adaptation gain
gamma_1 = 20;

%--------------- Second set of parameters ----------------

%P
Np_2 = [1];
Dp_2 = [1 7 12];
kp_2 = 2;
P_2 = tf(kp_2*Np_2,Dp_2);
gP_2 = length(Dp_2) - 1; % system order

%Pm
Nm_2 = [1];
Dm_2 = [1 5 6];
km_2 = 1;
Pm_2 = tf(km_2*Nm_2,Dm_2);
gPm_2 = length(Dm_2)-1;

%Initial conditions
y0_2  = 10;

%Adaptation gain
gamma_2 = 1;

run sim_mrac.m;

%% 3rd order system

%Simulation time
tfinal = 500;

%Reference
w = [2 3 5];
A = zeros(1,length(w));
A(1) = 10;
A(2) = 25;
A(3) = 10;

%Observer
A0 = tf([1 1],[1]);

%Monopoli multiplier
Li = [1 1.2];

%--------------- First set of parameters ----------------

%P
Np_1 = [1 1];
Dp_1 = [1 -6 12 -8];
kp_1 = 1;
P_1 = tf(kp_1*Np_1,Dp_1);
gP_1 = length(Dp_1) - 1; % system order

%Pm
Nm_1 = [1 0.8];
Dm_1 = [1 3 2.75 0.75];
km_1 = 1;
Pm_1 = tf(km_1*Nm_1,Dm_1);
gPm_1 = length(Dm_1)-1;

%Initial conditions
y0_1  = 0;

%Adaptation gain
gamma_1 = 20;

%--------------- Second set of parameters ----------------

%P
Np_2 = [1 2];
Dp_2 = [1 6 11 6];
kp_2 = 2;
P_2 = tf(kp_2*Np_2,Dp_2);
gP_2 = length(Dp_2) - 1; % system order

%Pm (not SPR!)
Nm_2 = [1 1.5];
Dm_2 = [1 5.5 8.5 3];
km_2 = 1;
Pm_2 = tf(km_2*Nm_2,Dm_2);
gPm_2 = length(Dm_2)-1;

%Initial conditions
y0_2  = 10;

%Adaptation gain
gamma_2 = 1;

run sim_mrac.m;
