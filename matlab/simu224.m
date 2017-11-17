%----------------------------------------------------------------------
%
%  COE-835  Controle adaptativo
%
%  Script para simular exemplo 
%
%  MRAC  : n  = 1     First order plant
%          n* = 1     Relative degree
%          np = 2     Adaptive parameters
%
%                                                        Ramon R. Costa
%                                                        30/abr/13, Rio
%----------------------------------------------------------------------
clear;
clc;

disp('-------------------------------')
disp('Script para simular o trabalho 6')
disp(' ')
disp('Caso: Planta ............. n = 2')
disp('      Grau relativo ..... n* = 2')
disp('      Parâmetros ........ np = 4')
disp(' ')
disp('Algoritmo: MRAC direto')
disp(' ')
disp('-------------------------------')

global Ay By Aym Bym Auf Buf Ayf Byf Ao Bo kp gamma w A Pg Pmg Lg thetag;


%% r
w1 = 2;
w2 = 3;
w = [w1 w2];

a1 = 10;
a2 = 25;
A = [a1 a2];

%% H, Hm
Z = [1];
P = [1 -2 1];
kp = 1;

Zm = [1];
Pm = [1 3 2];
km = 1;

H = tf(kp*Z,P);
Hm = tf(km*Zm,Pm);
A0 = tf([1 1],[1]);
Li = [1 1.5];

%% Init
[theta_1, theta_n, theta_2, theta_2n, L] = find2DOFparameters(H,Hm,A0); 
thetag = length(theta_1)+length(theta_2)+2;
gamma = 100*eye(thetag);

%% Sistemas
% Planta
ss_H = canon(ss(H), 'companion');
Ay = ss_H.A';
By = ss_H.C'; % Forma canônica observável

% Modelo
ss_Hm = canon(ss(Hm), 'companion');
Aym = ss_Hm.A';
Bym = ss_Hm.C'; % Forma canônica observável

% Filtro u
ss_uf = canon(ss(tf(1,L)), 'companion');
Auf = ss_uf.A';
Buf = ss_uf.C'; % Forma canônica observável

% Filtro y
ss_yf = canon(ss(tf(1,L)), 'companion');
Ayf = ss_yf.A';
Byf = ss_yf.C'; % Forma canônica observável

% Filtro omega
ss_o = canon(ss(tf(1,Li)), 'companion');
Ao = ss_o.A';
Bo = ss_o.C'; % Forma canônica observável

%% Init
Pg = length(P)-1;
Pmg = length(Pm)-1;
Lg = length(L)-1;

y0  = zeros(Pg,1);
ym0 = zeros(Pmg,1);
uf0 = zeros(Lg,1);
yf0 = zeros(Lg,1);
theta0 = zeros(thetag,1);
zeta0 = [uf0' y0(1) yf0' 0]';

clf;
tfinal = 200;

init = [y0' ym0' uf0' yf0' theta0' zeta0']';


%% Plots
[T,X] = ode23s('mrac224',tfinal,init,'');

y      = X(:,1);
ym     = X(:,Pg+1);
theta =  X(:,Pg+Pmg+Lg+Lg+1:Pg+Pmg+Lg+Lg+thetag);

e =  y - ym;
r = A(1)*sin(w(1).*T) + A(2)*sin(w(2).*T);

%Set matlab interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

figure(1)
clf
plot(T,e);grid;shg
legend('e','Location','SouthEast')
title('$\theta(0)=0$')
print -depsc2 en2t0

figure(2)
clf
plot(T,theta(:,1)-theta_1,T,theta(:,2)-theta_n,T,theta(:,3)-theta_2,T,theta(:,4)-theta_2n);grid;shg
legend('$\tilde{\theta_1}$','$\tilde{\theta_n}$','$\tilde{\theta_2}$','$\tilde{\theta_{2n}}$','Location','SouthEast')
title('$\theta(0)=\textrm{rand}$')
print -depsc2 tiln2t0
%---------------------------------------------------------------------

