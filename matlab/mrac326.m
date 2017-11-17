%======================================================================
%
%  COE-835  Controle adaptativo
%
%  Script para simular o exemplo
%
%  MRAC  : n  = 2     First order plant
%          n* = 1     Relative degree
%          np = 4     Adaptive parameters
%
%======================================================================

function dx=mrac326(t,x)

global Ay By Aym Bym Auf Buf Ayf Byf Ao Bo kp gamma w A Pg Pmg Lg thetag;

y      = x(1:Pg);
ym     = x(Pg+1:Pg+Pmg);
uf     = x(Pg+Pmg+1:Pg+Pmg+Lg);
yf     = x(Pg+Pmg+Lg+1:Pg+Pmg+Lg+Lg);
theta  = x(Pg+Pmg+Lg+Lg+1:Pg+Pmg+Lg+Lg+thetag);
zeta = x(Pg+Pmg+Lg+Lg+thetag+1:end);

%--------------------------
r = A(1)*sin(w(1)*t) + A(2)*sin(w(2)*t) + A(3)*sin(w(3)*t);

omega = [uf' y(1) yf' r]';
e  = y(1) - ym(1);
dtheta = -sign(kp)*gamma*zeta*e;

u = theta'*omega + dtheta'*zeta;

%------- Cálculo de y --------
dy = Ay*y + By*u;

%------- Cálculo de ym --------
dym = Aym*ym + Bym*r;

%------- Cálculo de uf --------
duf = Auf*uf + Buf*u;

%------- Cálculo de yf --------
dyf = Ayf*yf + Byf*y(1);

%------- Cálculo de zeta --------
dzeta = Ao*zeta + Bo*omega;

%--------------------------
dx = [dy' dym' duf' dyf' dtheta' dzeta']';    %Translation

%---------------------------
