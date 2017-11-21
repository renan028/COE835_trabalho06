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
function dx=mrac(t,x)

global Ay By Aym Bym Auf Buf Ayf Byf Ao Bo kp gamma w A gP gPm;

y      = x(1:gP);
ym     = x(gP+1:gP+gPm);
uf     = x(gP+gPm+1:2*gP+gPm-1);
yf     = x(2*gP+gPm:3*gP+gPm-2);
theta  = x(3*gP+gPm-1:5*gP+gPm-2);
zeta   = x(5*gP+gPm-1:end);

%--------------------------
r = 0;
for i=1:gP
    r = r + A(i)*sin(w(i)*t);
end

omega = [uf' y(1) yf' r]';
e  = y(1) - ym(1);
dtheta = -sign(kp)*gamma*zeta*e;

u = theta'*omega + dtheta'*zeta;

%------- Calculo de y --------
dy = Ay*y + By*u;

%------- Calculo de ym --------
dym = Aym*ym + Bym*r;

%------- Calculo de uf --------
duf = Auf*uf + Buf*u;

%------- Calculo de yf --------
dyf = Ayf*yf + Byf*y(1);

%------- Calculo de zeta --------
dzeta = Ao*zeta + Bo*omega;

%--------------------------
dx = [dy' dym' duf' dyf' dtheta' dzeta']';    %Translation

%---------------------------