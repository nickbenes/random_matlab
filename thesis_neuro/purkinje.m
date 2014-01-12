%% Purkinje ODE
function xdot=purkinje(t,x)

% Purkinje cell model from Kramer, Traub, and Kopell 2008
% V' = -J - gNa*minf[V]^3*h*(V-50) - gK*n^4*(V+95) - gCa*c^2*(V-125) 
%         - gM*M*(V+95) - gL*(V+70)
% x' = (xinf[V]-x)/taux[V] = ax[V]*(1-x)-bx[V]*x
% x(1)=V; x(2)=h; x(3)=n; x(4)=c; x(5)=M;

% Parameters
global J
gNa=125; gK=10; gCa=1; gM=0.75; gL=2;

% Voltage-dependent steady state and gating values
V=x(1);
minf=1.0/(1.0 + Exp((-V-34.5)/10.0));
ah=(1.0/(1.0 + Exp( (V+59.4)/10.7 ))) / (0.15 + 1.15 / (1.0 + Exp( (V+33.5)/15.0 )));
bh=(1.0 - 1.0 / (1.0 + Exp( (V+59.4)/10.7 ))) / (0.15 + 1.15 / (1.0 + Exp( (V+33.5)/15.0 )));
an=(1.0 / (1.0 + Exp( (-V-29.5)/10.0 ))) / (0.25 + 4.35*Exp(-abs(V+10.0)/10.0));
bn=(1.0 - 1.0 / (1.0 + Exp( (-V-29.5)/10.0 ))) / (0.25 + 4.35*Exp(-abs(V+10.0)/10.0));
aM = 0.02/(1.0 + Exp((-20 - V)/5.));
bM = 0.01*Exp((-43 - V)/18.);

% ODEs
xdot(1) = -J -gNa*minf^3*x(2)*(-50.0 + x(1)) -gK*x(3)^4.0*(95.0 + x(1))...
                -gCa*x(4)^2*(-125 + x(1)) -gM*x(5)*(95 + x(1)) -gL*(70 + x(1));
xdot(2) = ah*(1-x(2)) - bh*x(2);
xdot(3) = an*(1-x(3)) - bn*x(3);
xdot(4) = ((1.6*(1 - x(4)))/(1 + Exp(-0.072*(-5 + x(1)))) - (0.02*x(4)*(8.9 + x(1)))/(-1 + Exp((8.9 + x(1))/5.)));
xdot(5) = aM*(1-x(5)) - bM*x(5);

xdot=xdot';