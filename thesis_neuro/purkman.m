function purkman(Jin,N)

% Plot the intersection of the attracting and repelling manifolds of the
% purkinje model as they cross a cylinder.
% The cylinder is determined by the V and h values of the system's limit
% cycle at the saddle-node of limit-cycle bifurcation value of M (as a
% parameter in the fast subsystem).  I.e. the V and h values of the
% non-hyperbolic limit cycle.

%% Initialization
global J
J=Jin;

% x(1)=V; x(2)=h; x(3)=n; x(4)=c; x(5)=M;
N=N+1;
angles=linspace(-pi,pi,50);
tspan=0:.1:5;
options=odeset('RelTol',1e-5);
att=[]; rep=[];
attxing=zeros(N,2);
repxing=zeros(N,2);
outcolors=colormap(winter(N));
incolors=colormap(autumn(N));

%% Initial conditions

% Load data for attracting and repelling rings from XPP data file
% Data is an m by 5 matrix with the following columns:
% 1=time, 2=V, 3=h, 4=n, 5=c
attring=load('stablepurk.dat');
%repring=load('unstablepurk.dat');

% Choose N representative initial conditions to follow
atttimes=linspace(attring(1,1),attring(end,1),N);
attic=interp1q(attring(:,1),attring(:,2:5),atttimes);

%reptimes=linspace(repring(1,1),repring(end,1),N);
%repic=interp1q(repring(:,1),repring(:,2:5),reptimes);

%% Convert non-hyperbolic l.c. data to polar coordinates for comparison
% Load data for non-hyperbolic limit cycle
snring=load('snpurk.dat');

% Remember that we only care about the V-h ring

% Find a point in the middle of the l.c.
mid(1,1:2)=(max(snring(:,2:3))+min(snring(:,2:3)))/2;

% Finda radial coordinates of the l.c.
snr=sqrt((snring(:,2)-mid(1)).^2+(snring(:,3)-mid(2)).^2);
sntheta=atan2(snring(:,2)-mid(1),snring(:,3)-mid(2));
sn=[snr sntheta];    

%% Integration

for i=1:N %integrate for each init cond

% Attracting limit cycle
while cross(mid,sn,attic(i,1:2)) < 0
    % Integrate for a while
    [t,x]=ode45(@purkinje,tspan,attic(i,:),options);
    att=[att; x];
    
    % Check to make sure V hasn't dropped to rest
    if att(end,1)<-50
        error('Attracting manifold dropped to rest');
    end
end
% x now contains last 5s of data, including time of crossing

% Find M value and angle at crossing
while length(x(:,1))>2 %perform binary search
    xmid=floor(length(x(:,1))/2);
    if cross(mid,sn,x(xmid,1:2))>0 %crossing is in first half
        x=x(1:xmid,:);
    else %crossing is in second half
        x=x(xmid:end,:);
    end
end
% x now contains data just before and just after crossing
% take weighted average of M and theta
attxing(i,:)=interp1q([cross(mid,sn,x(1,1:2));cross(mid,sn,x(1,1:2))],...
             [x(1,5) atan2(x(1,1),x(1,2));x(2,5) atan2(x(2,1),x(2,2))],0);
  
% Repelling limit cycle


end

%% Plots

% Plot M vs angle at crossing for each manifold

% Plot trajectories in 3D


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

%% Function for determining if cylinder has been crossed
function dist = cross(cent,cyl,pt)

% Determine angle spacing from size of cyl

% Determine angle of pt

% Determine distance of pt from cent and compare with interp value from cyl
