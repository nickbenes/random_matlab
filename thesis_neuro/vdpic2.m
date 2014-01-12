% vdpic.m
% Follows N initial conditions around a ring on the attracting branch of a
% singular manifold.  ICs evolve according to a modified version of the van
% der Pol oscillator.  Two figures are plotted: Figure 10 plots radius vs.
% height of the trajectories; Figure 11 plots initial condition vs. max
% height, along with a reference height of 1.
%
% a,b, and epsilon are as in the rotating van der pol model, and N is the
% number of initial conditions to follow (evenly distributed around a 
% circle). "out" is an N by T by 3 array, where T is the length of the time
% vector; e.g. "out(i,:,1)" is the time series for the x coordinate of the 
% i-th initial condition, "out(i,:,2)" is the corresponding y coordinate, 
% etc. "r" is an N by T array of the radial distances, and "sup" is a 
% vector of length N storing the maximum of each initial condition.
% 
% Nick Benes, 18 Dec 2009

function [out r sup]=vdpic2(a,b,eps,N)

%clear all; 
%close all;

N=N+1;
%a=0.99403; b=.001; eps=.1;
%N=50;
outcolors=colormap(winter(N));
incolors=colormap(autumn(N));

tf=40; dt=.005;

tspan=0:dt:tf;
angles=linspace(0,2*pi,N);
out=zeros(N,length(tspan),3);
in=zeros(N,length(tspan),3);
r=zeros(N,length(tspan));
rin=zeros(N,length(tspan));
theta=zeros(N,length(tspan));
thetain=zeros(N,length(tspan));
sup=zeros(N,1);

%attracting branch
outer=(1+sqrt(3))/2;
shootx=outer*cos(angles);
shooty=outer*sin(angles);
shootz=linspace(.5,.5,N);
shoot=[shootx; shooty; shootz];

%repelling branch
shootinx=.5*cos(angles);
shootiny=.5*sin(angles);
shootinz=linspace(.5,.5,N);
shootin=[shootinx; shootiny; shootinz];


for i=1:N
    [t,x]=ode45(@(t,x) rotvdp(t,x,eps,a,b),tspan,shoot(:,i)');
    out(i,:,:)=x;
    [inf ind]=min(x(:,3));
    sup(i)=max(x(ind:end,3));
end

for i=1:N
    [t,x]=ode45(@(t,x) backrotvdp(t,x,eps,a,b),tspan,shootin(:,i)');
    in(i,:,:)=x;
%    [inf ind]=min(x(:,3));
%    sup(i)=max(x(ind:end,3));
end

% determine polar coordinates and plot r vs. z for each trajectory
for j=1:N
    figure(10)
    r(j,:)=(out(j,:,1).^2+out(j,:,2).^2).^(1/2);
    rin(j,:)=(in(j,:,1).^2+in(j,:,2).^2).^(1/2);
    theta(j,:)=atan2(out(j,:,2),out(j,:,1));
    thetain(j,:)=atan2(in(j,:,2),in(j,:,1));
    plot(r(j,:),out(j,:,3),'Color',outcolors(j,:));
    hold on;
    plot(rin(j,:),in(j,:,3),'Color',incolors(j,:));
    hold on;
end

% 3d plot of i.c. vs. max height
figure(11)
plot3(cos(angles),sin(angles),sup);
hold on;
plot3(cos(angles),sin(angles),ones(1,N),'Color','k');

% unfolded plot of i.c. vs. max height
figure(12)
plot(angles,sup);

%% plot in cylindrical coordinates of attracting & repelling manifolds

% only look at portion of manifold with r in intersection region and on
% first pass
cut=zeros(N,1);
cutin=zeros(N,1);
for i=1:N % find first pass
    cut(i)=find(r(i,:)<1-3*eps/2,1);
    cutin(i)=find(rin(i,:)>1-eps/2,1);
end
rcut=r(:,1:max(cut)); % limit to first pass
rcutin=rin(:,1:max(cutin));

sec=find(rcut>1-3*eps/2 & rcut<1-eps/2); % intersection region
secin=find(rcutin>1-3*eps/2 & rcutin<1-eps/2); % intersection region


% plot data on fixed grid
zout=out(:,:,3);
zin=in(:,:,3);
xlin=linspace(1-3*eps/2,1-eps/2,35);
ylin=linspace(-pi,pi,35);
[X,Y]=meshgrid(xlin,ylin);
Zout=griddata(r(sec),theta(sec),zout(sec),X,Y,'cubic');
Zin=griddata(rin(secin),thetain(secin),zin(secin),X,Y,'cubic');

figure(13)
surf(X,Y,Zout,Zout-1e-2)
hold on
surf(X,Y,Zin,Zin+1e-2)
colormap(jet)

%% rotvdp
function xdot=rotvdp(t,x,eps,a,b)
% Differential equations for rotating van der Pol system
% Used by Barry, Benes, Burke, Kaper, and Kramer
% Based on rotating the van der Pol oscillator about the z-axis
% x(1)=x, x(2)=y, x(3)=z
% good parameters: a=0.99403 +/- 1e-5 (maybe 1.5?); b=.001; eps=.1;

xdot(1) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(1) - x(3).*x(2); 
xdot(2) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(2) + x(3).*x(1); 
xdot(3) = eps.*(a - ((x(1) - b).^2 + x(2).^2).^(1/2));

xdot=xdot';

%% backward rotvdp
function xdot=backrotvdp(t,x,eps,a,b)
% Differential equations for rotating van der Pol system
% Used by Barry, Benes, Burke, Kaper, and Kramer
% Based on rotating the van der Pol oscillator about the z-axis
% x(1)=x, x(2)=y, x(3)=z
% good parameters: a=0.99403 +/- 1e-5 (maybe 1.5?); b=.001; eps=.1;

xdot(1) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(1) - x(3).*x(2); 
xdot(2) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(2) + x(3).*x(1); 
xdot(3) = eps.*(a - ((x(1) - b).^2 + x(2).^2).^(1/2));
xdot=-xdot;

xdot=xdot';