% manisect.m
% plot intersecting attracting and repelling manifolds
function [out r]=manisecttest(a,b,eps,tatt,trep)

%clear all; 
close all;

%a=0.99403; b=.001; eps=.1;
N=5;
%tf=75; 
dt=.03;
dx=500; %Number of x points
%tolin=0.5; tolout=0.1;
tspanatt=0:dt:tatt;
tspanrep=0:dt:trep;
angles=linspace(2*pi/N,2*pi,N);
out=zeros(N,length(tspanatt),3);
r=zeros(N,length(tspanrep));

attracting=colormap(winter(N));
repelling=colormap(autumn(N));

outer=(1+sqrt(3))/2;
shootxatt=outer*cos(angles);
shootyatt=outer*sin(angles);
shootzatt=linspace(.5,.5,N);
shootatt=[shootxatt; shootyatt; shootzatt];
shootxrep=.5*cos(angles);
shootyrep=.5*sin(angles);
shootzrep=linspace(.5,.5,N);
shootrep=[shootxrep; shootyrep; shootzrep];

att=[];
rep=[];

%attracting manifold
for i=1:N
    [tatt,xatt]=ode45(@(tatt,xatt) fwdrotvdp(tatt,xatt,eps,a,b),tspanatt,shootatt(:,i)');
    att=[att; xatt];
    outa(i,:,:)=xatt;
end

rada=(att(:,1).^2+att(:,2).^2).^(1/2);
slicea=find(att(:,3)<.01 & att(:,3)>0 & rada(:)<1 & att(:,1)>0 & att(:,2)>0);

xlin=linspace(min(att(slicea,1)),max(att(slicea,1)),dx);
ylin=linspace(min(att(slicea,2)),max(att(slicea,2)),dx);
[Xatt,Yatt] = meshgrid(xlin,ylin);
Zatt = griddata(att(slicea,1),att(slicea,2),att(slicea,3),Xatt,Yatt,'cubic');

%repelling manifold
for i=1:N
    [trep,xrep]=ode45(@(trep,xrep) backrotvdp(trep,xrep,eps,a,b),tspanrep,shootrep(:,i)');
    rep=[rep; xrep];
    outr(i,:,:)=xrep;
end

radr=(rep(:,1).^2+rep(:,2).^2).^(1/2);
slicer=find(rep(:,3)<.01 & rep(:,3)>0 & radr(:)<1 & rep(:,1)>0 & rep(:,2)>0);

xlin=linspace(min(rep(slicer,1)),max(rep(slicer,1)),dx);
ylin=linspace(min(rep(slicer,2)),max(rep(slicer,2)),dx);
[Xrep,Yrep] = meshgrid(xlin,ylin);
Zrep = griddata(rep(slicer,1),rep(slicer,2),rep(slicer,3),Xrep,Yrep,'cubic');

for j=1:N
    figure(10)
    ra(j,:)=(outa(j,:,1).^2+outa(j,:,2).^2).^(1/2);
    rr(j,:)=(outr(j,:,1).^2+outr(j,:,2).^2).^(1/2);
    plot(ra(j,:),outa(j,:,3),'Color',attracting(j,:));
    hold on;
    plot(rr(j,:),outr(j,:,3),'Color',repelling(j,:));
    hold on;
%   figure(11)
%   plot3(out(j,:,1),out(j,:,2),out(j,:,3),'Color',colors(j,:));
%   hold on;
end

figure(50)

surf(Xatt,Yatt,Zatt,-.0025-Zatt,'EdgeColor','none','FaceColor','interp')
hold on
surf(Xrep,Yrep,Zrep,.0025+Zrep,'EdgeColor','none','FaceColor','interp')
colormap jet

%% forward rotvdp
function xdot=fwdrotvdp(t,x,eps,a,b)
% Differential equations for rotating van der Pol system
% Used by Barry, Benes, Burke, Kaper, and Kramer
% Based on rotating the van der Pol oscillator about the z-axis
% x(1)=x, x(2)=y, x(3)=z
% good parameters: a=0.99403 +/- 1e-5 (maybe 1.5?); b=.001; eps=.1;

xdot(1) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(1) - 5*x(2); 
xdot(2) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(2) + 5*x(1); 
xdot(3) = eps.*(a - ((x(1) - b).^2 + x(2).^2).^(1/2));

xdot=xdot';

%% backward rotvdp
function xdot=backrotvdp(t,x,eps,a,b)
% Differential equations for rotating van der Pol system
% Used by Barry, Benes, Burke, Kaper, and Kramer
% Based on rotating the van der Pol oscillator about the z-axis
% x(1)=x, x(2)=y, x(3)=z
% good parameters: a=0.99403 +/- 1e-5 (maybe 1.5?); b=.001; eps=.1;

xdot(1) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(1) - 5*x(2); 
xdot(2) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(2) + 5*x(1); 
xdot(3) = eps.*(a - ((x(1) - b).^2 + x(2).^2).^(1/2));
xdot=-xdot;

xdot=xdot';