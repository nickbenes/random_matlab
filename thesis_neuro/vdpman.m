% vdpman.m
function [out r]=vdpman(a,b,eps,N)

%clear all; 
close all;

%a=0.99403; b=.001; eps=.1;
N=50;
colors=colormap(hsv(N));
tf=75; dt=.01;
tolin=0.5; tolout=0.1;

tspan=0:dt:tf;
%angles=linspace(2*pi/N,2*pi,N);
angles=linspace(0,pi/4,N);
out=zeros(N,length(tspan),3);
r=zeros(N,length(tspan));

shootx=1.5*cos(angles);
shooty=1.5*sin(angles);
shootz=linspace(1,1,N);
shoot=[shootx; shooty; shootz];

for i=1:N
    [t,x]=ode45(@(t,x) rotvdp(t,x,eps,a,b),tspan,shoot(:,i)');
    out(i,:,:)=x;
end


for j=1:N
    figure(10)
    r(j,:)=(out(j,:,1).^2+out(j,:,2).^2).^(1/2);
    plot(r(j,:),out(j,:,3),'Color',colors(j,:));
    hold on;
    figure(11)
    plot3(out(j,:,1),out(j,:,2),out(j,:,3),'Color',colors(j,:));
    hold on;
    figure(12)
    slice=find(r(j,:)<1+tolout & r(j,:)>1-tolin & out(j,:,3)<1 & out(j,:,1)>0 & out(j,:,2)>0);
    %firsts=find(out(j,slice(:),1)
    plot3(out(j,slice,1),out(j,slice,2),out(j,slice,3),'s',...
        'MarkerEdgeColor',colors(j,:),'MarkerFaceColor',colors(j,:),...
        'MarkerSize',3.5);
    grid on; axis equal;
    hold on;
end

%% rotvdp
function xdot=rotvdp(t,x,eps,a,b)
% Differential equations for rotating van der Pol system
% Used by Barry, Benes, Burke, Kaper, and Kramer
% Based on rotating the van der Pol oscillator about the z-axis
% x(1)=x, x(2)=y, x(3)=z
% good parameters: a=0.99403 +/- 1e-5 (maybe 1.5?); b=.001; eps=.1;

xdot(1) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(1) - 5*x(2); 
xdot(2) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(2) + 5*x(1); 
xdot(3) = eps.*(a - ((x(1) - b).^2 + x(2).^2).^(1/2));

xdot=xdot';