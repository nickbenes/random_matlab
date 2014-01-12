% manisect.m
% plot intersecting attracting and repelling manifolds
function [out r]=manisect(a,b,eps,tatt,trep)

%clear all; 
close all;

%a=0.99403; b=.001; eps=.1;
N=3;
%tf=75; 
dt=.03;
%tolin=0.5; tolout=0.1;
tspanatt=0:dt:tatt;
tspanrep=0:dt:trep;
angles=linspace(2*pi/N,2*pi,N);
out=zeros(N,length(tspanatt),3);
r=zeros(N,length(tspanrep));

attracting=colormap(winter(length(tspanatt)));
repelling=colormap(autumn(length(tspanrep)));

shootxatt=1.5*cos(angles);
shootyatt=1.5*sin(angles);
shootzatt=linspace(1,1,N);
shootatt=[shootxatt; shootyatt; shootzatt];
shootxrep=.5*cos(angles);
shootyrep=.5*sin(angles);
shootzrep=linspace(.5,.5,N);
shootrep=[shootxrep; shootyrep; shootzrep];

%attracting manifold
for i=1:N
    [tatt,xatt]=ode45(@(tatt,xatt) fwdrotvdp(tatt,xatt,eps,a,b),tspanatt,shootatt(:,i)');
    att(i,:,:)=xatt;
end

%repelling manifold
for i=1:N
    [trep,xrep]=ode45(@(trep,xrep) backrotvdp(trep,xrep,eps,a,b),tspanrep,shootrep(:,i)');
    rep(i,:,:)=xrep;
end

for j=1:N
    %figure(10)
    %r(j,:)=(out(j,:,1).^2+out(j,:,2).^2).^(1/2);
    %plot(r(j,:),out(j,:,3),'Color',colors(j,:));
    %hold on;
    %figure(11)
    for k=1:(length(tatt)-1)
        plot3([att(j,k,1);att(j,k+1,1)],[att(j,k,2);att(j,k+1,2)],...
            [att(j,k,3);att(j,k+1,3)],'Color',attracting(k,:));
        hold on;
    end
    for k=1:(length(trep)-1)
        plot3([rep(j,k,1);rep(j,k+1,1)],[rep(j,k,2);rep(j,k+1,2)],...
            [rep(j,k,3);rep(j,k+1,3)],'Color',repelling(k,:));
        hold on;
    end
    %figure(12)
    %slice=find(r(j,:)<1+tolout & r(j,:)>1-tolin & out(j,:,3)<1 & out(j,:,1)>0 & out(j,:,2)>0);
    %firsts=find(out(j,slice(:),1)
    %plot3(out(j,slice,1),out(j,slice,2),out(j,slice,3),'s',...
    %    'MarkerEdgeColor',colors(j,:),'MarkerFaceColor',colors(j,:),...
    %    'MarkerSize',3.5);
    %grid on; axis equal;
    %hold on;
end

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