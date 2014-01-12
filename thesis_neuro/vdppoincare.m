% Follows N initial conditions around a ring on the attracting and 
% repelling branches of a singular manifold.  ICs evolve according to a 
% modified version of the van der Pol oscillator.
% 
% We plot the intersection of the manifolds with a cylinder of radius 1
%
% a,b, and epsilon are as in the rotating van der pol model, and N is the
% number of initial conditions to follow (evenly distributed around a 
% circle). 
% 
% Nick Benes, 13 Jan 2010
% Revised with improved error tolerances 20 Jan 2010

function vdppoincare(a,b,eps,N)
%figure(20); figure(30);
%close([20 30]);

N=N+1;
% Typical values: a=0.99403; b=.001; eps=.1; N=50;

tf=25; dt=.001;
tspan=0:dt:tf;
tspanplus=0:dt:5;
options=odeset('RelTol',6e-6);
outcolors=colormap(winter(N));
incolors=colormap(autumn(N));


% Establishing initial conditions and preallocating
angles=linspace(-pi,pi,N);
att=zeros(N,length(tspan),3);
%rep=zeros(N,2*length(tspan),3);
rep=zeros(N,length(tspan),3);
attplot=zeros(N,2);
repplot=zeros(N,2);
r=zeros(N,length(tspan));
%rin=zeros(N,2*length(tspan));
rin=zeros(N,length(tspan));

outer=(1+sqrt(3))/2;
shootxatt=outer*cos(angles);
shootyatt=outer*sin(angles);
shootzatt=linspace(.5,.5,N);
shootatt=[shootxatt; shootyatt; shootzatt];
shootxrep=.5*cos(angles);
shootyrep=.5*sin(angles);
shootzrep=linspace(.5,.5,N);
shootrep=[shootxrep; shootyrep; shootzrep];

% Attracting manifold
for i=1:N
    [t,x]=ode45(@(t,x) rotvdp(t,x,eps,a,b),tspan,shootatt(:,i)',options); %Flow fwd
    r(i,1:size(x,1))=(x(:,1).^2+x(:,2).^2).^(1/2); %compute radial coord
    cross=find(r(i,:)<1,1); %find first crossing of r=1 cylinder
    while isempty(cross)
        [t,x2]=ode45(@(t,x) rotvdp(t,x,eps,a,b),tspanplus,x(end,:),options); %Flow fwd
        while size(r,2) < size(x,1)+size(x2,1)
            r=[r ones(N,length(tspanplus))];
            att=[att zeros(N,length(tspanplus),3)];
        end
        r(i,1+size(x,1):size(x,1)+size(x2,1))=(x2(:,1).^2+x2(:,2).^2).^(1/2); %compute radial coord
        x=[x; x2];
        cross=find(r(i,:)<1,1); %find first crossing of r=1 cylinder
    end
    att(i,1:size(x,1),:)=x; 
    pt=x(cross,:).*(1-r(i,cross))./(r(i,cross-1)-r(i,cross))...
        +x(cross-1,:).*(r(i,cross-1)-1)./(r(i,cross-1)-r(i,cross));
    % pt is the linear interpolation of where the trajectory crosses r=1
    % It is an average of the x coordinates weighted by their dist from r=1
    attplot(i,1)=atan2(pt(2),pt(1)); %for plotting
    attplot(i,2)=pt(3);
end

clear x
clear cross

% Repelling manifold
for i=1:N
    % Flow back
    [t,x]=ode45(@(t,x) backrotvdp(t,x,eps,a,b),tspan,shootrep(:,i)',options); 
    rin(i,1:size(x,1))=(x(:,1).^2+x(:,2).^2).^(1/2); %compute radial coord
    cross=find(rin(i,:)>1,1); %find first crossing of r=1 cylinder
    while isempty(cross)
        [t,x2]=ode45(@(t,x) backrotvdp(t,x,eps,a,b),tspanplus,x(end,:),options); %Flow fwd
        while size(rin,2) < size(x,1)+size(x2,1)
            rin=[rin ones(N,length(tspanplus))];
            rep=[rep zeros(N,length(tspanplus),3)];
        end
        rin(i,size(x,1)+1:size(x,1)+size(x2,1))=(x2(:,1).^2+x2(:,2).^2).^(1/2); %compute radial coord
        x=[x; x2];
        cross=find(rin(i,:)>1,1); %find first crossing of r=1 cylinder
    end
    rep(i,1:size(x,1),:)=x; 
    pt=x(cross,:).*(1-rin(i,cross))./(rin(i,cross-1)-rin(i,cross))...
        +x(cross-1,:).*(rin(i,cross-1)-1)./(rin(i,cross-1)-rin(i,cross));
    % pt is the linear interpolation of where the trajectory crosses r=1
    % It is an average of the x coordinates weighted by their dist from r=1
    repplot(i,1)=atan2(pt(2),pt(1)); %for plotting
    repplot(i,2)=pt(3);
end

%sort to eliminate wrap-around artifact
attplot=sortrows(attplot);
repplot=sortrows(repplot);

%use periodicity to add an extra point
attplot=[attplot(end,:)-[2*pi 0]; attplot; attplot(1,:)+[2*pi 0]];
repplot=[repplot(end,:)-[2*pi 0]; repplot; repplot(1,:)+[2*pi 0]];

%for j=1:N
%    figure(20)
%    plot(r(j,:),att(j,:,3),'Color',outcolors(j,:));
%    hold on;
%    plot(rin(j,:),rep(j,:,3),'Color',incolors(j,:));
%    hold on;
%end


figure(30)
plot1=plot([attplot(:,1) repplot(:,1)],[attplot(:,2) repplot(:,2)]);
set(plot1(1),'Color','b','LineStyle',':','DisplayName','Attracting Manifold');
set(plot1(2),'Color','r','LineStyle','--','DisplayName','Repelling Manifold');
title({'Intersections of attracting and repelling branches with the cylinder r=1';...
    ['(a = ',num2str(a),', b = ',num2str(b),', \epsilon = ',num2str(eps),')']});
legend('Location','Best');
axis([-pi pi -.0195 -.0145]);

figure(40)
%attplot2=interp1([attplot(end-1:end,1)-2*pi; attplot(:,1); attplot(1:2,1)+2*pi],...
%                [attplot(end-1:end,2); attplot(:,2); attplot(1:2,2)],...
%                angles,'cubic');
%repplot2=interp1([repplot(end-1:end,1)-2*pi; repplot(:,1); repplot(1:2,1)+2*pi],...
%                [repplot(end-1:end,2); repplot(:,2); repplot(1:2,2)],...
%                angles,'cubic');
attplot2=interp1(attplot(:,1),attplot(:,2),angles','cubic');
repplot2=interp1(repplot(:,1),repplot(:,2),angles','cubic');
dif=attplot2-repplot2;
plot(angles(2:end-1),dif(2:end-1),'Color','k');
title('Height of attracting branch above repelling branch at r=1');


%% rotvdp
function xdot=rotvdp(t,x,eps,a,b)
% Differential equations for rotating van der Pol system
% Used by Barry, Benes, Burke, Kaper, and Kramer
% Based on rotating the van der Pol oscillator about the z-axis
% x(1)=x, x(2)=y, x(3)=z
% good parameters: a=0.99403 +/- 1e-5 (maybe 1.5?); b=.001; eps=.1;

xdot(1) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(1) - x(2); 
xdot(2) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(2) + x(1); 
xdot(3) = eps.*(a - ((x(1) - b).^2 + x(2).^2).^(1/2));

xdot=xdot';

%% backward rotvdp
function xdot=backrotvdp(t,x,eps,a,b)
% Differential equations for rotating van der Pol system
% Used by Barry, Benes, Burke, Kaper, and Kramer
% Based on rotating the van der Pol oscillator about the z-axis
% x(1)=x, x(2)=y, x(3)=z
% good parameters: a=0.99403 +/- 1e-5 (maybe 1.5?); b=.001; eps=.1;

xdot(1) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(1) - x(2); 
xdot(2) = (x(3) - (2*(x(1).^2 + x(2).^2).^(3/2) - 3*(x(1).^2 + x(2).^2) + 1)).*x(2) + x(1); 
xdot(3) = eps.*(a - ((x(1) - b).^2 + x(2).^2).^(1/2));
xdot=-xdot;

xdot=xdot';