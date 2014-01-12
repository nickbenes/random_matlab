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