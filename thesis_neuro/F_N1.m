function vdot=F_N1(t,V)

%Parameters
a=1.5; b=1; p=.08; I=1.5;
%dv/dt
vdot(1)=10*(V(1)-V(1).^3/3-V(2)+I);
%dr/dt
vdot(2)=p*(1.25*V(1)+a-b*V(2));
vdot=vdot';