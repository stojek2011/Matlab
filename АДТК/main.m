global eps  w  a  b g
eps=0.05;
w=6.5;
a=-0.15;
b=-0.01;
g=0.01;
T=10;
h1=0.05;
h2=0;
n=5;
i=1
A1=(-16*eps/(b*w^2+5*g*w^4))^(1/4)
A0=sqrt(w^2/(abs(a)))
%%%%
A=16.8
Y=0
i=1;
j=1;
while i<=n
while j<=n    
A=A+h2;
A=(A-0.01);
dA=2*eps*w-25*(A^4*g*w^5)/8-(5*A^4*b*w^3)/8;
tspan = [0 T];
y0=[A Y];
[tout1,yout1] = ode113(@Du, tspan, y0);
r1=yout1(:,1);
r21=yout1(:,2);
 figure(1)
hold on;grid on
title('Фазовый портрет')
ylabel('dx');
xlabel('x');
plot(r1,-r21,'b')
plot(r1,r21,'b')
%%%%
y0=[-A -Y];
[tout2,yout2] = ode113(@Du, tspan, y0);
r2=yout2(:,1);
r22=yout2(:,2);
plot(r2,r22,'b')
plot(r2,-r22,'b')
g=g+h;
j=j+1;
end
j=1;
b=b+1;
i=i+1;
end




