clc 
close all
global beta l wch D m q dd Nx0 E h R nu
m=4*10^3
q=10*10^3
l=1.5;
nu=0.3;
dd=-0.2/2*10^-3
R=0.6;
h=0.012
po=1.5*10^6;
l0=2.5*sqrt(h*R)
m=4*10^3;
E=2*10^11;
D=(h^3)*E/(12*(1-nu)^2)
Nx0=-q*l/(2*pi*R)+po*R/4
wch=R^2/E*(po/2-nu*Nx0/R);
beta=((3*(1-nu^2)/(R^2*h^2)))^0.25;
xx1=[0:0.01:l/2];
xx2=[l/2:0.01:l];
xx3=[l:0.01:3*l/2];
xx4=[3*l/2:0.01:2*l];
[WW ww1]=ww();
kk=-SM11(ww1)+Nx0;
KK=KeK(kk)
WW=KK;
WW=WW;
hold on; grid on;
xlabel('x m')
ylabel('\sigma 22 Pa')
WW1=WW(1,:)
plot(xx1,WW1);
WW2=WW(2,:)
plot(xx2,WW2);
WW3=WW(3,:)
plot(xx3,WW3);
WW4=WW(4,:)
plot(xx4,WW4);
function [L1] = S(w)
global DD 
syms s11 s22 s12 s13  s33 w
SS=[s11 0 0
    0 s22 0
    0   0  s33];
ll=eig(SS);
l=subs(ll,s11,(SM11(w)+Nx0/h));
l=subs(l,s22,SM11(w)*nu+SN22(w));
l=subs(l,s33,p0/2);
L1=l;
end
function [WWW www] = ww()
global l wch
syms C11 C22 x
xx=[0:0.01:l/2];
xx2=[-l/2:0.01:0];
ww1=W1(wch);
ww2=W2(wch);
WW1=subs(ww1,[C11 C22],Kr1(ww1));
ww11=WW1;
WW1=subs(WW1,x,xx)
WW2=subs(ww2,[C11 C22],Kr2(ww2));
ww22=WW2;
WW2=subs(WW2,x,xx2);
WW3=subs(ww1,[C11 C22],Kr3(ww1))
ww33=WW3;
WW3=subs(WW3,x,xx);
WW4=subs(ww2,[C11 C22],Kr4(ww2));
ww44=WW4;
WW4=subs(WW4,x,xx2);
WWW=[WW1
     WW2
     WW3
     WW4];
 www=[ww11
      ww22
      ww33
      ww44];
end
function [WWW1] = KeK(WWW)
global  l
syms C11 C22 x
xx=[0:0.01:l/2];
xx2=[-l/2:0.01:0];
ww11=WWW(1,:);
ww22=WWW(2,:);
ww33=WWW(3,:);
ww44=WWW(4,:);
ww11=subs(ww11,x,xx);
ww22=subs(ww22,x,xx2);
ww33=subs(ww33,x,xx);
ww44=subs(ww44,x,xx2);
WWW1=[ww11
     ww22
     ww33
     ww44];

end
function [C] = Kr1(W)
global m
syms C11 C22 x
k1=subs(W,x,0);
k2=subs(Mx(W),x,0);
K=[k1
    k2-m];
[C1 C2]=solve(K,[C11,C22])
C=[C1 C2]
end
function [C] = Kr2(W)
global m q dd
syms C11 C22 x
k1=subs(W,x,0);
k2=subs(Mx(W),x,0);
K=[k1-dd
    k2+m];
[C1 C2]=solve(K,[C11 C22]);
C=[C1 C2]
end
function [C] = Kr3(W)
global m q dd
syms C11 C22 x
k1=subs(W,x,0);
k2=subs(Mx(W),x,0);
K=[k1-dd
    k2-m];
[C1 C2]=solve(K,[C11 C22]);
C=[C1 C2]
end
function [C] = Kr4(W)
global m dd q 
syms C11 C22 x
k1=subs(Q(W),x,0);
k2=subs(Mx(W),x,0);
K=[k1+q
    k2];
[C1 C2]=solve(K,[C11 C22]);
C=[C1 C2];
end
function [W] = W1(wch)
global beta
syms C11 C22 x
W=exp(-beta*x)*(C11*cos(beta*x)+C22*sin(beta*x))+wch;
end
function [W] = W2(wch)
global beta
syms C11 C22 x
W=exp(+beta*x)*(C11*cos(beta*x)+C22*sin(beta*x))+wch;
end
function [ddw] = dw(W)
syms C11 C22 x
ddw=diff(W,x);
end
function [Mxx] = Mx(W)
global D
syms C11 C22 x
Mxx=D*diff(W,x,2);
end
function [QQ] = Q(W)
global D
syms C11 C22 x
QQ=D*diff(W,x,3);
end
function [nx2] = N22(w)
global Nx0 h R E nu
nx2=E*h*w/R+Nx0*nu
end
function [sm11] = SM11(w)
global h
sm11=6*Mx(w)/h^2;
end
function [sm22] = SN22(w)
global h
sm22=N22(w)/h;
end

