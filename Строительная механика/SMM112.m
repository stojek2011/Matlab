close all
clear all
clc
clear
napr(1)
function [u1 u2 u3 u4 v1 v2 v3 v4  Ju Jv Jw W1 W2 W3 W4 F] = otkr()
global alfa 
n=4;
E=0.7*10^9;
t1=1/2;
t2=1/2;
t3=1/2;
t4=2/2;
x1=[t1 10.7];
y1=[t1 t1];
x2=[t2 t2];
y2=[t2 13-t2];
x3=[t2 20.3];
y3=[13-t3 13-t3];
x4=[t2 20.3];
y4=[6.5 6.5];
%%
figure(1)
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
s1=10.7*1;
s2=(13-1*2)*1;
s3=20.3*1;
s4=(20.3-1)*2;
F=s1+s2+s3+s4;
Xc1=10.7/2;
Xc2=0.5;
Xc3=20.3/2;
Xc4=(20.3-1)/2+1;
Xc=(Xc1*s1+Xc2*s2+Xc3*s3+Xc4*s4)/(s1+s2+s3+s4);
Yc1=0.5;
Yc2=6.5;
Yc3=13-0.5;
Yc4=6.5;
Yc=(Yc1*s1+Yc2*s2+Yc3*s3+Yc4*s4)/(s1+s2+s3+s4);
Jx1=(10.7*1^3)/(12);
Jx2=(1*((13-2)^3))/12;
Jx3=(20.3*1^3)/12;
Jx4=((20.3-1)*2^3)/12;
Jx=(Jx1+s1*(Yc1-Yc)^2)+(Jx2+s2*(Yc2-Yc)^2)+(Jx3+s3*(Yc3-Yc)^2)+(Jx4+s4*(Yc4-Yc)^2);
Jy1=(10.7^3*1)/(12);
Jy2=(1^3*((13-2)))/12;
Jy3=(20.3^3*1)/12;
Jy4=((20.3-1)^3*2)/12;
Jy=(Jy1+s1*(Xc1-Xc)^2)+(Jy2+s2*(Xc2-Xc)^2)+(Jy3+s3*(Xc3-Xc)^2)+(Jy4+s4*(Xc4-Xc)^2);
Jxy=s1*(Xc1-Xc)*(Yc1-Yc)+s2*(Xc2-Xc)*(Yc2-Yc)+s3*(Xc3-Xc)*(Yc3-Yc)+s4*(Xc4-Xc)*(Yc4-Yc);
alfa=atan(2*Jxy/(Jy-Jx))/2;
Ju=(Jx+Jy)/2-sqrt(((Jx-Jy)^2)/4+Jxy^2);
Jv=(Jx+Jy)/2+sqrt(((Jx-Jy)^2)/4+Jxy^2);
u1=x1;
u2=x2;
u3=x3;
u4=x4;
v1=y1;
v2=y2;
v3=y3;
v4=y4;
[u1(1),v1(1)]=uv(x1(1),y1(1),Xc,Yc,alfa);
[u1(2),v1(2)]=uv(x1(2),y1(2),Xc,Yc,alfa);
[u2(1),v2(1)]=uv(x2(1),y2(1),Xc,Yc,alfa);
[u2(2),v2(2)]=uv(x2(2),y2(2),Xc,Yc,alfa);
[u3(1),v3(1)]=uv(x3(1),y3(1),Xc,Yc,alfa);
[u3(2),v3(2)]=uv(x3(2),y3(2),Xc,Yc,alfa);
[u4(1),v4(1)]=uv(x4(1),y4(1),Xc,Yc,alfa);
[u4(2),v4(2)]=uv(x4(2),y4(2),Xc,Yc,alfa);
%%
figure(2)
n=20;
m=0.25;
z2=num2str(abs(u1(1)));
z1=num2str(abs(u1(2)));
gr2(x1,y1,u1,0,1,z1,z2,m,n);
z2=' ';
z1=' ';
gr2(x2,y2,u2,1,1,z1,z2,m,n);
z2=num2str(abs(u3(1)));
z1=num2str(abs(u3(2)));
gr2(x3,y3,u3,0,0,z1,z2,m,n);
z2=num2str(abs(u4(1)));
z1=num2str(abs(u4(2)));
gr2(x4,y4,u4,0,1,z1,z2,m,n);
figure(32)
z2=num2str(abs(v1(1)));
z1=num2str(abs(v1(2)));
gr2(x1,y1,v1,0,1,z1,z2,m,n);
z2=' ';
z1=' ';
gr2(x2,y2,v2,1,1,z1,z2,m,n);
z2=num2str(abs(v3(1)));
z1=num2str(abs(v3(2)));
gr2(x3,y3,v3,0,0,z1,z2,m,n);
z2=num2str(abs(v4(1)));
z1=num2str(abs(v4(2)));
gr2(x4,y4,v4,0 ,1,z1,z2,m,n);

%%
figure(4)
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
text(Xc,Yc,' C','FontSize', 15)
S1=[(x2(1)+x2(2))/2,(y2(1)+y2(2))/2];
O1=S1;
plot(Xc,Yc,'h','linewidth',1);

%%
w11=0;
w12=(x1(2)-x1(1))*(y4(1)-y1(1));
w21=0;
w22=0;
w31=0;
w32=-(x3(2)-x3(1))*(y3(1)-y4(1));
%%
figure(5)
m=0.05;
W1=[w11 w12];
W2=[0 0];
W3=[w31 w32];
W4=[0  0];
z1='-61.2';
z2='0';
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1='0';
z2='0';
gr2(x2,y2,W2,1,1,z1,z2,m,n);
z1='118';
z2='0';
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z1='0';
z2='0';
gr2(x4,y4,W4,0 ,0,z1,z2,m,n);
%%
Sv1=(x1(2)-x1(1))*simp(W1,u1,2);
Sv2=(y2(2)-y2(1))*simp(W2,u2,1);
Sv3=(x3(2)-x3(1))*simp(W3,u3,2);
Sv4=(x4(2)-x4(1))*2*simp(W4,u4,2);
Svw1=Sv1+Sv2+Sv3+Sv4;
%%
Su1=(x1(2)-x1(1))*simp(W1,v1,1);
Su2=(y2(2)-y1(1))*simp(W2,v2,2);
Su3=(x3(2)-x3(1))*simp(W3,v3,1);
Su4=(x4(2)-x4(1))*2*simp(W4,v4,1);
Suw1=Su1+Su2+Su3+Su4;
%%
au1=Suw1/Ju;
av1=-Svw1/Jv;
ax1=au1*cos(alfa)-av1*sin(alfa);
ay1=au1*sin(alfa)+av1*cos(alfa);
%%
figure(6)
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
text(Xc,Yc,'C','FontSize', 15)
plot(Xc,Yc,'h','linewidth',1);
A=[S1(1)+ax1 S1(2)+ay1]
O=[S1(1)+ax1 S1(2)];
w42=(x4(2)-A(1))*(ay1);
w41=(x4(1)-A(1))*ay1;
w21=w41+(+A(1)-x2(1))*(y4(1)-y2(1));
w22=w41+(y2(2)-y4(1))*(x2(2)-A(1));
w31=w22;
w32=w22-(x3(2)-x3(1))*(y3(2)-A(2));
w12=w21;
w11=w12+(x1(2)-x1(1))*(A(2)-y1(1));
text(A(1),A(2),'A','FontSize', 15)
plot(A(1),A(2),'h','linewidth',1);
text(O(1),O(2),'.S','FontSize', 15)
plot(O(1),O(2),'h','linewidth',1);
%%
figure(7)
m=0.025;
W1=[w12 w11];
W2=[w21 w22];
W3=[w31 w32];
W4=[w41 w42];
z1='-22.90';
z2='20.50';
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1='-37.30';
z2='20.50';
gr2(x2,y2,W2,1,1,z1,z2,m,n);
z1='20.5';
z2='-37.30';
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z1='-37.30';
z2='116.02';
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z2='20.50';
z1='-63.77';
gr2(x4,y4,W4,0 ,0,z1,z2,m,n);
%%
figure(8)
m=0.01
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
S2=[x1(1) y1(1)];
text(S2(1),S2(2),'S12','FontSize', 15)
plot(S2(1),S2(2),'h','linewidth',1);
%%
figure(9)
w11=0;
w12=0;
w21=0;
w22=0;
w31=0;
w32=-(y3(1)-y1(1))*(x3(2)-x3(1));
w41=0
w42=-(x4(2)-x4(1))*(y4(1)-y1(1));
W1=[w12 w11];
W2=[w21 w22];
W3=[w31 w32];
W4=[w41 w42];
%%
Sv1=(x1(2)-x1(1))*simp(W1,u1,2);
Sv2=(y2(2)-y2(1))*simp(W2,u2,1);
Sv3=(x3(2)-x3(1))*simp(W3,u3,2);
Sv4=(x4(2)-x4(1))*2*simp(W4,u4,2);
Svw1=Sv1+Sv2+Sv3+Sv4;
%%
Su1=(x1(2)-x1(1))*simp(W1,v1,1);
Su2=(y2(2)-y1(1))*simp(W2,v2,2);
Su3=(x3(2)-x3(1))*simp(W3,v3,1);
Su4=(x4(2)-x4(1))*2*simp(W4,v4,1);
Suw1=Su1+Su2+Su3+Su4;
%%
au2=Suw1/Ju;
av2=-Svw1/Jv;
ax2=au2*cos(alfa)-av2*sin(alfa);
ay2=au2*sin(alfa)+av2*cos(alfa);
z1='';
z2='';
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1='';
z2='';
gr2(x2,y2,W2,1,1,z1,z2,m,n);
z1='';
z2='';
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z1='';
z2='';
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z2='';
z1='';
gr2(x4,y4,W4,0 ,0,z1,z2,m,n);
%%
figure(10)
S2=[x1(1) y1(1)];
A2=[S2(1)+ax2 S2(2)+ay2];
plot(A2(1),A2(2),'h','linewidth',1);
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
S3=[x3(1) y3(1)];
text(A2(1),A2(2),'A3','FontSize', 15)
figure(55)
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
hold on; grid on
text(S2(1),S2(2),'S12','FontSize', 15)
plot(S2(1),S2(2),'h','linewidth',1);

%%
figure(11)
S3=[x3(1) y3(1)];
w12=(x1(2)-x1(1))*(y3(2)-y1(1));
w11=0;
w21=0;
w22=0;
w31=0;
w32=0;
w41=0;
w42=(x4(2)-x4(1))*(y3(1)-y4(1));
W1=[w11 w12];
W2=[w21 w22];
W3=[w31 w32];
W4=[w41 w42];
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1='';
z2='';
gr2(x2,y2,W2,1,1,z1,z2,m,n);
z1='';
z2='';
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z1='';
z2='';
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z2='';
z1='';
gr2(x4,y4,W4,0 ,0,z1,z2,m,n); 
%%
Sv1=(x1(2)-x1(1))*simp(W1,u1,2);
Sv2=(y2(2)-y2(1))*simp(W2,u2,1);
Sv3=(x3(2)-x3(1))*simp(W3,u3,2);
Sv4=(x4(2)-x4(1))*2*simp(W4,u4,2);
Svw1=Sv1+Sv2+Sv3+Sv4;
%%
Su1=(x1(2)-x1(1))*simp(W1,v1,1);
Su2=(y2(2)-y1(1))*simp(W2,v2,2);
Su3=(x3(2)-x3(1))*simp(W3,v3,1);
Su4=(x4(2)-x4(1))*2*simp(W4,v4,1);
Suw1=Su1+Su2+Su3+Su4;
%%
au3=Suw1/Ju;
av3=-Svw1/Jv;
ax3=au3*cos(alfa)-av3*sin(alfa);
ay3=au3*sin(alfa)+av3*cos(alfa);
A3=[S3(1)+ax3 S3(2)+ay3];
%%
figure(12)
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
axis([-8 22.5 -2 15])
text(A3(1),A3(2),'A3','FontSize', 10)
text(A2(1),A2(2),'A3','FontSize', 10)
text(A(1),A(2),'A3','FontSize', 10)
%%
S=[x1(1) A(2)];
figure(13)
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
axis([-8 22.5 -2 15])
text(A(1),A(2),'A','FontSize', 10)
text(S(1),S(2),' S1','FontSize', 10)
plot(S(1),S(2),'h','linewidth',1);
plot(A(1),A(2),'h','linewidth',1);
%%
figure(14)
w22=(S(1)-A(1))*(y2(2)-S(2));
w31=w22;
w32=w22-(x3(2)-x3(1))*(y3(2)-S(2));
w11=-(S(1)-A(1))*(S(2)-y2(1));
w21=w11;
w12=w11+(x1(2)-x1(1))*((S(2)-y1(2)));
w41=((w22-w21)/2+w21);
w42=w41+(x4(2)-x4(1))*(S(2)-y4(2));
W1=[w11 w12];
W2=[w21 w22];
W3=[w31 w32];
W4=[w41 w42];
m=0.028
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1='';
z2='';
gr2(x2,y2,W2,1,1,z1,z2,m,n);
z1='';
z2='';
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z1='';
z2='';
gr2(x4,y4,W4,0,0,z1,z2,m,n);
z2='';
z1='';
%
cc=[1 1];
figure(15)
Sw1=(x1(2)-x1(1))*simp(W1,cc,1);
Sw2=(y2(2)-y1(1))*simp(W2,cc,1);
Sw3=(x3(2)-x3(1))*simp(W3,cc,1);
Sw4=(x4(2)-x4(1))*2*simp(W4,cc,1);
Sw=Sw1+Sw2+Sw3+Sw4;
C=-Sw/(s1+s2+s3+s4);
C=[C C];
W1=W1+C;
W2=W2+C;
W3=W3+C;
W4=W4+C;
m=0.07;
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1='';
z2='';
gr2(x2,y2,W2,1,1,z1,z2,m,n);
z1='';
z2='';
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z1='';
z2='';
gr2(x4,y4,W4,0,0,z1,z2,m,n);
z2='';
z1='';
Sw1=(x1(2)-x1(1))*simp(W1,cc,1);
Sw2=(y2(2)-y1(1))*simp(W2,cc,1);
Sw3=(x3(2)-x3(1))*simp(W3,cc,1);
Sw4=(x4(2)-x4(1))*2*simp(W4,cc,1);
Sw=Sw1+Sw2+Sw3+Sw4

%%
l1=(x1(2)-x1(1));
l2=(y2(2)-y1(1));
l3=(x3(2)-x3(1));
l4=(x4(2)-x4(1));
Sw1=l1*simp(W1,W1,1);
Sw2=l2*simp(W2,W2,1);
Sw3=l3*simp(W3,W3,1);
Sw4=l4*2*simp(W4,W4,1);
Sww1=Sw1+Sw2+Sw3+Sw4;
Jw=Sww1
Ww=Jw/58.47
Jk=((y2(2)-y2(1))*1^3+(x4(2)-x4(1))*2^3+(x1(2)-x1(1))*1^3+(x3(2)-x3(1))*1^3)/3

end
function [Mu Mv] = mom()
otkr;
close all
clc
global alfa P q 
P=1.7
l=1.4
q=3.5
h=0.007
M1=2*P*h+P*h;
alfa
%%
k=0.05
x=[0:k:l/2];
n=length(x)
x1=[0:k:l/2]
x2=[l/2:k:l];
x3=[l:k:l*3/2];
x4=[l*3/2:k:2*l];
m1=M1+(q*x.^2)/2-2*P*x;
m2=m1(n)-3/4*P*x-2*P*x+q*(l/2)*x;
m3=m2(n)-(q*x.^2)/2-3/4*P*x-2*P*x+q*(l/2)*x;
m4=m3(n)+(2*P/3)*x-(q*l/2)*x-3/4*P*x-2*P*x+q*(l/2)*x;
%%
figure(16)
Mx=[m1 m2 m3 m4];
X=[x1 x2 x3 x4];
hold on;grid on;
[G,Mx]=inv(X,Mx);
m=1;
mm(X,Mx,m);
marcMx=[m1(1) m2(1) m3(1) m4(1) m4(n)]
marcx=[0 l/2 l 3*l/2 l*2]
y=[0 -10]
ras(marcx,y)
xx=[x1 x2];
M11=xx*1
[G,xx]=inv(M11,xx)
%%
figure(17)
hold on; grid on;
n=length(M11)
n1=length(m3)
line([0 2*l],[0 0]);
mm(xx,M11,m);
marcx=[0 l/2 l 3*l/2 l*2]
y=[0 1.5]
ras(marcx,y)
d1=l*(M11(1)*m3(1)+4*(M11(n)/4)*m3( round(n1/2))+(m3(n1))*M11(n)/2)/(2*6)+l*(M11(n)*m4(1)/2+4*(3*M11(n)/4)*m4( round(n1/2))+(m4(n1))*M11(n))/(2*6)
dd1=l*(M11(1)^2+4*((M11(1)+M11(n))/2)^2+M11(n)^2)/(6)
R1=-d1/dd1
%%
m1=M1+(q*x.^2)/2-2*P*x;
m2=m1(n1)-3/4*P*x-2*P*x+q*(l/2)*x;
m3=m2(n1)-(q*x.^2)/2-3/4*P*x-2*P*x+q*(l/2)*x+R1*x;
m4=m3(n1)+(2*P/3)*x-(q*l/2)*x-3/4*P*x-2*P*x+q*(l/2)*x+R1*x;
Mx=[m1 m2 m3 m4];
X=[x1 x2 x3 x4];
[G,Mx]=inv(X,Mx);
m=1;

ff=l*(M11(1)*m3(1)+4*(M11(n)/4)*m3( round(n1/2))+(m3(n1))*M11(n)/2)/(2*6)+l*(M11(n)*m4(1)/2+4*(3*M11(n)/4)*m4( round(n1/2))+(m4(n1))*M11(n))/(2*6)
figure(18)
hold on; grid on;
mm(X,Mx,m);
marcMx=[m1(1) m2(1) m3(1) m4(1) m4(n1)]
marcx=[0 l/2 l 3*l/2 l*2]
y=[0 -3]
ras(marcx,y);
%%
h3=0.115
h2=0.07969
x=[0:k:l/2];
M2=2*P*h2-2*P*h3
m1=M2+x*0;
m2=m1(n1)+x*0;
m3=m2(n1)-(q*x.^2)/2
m4=m3(n1)+(q*x.^2)/2-((q*l/2))*(x)
My=[m1 m2 m3 m4]
X=[x1 x2 x3 x4];
%%
figure(19)
hold on;grid on;
[G,My]=inv(X,My);
m=1;
xlabel('z')
ylabel('Mp')
mm(X,My,m);
marcMx=[m1(1) m2(1) m3(1) m4(1) m4(n1)]
marcx=[0 l/2 l 3*l/2 l*2]
y=[0 -2]
ras(marcx,y)
%%
d2=l*(M11(1)*m3(1)+4*(M11(n)/4)*m3( round(n1/2))+(m3(n1))*M11(n)/2)/(2*6)+l*(M11(n)*m4(1)/2+4*(3*M11(n)/4)*m4( round(n1/2))+(m4(n1))*M11(n))/(2*6)
dd2=l*(M11(1)^2+4*((M11(1)+M11(n))/2)^2+M11(n)^2)/(6)
R2=-d2/dd2
%%
figure(20)
m2=m1+x*0;
m3=m2(n1)-(q*x.^2)/2+R2*x
m4=m3(n1)+(q*x.^2)/2-((q*l/2))*(x)+R2*x
My=[m1 m2 m3 m4];
hold on;grid on;
xlabel('z')
ylabel('My')
[G,My]=inv(X,My);
m=1;
mm(X,My,m);
marcMx=[m1(1) m2(n1) m3(n1) m4(n1) m4(n1)];
marcx=[0 l/2 l 3*l/2 l*2]
y=[0.6 -0.6]
ras(marcx,y)
f=l*(M11(1)*m3(1)+4*(M11(n)/4)*m3( round(n1/2))+(m3(n1))*M11(n)/2)/(2*6)+l*(M11(n)*m4(1)/2+4*(3*M11(n)/4)*m4( round(n1/2))+(m4(n1))*M11(n))/(2*6)

%%
 Mu=Mx.*cos(alfa)+My.*sin(alfa);
 Mv=-Mx.*sin(alfa)+My.*cos(alfa);
 n=length(Mu)
 figure(21)
hold on;grid on;
xlabel('z')
ylabel('Mu')
mm(X,Mu,m);
marcMu=[Mu(1) Mu(round(n*2/5)) Mu(round(n*3/5)) Mu(round(n*4/5)) Mu(n)];
marcx=[0 l/2 l 3*l/2 l*2];
y=[0 -2];
marcMv=[Mv(1) Mv(round(n*2/5)) Mv(round(n*3/5)) Mv(round(n*4/5)) Mv(n)];
ras(marcx,y)
figure(22)
hold on;grid on;
xlabel('z')
ylabel('Mv')
mm(X,Mv,m);
ras(marcx,y)
end
function [Bw Mk] = MNP()
global k Jw G Jk E P q W4 Bw Mk
syms M0 Bw0 Md 
Jk=66.8;
Jk=Jk*10^-8;
Jw=3.3550*10^5*10^-12;
E=70*10^9;
P=1.7;
l=1.4;
q=3.5;
nu=0.3;
G=E/(2*(1+nu));
k=sqrt(G*Jk/(E*Jw))
ksi1=k*l/2
ksi2=k*l
ksi3=3*k*l/2
ksi4=2*k*l
%%
Bwk=2*P*(-4.01)*10^-4+P*(21.42)*10^-4
m1=0.003*q;
m2=-(0.263+0.003+0.002)*q;
m3=0;
m4=(0.06)*q;
M1=P*(2/3)*(0.263);
M2=-P*(3/4)*(0.167);
M3=-2*P*(0.236);
V0=[0
    0
    Bw0
    M0]
dV1=[0
     0
     0
     M1]
dV2=[0
    0
     0
     Md]
 
dV3=[0
    0
    0
    M2]
%% ћатрица U при ksi=2  v
U2=c(ksi2)*V0+uh(ksi2,m1)+c((ksi2-ksi1))*dV1+uh((ksi2-ksi1),(m2-m1))
%% ћатрица U при ksi=4
U4=c(ksi4)*V0+uh(ksi4,m1)+c((ksi4-ksi1))*dV1+uh((ksi4-ksi1),(m2-m1))+c((ksi4-ksi2))*dV2+uh((ksi4-ksi2),(m3-m2))+c((ksi4-ksi3))*dV3+uh((ksi4-ksi3),(m4-m3))
R=[U2(1)
   U4(3)
   U4(4)+G*Jk*U4(2)]
RR=[0
    Bwk
    M3]
A=R-RR
%% –ешим полученную систему при помощи функции solve
[M0,Bw0,Md]=solve(A==0,M0,Bw0,Md)
V0=[0
    0
    Bw0
    M0]
dV1=[0
     0
     0
     M1]
dV2=[0
     0
     0
     Md]
 
dV3=[0
    0
    0
    M2]
ksi=[0:0.005:ksi4];
Bw=ksi;
fi=ksi;
M=ksi;
Mk=ksi;
Mw=ksi;
dfi=ksi;
nn=length(ksi);
i=1;

while i<=nn
 if((ksi(i)<ksi1))
  U=c(ksi(i))*V0+uh(ksi(i),m1);
  fi(i)=U(1);
  dfi(i)=U(2);
  Bw(i)=U(3);
  Mw(i)=U(4);
  Mk(i)=k*G*Jk*dfi(i);
  M(i)=Mw(i)+Mk(i);
 end
 if((ksi(i)>ksi1)&&(ksi(i)<ksi2))
  U=c(ksi(i))*V0+uh(ksi(i),m1)+c(ksi(i)-ksi1)*dV1+uh((ksi(i)-ksi1),(m2-m1));
  fi(i)=U(1);
  dfi(i)=U(2);
  Bw(i)=U(3);
  Mw(i)=U(4);
  Mk(i)=k*G*Jk*dfi(i);
  M(i)=Mw(i)+Mk(i);
 end
 if((ksi(i)>ksi2)&&(ksi(i)<ksi3))
  U=c(ksi(i))*V0+uh(ksi(i),m1)+c(ksi(i)-ksi1)*dV1+uh((ksi(i)-ksi1),(m2-m1))+c(ksi(i)-ksi2)*dV2+uh((ksi(i)-ksi2),(-m2));
  fi(i)=U(1);
  dfi(i)=U(2);
  Bw(i)=U(3);
  Mw(i)=U(4);
  Mk(i)=k*G*Jk*dfi(i);
  M(i)=Mw(i)+Mk(i);
 end
 if((ksi(i)>ksi3)&&(ksi(i)<=ksi4))
  U=c(ksi(i))*V0+uh(ksi(i),m1)+c(ksi(i)-ksi1)*dV1+uh((ksi(i)-ksi1),(m2-m1))+c(ksi(i)-ksi2)*dV2+uh((ksi(i)-ksi2),(m3-m2))+c(ksi(i)-ksi3)*dV3+uh((ksi(i)-ksi3),(m4-m3));
 fi(i)=U(1);
  dfi(i)=U(2);
  Bw(i)=U(3);
  Mw(i)=U(4);
  Mk(i)=k*G*Jk*dfi(i);
  M(i)=Mw(i)+Mk(i);
 end
 i=i+1;
end
%%
m=1
figure(8)
hold on; grid on;
mm(ksi,fi,m);
xlabel('ksi')
ylabel('fi')
figure(9)
hold on; grid on;
mm(ksi,dfi,m);
xlabel('ksi')
ylabel('dfi/dz')
figure(10)
hold on; grid on;
mm(ksi,Bw,m);
xlabel('ksi')
ylabel('Bw кЌм^2')
figure(11)
hold on; grid on;
mm(ksi,Mk,m);
xlabel('ksi')
ylabel('Mk кЌм')
figure(12)
hold on; grid on;
mm(ksi,Mw,m)
xlabel('ksi')
ylabel('Mw кЌм')
figure(13);
hold on; grid on;
xlabel('ksi')
ylabel('M кЌм')
mm(ksi,M,m);

figure(21)
hold on; grid on;
plot(ksi,M)
end
function [u1 u2 u3 u4 u5 v1 v2 v3 v4 v5 Ju Jv Jw W1 W21 W22 W3 W4 W5 F Jk om] = zakr()
t1=1/2;
t2=1/2;
t3=1/2;
t4=2/2;
t5=1/1;
x1=[t1 10.7]
y1=[t1 t1]
x2=[t2 t2];
y2=[t2 13-t2];
x3=[t2 20.3-0.5];
y3=[13-t3 13-t3];
x4=[t2 20.3-0.5];
y4=[6.5 6.5];
x5=[x4(2) x4(2)];
y5=[y3(2) y4(2)];
%%
figure(29)
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
gr1(x5,y5);
%%
s1=10.7*1;
s2=(13-1*2)*1;
s3=20.3*1;
s4=(20.3-1)*2;
s5=(4.5)*1;
F=s1+s2+s3+s4+s5
Xc1=10.7/2;
Xc2=0.5;
Xc3=20.3/2;
Xc4=(20.3-1)/2+1;
Xc5=(19.8);
Xc=(Xc1*s1+Xc2*s2+Xc3*s3+Xc4*s4+Xc5*s5)/(s1+s2+s3+s4+s5);
Yc1=0.5;
Yc2=6.5;
Yc3=13-0.5;
Yc4=6.5;
Yc5=9.75;
Yc=(Yc1*s1+Yc2*s2+Yc3*s3+Yc4*s4+Yc5*s5)/(s1+s2+s3+s4+s5);
Jx1=(10.7*1^3)/(12);
Jx2=(1*((13-2)^3))/12;
Jx3=(20.3*1^3)/12;
Jx4=((20.3-1)*2^3)/12;
Jx5=((1)*4.5^3)/12;
Jx=(Jx1+s1*(Yc1-Yc)^2)+(Jx2+s2*(Yc2-Yc)^2)+(Jx3+s3*(Yc3-Yc)^2)+(Jx4+s4*(Yc4-Yc)^2)+(Jx5+s5*(Yc5-Yc)^2);
Jy1=(10.7^3*1)/(12);
Jy2=(1^3*((13-2)))/12;
Jy3=(20.3^3*1)/12;
Jy4=((20.3-1)^3*2)/12;
Jy5=(1^3*((4.5)))/12;
Jy=(Jy1+s1*(Xc1-Xc)^2)+(Jy2+s2*(Xc2-Xc)^2)+(Jy3+s3*(Xc3-Xc)^2)+(Jy4+s4*(Xc4-Xc)^2)+(Jy5+s5*(Xc5-Xc)^2);
Jxy=s1*(Xc1-Xc)*(Yc1-Yc)+s2*(Xc2-Xc)*(Yc2-Yc)+s3*(Xc3-Xc)*(Yc3-Yc)+s4*(Xc4-Xc)*(Yc4-Yc)+s5*(Xc5-Xc)*(Yc5-Yc);
alfa=atan(2*Jxy/(Jy-Jx))/2
agrad=180*alfa/pi
Ju=(Jx+Jy)/2-sqrt(((Jx-Jy)^2)/4+Jxy^2);
Jv=(Jx+Jy)/2+sqrt(((Jx-Jy)^2)/4+Jxy^2);
%%
figure(30)
hold on;grid on
text(Xc,Yc,'C','linewidth',15)
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
gr1(x5,y5);

%%
u1=x1;
u2=x2;
u3=x3;
u4=x4;
u5=x5;
v1=y1;
v2=y2;
v3=y3;
v4=y4;
v5=y5;
[u1(1),v1(1)]=uv(x1(1),y1(1),Xc,Yc,alfa);
[u1(2),v1(2)]=uv(x1(2),y1(2),Xc,Yc,alfa);
[u2(1),v2(1)]=uv(x2(1),y2(1),Xc,Yc,alfa);
[u2(2),v2(2)]=uv(x2(2),y2(2),Xc,Yc,alfa);
[u3(1),v3(1)]=uv(x3(1),y3(1),Xc,Yc,alfa);
[u3(2),v3(2)]=uv(x3(2),y3(2),Xc,Yc,alfa);
[u4(1),v4(1)]=uv(x4(1),y4(1),Xc,Yc,alfa);
[u4(2),v4(2)]=uv(x4(2),y4(2),Xc,Yc,alfa);
[u5(1),v5(1)]=uv(x5(1),y5(1),Xc,Yc,alfa);
[u5(2),v5(2)]=uv(x5(2),y5(2),Xc,Yc,alfa);
%%
figure(31)
n=20;
m=0.25;
z2=' ';
z1=' ';
gr2(x2,y2,u2,1,1,z1,z2,m,n);
z2=num2str(abs(u3(1)));
z1=num2str(abs(u3(2)));
gr2(x3,y3,u3,0,0,z1,z2,m,n);
z2=num2str(abs(u4(1)));
z1=num2str(abs(u4(2)));
gr2(x4,y4,u4,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x5,y5,u5,1 ,0,z1,z2,m,n);
z2=num2str(abs(u1(1)));
z1=num2str(abs(u1(2)));
gr2(x1,y1,u1,0,1,z1,z2,m,n);
figure(32)
%%
z2=' ';
z1=' ';
gr2(x2,y2,v2,1,1,z1,z2,m,n);
z2=num2str(abs(v3(1)));
z1=num2str(abs(v3(2)));
gr2(x3,y3,v3,0,0,z1,z2,m,n);
z2=num2str(abs(v4(1)));
z1=num2str(abs(v4(2)));
gr2(x4,y4,v4,0 ,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x5,y5,v5,1 ,0,z1,z2,m,n);
z2=num2str(abs(v1(1)));
z1=num2str(abs(v1(2)));
gr2(x1,y1,v1,0,1,z1,z2,m,n);

%%
figure(33)
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
gr1(x5,y5);
A1= zA1(u1,x1,u2,x2,u3,x3,u4,x4,u5,x5,v1,y1,v2,y2,v3,y3,v4,y4,v5,y5,Ju,Jv,alfa)
A2= zA2(u1,x1,u2,x2,u3,x3,u4,x4,u5,x5,v1,y1,v2,y2,v3,y3,v4,y4,v5,y5,Ju,Jv,alfa)
A3= zA3(u1,x1,u2,x2,u3,x3,u4,x4,u5,x5,v1,y1,v2,y2,v3,y3,v4,y4,v5,y5,Ju,Jv,alfa)
figure(60)
hold on; grid on
text(A1(1),A1(2),'A1')
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
gr1(x5,y5);
[W1 W22 W21 W3 W4 W5 Jw]=zw(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,A1,F);
S=(x3(2)-x3(1))*3/2+(y3(2)-y4(2))*2
om=(x3(2)-x3(1))*(y3(2)-y4(2))*2
(om^2)/20
Jk=om^2/S+((y4(1)-y2(1))*1^3+(x1(2)-x1(1))*1^3)/3

end
function A = zA1(u1,x1,u2,x2,u3,x3,u4,x4,u5,x5,v1,y1,v2,y2,v3,y3,v4,y4,v5,y5,Ju,Jv,alfa)
%%
n=20;
S1=[x4(1) y4(1)];
w41=0;
w42=0;
w52=0;
w51=(x4(2)-x4(1))*(y3(2)-y4(1));
w32=w51;
w31=w32+(x3(2)-x3(1))*(y3(1)-y4(2));
w21=w31;
w22=w21;
w11=w22;
w12=w22+(x1(2)-x1(1))*(y4(1)-y1(1));
%%
figure(34)
m=0.005
W1=[w11 w12];
W2=[w11 w22];
W3=[w31 w32];
W4=[0  0];
W5=[w51 w52];
z2=num2str(abs(W1(1)));
z1=num2str(abs(W1(2)));
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x2,y2,W2,1,1,z1,z2,m,n);
z2=num2str(abs(W3(1)));
z1=num2str(abs(W3(2)));
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x4,y4,W4,0 ,0,z1,z2,m,n);
gr2(x5,y5,W5,1 ,0,z1,z2,m,n);
%%
S=(x3(2)-x3(1))+(x3(2)-x3(1))/2+(y3(2)-y4(2))*2;
om=(x3(2)-x3(1))*(y3(2)-y4(2))*2;
s1=0;
s2=(x4(2)-x4(1))/2;
s3=s2+(y3(2)-y4(2));
s4=s3+x4(2)-x4(1);
s5=s4+(y3(2)-y4(2));
W4=[wp(W4(1),s1,S,om) wp(W4(2),s2,S,om) ];
W5=[wp(W5(1),s3,S,om) wp(W5(2),s2,S,om)];
W3=[wp(W3(1),s4,S,om) wp(W3(2),s3,S,om) ];
W22=[wp(W2(1),s5,S,om) wp(W2(2),s4,S,om)];
W12=[W2(2) W2(1)];
W1=W1-W12;
W12=W12-W12;

%%
figure(35)
z2=num2str(abs(W1(1)));
z1=num2str(abs(W1(2)));
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
x21=[0.5 0.5];
y21=[0.5 y4(1)];
x22=[0.5 0.5];
y22=[y4(1) y2(2)];
z2='';
z1='';
gr2(x21,y21,W12,1,1,z1,z2,m,n);
gr2(x22,y22,W22,1,1,z1,z2,m,n);
z2=num2str(abs(W3(1)));
z1=num2str(abs(W3(2)));
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z2=num2str(abs(W4(1)));
z1=num2str(abs(W4(2)));
gr2(x4,y4,W4,0 ,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x5,y5,W5,1 ,0,z1,z2,m,n);
%
u22=[u4(1) u2(2)];
v22=[v4(1) v2(2)];
%%
Sv1=(x1(2)-x1(1))*simp(W1,u1,1);
Sv2=(y2(2)-y4(1))*simp(W22,u22,1);
Sv3=(x3(2)-x3(1))*simp(W3,u3,1);
Sv4=(x4(2)-x4(1))*2*simp(W4,u4,1);
Sv5=(y3(2)-y4(1))*simp(W5,u5,1);
Svw1=Sv1+Sv2+Sv3+Sv4+Sv5;
%%
Su1=(x1(2)-x1(1))*simp(W1,v1,1);
Su2=(y2(2)-y4(1))*simp(W22,v22,1);
Su3=(x3(2)-x3(1))*simp(W3,v3,1);
Su4=(x4(2)-x4(1))*2*simp(W4,v4,1);
Su5=(y3(2)-y4(1))*simp(W5,v5,1);
Suw1=Su1+Su2+Su3+Su4+Su5;
%%
au1=Suw1/Ju
av1=-Svw1/Jv
ax1=au1*cos(alfa)-av1*sin(alfa)
ay1=au1*sin(alfa)+av1*cos(alfa)
A=[ax1+S1(1) ay1+S1(2)];
end
function A = zA2(u1,x1,u2,x2,u3,x3,u4,x4,u5,x5,v1,y1,v2,y2,v3,y3,v4,y4,v5,y5,Ju,Jv,alfa)
n=20 ;
w22=0;
w21=0;
w41=0;
w42=(y3(2)-y4(2))*(x4(2)-x4(1));
w52=w42;
w51=w52+(x4(2)-x4(1))*(y3(2)-y4(1));
w32=w51;
w31=w32;
w11=0;
w12=w22+(x1(2)-x1(1))*(y3(1)-y1(1));
%%
figure(36)
m=0.025;
W1=[w11 w12];
W2=[w11 w22];
W3=[w31 w32];
W4=[w41 w42];
W5=[w51 w52];
z2=num2str(abs(W1(1)));
z1=num2str(abs(W1(2)));
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x2,y2,W2,1,1,z1,z2,m,n);
z2=num2str(abs(W3(1)));
z1=num2str(abs(W3(2)));
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z2=num2str(abs(W4(1)));
z1=num2str(abs(W4(2)));
gr2(x4,y4,W4,0 ,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x5,y5,W5,1 ,0,z1,z2,m,n);
%%
S=(x3(2)-x3(1))+(x3(2)-x3(1))/2+(y3(2)-y4(2))*2;
om=(x3(2)-x3(1))*(y3(2)-y4(2))*2;
s1=0;
s2=y3(1)-y4(1);
s3=s2+(x4(2)-x4(1))/2;
s4=s3+(y3(2)-y4(2));
s5=s4+x4(2)-x4(1);
%%
W4=[wp(W4(1),s2,S,om) wp(W4(2),s3,S,om) ];
W5=[wp(W5(1),s4,S,om) wp(W5(2),s3,S,om)];
W3=[wp(W3(1),s5,S,om) wp(W3(2),s4,S,om) ];
W22=[wp(W2(1),s2,S,om) wp(W2(2),s1,S,om)];
W12=[W4(1) W4(1)];
W1=W1+W12;
figure(37)
z2=num2str(abs(W1(1)));
z1=num2str(abs(W1(2)));
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
x21=[0.5 0.5];
y21=[0.5 y4(1)];
x22=[0.5 0.5];
y22=[y4(1) y2(2)];
z2='';
z1=num2str(abs(W12(2)));
gr2(x21,y21,W12,1,1,z1,z2,m,n);
gr2(x22,y22,W22,1,1,z1,z2,m,n);
z2=num2str(abs(W3(1)));
z1=num2str(abs(W3(2)));
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z2=num2str(abs(W4(1)));
z1=num2str(abs(W4(2)));
gr2(x4,y4,W4,0 ,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x5,y5,W5,1 ,0,z1,z2,m,n);
%%
u22=[u4(1) u2(2)];
v22=[v4(1) v2(2)];
u12=[u2(1) u4(1)];
v12=[v2(1) v4(1)];
Sv1=(x1(2)-x1(1))*simp(W1,u1,1);
Sv2=(y3(2)-y4(1))*simp(W22,u22,1);
Sv3=(x3(2)-x3(1))*simp(W3,u3,1);
Sv4=(x4(2)-x4(1))*2*simp(W4,u4,1);
Sv21=(y4(1)-y1(1))*simp(W12,u12,1);
Sv5=(y3(2)-y4(1))*simp(W5,u5,1);
Svw1=Sv1+Sv2+Sv3+Sv4+Sv21+Sv5;
%%
Su1=(x1(2)-x1(1))*simp(W1,v1,1);
Su2=(y2(2)-y4(1))*simp(W22,v22,1);
Su3=(x3(2)-x3(1))*simp(W3,v3,1);
Su4=(x4(2)-x4(1))*2*simp(W4,v4,1);
Su21=(y4(1)-y2(1))*simp(W12,v12,1);
Su5=(y3(2)-y4(1))*simp(W5,v5,1);
Suw1=Su1+Su2+Su3+Su4+Su21+Su5;
%%
au1=Suw1/Ju;
av1=-Svw1/Jv;
ax1=au1*cos(alfa)-av1*sin(alfa);
ay1=au1*sin(alfa)+av1*cos(alfa);
S1=[x3(1) y3(1)];
A=[ax1+S1(1) ay1+S1(2)];
end
function A = zA3(u1,x1,u2,x2,u3,x3,u4,x4,u5,x5,v1,y1,v2,y2,v3,y3,v4,y4,v5,y5,Ju,Jv,alfa)
n=20 
w22=0;
w21=0;
w41=0;
w42=(y3(2)-y4(2))*(x4(2)-x4(1));
w52=w42;
w51=w52+(x4(2)-x4(1))*(y3(2)-y4(1));
w32=w51;
w31=w32;
w11=0;
w12=w22+(x1(2)-x1(1))*(y3(1)-y1(1));
%%
figure(36)
m=0.025;
W1=[w11 w12];
W2=[w11 w22];
W3=[w31 w32];
W4=[w41 w42];
W5=[w51 w52];
z2=num2str(abs(W1(1)));
z1=num2str(abs(W1(2)));
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x2,y2,W2,1,1,z1,z2,m,n);
z2=num2str(abs(W3(1)));
z1=num2str(abs(W3(2)));
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z2=num2str(abs(W4(1)));
z1=num2str(abs(W4(2)));
gr2(x4,y4,W4,0 ,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x5,y5,W5,1 ,0,z1,z2,m,n);
%%
S=(x3(2)-x3(1))+(x3(2)-x3(1))/2+(y3(2)-y4(2))*2;
om=(x3(2)-x3(1))*(y3(2)-y4(2))*2;
s1=0;
s2=y3(1)-y4(1);
s3=s2+(x4(2)-x4(1))/2;
s4=s3+(y3(2)-y4(2));
s5=s4+x4(2)-x4(1);
%%
W4=[wp(W4(1),s2,S,om) wp(W4(2),s3,S,om) ];
W5=[wp(W5(1),s4,S,om) wp(W5(2),s3,S,om)];
W3=[wp(W3(1),s5,S,om) wp(W3(2),s4,S,om) ];
W22=[wp(W2(1),s2,S,om) wp(W2(2),s1,S,om)];
W12=[W4(1) W4(1)];
W1=W1+W12;
figure(37)
z2=num2str(abs(W1(1)));
z1=num2str(abs(W1(2)));
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
x21=[0.5 0.5];
y21=[0.5 y4(1)];
x22=[0.5 0.5];
y22=[y4(1) y2(2)];
z2='';
z1=num2str(abs(W12(2)));
gr2(x21,y21,W12,1,1,z1,z2,m,n);
gr2(x22,y22,W22,1,1,z1,z2,m,n);
z2=num2str(abs(W3(1)));
z1=num2str(abs(W3(2)));
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z2=num2str(abs(W4(1)));
z1=num2str(abs(W4(2)));
gr2(x4,y4,W4,0 ,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x5,y5,W5,1 ,0,z1,z2,m,n);
%%
u22=[u2(2) u4(1)];
v22=[v2(2) v4(1)];
u12=[u2(1) u4(1)];
v12=[v2(1) v4(1)];
Sv1=(x1(2)-x1(1))*simp(W1,u1,1);
Sv2=(y3(2)-y4(1))*simp(W22,u22,1);
Sv3=(x3(2)-x3(1))*simp(W3,u3,1);
Sv4=(x4(2)-x4(1))*2*simp(W4,u4,1);
Sv21=(y4(1)-y1(1))*simp(W12,u12,1);
Sv5=(y3(2)-y4(1))*simp(W5,u5,1);
Svw1=Sv1+Sv2+Sv3+Sv4+Sv21+Sv5;
%%
Su1=(x1(2)-x1(1))*simp(W1,v1,1);
Su2=(y2(2)-y4(1))*simp(W22,v22,1);
Su3=(x3(2)-x3(1))*simp(W3,v3,1);
Su4=(x4(2)-x4(1))*2*simp(W4,v4,1);
Su21=(y4(1)-y2(1))*simp(W12,v12,1);
Su5=(y3(2)-y4(1))*simp(W5,v5,1);
Suw1=Su1+Su2+Su3+Su4+Su21+Su5;
%%
au1=Suw1/Ju;
av1=-Svw1/Jv;
ax1=au1*cos(alfa)-av1*sin(alfa);
ay1=au1*sin(alfa)+av1*cos(alfa);
S1=[x3(1) y3(1)];
A=[ax1+S1(1) ay1+S1(2)];
end
function [W1 W22 W21 W3 W4 W5 Jw] = zw(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,A,F)
S2=[x4(1) A(2)];
w41=0;
w42=(x4(2)-x4(1))*(A(2)-y4(2));
w52=w42;
w51=w52+(y3(2)-y4(2))*(x4(2)-A(1));
w32=w51;
w31=w32+(x3(2)-x3(1))*(y3(1)-A(2));
w222=w31;
w221=w222+(y2(2)-y4(1))*(A(1)-x2(2));
w212=w221;
w211=w212+(y4(1)-y2(1))*(A(1)-x2(1));
w11=w211;
w12=w11+(x1(2)-x1(1))*(A(2)-y1(1));
x21=[x1(1) x1(1)];
y21=[y2(1) y4(1)];
x22=[x1(1) x1(1)];
y22=[y4(1) y2(2)];
%%
figure(110)
m=0.005;
n=20;
W1=[w11 w12];
W21=[w211 w212]
W22=[w221 w222];
W3=[w31 w32];
W4=[w41 w42];
W5=[w51 w52];
z2=num2str(abs(W1(1)));
z1=num2str(abs(W1(2)));
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x21,y21,W21,1,1,z1,z2,m,n);
gr2(x22,y22,W22,1,1,z1,z2,m,n);
z2=num2str(abs(W3(1)));
z1=num2str(abs(W3(2)));
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x4,y4,W4,0 ,1,z1,z2,m,n);
gr2(x5,y5,W5,1 ,0,z1,z2,m,n);
text(A(1),A(2),'A')
plot(A(1),A(2),'h','linewidth',1);
%%
S=(x3(2)-x3(1))+(x3(2)-x3(1))/2+(y3(2)-y4(2))*2;
om=(x3(2)-x3(1))*(y3(2)-y4(2))*2;
s1=0;
s2=(x4(2)-x4(1))/2;
s3=s2+(y3(2)-y4(2));
s4=s3+x4(2)-x4(1);
s5=s4+(y3(2)-y4(2));
W4=[wp(W4(1),s1,S,om) wp(W4(2),s2,S,om) ];
W5=[wp(W5(1),s3,S,om) wp(W5(2),s2,S,om)];
W3=[wp(W3(1),s4,S,om) wp(W3(2),s3,S,om) ];
W22=[wp(W22(1),s5,S,om) wp(W22(2),s4,S,om)];
H=W21(2)-W22(1);
H=[H H];
W21=W21-H
W1=W1-H;
%%
figure(111)
m=0.05
z2=num2str(abs(W1(1)));
z1=num2str(abs(W1(2)));
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
x21=[0.5 0.5];
y21=[0.5 y4(1)];
x22=[0.5 0.5];
y22=[y4(1) y2(2)];
z2='';
z1='';
gr2(x21,y21,W21,1,1,z1,z2,m,n);
gr2(x22,y22,W22,1,1,z1,z2,m,n);
z2=num2str(abs(W3(1)));
z1=num2str(abs(W3(2)));
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z2=num2str(abs(W4(1)));
z1=num2str(abs(W4(2)));
gr2(x4,y4,W4,0 ,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x5,y5,W5,1 ,0,z1,z2,m,n);
plot(A(1),A(2),'h','linewidth',1);
text(A(1),A(2),'A')
%%
k=[1 1]
Sw1=(x1(2)-x1(1))*simp(W1,k,1);
Sw2=(y2(2)-y4(1))*simp(W22,k,1);
Sw22=(y4(1)-y2(1))*simp(W21,k,1);
Sw3=(x3(2)-x3(1))*simp(W3,k,1);
Sw4=(x4(2)-x4(1))*2*simp(W4,k,1);
Sw5=(y3(2)-y4(1))*simp(W5,k,1);
Sw=Sw1+Sw2+Sw22+Sw3+Sw4+Sw5;
C=-Sw/(F)
C=[C C]
W1=W1+C;
W3=W3+C;
W4=W4+C;
W21=W21+C;
W22=W22+C;
%%
figure(112)
m=0.05
z2=num2str(abs(W1(1)));
z1=num2str(abs(W1(2)));
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
x21=[0.5 0.5];
y21=[0.5 y4(1)];
x22=[0.5 0.5];
y22=[y4(1) y2(2)];
z2='';
z1='';
gr2(x21,y21,W21,1,1,z1,z2,m,n);
gr2(x22,y22,W22,1,1,z1,z2,m,n);
z2=num2str(abs(W3(1)));
z1=num2str(abs(W3(2)));
gr2(x3,y3,W3,0,0,z1,z2,m,n);
z2=num2str(abs(W4(1)));
z1=num2str(abs(W4(2)));
gr2(x4,y4,W4,0 ,1,z1,z2,m,n);
z1=' ';
z2=' ';
2*2
gr2(x5,y5,W5,1 ,0,z1,z2,m,n);
Sw1=(x1(2)-x1(1))*simp(W1,k,1);
Sw2=(y2(2)-y4(1))*simp(W22,k,1);
Sw22=(y4(1)-y2(1))*simp(W21,k,1);
Sw3=(x3(2)-x3(1))*simp(W3,k,1);
Sw4=(x4(2)-x4(1))*2*simp(W4,k,1);
Sw5=(y3(2)-y4(1))*simp(W5,k,1);
Sw=Sw1+Sw2+Sw22+Sw3+Sw4+Sw5;
C=-Sw/(F)
Sww1=(x1(2)-x1(1))*simp(W1,W1,1);
Sww2=(y3(2)-y4(1))*simp(W22,W22,1);
Sww22=(y4(1)-y2(1))*simp(W21,W21,1);
Sww3=(x3(2)-x3(1))*simp(W3,W3,1);
Sww4=(x4(2)-x4(1))*2*simp(W4,W4,1);
Sww5=(y3(2)-y4(1))*simp(W5,W5,1);
Sww=Sww1+Sww2+Sww22+Sww3+Sww4+Sww5;
Jw=Sww
Ww=Jw/95.9

end
function  []= napr(k)
clc
global alfa 
[u1 u2 u3 u4 v1 v2 v3 v4  Ju Jv Jw W1 W2 W3 W4 F]=otkr;
[Mu Mv] = mom;
[Bw Mk]=MNP;
close all
clc
n1=length(Mu)
x1=[0 10.7]
y1=[0 0]
x2=[0 0];
y2=[0 13];
x3=[0 20.3];
y3=[13 13];
x4=[0 20.3];
y4=[6.5 6.5];

mu=Mu(n1/k)*10^-3
mv=Mv(n1/k)*10^-3
%%
sigu1=-mu/(Ju*10^-8)*(v1*10^-2)
sigu2=-mu/(Ju*10^-8)*(v2*10^-2)
sigu3=-mu/(Ju*10^-8)*(v3*10^-2)
sigu4=-mu/(Ju*10^-8)*(v4*10^-2)
%%
sigv1=-mv/(Jv*10^-8)*(u1*10^-2)
sigv2=-mv/(Jv*10^-8)*(u2*10^-2)
sigv3=-mv/(Jv*10^-8)*(u3*10^-2)
sigv4=-mv/(Jv*10^-8)*(u4*10^-2)
%графих u дл€ нагл€дности
figure(50)
n=20;
m=0.25;
z2=num2str(abs(v1(1)));
z1=num2str(abs(v1(2)));
gr2(x1,y1,v1,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x2,y2,v2,1,1,z1,z2,m,n);
z2=num2str(abs(v3(1)));
z1=num2str(abs(v3(2)));
gr2(x3,y3,v3,0,0,z1,z2,m,n);
z2=num2str(abs(v4(1)));
z1=num2str(abs(v4(2)));
gr2(x4,y4,v4,0,1,z1,z2,m,n);
figure(66)
n=20;
m=0.25;
z2=num2str(abs(u1(1)));
z1=num2str(abs(u1(2)));
gr2(x1,y1,u1,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x2,y2,u2,1,1,z1,z2,m,n);
z2=num2str(abs(u3(1)));
z1=num2str(abs(u3(2)));
gr2(x3,y3,u3,0,0,z1,z2,m,n);
z2=num2str(abs(u4(1)));
z1=num2str(abs(u4(2)));
gr2(x4,y4,u4,0,1,z1,z2,m,n);
%графих напр€жений Mu
figure(51)
m=abs(y3(1)/(sigu3(1)*5))
z2=num2str(abs(sigu1(1)));
z1=num2str(abs(sigu1(2)));
gr2(x1,y1,sigu1,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x2,y2,sigu2,1,1,z1,z2,m,n);
z2=num2str(abs(sigu3(1)));
z1=num2str(abs(sigu3(2)));
gr2(x3,y3,sigu3,0,0,z1,z2,m,n);
z2=num2str(abs(sigu4(1)));
z1=num2str(abs(sigu4(2)));
gr2(x4,y4,sigu4,0 ,1,z1,z2,m,n);
%графих напр€жений Mv
figure(52)
m=abs(y3(1)/(sigv3(2)*5))
z2=' '
z1=' '
gr2(x2,y2,sigv2,1,1,z1,z2,m,n);
z2=num2str(abs(sigv3(1)));
z1=num2str(abs(sigv3(2)));
gr2(x3,y3,sigv3,0,0,z1,z2,m,n);
z2=num2str(abs(sigv4(1)));
z1=num2str(abs(sigv4(2)));
gr2(x4,y4,sigv4,0 ,1,z1,z2,m,n);
z2=num2str(abs(sigv1(1)));
z1=num2str(abs(sigv1(2)));
gr2(x1,y1,sigv1,0,1,z1,z2,m,n);
%графих напр€жений Bw
n2=length(Bw)
bw=Bw(round(n2/k))*10^-3
Jw=Jw*10;
sigw1=bw/(Jw)*(W1*10^-4)
sigw2=bw/(Jw)*(W2*10^-4)
sigw3=bw/(Jw)*(W3*10^-4)
sigw4=bw/(Jw)*(W4*10^-4)
%%
figure(53)
m=y3(1)/(sigw3(1)*5)
z2=' '
z1=' '
gr2(x2,y2,sigw2,1,1,z1,z2,m,n);
z2=num2str(abs(sigw3(1)*10^6));
z1=num2str(abs(sigw3(2)*10^6));
gr2(x3,y3,sigw3,0,0,z1,z2,m,n);
z2=num2str(abs(sigw4(1)*10^6));
z1=num2str(abs(sigw4(2)*10^6));
gr2(x4,y4,sigw4,0 ,1,z1,z2,m,n);
z2=num2str(abs(sigw1(1)*10^6));
z1=num2str(abs(sigw1(2)*10^6));
gr2(x1,y1,sigw1,0,1,z1,z2,m,n);
%%
sigmaz=10^3*5.1/(F*10^-4)*10^-6;
sig11=sigu1+sigv1+sigw1+sigmaz;
sig12=sigu2+sigv2+sigw2+sigmaz;
sig13=sigu3+sigv2+sigw2+sigmaz;
sig14=sigu4+sigv2+sigw2+sigmaz;
%%
figure(54)
m=y3(1)/(sig13(1)*5)
z2=num2str(abs(sig11(1)));
z1=num2str(abs(sig11(2)));
gr2(x1,y1,sig11,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x2,y2,sig12,1,1,z1,z2,m,n);
z2=num2str(abs(sig13(1)));
z1=num2str(abs(sig13(2)));
gr2(x3,y3,sig13,0,0,z1,z2,m,n);
z2=num2str(abs(sig14(1)));
z1=num2str(abs(sig14(2)));
gr2(x4,y4,sig14,0 ,1,z1,z2,m,n);
%
figure(55)
mk=Mk(round(n/k))*10^-3;
Jk=66.8*10^-8;
Wk=10^-6*((y2(2)-y2(1))*1^2+(x4(2)-x4(1))*2^2+(x1(2)-x1(1))*1^2+(x3(2)-x3(1))*1^2)/3;
Jk1=(10^-8)*((x1(2)-x1(1)*1^3))/3;
Jk2=(10^-8)*((y2(2)-y2(1)*1^3))/3;
Jk3=(10^-8)*((x3(2)-x3(1)*1^3))/3;
Jk4=(10^-8)*((x4(2)-x1(1)*2^3))/3;
tau1=Jk1*mk/(Jk*Wk)
T1=[tau1 tau1]
tau2=Jk2*mk/(Jk*Wk)
T2=[tau2 tau2]
tau3=Jk3*mk/(Jk*Wk)
T3=[tau3 tau3]
tau4=Jk4*mk/(Jk*Wk)
T4=[tau4 tau4];
m=y3(1)/(T3(2)*5)
z2=num2str(abs(T1(1)));
z1=num2str(abs(T1(2)));
gr2(x1,y1,T1,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x2,y2,T2,1,1,z1,z2,m,n);
z2=num2str(abs(T3(1)));
z1=num2str(abs(T3(2)));
gr2(x3,y3,T3,0,0,z1,z2,m,n);
z2=num2str(abs(T4(1)));
z1=num2str(abs(T4(2)));
gr2(x4,y4,T4,0 ,1,z1,z2,m,n);
%%
se1=sqrt(sig11.^2+4*T1.^2)
se2=sqrt(sig12.^2+4*T2.^2)
se3=sqrt(sig13.^2+4*T3.^2)
se4=sqrt(sig14.^2+4*T4.^2)
figure(56)
m=y3(1)/(se3(1)*5);
z1=num2str(abs(se1(2)));
z2=num2str(abs(se1(1)));
gr2(x1,y1,se1,0,1,z1,z2,m,n);
z1=num2str(abs(se2(2)));
z2=num2str(abs(se2(1)));
gr2(x2,y2,se2,1,1,z1,z2,m,n);
z2=num2str(abs(se3(1)));
z1=num2str(abs(se3(2)));
gr2(x3,y3,se3,0,0,z1,z2,m,n);
z2=num2str(abs(se4(1)));
z1=num2str(abs(se4(2)));
gr2(x4,y4,se4,0 ,1,z1,z2,m,n);
end
function  []= naprz(k)
clc
global alfa u1 u2 u3 u4 v1 v2 v3 v4  Ju Jv Jw W1 W2 W3 W4 Mu Mv Mk Bw
[u1 u2 u3 u4 u5 v1 v2 v3 v4 v5 Ju Jv Jw W1 W21 W22 W3 W4 W5 F]=zakr;
[Mu Mv] = mom;
[Bw Mk]=MNP;
close all
clc
n1=length(Mu)
x1=[0 10.7]
y1=[0 0]
x2=[0 0];
y2=[0 13];
x3=[0 20.3];
y3=[13 13];
x4=[0 20.3];
y4=[6.5 6.5];
y5=[y3(2) y4(2)]
x5=[x4(2) x4(2)]
mu=Mu(n1/k)*10^-3
mv=Mv(n1/k)*10^-3
Ju
Jv
%%
sigu1=-mu/(Ju*10^-8)*(v1*10^-2)
sigu2=-mu/(Ju*10^-8)*(v2*10^-2)
sigu3=-mu/(Ju*10^-8)*(v3*10^-2)
sigu4=-mu/(Ju*10^-8)*(v4*10^-2)
sigu5=-mu/(Ju*10^-8)*(v5*10^-2)
%%
sigv1=-mv/(Jv*10^-8)*(u1*10^-2)
sigv2=-mv/(Jv*10^-8)*(u2*10^-2)
sigv3=-mv/(Jv*10^-8)*(u3*10^-2)
sigv4=-mv/(Jv*10^-8)*(u4*10^-2)
sigv5=-mv/(Jv*10^-8)*(u5*10^-2)
%графих u дл€ нагл€дности
figure(50)
n=20;
m=0.25;
z2=num2str(abs(u1(1)));
z1=num2str(abs(u1(2)));
gr2(x1,y1,u1,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x2,y2,u2,1,1,z1,z2,m,n);
z2=num2str(abs(u3(1)));
z1=num2str(abs(u3(2)));
gr2(x3,y3,u3,0,0,z1,z2,m,n);
z2=num2str(abs(u4(1)));
z1=num2str(abs(u4(2)));
gr2(x4,y4,u4,0,1,z1,z2,m,n);
%%
figure(66)
n=20;
m=0.25;
z2=num2str(abs(v1(1)));
z1=num2str(abs(v1(2)));
gr2(x1,y1,v1,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x2,y2,v2,1,1,z1,z2,m,n);
z2=num2str(abs(v3(1)));
z1=num2str(abs(v3(2)));
gr2(x3,y3,v3,0,0,z1,z2,m,n);
z2=num2str(abs(v4(1)));
z1=num2str(abs(v4(2)));
gr2(x4,y4,v4,0,1,z1,z2,m,n);

%графих напр€жений Mu
figure(51)
m=y3(1)/(sigu3(1)*5)
z2=num2str(abs(sigu1(1)));
z1=num2str(abs(sigu1(2)));
gr2(x1,y1,sigu1,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x2,y2,sigu2,1,1,z1,z2,m,n);
gr2(x5,y5,sigu5,1,0,z1,z2,m,n);
z2=num2str(abs(sigu3(1)));
z1=num2str(abs(sigu3(2)));
gr2(x3,y3,sigu3,0,0,z1,z2,m,n);
z2=num2str(abs(sigu4(1)));
z1=num2str(abs(sigu4(2)));
gr2(x4,y4,sigu4,0 ,1,z1,z2,m,n);
%графих напр€жений Mv
figure(52)
m=y3(1)/(sigv3(2)*5)
z2=num2str(abs(sigv1(1)));
z1=num2str(abs(sigv1(2)));
gr2(x1,y1,sigv1,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x5,y5,sigv5,1,0,z1,z2,m,n);
gr2(x2,y2,sigv2,1,1,z1,z2,m,n);
z2=num2str(abs(sigv3(1)));
z1=num2str(abs(sigv3(2)));
gr2(x3,y3,sigv3,0,0,z1,z2,m,n);
z2=num2str(abs(sigv4(1)));
z1=num2str(abs(sigv4(2)));
gr2(x4,y4,sigv4,0 ,1,z1,z2,m,n);
%графих напр€жений Bw
figure(53)
z2=num2str(abs(W1(1)));
z1=num2str(abs(W1(2)));
gr2(x1,y1,W1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
x21=[0.5 0.5];
y21=[0.5 y4(1)];
x22=[0.5 0.5];
y22=[y4(1) y2(2)];

%%
figure(54)
m=y3(1)/(sigw3(2)*5)
z2=' '
z1=' '
gr2(x2,y2,sigw2,1,1,z1,z2,m,n);
z2=num2str(abs(sigw3(1)));
z1=num2str(abs(sigw3(2)));
gr2(x3,y3,sigw3,0,0,z1,z2,m,n);
z2=num2str(abs(sigw4(1)));
z1=num2str(abs(sigw4(2)));
gr2(x4,y4,sigw4,0 ,1,z1,z2,m,n);
z2=num2str(abs(sigw1(1)));
z1=num2str(abs(sigw1(2)));
gr2(x1,y1,sigw1,0,1,z1,z2,m,n);
%%
F
sigmaz=10^3*5.1/(F*10^-4)*10^-6;
sig11=sigu1+sigv1+sigw1+sigmaz;
sig12=sigu2+sigv2+sigw2+sigmaz;
sig13=sigu3+sigv2+sigw2+sigmaz;
sig14=sigu4+sigv2+sigw2+sigmaz;
%%
figure(55)
m=15
z2=num2str(abs(sig11(1)));
z1=num2str(abs(sig11(2)));
gr2(x1,y1,sig11,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x2,y2,sig12,1,1,z1,z2,m,n);
z2=num2str(abs(sig13(1)));
z1=num2str(abs(sig13(2)));
gr2(x3,y3,sig13,0,0,z1,z2,m,n);
z2=num2str(abs(sig14(1)));
z1=num2str(abs(sig14(2)));
gr2(x4,y4,sig14,0 ,1,z1,z2,m,n);
%
figure(56)
m=40;
mk=Mk(round(n/k))*10^-3;
Jk=66.8*10^-8;
Wk=10^-6*((y2(2)-y2(1))*1^2+(x4(2)-x4(1))*2^2+(x1(2)-x1(1))*1^2+(x3(2)-x3(1))*1^2)/3;
Jk1=(10^-8)*((x1(2)-x1(1)*1^3))/3;
Jk2=(10^-8)*((y2(2)-y2(1)*1^3))/3;
Jk3=(10^-8)*((x3(2)-x3(1)*1^3))/3;
Jk4=(10^-8)*((x4(2)-x1(1)*2^3))/3;
tau1=Jk1*mk/(Jk*Wk)
T1=[tau1 tau1]
tau2=Jk2*mk/(Jk*Wk)
T2=[tau2 tau2]
tau3=Jk3*mk/(Jk*Wk)
T3=[tau3 tau3]
tau4=Jk4*mk/(Jk*Wk)
T4=[tau4 tau4];
z2=num2str(abs(T1(1)));
z1=num2str(abs(T1(2)));
gr2(x1,y1,T1,0,1,z1,z2,m,n);
z2=' '
z1=' '
gr2(x2,y2,T2,1,1,z1,z2,m,n);
z2=num2str(abs(T3(1)));
z1=num2str(abs(T3(2)));
gr2(x3,y3,T3,0,0,z1,z2,m,n);
z2=num2str(abs(T4(1)));
z1=num2str(abs(T4(2)));
gr2(x4,y4,T4,0 ,1,z1,z2,m,n);
%%
se1=sqrt(sig11.^2+4*T1.^2)
se2=sqrt(sig12.^2+4*T2.^2)
se3=sqrt(sig13.^2+4*T3.^2)
se4=sqrt(sig14.^2+4*T4.^2)
figure(57)
m=10;
z1=num2str(abs(se1(2)));
z2=num2str(abs(se1(1)));
gr2(x1,y1,se1,0,1,z1,z2,m,n);
z1=num2str(abs(se2(1)));
z2=num2str(abs(se2(2)));
gr2(x2,y2,se2,1,1,z1,z2,m,n);
z2=num2str(abs(se3(1)));
z1=num2str(abs(se3(2)));
gr2(x3,y3,se3,0,0,z1,z2,m,n);
z2=num2str(abs(se4(1)));
z1=num2str(abs(se4(2)));
gr2(x4,y4,se4,0 ,1,z1,z2,m,n);
end
function [] = npwk(k)
clc
[u1 u2 u3 u4 u5 v1 v2 v3 v4 v5 Ju Jv Jw W1 W21 W22 W3 W4 W5 F Jk om]=zakr;
[Bw Mk]=MNP;
close all
clc 
x1=[0 10.7];
y1=[0 0];
x2=[0 0];
y2=[0 13];
x3=[0 20.3];
y3=[13 13];
x4=[0 20.3];
y4=[6.5 6.5];
y5=[y3(2) y4(2)];
x5=[x4(2) x4(2)];
%%
figure(1)
n=20;
x21=[0 0];
y21=[0 y4(1)];
x22=[0 0];
y22=[y4(1) y2(2)];
n2=length(Bw)
bw=Bw(round(n2/k))*10^3
Jw=Jw*10*10^-8
sigw1=bw/(Jw)*(W1*10^-4)
sigw21=bw/(Jw)*(W21*10^-4)
sigw22=bw/(Jw)*(W22*10^-4)
sigw3=bw/(Jw)*(W3*10^-4)
sigw4=bw/(Jw)*(W4*10^-4)
sigw5=bw/(Jw)*(W5*10^-4)
%%
m=x3(2)/(sigw1(2)*5)
z1=' ';
z2=' ';
gr2(x21,y21,sigw21,1,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x22,y22,sigw22,1,1,z1,z2,m,n);
z2=num2str(abs(sigw3(1)));
z1=num2str(abs(sigw3(2)));
gr2(x3,y3,sigw3,0,0,z1,z2,m,n);
z2=num2str(abs(sigw4(1)));
z1=num2str(abs(sigw4(2)));
gr2(x4,y4,sigw4,0 ,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x5,y5,sigw5,1 ,0,z1,z2,m,n);
z2=num2str(abs(sigw1(1)));
z1=num2str(abs(sigw1(2)));
gr2(x1,y1,sigw1,0,1,z1,z2,m,n);
%%
figure(2)
n=20;
mk=Mk(round(n/k))*10^3;
om=om*10^-4
Jk1=(10^-8)*((x1(2)-x1(1)*1^3))/3;
Jk21=(10^-8)*((y4(1)-y2(1)*1^3))/3;
Jkk=Jk-Jk1-Jk21
Wk=(Jk1)
tau1=Jk1*mk/(Jk*Wk)
T1=[tau1 tau1]
Wk=(Jk21)
tau21=Jk21*mk/(Jk*Wk)
T21=[tau21 tau21]
mkk=mk*Jkk/Jk;
tau22=mkk/(om*10^-2)
T22=[tau22 tau22]
tau3=mkk/(om*10^-2)
T3=[tau3 tau3]
tau4=mkk/(om*2*10^-2)
T4=[tau4 tau4]
tau5=mkk/(om*10^-2)
T5=[tau5 tau5]

%%
m=y4(2)/(T4(2)*5)
%%
z2=num2str(abs(T1(1)));
z1=num2str(abs(T1(2)));
gr2(x1,y1,T1,0,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x21,y21,T21,1,1,z1,z2,m,n);
z1=num2str(abs(T22(2)));
gr2(x22,y22,T22,1,1,z1,z2,m,n);
z2=num2str(abs(T3(1)));
z1=num2str(abs(T3(2)));
gr2(x3,y3,T3,0,0,z1,z2,m,n);
z2=num2str(abs(T4(1)));
z1=num2str(abs(T4(2)));
gr2(x4,y4,T4,0 ,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x5,y5,T5,1 ,0,z1,z2,m,n);
T3*10^-6
%%
end
function [C] = c(ksi)
global G Jk Jw k E
C=[ 1 sinh(ksi) -(cosh(ksi)-1)/(k^2*E*Jw) (-sinh(ksi)+ksi)/(k*G*Jk)
    0 cosh(ksi) -(sinh(ksi))/(k^2*E*Jw) (-cosh(ksi)+1)/(k*G*Jk)
    0 -G*Jk*sinh(ksi) cosh(ksi) sinh(ksi)/k
    0 -G*Jk*k*cosh(ksi) k*sinh(ksi) cosh(ksi)];
end
function [sr] = tr(u)
sr=(u(2)-u(1))/2+u(1);
end
function [sr] = os(u)
sr=(u(1)+u(2))/2;
end
function [Uh] = uh(ksi,m)
  global G Jk k
   Uh=[m*(1-cosh(ksi)+((ksi)^2)/2)/(G*Jk*k^2)
      m*(-sinh(ksi)+ksi)/(G*Jk*k^2)
      -m*(-cosh(ksi)+1)/k^2
      m*sinh(ksi)/k];
end
function [] = gr1(x,y)
  line([x(1) x(2)],[y(1) y(2)]);
  hold on;grid on
end
function [] = gr2(x,y,u,i,j,g1,g2,m,n)
n1=2;
ff=[0.3010 0.7450 0.9330];
if(g1~=' ')
 g1=str2double(g1);
 g1=round(g1,n1);
 g1=num2str(g1);
end
if(g2~=' ')
g2=str2double(g2);
g2=round(g2,n1);
g2=num2str(g2);
end
u=u*m;
ll=0.5;
kk=0.25;
line([x(1) x(2)],[y(1) y(2)],'LineWidth',4,'Color',ff);
if(j==1)
u=-u;
end
if(i==1)  
 line([(x(1)+u(1)) (x(2)+u(2))],[(y(1)) (y(2))],'Color',ff);
 k=1;
 l=1;
 z1=y(1);
 z2=x(1)+u(1);
 z11=z1+k*u(1);
 z21=z2-j*u(1);

 while k<=n;
 z1=x(1)+u(1)+k*(((u(2))-(u(1))))/n;
 z2=y(1)+(k*((y(2))-(y(1)))/(n));
 line([z1 x(1)],[z2 z2],'Color',ff);
 k=k+1;
 end
 
 z12=z1-kk*u(2);
 z22=z2+ll*u(2);
text(z11,z21,g2,'FontSize', 15);
text(z12,z22,g1,'FontSize', 15);
end
if(i==0)
 line([(x(1)) x(2)],[(y(1)+u(1)) (y(2))+u(2)],'Color',ff);
 k=1;
 z1=x(1);
 z2=y(1)+ll*u(1);
 line([z1 z1],[y(1) z2]);
 z11=z1-kk*u(1);;
 z21=z2+ll*u(1);
 l=1;
 
 while k<=n;
 z1=x(1)+abs(k*((x(2)-x(1))/(n)));
 
 z2=y(1)+u(1)+k*(u(2)-u(1))/n;
 if(l==0)
     z2=-(y(1)+u(1)+k*abs(((u(1))-(u(2))))/n);
 end
 line([z1 z1],[y(1) z2],'Color',ff);
 k=k+1;
 end
 z12=z1+kk*u(2);
 z22=z2;
text(z11,z21,g2,'FontSize', 15);
text(z12,z22,g1,'FontSize', 15);
end    
hold on;grid on
end
function [U,V] = uv(x,y,Xc,Yc,alfa)
  U=(x-Xc)*cos(alfa)+(y-Yc)*sin(alfa);
  V=-(x-Xc)*sin(alfa)+(y-Yc)*cos(alfa);
end
function w = W(S,O,r)
 w=(S(1)-r(1))*(S(2)-r(2));
end
function S = simp(W,v,l)
if(l==1)
S=(v(1)*W(1)+v(2)*W(2)+4*(os(W))*(os(v)))/6;   
end
if(l==2)
S=(v(1)*W(1)+v(2)*W(2)+4*(os(W))*(tr(v)))/6;  
end
end
function S = simp2(W,v,l)
if(l==1)
S=(v(1)*W(1)+v(2)*W(2)+4*(os(W))*(os(v)))/6;   
end
if(l==2)
S=(v(1)*W(1)+v(2)*W(2)+4*(tr(W))*(tr(v)))/6;  
end
end
function [x,y] = inv(x,y)
n=length(x)
i=1;
x1=x;
y1=y;
while i<=n/2
 x(i)=x(n-i+1);
 x(n-i+1)=x1(i);
 y(i)=y(n-i+1);  
 y(n-i+1)=y1(i);
 i=i+1;
end
end
function [] = mm(x,M,m)
n=length(x);
i=1;
while i<n
z1=[x(i) x(i+1)];
z2=-[M(i) M(i+1)];
z3=[0 0];
b='';
m=1;
l=1;
gr2(z1,z3,z2,0,1,b,b,m,l);
i=i+1;
end
end
function  x= ras(x,y)
n=length(x);
i=1;
while i<=n
line([(x(i)) x(i)],[y(1) y(2)],'Color','k','LineWidth',2);
   i=i+1;
end
end
function w = wp(w,s,S,o)
 w=w-o*s/S;
end

