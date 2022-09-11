close all
clear all
clc
zakr
function [u1 u2 u3 u4 u5 v1 v2 v3 v4 v5 Ju Jv Jw W1 W2 W3 W4 W5 F] = otkr()
%%
x1=[0 0];
y1=[18.5 0];
x2=[0 21.5];
y2=[0 0];
x3=[0 21.5];
y3=[9.5 9.5];
x4=[11.5 11.5];
y4=[9.5 23.75];
x5=[0.5 11.5 ];
y5=[23.75 23.75];
alfa=0.52;
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
Xc=7.82;
Yc=9.89;
Ju=5.78*10^3;
Jv=4.44*10^3;
Jw=4.6395*10^5*10^(-12)
F=19.5*1+18.5*2+19.5*2+9.5*1+12.25*1;
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
W1=[-53.56 12.036];
W2=[12.036 139.63];
W3=[-20.76 -77.6688];
W4=[-49.21 38.86];
W5=[204.12 38.86 ];
end
function [Mu Mv] = mom()
P=2
l=1.4
q=3.5
h1=(9.896-0.5)*10^-2
h2=(9.896-9.5)*10^-2
M1=-5*P*h1/8+P*h2*3/4
alfa=(30/180)*pi
%%
k=0.01
x=[0:k:l/2];
n=length(x)
x1=[0:k:l/2];
x2=[l/2:k:l];
x3=[l:k:l*3/2];
x4=[l*3/2:k:2*l];
m1=M1+x*0+P*x;
m2=m1(n)-3*P/4*x+P*x;
m3=m2(n)-3*P/4*x-q*(x.^2)/2+P*x;
m4=m3(n)-q*(l/2)*x-3*P/4*x-q*(x.^2)/2+P*x;
%%
figure(1)
Mx=[m1 m2 m3 m4];
X=[x1 x2 x3 x4];
hold on;
ylabel('Mp кЌм')
[G,Mx]=inv(X,Mx);
m=1;
mm(X,Mx,m);
marcMx=[m1(1) m2(1) m3(1) m4(1) m4(n)]
marcx=[0 l/2 l 3*l/2 l*2]
y=[-1 1]
ras(marcx,y)
xx=[x1 x2];
M11=xx*1;
[G,xx]=inv(M11,xx);
%%
figure(2)
hold on; 
ylabel('M1')
n=length(M11)
n1=length(m3)
line([0 2*l],[0 0]);
mm(xx,M11,m);
marcx=[0 l/2 l 3*l/2 l*2]
y=[-1 1]
ras(marcx,y)
d1=l*(M11(1)*m3(1)+4*(M11(n)/4)*m3( round(n1/2))+(m3(n1))*M11(n)/2)/(2*6)+l*(M11(n)*m4(1)/2+4*(3*M11(n)/4)*m4( round(n1/2))+(m4(n1))*M11(n))/(2*6)
dd1=l*(M11(1)^2+4*((M11(1)+M11(n))/2)^2+M11(n)^2)/(6)
R1=-d1/dd1
%%
n=length(x)

m1=M1+x*0+P*x;
m2=m1(n)-3*P/4*x+P*x;
m3=m2(n)-3*P/4*x-q*(x.^2)/2+R1*x+P*x;
m4=m3(n)-q*(l/2)*x-3*P/4*x-q*(x.^2)/2+R1*x+P*x;
Mx=[m1 m2 m3 m4];
X=[x1 x2 x3 x4];
[G,Mx]=inv(X,Mx);
m=1;
round(n1/2)
ff=l*(M11(1)*m3(1)+4*(M11(n)/4)*m3( round(n1/2))+(m3(n1))*M11(n)/2)/(2*6)+l*(M11(n)*m4(1)/2+4*(3*M11(n)/4)*m4( round(n1/2))+(m4(n1))*M11(n))/(2*6)
figure(3)
hold on; 
ylabel('Mx кЌм')

mm(X,Mx,m);
marcMx=[m1(1) m2(1) m3(1) m4(1) m4(n1)]
marcx=[0 l/2 l 3*l/2 l*2]
y=[-1 1]
ras(marcx,y);
%%
h1=(7.8-1)*10^-2
h2=(10.5-7.8)*10^-2
M1=-5*P*h1/8-P*h2*3/4
x=[0:k:l/2];
n=length(x)
m1=M1+0*x;
m2=m1(n)+0*x;
m3=m2(n)+0*x;
m4=m3(n)-3*P/2*x;
My=[m1 m2 m3 m4];
X=[x1 x2 x3 x4];
%%
figure(4)
hold on;
ylabel('Mp кЌм')

[G,My]=inv(X,My);
m=1;
mm(X,My,m);
marcMx=[m1(1) m2(1) m3(1) m4(1) m4(n1)]
marcx=[0 l/2 l 3*l/2 l*2]
y=[-1 1]
ras(marcx,y)
%%
n=length(M11)
n1=length(m3)
d2=l*(M11(1)*m3(1)+4*(M11(n)/4)*m3( round(n1/2))+(m3(n1))*M11(n)/2)/(2*6)+l*(M11(n)*m4(1)/2+4*(3*M11(n)/4)*m4( round(n1/2))+(m4(n1))*M11(n))/(2*6)
dd2=l*(M11(1)^2+4*((M11(1)+M11(n))/2)^2+M11(n)^2)/(6)
R2=-d2/dd2
figure(5)
h1=(6.8)*10^-2
h2=(3.2)*10^-2
m1=M1+0*x;
m2=m1(n1)+0*x;
m3=m2(n1)+0*x+R2*x;
m4=m3(n1)-3*P/2*x+R2*x;
My=[m1 m2 m3 m4];
hold on;
ylabel('My кЌм')

[G,My]=inv(X,My);
m=1;
mm(X,My,m);
marcMx=[m1(1) m2(n1) m3(n1) m4(n1) m4(n1)]
marcx=[0 l/2 l 3*l/2 l*2]
y=[-1 1]
ras(marcx,y)
f=l*(M11(1)*m3(1)+4*(M11(n)/4)*m3( round(n1/2))+(m3(n1))*M11(n)/2)/(2*6)+l*(M11(n)*m4(1)/2+4*(3*M11(n)/4)*m4( round(n1/2))+(m4(n1))*M11(n))/(2*6)
%%
Mu=Mx*cos(alfa)+My*sin(alfa);
Mv=-Mx*sin(alfa)+My*cos(alfa);
 n=length(Mu)
 figure(6)
ylabel('Mu кЌм')
hold on;
mm(X,Mu,m);
marcMu=[Mu(1) Mu(round(n*2/5)) Mu(round(n*3/5)) Mu(round(n*4/5)) Mu(n)]
marcx=[0 l/2 l 3*l/2 l*2]
y=[-1 1]
marcMv=[Mv(1) Mv(round(n*2/5)) Mv(round(n*3/5)) Mv(round(n*4/5)) Mv(n)]
ras(marcx,y)
figure(7)
ylabel('Mv кЌм')
hold on;
mm(X,Mv,m);
ras(marcx,y)
end
function [Bw Mk] = MNP()
global k Jw G Jk E
 syms M0 Bw0 Md 
Jk=(1/3)*(19.5*(2)^3+20.3*(1)^3+18.5*(2)^3+19.5*(1)^3+12.25*(1)^3)
Jk=Jk*10^-8
Jw=4.6395*10^5*10^(-12)
E=70*10^9
nu=0.3
P=2
l=1.4
q=3.5
G=E/(2*(1+nu))
k=sqrt(G*Jk/(E*Jw))
l=1.4
ksi1=k*l/2
ksi2=k*l
ksi3=3*k*l/2
ksi4=2*k*l
%%
%%
Bwk=P*(5/8)*(12*10^-4)+3/4*49.21*10^-4
m1=-(0.168)*q
m2=m1
m3=0
m4=0
M1=-P*(3/2)*(0.059)
M2=P*(3/4)*(0.0364)
M3=-P*(1)*(0.046)

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
U4=c(ksi4)*V0+uh(ksi4,m1)+c((ksi4-ksi1))*dV1+uh((ksi4-ksi1),(m2-m1))+c((ksi4-ksi2))*dV2+uh((ksi4-ksi2),(m3-m2))+c((ksi4-ksi3))*dV3+uh((ksi4-ksi3),(m3))
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
  Mk(i)=G*Jk*dfi(i);
  M(i)=Mw(i)+Mk(i);
 end
 if((ksi(i)>ksi1)&&(ksi(i)<ksi2))
  U=c(ksi(i))*V0+uh(ksi(i),m1)+c(ksi(i)-ksi1)*dV1+uh((ksi(i)-ksi1),(m2-m1));
  fi(i)=U(1);
  dfi(i)=U(2);
  Bw(i)=U(3);
  Mw(i)=U(4);
  Mk(i)=G*Jk*dfi(i);
  M(i)=Mw(i)+Mk(i);
 end
 if((ksi(i)>ksi2)&&(ksi(i)<ksi3))
  U=c(ksi(i))*V0+uh(ksi(i),m1)+c(ksi(i)-ksi1)*dV1+uh((ksi(i)-ksi1),(m2-m1))+c(ksi(i)-ksi2)*dV2+uh((ksi(i)-ksi2),(-m2));
  fi(i)=U(1);
  dfi(i)=U(2);
  Bw(i)=U(3);
  Mw(i)=U(4);
  Mk(i)=G*Jk*dfi(i);
  M(i)=Mw(i)+Mk(i);
 end
 if((ksi(i)>ksi3)&&(ksi(i)<=ksi4))
  U=c(ksi(i))*V0+uh(ksi(i),m1)+c(ksi(i)-ksi1)*dV1+uh((ksi(i)-ksi1),(m2-m1))+c(ksi(i)-ksi2)*dV2+uh((ksi(i)-ksi2),(-m2))+c(ksi(i)-ksi3)*dV3+uh((ksi(i)-ksi3),(m3));
 fi(i)=U(1);
  dfi(i)=U(2);
  Bw(i)=U(3);
  Mw(i)=U(4);
  Mk(i)=G*Jk*dfi(i);
  M(i)=Mw(i)+Mk(i);
 end
 i=i+1;
end
%%
m=1
figure(8)
hold on; grid on
mm(ksi,fi,m);
marcx=[0 ksi1 ksi2 ksi3 ksi4]
y=10^-6*[-3.5 0.7]
ras(marcx,y)
figure(9)
hold on; grid on
mm(ksi,dfi,m);
marcx=[0 ksi1 ksi2 ksi3 ksi4]
y=10^-6*[-2.6 1.2]
ras(marcx,y)
figure(10)
hold on;grid on
mm(ksi,Bw,m);
marcx=[0 ksi1 ksi2 ksi3 ksi4]
y=[-0.2 0.125]
ras(marcx,y)
figure(11)
hold on; grid on
ras(marcx,y)
mm(ksi,Mk,m);
marcx=[0 ksi1 ksi2 ksi3 ksi4]
y=[-0.1 0.05]
ras(marcx,y)
figure(12)
hold on; grid on
mm(ksi,Mw,m)
ras(marcx,y)
marcx=[0 ksi1 ksi2 ksi3 ksi4]
y=[-0.4 0.6]
ras(marcx,y)
figure(13);
hold on; grid on
ras(marcx,y)
mm(ksi,M,m);
marcx=[0 ksi1 ksi2 ksi3 ksi4]
y=[-0.4 0.6]
ras(marcx,y)

end
function  []= napr(k,k1)
clc
global alfa 
[Mu Mv] = mom;
[Bw Mk]=MNP;
[u1 u2 u3 u4 u5 v1 v2 v3 v4 v5 Ju Jv Jw W1 W2 W3 W4 W5 F] = otkr()
close all
clc
n1=length(Mu)
x1=[0 0]
y1=[18.5 0]
x2=[0 21.5];
y2=[0 0];
x3=[0 21.5];
y3=[9.5 9.5];
x4=[11.5 11.5];
y4=[9.5 23.75];
x5=[0.5 11.5 ];
y5=[23.75 23.75];
%%
P=2;
figure(31)
n=20;
m=0.25;
z2=' ';
z1=' ';
u1
u2
u3
u4 
u5
gr2(x2,y2,u2,0,1,z1,z2,m,n);
z2=num2str(abs(u3(1)));
z1=num2str(abs(u3(2)));
gr2(x3,y3,u3,0,1,z1,z2,m,n);
z2=num2str(abs(u4(1)));
z1='';
gr2(x4,y4,u4,1,0,z1,z2,m,n);
z2=num2str(abs(u5(1)));
z1=num2str(abs(u5(2)));
gr2(x5,y5,u5,0 ,0,z1,z2,m,n);
z2=num2str(abs(u1(1)));
z1=num2str(abs(u1(2)));
gr2(x1,y1,u1,1,1,z1,z2,m,n);
mu=Mu(n1/k)*10^-3
mv=Mv(n1/k)*10^-3
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
%графих sigu 
figure(50)
n=20;
m=abs(x2(2)/(sigu2(2)*5))

z2=num2str(abs(sigu4(1)));
z1='';
gr2(x4,y4,sigu4,1,0,z1,z2,m,n);
z2=num2str(abs(sigu5(1)));
z1=num2str(abs(sigu5(2)));
gr2(x5,y5,sigu5,0 ,0,z1,z2,m,n);
z1='';
z2='';
gr2(x1,y1,sigu1,1,1,z1,z2,m,n);
z2=num2str(abs(sigu3(1)));
z1=num2str(abs(sigu3(2)));
gr2(x3,y3,sigu3,0,1,z1,z2,m,n);
z2=num2str(abs(sigu2(1)));
z1=num2str(abs(sigu2(2)));
gr2(x2,y2,sigu2,0,1,z1,z2,m,n);
%графих напр€жений Mv
figure(52)
m=abs(y3(2)/(sigv3(2)*5))
z2=num2str(abs(sigv4(1)));
z1='';
gr2(x4,y4,sigv4,1,0,z1,z2,m,n);
z2=num2str(abs(sigv5(1)));
z1=num2str(abs(sigv5(2)));
gr2(x5,y5,sigv5,0 ,0,z1,z2,m,n);
z1='';
z2='';
gr2(x1,y1,sigv1,1,1,z1,z2,m,n);
z2=num2str(abs(sigv3(1)));
z1=num2str(abs(sigv3(2)));
gr2(x3,y3,sigv3,0,1,z1,z2,m,n);
z2=num2str(abs(sigv2(1)));
z1=num2str(abs(sigv2(2)));
gr2(x2,y2,sigv2,0,1,z1,z2,m,n);
%графих напр€жений Bw
n2=length(Bw)
bw=Bw(round(n2/k))*10^-3
Jw=Jw;
sigw1=bw/(Jw)*(W1*10^-4)
sigw2=bw/(Jw)*(W2*10^-4)
sigw3=bw/(Jw)*(W3*10^-4)
sigw4=bw/(Jw)*(W4*10^-4)
sigw5=bw/(Jw)*(W5*10^-4)
%%
figure(53)
m=abs(y5(1)/(sigw5(1)*5))
z2=num2str(abs(sigw4(1)));
z1='';
gr2(x4,y4,sigw4,1,0,z1,z2,m,n);
z2=num2str(abs(sigw5(1)));
z1=num2str(abs(sigw5(2)));
gr2(x5,y5,sigw5,0 ,0,z1,z2,m,n);
z1='';
z2='';
gr2(x1,y1,sigw1,1,1,z1,z2,m,n);
z2=num2str(abs(sigw2(1)));
z1=num2str(abs(sigw2(2)));
gr2(x2,y2,sigw2,0,1,z1,z2,m,n);
z2=num2str(abs(sigw3(1)));
z1=num2str(abs(sigw3(2)));
gr2(x3,y3,sigw3,0,1,z1,z2,m,n);
%%
sigmaz=-(P/8)*10^3*5.1/(F*10^-4)*10^-6
sig11=sigu1+sigv1+sigw1+sigmaz;
sig12=sigu2+sigv2+sigw2+sigmaz;
sig13=sigu3+sigv3+sigw3+sigmaz;
sig14=sigu4+sigv4+sigw4+sigmaz;
sig15=sigu5+sigv5+sigw5+sigmaz;
%%
figure(54)
m=x2(2)/(sig12(2)*5)
z2=num2str(abs(sig14(1)));
z1='';
gr2(x4,y4,sig14,1,1,z1,z2,m,n);
z2=num2str(abs(sig15(1)));
z1=num2str(abs(sig15(2)));
gr2(x5,y5,sig15,0 ,1,z1,z2,m,n);
z1='';
z2='';
gr2(x1,y1,sig11,1,0,z1,z2,m,n);
z2=num2str(abs(sig13(1)));
z1=num2str(abs(sig13(2)));
gr2(x3,y3,sig13,0,0,z1,z2,m,n);
z2=num2str(abs(sig12(1)));
z1=num2str(abs(sig12(2)));
gr2(x2,y2,sig12,0,0,z1,z2,m,n);
%
figure(55)
n2=length(Mk)
mk=Mk(round(n/k))*10^-3;
Jk=115.41*10^-8
Wk=64.75*10^-6
Jk1=(10^-8)*((y1(1)-y1(2))*2^3)/3
Jk2=(10^-8)*((x2(2)-x2(1))*1^3)/3
Jk3=(10^-8)*((x3(2)-x3(1))*2^3)/3
Jk4=(10^-8)*((y4(2)-y4(1))*1^3)/3
Jk5=(10^-8)*((x5(2)-x5(1)*1^3))/3
J=Jk1+Jk2+Jk3+Jk4+Jk5
tau1=Jk1*mk/(Jk*Wk)
T1=[tau1 tau1]
tau2=Jk2*mk/(Jk*Wk)
T2=[tau2 tau2]
tau3=Jk3*mk/(Jk*Wk)
T3=[tau3 tau3]
tau4=Jk4*mk/(Jk*Wk)
T4=[tau4 tau4];
tau5=Jk5*mk/(Jk*Wk)
T5=[tau5 tau5];
m=y3(1)/(T3(2)*2)
z2=num2str(abs(T4(1)));
z1='';
gr2(x4,y4,T4,1,0,z1,z2,m,n);
z2=num2str(abs(T5(1)));
z1=num2str(abs(T5(2)));
gr2(x5,y5,T5,0 ,0,z1,z2,m,n);
z1='';
z2='';
gr2(x1,y1,T1,1,1,z1,z2,m,n);
z2=num2str(abs(T3(1)));
z1=num2str(abs(T3(2)));
gr2(x3,y3,T3,0,1,z1,z2,m,n);
z2=num2str(abs(T2(1)));
z1=num2str(abs(T2(2)));
gr2(x2,y2,T2,0,1,z1,z2,m,n);
%%
se1=sqrt(sig11.^2+4*T1.^2)
se2=sqrt(sig12.^2+4*T2.^2)
se3=sqrt(sig13.^2+4*T3.^2)
se4=sqrt(sig14.^2+4*T4.^2)
se5=sqrt(sig15.^2+4*T5.^2)
figure(56)
m=abs(y3(1)/(se3(1)*2));
z2=num2str(abs(se4(1)));
z1='';
gr2(x4,y4,se4,1,0,z1,z2,m,n);
z2=num2str(abs(se5(1)));
z1=num2str(abs(se5(2)));
gr2(x5,y5,T5,0 ,0,z1,z2,m,n);
z1='';
z2='';
gr2(x1,y1,T1,1,1,z1,z2,m,n);
z2=num2str(abs(se3(1)));
z1=num2str(abs(se3(2)));
gr2(x3,y3,se3,0,1,z1,z2,m,n);
z2=num2str(abs(se2(1)));
z1=num2str(abs(se2(2)));
gr2(x2,y2,T2,0,1,z1,z2,m,n);
end
function [] = npwk(k)
clc
%[u1 u2 u3 u4 u5 v1 v2 v3 v4 v5 Ju Jv Jw W1 W21 W22 W3 W4 W5 F Jk om]=zakr;
[Bw Mk]=MNP;
close all
clc 
x11=[1 1]
y11=[18.5 9.5]
x12=[1 1]
y12=[9.5 0.5]
x2=[1 21];
y2=[0.5 0.5];
x3=[1 21];
y3=[9.5 9.5];
x4=[11.5 11.5];
y4=[9.5 23.02];
x5=[2 11.5 ];
y5=[23.02 23.02];
x6=[21 21];
y6=[0.5 9.5];
%%

%%
figure(2)
n=20;
n1=length(Mk);
mk=Mk(round(n1/k))*10^3
Jk=3*10^3*10^-8
om=360*10^-4
Jk5=(10^-8)*((12.25)*1^3)/3;
Jk4=(10^-8)*((10.5)*1^3)/3;
Jk1=(10^-8)*((8)*2^3)/3;
Jkk=Jk-Jk1-Jk4-Jk5
Wk=(Jk1*10^2)
tau11=Jk1*mk/(Jk*Wk)
T11=[tau11 tau11]
Wk=(Jk4*10^2)
tau4=Jk4*mk/(Jk*Wk)
T4=[tau4 tau4]
Wk=(Jk5*10^2)
tau5=Jk5*mk/(Jk*Wk)
T5=[tau5 tau5]
%%
mkk=mk*Jkk/Jk
tau2=mkk/(om*10^-2)
T2=[tau2 tau2]
tau3=mkk/(om*2*10^-2)
T3=[tau3 tau3]
tau6=mkk/(om*10^-2)
T6=[tau6 tau6]
tau12=mkk/(om*2*10^-2)
T12=[tau12 tau12]
%%
m=y4(2)/(T3(2)*8)
z1='';
z2='';
gr2(x4,y4,T4,1,0,z1,z2,m,n);
gr2(x5,y5,T5,0 ,0,z1,z2,m,n);
gr2(x6,y6,T6,1,0,z1,z2,m,n);
gr2(x11,y11,T11,1,1,z1,z2,m,n);
gr2(x12,y12,T12,1,1,z1,z2,m,n);
gr2(x3,y3,T3,0,1,z1,z2,m,n);
gr2(x2,y2,T2,0,1,z1,z2,m,n);

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
function x = gr1(x,y)
  line([x(1) x(2)],[y(1) y(2)]);
  hold on;
end
function [] = gr2(x,y,u,i,j,g1,g2,m,n)
n1=2;
ff='g';
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
ll=0.25;
kk=0.25;
line([x(1) x(2)],[y(1) y(2)],[0 1],'LineWidth',4,'Color',ff);
if(j==1)
u=-u;
end
if(i==1)  
 line([(x(1)+u(1)) (x(2)+u(2))],[(y(1)) (y(2))],[0 1],'Color',ff);
 k=1;
 l=1;
 z1=y(1);
 z2=x(1)+u(1);
 z11=z1+k*u(1);
 z21=z2-j*u(1);

 while k<=n;
 z1=x(1)+u(1)+k*(((u(2))-(u(1))))/n;
 z2=y(1)+(k*((y(2))-(y(1)))/(n));
 line([z1 x(1)],[z2 z2],[0 1],'Color',ff);
 k=k+1;
 end
 
 z12=z1-kk*u(2);
 z22=z2+ll*u(2);
text(z11,z21,g2,'FontSize', 15);
text(z12,z22,g1,'FontSize', 15);
end
if(i==0)
 line([(x(1)) x(2)],[(y(1)+u(1)) (y(2))+u(2)],[0 1],'Color',ff);
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
 line([z1 z1],[y(1) z2],[0 1],'Color',ff);
 k=k+1;
 end
 z12=z1+kk*u(2);
 z22=z2;
text(z11,z21,g2,'FontSize', 15);
text(z12,z22,g1,'FontSize', 15);
end    
hold on;
end
function [U,V] = uv(x,y,Xc,Yc,alfa)
  U=(x-Xc)*cos(alfa)+(y-Yc)*sin(alfa);
  V=-(x-Xc)*sin(alfa)+(y-Yc)*cos(alfa);
end
function w = W(S,O,r)

 w=(S(1)-r(1))*(S(2)-r(2));
end
function S = simp(w,u,l)
S=l*(w(1)*u(1)+4*((u(1)+u(2))/2-u(2))*((w(1)+w(2))/2-w(1))+w(2)*u(2));
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
function x = mm(x,M,m)
n=length(x);
i=1;
x;
M;
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
function [u1 u2 u3 u4 u5 v1 v2 v3 v4 v5 Ju Jv Jw W1 W21 W22 W3 W4 W5 F Jk om] = zakr()
x1=[1 1]
y1=[18.5 0.5]
x2=[1 21];
y2=[0.5 0.5];
x3=[1 21];
y3=[9.5 9.5];
x4=[11.5 11.5];
y4=[9.5 23.02];
x5=[2 11.5 ];
y5=[23.02 23.02];
x6=[21 21];
y6=[0.5 9.5];
%%
figure(29)
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
gr1(x5,y5);
gr1(x6,y6);
%%
om=2*((x2(2)-x2(1)))*(y6(2)-y6(1))
s1=2*(18.5)
s2=1*(19.5)
s3=s2*2
s4=1*(12.25)
s5=1*(10.5)
s6=1*(7.5)
F=s1+s2+s3+s4+s5+s6
Xc1=1;
Xc2=11.75
Xc3=11.75;
Xc4=11;
Xc5=6.25;
Xc6=21;
Xc=(Xc1*s1+Xc2*s2+Xc3*s3+Xc4*s4+Xc5*s5+Xc6*s6)/(s1+s2+s3+s4+s5+s6)
Yc1=9.25;
Yc2=0.5;
Yc3=9.5;
Yc4=16.625;
Yc5=23.25;
Yc6=4.75;
Yc=(Yc1*s1+Yc2*s2+Yc3*s3+Yc4*s4+Yc5*s5+Yc6*s6)/(s1+s2+s3+s4+s5+s6)
Jx1=(2*(18.5)^3)/12
Jx2=1^3*(19.5)/12
Jx3=2^3*(19.5)/12
Jx4=1*((12.25)^3)/12
Jx5=1^3*(10.5)/12
Jx6=2*((7.5)^3)/12
Jx=(Jx1+s1*(Yc1-Yc)^2)+(Jx2+s2*(Yc2-Yc)^2)+(Jx3+s3*(Yc3-Yc)^2)+(Jx4+s4*(Yc4-Yc)^2)+(Jx5+s5*(Yc5-Yc)^2)+(Jx6+s6*(Yc6-Yc)^2)
Jy1=2^3*(18.5)/12
Jy2=1*(19.5)^3/12
Jy3=2*(19.5)^3/12
Jy4=1^3*(12.25)/12
Jy5=1*(10.5)^3/12
Jy6=2^3*(7.5)/12
Jy=(Jy1+s1*(Xc1-Xc)^2)+(Jy2+s2*(Xc2-Xc)^2)+(Jy3+s3*(Xc3-Xc)^2)+(Jy4+s4*(Xc4-Xc)^2)+(Jy5+s5*(Xc5-Xc)^2)+(Jy6+s6*(Xc6-Xc)^2)
Jxy=s1*(Xc1-Xc)*(Yc1-Yc)+s2*(Xc2-Xc)*(Yc2-Yc)+s3*(Xc3-Xc)*(Yc3-Yc)+s4*(Xc4-Xc)*(Yc4-Yc)+s5*(Xc5-Xc)*(Yc5-Yc)+s6*(Xc6-Xc)*(Yc6-Yc)
alfa=atan(2*Jxy/(Jy-Jx))/2
agrad=180*alfa/pi
Ju=(Jx+Jy)/2-sqrt(((Jx-Jy)^2)/4+Jxy^2)
Jv=(Jx+Jy)/2+sqrt(((Jx-Jy)^2)/4+Jxy^2)
Juv=(Jx-Jy)/2*sin(2*alfa)+Jxy*cos(2*alfa)
%%

%%
u1=x1;
u2=x2;
u3=x3;
u4=x4;
u5=x5;
u6=x6;
v1=y1;
v2=y2;
v3=y3;
v4=y4;
v5=y5;
v6=y6;
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
[u6(1),v6(1)]=uv(x6(1),y6(1),Xc,Yc,alfa);
[u6(2),v6(2)]=uv(x6(2),y6(2),Xc,Yc,alfa);
u1
u2
u3
u4
u5
u6
v1
v2
v3
v4
v5
v6
%%
figure(31)
n=20;
m=abs(y3(2)/(v3(2)*5))
z2=num2str(abs(v4(1)));
z1='';
gr2(x4,y4,v4,1,0,z1,z2,m,n);
z2=num2str(abs(v5(1)));
z1=num2str(abs(v5(2)));
gr2(x5,y5,v5,0 ,0,z1,z2,m,n);
z1='';
z2='';
gr2(x6,y6,v6,1,0,z1,z2,m,n);
gr2(x1,y1,v1,1,1,z1,z2,m,n);
z2=num2str(abs(v3(1)));
z1=num2str(abs(v3(2)));
gr2(x3,y3,v3,0,0,z1,z2,m,n);
z2=num2str(abs(v2(1)));
z1=num2str(abs(v2(2)));
gr2(x2,y2,v2,0,1,z1,z2,m,n);
figure(32)
%%
m=abs(y3(2)/(u3(2)*5))
z2=num2str(abs(u4(1)));
z1='';
gr2(x4,y4,u4,1,0,z1,z2,m,n);
z2=num2str(abs(u5(1)));
z1=num2str(abs(u5(2)));
gr2(x5,y5,u5,0 ,0,z1,z2,m,n);
z1='';
z2='';
gr2(x6,y6,u6,1,0,z1,z2,m,n);
gr2(x1,y1,u1,1,1,z1,z2,m,n);
z2=num2str(abs(u3(1)));
z1=num2str(abs(u3(2)));
gr2(x3,y3,u3,0,0,z1,z2,m,n);
z2=num2str(abs(u2(1)));
z1=num2str(abs(u2(2)));
gr2(x2,y2,u2,0,1,z1,z2,m,n);

%%
A1= zA1(u1,x1,u2,x2,u3,x3,u4,x4,u5,x5,u6,x6,v1,y1,v2,y2,v3,y3,v4,y4,v5,y5,v6,y6,Ju,Jv,alfa)
%A2= zA2(u1,x1,u2,x2,u3,x3,u4,x4,u5,x5,v1,y1,v2,y2,v3,y3,v4,y4,v5,y5,Ju,Jv,alfa)
%A3= zA3(u1,x1,u2,x2,u3,x3,u4,x4,u5,x5,v1,y1,v2,y2,v3,y3,v4,y4,v5,y5,Ju,Jv,alfa)
figure(60)
hold on; grid on
text(A1(1),A1(2),'A1')
gr1(x1,y1);
gr1(x2,y2);
gr1(x3,y3);
gr1(x4,y4);
gr1(x5,y5);
S=(x3(2)-x3(1))/2+(x3(2)-x3(1))/1+(y6(2)-y6(1))/2+(y6(2)-y6(1))/1
om=(x3(2)-x3(1))*(y6(2)-y6(1))*2
Jk=om^2/S+(7.5*2^3)/3+(10.5*1^3)/3+(12.25*1^3)/3
[W1 W22 W21 W3 W4 W5 Jw]=zw(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,A1,F);
%Jk=om^2/S+((y4(1)-y2(1))*1^3+(x1(2)-x1(1))*1^3)/3

end
function A = zA1(u1,x1,u2,x2,u3,x3,u4,x4,u5,x5,u6,x6,v1,y1,v2,y2,v3,y3,v4,y4,v5,y5,v6,y6,Ju,Jv,alfa)
%%
S1=[x4(1) y4(1)];
%%
n=20;
x11=[x1(1) x1(1)];
y11=[y1(1) y3(1)];
x12=[x1(1) x1(1)];
y12=[y3(1) y1(2)];
x31=[x3(1) x4(1)];
y31=[y3(1) y3(1)];
x32=[x4(2) x6(1)];
y32=[y3(1) y3(1)];
w311=0;
w312=0;
w112=0;
w121=0;
w111=-(x4(1)-x1(1))*(y1(1)-y4(1));
w122=(x4(1)-x1(1))*(-y1(2)+y4(1));
w21=w122;
w22=w21+(x2(2)-x2(1))*(y4(1)-y2(1));
w61=w22;
w62=w61+(y6(2)-y6(1))*(x6(2)-x4(1))
w322=w62;
w321=w322;
w41=w321;
w42=w41+0;
w52=w42
w51=w52+(x4(2)-x5(1))*(y5(1)-y4(1))
%%
figure(34)
m=0.008
n=30
W11=[w111 w112]
W12=[w121 w122]
W2=[w21 w22]
W31=[w311 w312];
W32=[w321 w322];
W4=[w41  w42]
W5=[w51 w52];
W6=[w61 w62];
z1=' ';
z2=' ';
gr2(x11,y11,W11,1,1,z1,z2,m,n);
gr2(x12,y12,W12,1,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x31,y31,W31,0,0,z1,z2,m,n);
gr2(x32,y32,W32,0,0,z1,z2,m,n);
z2='';
z1=num2str(abs(W4(2)));
gr2(x4,y4,W4,1 ,0,z1,z2,m,n);
z1=num2str(abs(W5(1)));
z2=' ';
gr2(x5,y5,W5,0 ,0,z1,z2,m,n);
gr2(x6,y6,W6,1 ,0,z1,z2,m,n);
z2=num2str(abs(W2(1)));
z1=num2str(abs(W2(2)));
gr2(x2,y2,W2,0,1,z1,z2,m,n);
%%
figure(55)
z1='';
z2='';
%z2=num2str(abs(W4(1)));
z1='';
gr2(x4,y4,W4,1,0,z1,z2,m,n);
%z2=num2str(abs(W5(1)));
%z1=num2str(abs(W5(2)));
gr2(x5,y5,W5,0 ,0,z1,z2,m,n);
%z2=num2str(abs(W6(1)));
z2='';
gr2(x6,y6,W6,1,0,z1,z2,m,n);
z1='';
z2='';
gr2(x11,y11,W11,1,1,z1,z2,m,n);
gr2(x12,y12,W12,1,1,z1,z2,m,n);
gr2(x32,y32,W32,0,0,z1,z2,m,n);
gr2(x31,y31,W31,0,0,z1,z2,m,n);
%z2=num2str(abs(W2(1)));
%z1=num2str(abs(W2(2)));
gr2(x2,y2,W2,0,1,z1,z2,m,n);
%%
S=(x3(2)-x3(1))/2+(x3(2)-x3(1))/1+(y6(2)-y6(1))/2+(y6(2)-y6(1))/1
om=(x3(2)-x3(1))*(y6(2)-y6(1))*2
s1=0;
s2=(x31(2)-x31(1))/2
s3=s2+(y3(1)-y2(1))/2
s4=s3+(x2(2)-x2(1))
s5=s4+(y6(2)-y6(1))
s6=s5+(x3(2)-x4(1))/2
W31=[wp(W31(1),s2,S,om) wp(W31(2),s1,S,om) ]
W12=[wp(W12(1),s2,S,om) wp(W12(2),s3,S,om)]
W2=[wp(W2(1),s3,S,om) wp(W2(2),s4,S,om)]
W6=[wp(W6(1),s4,S,om) wp(W6(2),s5,S,om) ]
W32=[wp(W32(1),s6,S,om) wp(W32(2),s5,S,om)]
W11=[W11(1)+W12(1) W11(2)+W12(1)]
dw4=(W32(1)-W4(1));
Dw4=[dw4 dw4];
W4=W4+Dw4;
W5=W5+Dw4
%%
figure(35)
m=0.04
z1='';
z2='';
%z2=num2str(abs(v4(1)));
z1='';
gr2(x4,y4,W4,1,0,z1,z2,m,n);
%z2=num2str(abs(v5(1)));
%z1=num2str(abs(v5(2)));
gr2(x5,y5,W5,0 ,0,z1,z2,m,n);
z1='';
z2='';
gr2(x6,y6,W6,1,0,z1,z2,m,n);
gr2(x11,y11,W11,1,1,z1,z2,m,n);
gr2(x12,y12,W12,1,1,z1,z2,m,n);
%z2=num2str(abs(v3(1)));
%z1=num2str(abs(v3(2)));
gr2(x31,y31,W31,0,0,z1,z2,m,n);
gr2(x32,y32,W32,0,0,z1,z2,m,n);
%z2=num2str(abs(v2(1)));
%z1=num2str(abs(v2(2)));
gr2(x2,y2,W2,0,1,z1,z2,m,n);
u11=[u1(1) u3(1)];
v11=[v1(1) v3(1)];
u12=[u3(1) u1(2)];
v12=[v3(1) v1(2)];
u31=[u3(1) u4(1)];
v31=[v3(1) v4(1)];
u32=[u4(1) u3(2)];
v32=[v4(1) v3(2)];
%%
l11=(y1(1)-y3(1));
l12=(y3(1)-y1(2));
l2=(x2(2)-x2(1));
l31=(x4(1)-x3(1));
l32=(x3(2)-x4(1));
l4=(y4(2)-y4(1));
l5=(x5(2)-x5(1));
l6=(y6(2)-y6(1));
Sv11=l11*2*simp(W11,u11,1);
Sv12=l12*2*simp(W12,u12,1);
Sv2=l2*simp(W2,u2,1);
Sv31=l31*2*simp(W31,u31,1);
Sv32=l32*2*simp(W32,u32,1);
Sv4=l4*simp(W4,u4,1);
Sv5=l5*simp(W5,u5,1);
Sv6=l6*simp(W6,u6,1);
Svw1=Sv12+Sv11+Sv2+Sv31+Sv32+Sv4+Sv5+Sv6
%%
Su11=l11*2*simp(W11,v11,1);
Su12=l12*2*simp(W12,v12,1);
Su2=l2*simp(W2,v2,1);
Su31=l31*2*simp(W31,v31,1);
Su32=l32*2*simp(W32,v32,1);
Su4=l4*simp(W4,v4,1);
Su5=l5*simp(W5,v5,1);
Su6=l6*simp(W6,v6,1);
Suw1=Su12+Su11+Su2+Su31+Su32+Su4+Su5+Su6
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
function [W1 W22 W21 W3 W4 W5 Jw] = zw(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,A,F)
S2=[x4(1) A(2)];
x11=[x1(1) x1(1)];
y11=[y1(1) y3(1)];
x12=[x1(1) x1(1)];
y12=[y3(1) y1(2)];
x31=[x3(1) x4(1)];
y31=[y3(1) y3(1)];
x32=[x4(2) x6(1)];
y32=[y3(1) y3(1)];
w312=0;
w311=(x31(2)-x31(1))*(y3(1)-A(2));
w111=w311-(y1(1)-y3(1))*(A(1)-x1(1));
w112=w311;
w121=w311;
w122=w121+(y3(2)-y2(1))*(A(1)-x1(1));
w21=w122;
w22=w21+(x2(2)-x2(1))*(A(2)-y2(2));
w61=w22;
w62=w61+(y6(2)-y6(1))*(x6(2)-A(1));
w322=w62;
w321=w322+(x2(2)-x4(1))*(y3(1)-A(2));
w41=w321;
w42=w41-(y4(2)-y4(1))*(A(1)-x4(1));
w52=w42;
w51=w52+(x5(2)-x5(1))*(y5(1)-A(2));
%%
figure(110)
m=0.008
n=20
W11=[w111 w112]
W12=[w121 w122]
W2=[w21 w22]
W31=[w311 w312]
W32=[w321 w322]
W4=[w41  w42]
W5=[w51 w52]
W6=[w61 w62]
z1=' ';
z2=' ';
gr2(x11,y11,W11,1,1,z1,z2,m,n);
gr2(x12,y12,W12,1,1,z1,z2,m,n);
z1=' ';
z2=' ';
gr2(x31,y31,W31,0,0,z1,z2,m,n);
gr2(x32,y32,W32,0,0,z1,z2,m,n);
z2='';
%z1=num2str(abs(W4(2)));
gr2(x4,y4,W4,1 ,0,z1,z2,m,n);
%z1=num2str(abs(W5(1)));
z2=' ';
gr2(x5,y5,W5,0 ,0,z1,z2,m,n);
gr2(x6,y6,W6,1 ,0,z1,z2,m,n);
%z2=num2str(abs(W2(1)));
%z1=num2str(abs(W2(2)));
gr2(x2,y2,W2,0,1,z1,z2,m,n);
%%
S=(x3(2)-x3(1))/2+(x3(2)-x3(1))/1+(y6(2)-y6(1))/2+(y6(2)-y6(1))/1
om=(x3(2)-x3(1))*(y6(2)-y6(1))*2
s1=0;
s2=(x31(2)-x31(1))/2;
s3=s2+(y3(1)-y2(1))/2;
s4=s3+(x2(2)-x2(1))
s5=s4+(y6(2)-y6(1));
s6=s5+(x3(2)-x4(1))/2
hh=W12(1)
W31=[wp(W31(1),s2,S,om) wp(W31(2),s1,S,om) ];
W12=[wp(W12(1),s2,S,om) wp(W12(2),s3,S,om)];
W2=[wp(W2(1),s3,S,om) wp(W2(2),s4,S,om)];
W6=[wp(W6(1),s4,S,om) wp(W6(2),s5,S,om) ];
W32=[wp(W32(1),s6,S,om) wp(W32(2),s5,S,om)];
dw1=(W12(1)-hh)
Dw1=[dw1 dw1]
dw4=(W32(1)-W4(1));
Dw4=[dw4 dw4];
W4=W4+Dw4;
W5=W5+Dw4
W11=W11+Dw1
W31
W32
W11
W12
W2
W6

%%
figure(111)
m=0.04
z1='';
z2='';
%z2=num2str(abs(v4(1)));
z1='';
gr2(x4,y4,W4,1,0,z1,z2,m,n);
%z2=num2str(abs(v5(1)));
%z1=num2str(abs(v5(2)));
gr2(x5,y5,W5,0 ,0,z1,z2,m,n);
z1='';
z2='';
gr2(x6,y6,W6,1,0,z1,z2,m,n);
gr2(x11,y11,W11,1,1,z1,z2,m,n);
gr2(x12,y12,W12,1,1,z1,z2,m,n);
%z2=num2str(abs(v3(1)));
%z1=num2str(abs(v3(2)));
gr2(x31,y31,W31,0,0,z1,z2,m,n);
gr2(x32,y32,W32,0,0,z1,z2,m,n);
%z2=num2str(abs(v2(1)));
%z1=num2str(abs(v2(2)));
gr2(x2,y2,W2,0,1,z1,z2,m,n);
%plot(A(1),A(2),'h','linewidth',1);
%text(A(1),A(2),'A')
%%
k=[1 1]
l11=(y1(1)-y3(1));
l12=(y3(1)-y1(2));
l2=(x2(2)-x2(1));
l31=(x4(1)-x3(1));
l32=(x3(2)-x4(1));
l4=(y4(2)-y4(1));
l5=(x5(2)-x5(1));
l6=(y6(2)-y6(1));
Sw11=l11*2*simp(W11,k,1);
Sw12=l12*2*simp(W12,k,1);
Sw2=l2*simp(W2,k,1);
Sw31=l31*2*simp(W31,k,1);
Sw32=l32*2*simp(W32,k,1);
Sw4=l4*simp(W4,k,1);
Sw5=l5*simp(W5,k,1);
Sw6=l6*simp(W6,k,1);
Sw=Sw12+Sw11+Sw2+Sw31+Sw32+Sw4+Sw5+Sw6
k=[1 1]
C=-Sw/(F)
C=[C C]
W4=W4+C;
W11=W11+C
W12=W12+C
W2=W2+C
W31=W31+C
W32=W32+C
W5=W5+C
W6=W6+C

%%
figure(112)
m=0.04
z1='';
z2='';
%z2=num2str(abs(v4(1)));
z1='';
gr2(x4,y4,W4,1,0,z1,z2,m,n);
%z2=num2str(abs(v5(1)));
%z1=num2str(abs(v5(2)));
gr2(x5,y5,W5,0 ,0,z1,z2,m,n);
z1='';
z2='';
gr2(x6,y6,W6,1,0,z1,z2,m,n);
gr2(x11,y11,W11,1,1,z1,z2,m,n);
gr2(x12,y12,W12,1,1,z1,z2,m,n);
%z2=num2str(abs(v3(1)));
%z1=num2str(abs(v3(2)));
gr2(x31,y31,W31,0,0,z1,z2,m,n);
gr2(x32,y32,W32,0,0,z1,z2,m,n);
%z2=num2str(abs(v2(1)));
%z1=num2str(abs(v2(2)));
gr2(x2,y2,W2,0,1,z1,z2,m,n);
%plot(A(1),A(2),'h','linewidth',1);
%text(A(1),A(2),'A')
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
function w = wp(w,s,S,o)
 w=w-o*s/S;
end
