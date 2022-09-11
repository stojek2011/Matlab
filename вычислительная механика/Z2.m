
load tochki.m
load tochki1.m
load tochki2.m
x=tochki(:,2);
y=tochki(:,1);
figure(1)
hold on ;grid on
ylabel('P ךÍ')
xlabel('dl, לל')
plot(x,y)
x=tochki1(:,2);
y=tochki1(:,1);
figure(2)
hold on ;grid on
ylabel('P ךÍ')
xlabel('dl, לל')
x0=[0:0.001:46];
y0=spline(x,y,x0);
plot(x0,y0)
figure(3)
hold on ;grid on
ylabel('P ךÍ')
xlabel('dl, לל')
p=polyfit(x,y,3)
y1=polyval(p,x0);
plot(x0,y1)
figure(4)
hold on ;grid on
ylabel('C ךÍ/לל')
xlabel('dl, לל')
ddl=diff(x,1);
dp=diff(y,1);
dl=tochki2(:,2);
C=dp./ddl;
plot(dl,C)
figure(5)
hold on ;grid on
ylabel('C ךÍ/לל')
xlabel('dl, לל')
dx0=diff(x0);
dy0=diff(y0);
C0=dy0./dx0;
x01=[0:0.001:46-0.001];
plot(x01,C0)
figure(6)
hold on ;grid on
ylabel('C ךÍ/לל')
xlabel('dl, לל')
dy1=diff(y1);
C1=dy1./dx0;
plot(x01,C1)
figure(7)
hold on ;grid on
ylabel('C ךÍ/לל')
xlabel('dl, לל')
plot(dl,C)
plot(x01,C0)
plot(x01,C1)
figure(8)
hold on ;grid on
ylabel('C ךÍ/לל')
xlabel('dl, לל')
ylabel('Ep Äז')
xlabel('dl, לל')
Ep=E(C,dl);
plot(dl,Ep)
figure(9)
hold on ;grid on
ylabel('C ךÍ/לל')
xlabel('dl, לל')
ylabel('Ep Äז')
xlabel('dl, לל')
Ep0=E(C0,x01);
plot(x01,Ep0)
figure(10)
hold on ;grid on
ylabel('C ךÍ/לל')
xlabel('dl, לל')
ylabel('Ep Äז')
xlabel('dl, לל')
Ep1=E(C1,x01);
plot(x01,Ep1)
figure(11)
hold on ;grid on
ylabel('C ךÍ/לל')
xlabel('dl, לל')
ylabel('Ep Äז')
xlabel('dl, לל')
plot(dl,Ep)
plot(x01,Ep0)
plot(x01,Ep1)
function [ Ep ] = E(C,dl)

 n=length(dl);
 Ep=dl;
 i=1;
while i<=n
  Ep(i)=(C(i)*dl(i)^2)/2;
  i=i+1;

end
end
