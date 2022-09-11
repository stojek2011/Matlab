x01=[0 2];
x02=[2 4]; 
x03=[4 6]; 
x04=[6 10];
y01=3*x01/2;
y02=-5*x02/2+8;
y03=4+x03*0;
y04=x04-10;
figure(1)
hold on;grid on
xlabel('x')
ylabel('F(x)')
axis([0 10  -6 6])
plot(x01,y01,x02,y02,x03,y03,x04,y04)
figure(2)
hold on;grid on
xlabel('x')
ylabel('Ï(x)')
%%%%%%%%%%%%%%%%%%%%%%
x1=[-5:0.01:2];
y1=(3/4)*x1.^2;
x2=[2:0.01:4];
n1=length(y1);
y2=-5*(x2.^2)/4+8*x2;
C1=y2(1)-y1(n1)
y2=y2-C1;
x3=[4:0.01:6];
n2=length(y2)
y3=4*x3;
C2=y3(1)-y2(n2)
y3=y3-C2;
n3=length(y3)
x4=[6:0.01:15];
y4=(x4.^2)/2-10*x4;
C3=y4(1)-y3(n3)
y4=y4-C3;
plot(x1,y1,x2,y2,x3,y3,x4,y4)
%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on;grid on
xlabel('x')
ylabel('dx')
x=[-100:0.01:100];
h1=100
h2=-1
h3=100
h4=100
%dx1=sqrt((2*(h1-2*(x.^2)/6)));
%dx2=sqrt((2*(h2+(sqrt(2)*x-sqrt(47/2)).^2)));
%dx3=sqrt((2*(h3-4*x)));
%dx4=sqrt((2*(h4+10*x-(x.^2)/2)));
dx1=sqrt((2*(h1+(x.^2)/2)));
dx2=sqrt((2*(h1-5*(x.^2)/4-8*x)));
dx3=sqrt((2*h1-4*x));
dx4=sqrt((2*(h1-(x.^2)/2+10*x)));
plot(x,dx1,x,-dx1)
plot(x,dx2,x,-dx2)
plot(x,dx3,x,-dx3)
plot(x,dx4,x,-dx4)
